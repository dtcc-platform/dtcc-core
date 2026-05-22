"""
Remote dataset descriptor and service discovery.

Enables transparent access to datasets hosted on remote DTCC services
(e.g., dtcc-sim running in Docker) via the Remote Dataset Protocol.
"""

import os
import json
import time
import logging
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pydantic import BaseModel, ConfigDict, ValidationError as PydanticValidationError

from .dataset import DatasetDescriptor
from .registry import register

logger = logging.getLogger(__name__)


class SSEStreamError(ConnectionError):
    """Raised when the SSE stream ends unexpectedly without a terminal event.

    This is a transport-category error that triggers polling fallback,
    distinct from RuntimeError which indicates a definitive job outcome
    (failed/cancelled) and should propagate immediately.
    """
    pass


class RemoteValidationError(ValueError):
    """Raised when a remote service returns a 422 validation error.

    Wraps the error details from the remote service's response.
    Atlas can catch this alongside pydantic.ValidationError to return 422.
    """

    def __init__(self, detail, status_code=422):
        self.detail = detail
        self.status_code = status_code
        super().__init__(f"Remote validation error: {detail}")

# Shared results directory for reading large results from shared volume
SHARED_RESULTS_DIR = os.environ.get("SHARED_RESULTS_DIR", "./data/shared-results")


class _PassthroughArgs(BaseModel):
    """Permissive model that accepts any kwargs. Real validation is server-side."""

    model_config = ConfigDict(extra="allow")


class RemoteDatasetDescriptor(DatasetDescriptor, register=False):
    """A DatasetDescriptor that delegates execution to a remote DTCC service."""

    def __init__(
        self,
        name: str,
        description: str,
        args_schema: Dict[str, Any],
        base_url: str,
        result_kind: str,
        supported_formats: List[str],
        source_service: str,
        timeout_hint: Optional[int] = None,
        data_category: Optional[str] = None,
    ):
        self.name = name
        self.description = description
        self._args_schema = args_schema
        self.base_url = base_url.rstrip("/")
        self.data_category = _remote_data_category(data_category, source_service)
        self.result_kind = result_kind
        self.supported_formats = supported_formats
        self.source_service = source_service
        self.timeout_hint = timeout_hint
        self.python_return_type = "tuple[bytes, str, str]"
        self.ArgsModel = _PassthroughArgs

    def show_options(self):
        """Return the raw JSON schema from discovery, not from _PassthroughArgs."""
        return self._args_schema

    def validate(self, kwargs):
        """Skip Pydantic validation but serialize Bounds objects for HTTP."""
        if "bounds" in kwargs:
            bounds = kwargs["bounds"]
            if hasattr(bounds, "tuple"):
                kwargs["bounds"] = bounds.tuple
        return kwargs

    def __call__(self, **kwargs):
        """Inject default format if missing, then delegate to build()."""
        if (
            "format" not in kwargs or kwargs["format"] is None
        ) and self.supported_formats:
            kwargs["format"] = self.supported_formats[0]
        validated = self.validate(kwargs)
        return self.build(validated)

    def build(self, validated_args, progress_callback=None, remote_info_callback=None):
        """Submit job to remote service, stream progress, return result.

        Returns:
            Tuple of (data_bytes, extension, content_type)
        """
        import httpx

        # 1. Submit job
        submit_url = f"{self.base_url}/api/v1/datasets/{self.name}/submit"
        resp = httpx.post(submit_url, json=validated_args, timeout=30)
        if resp.status_code == 422:
            error_detail = resp.json().get("detail", resp.text)
            raise RemoteValidationError(error_detail, status_code=422)
        resp.raise_for_status()
        task_data = resp.json()
        task_id = task_data["task_id"]

        # 2. Send remote task info back for cancellation support
        if remote_info_callback:
            remote_info_callback(
                {
                    "remote_task_id": task_id,
                    "cancel_url": (
                        f"{self.base_url}/api/v1/datasets/{self.name}"
                        f"/cancel/{task_id}"
                    ),
                }
            )

        # 3. Stream status via SSE, fall back to polling
        result_file = self._stream_status(task_id, progress_callback)

        # 4. Read result from shared volume (with path containment)
        safe_name = os.path.basename(result_file)
        if not safe_name or safe_name != result_file:
            raise RuntimeError(
                f"Remote service returned unsafe result_file path: {result_file!r}"
            )
        result_path = os.path.join(SHARED_RESULTS_DIR, safe_name)
        result_path = os.path.realpath(result_path)
        shared_dir = os.path.realpath(SHARED_RESULTS_DIR)
        if not result_path.startswith(shared_dir + os.sep):
            raise RuntimeError(
                f"Result path {result_path} escapes shared directory {shared_dir}"
            )

        with open(result_path, "rb") as f:
            data = f.read()

        # 5. Clean up shared volume file
        try:
            os.unlink(result_path)
        except OSError:
            logger.warning(f"Failed to clean up shared result: {result_path}")

        # 6. Determine extension and content type from result filename
        if result_file.endswith(".tar.gz"):
            extension = "tar.gz"
            content_type = "application/gzip"
        else:
            extension = (
                result_file.rsplit(".", 1)[-1] if "." in result_file else "bin"
            )
            content_type = "application/octet-stream"

        return data, extension, content_type

    def _stream_status(self, task_id, progress_callback):
        """Connect to SSE status stream, return result_file on completion."""
        import httpx

        status_url = (
            f"{self.base_url}/api/v1/datasets/{self.name}/status/{task_id}"
        )

        try:
            return self._stream_sse(status_url, progress_callback)
        except (httpx.TransportError, ConnectionError, OSError) as e:
            # Only fall back to polling on transport/connection errors.
            # httpx.TransportError is the base for all transport-level
            # exceptions (ConnectError, ReadError, RemoteProtocolError, etc.)
            # Terminal job states (failed, cancelled) are RuntimeError
            # and should propagate immediately, not trigger a retry.
            logger.warning(f"SSE connection failed, falling back to polling: {e}")
            return self._poll_status(status_url, progress_callback)

    def _stream_sse(self, status_url, progress_callback):
        """Read SSE events from status endpoint."""
        import httpx

        with httpx.stream(
            "GET",
            status_url,
            headers={"Accept": "text/event-stream"},
            timeout=httpx.Timeout(connect=10, read=None, write=10, pool=10),
        ) as response:
            response.raise_for_status()
            for line in response.iter_lines():
                if not line.startswith("data: "):
                    continue
                event = json.loads(line[6:])
                status = event.get("status")

                if status == "running" and progress_callback:
                    progress_callback(
                        {
                            "percent": event.get("progress", 0) * 100,
                            "message": event.get("message", ""),
                        }
                    )
                elif status == "completed":
                    return event["result_file"]
                elif status == "failed":
                    raise RuntimeError(
                        f"Remote job failed: {event.get('error', 'Unknown error')}"
                    )
                elif status == "cancelled":
                    raise RuntimeError("Remote job was cancelled")

        raise SSEStreamError("SSE stream ended without completion event")

    def _poll_status(self, status_url, progress_callback):
        """Poll status endpoint as fallback when SSE fails."""
        import httpx

        while True:
            resp = httpx.get(status_url, timeout=10)
            resp.raise_for_status()
            event = resp.json()
            status = event.get("status")

            if status == "running" and progress_callback:
                progress_callback(
                    {
                        "percent": event.get("progress", 0) * 100,
                        "message": event.get("message", ""),
                    }
                )
            elif status == "completed":
                return event["result_file"]
            elif status == "failed":
                raise RuntimeError(
                    f"Remote job failed: {event.get('error', 'Unknown error')}"
                )
            elif status == "cancelled":
                raise RuntimeError("Remote job was cancelled")

            time.sleep(2)


# Module-level cache of discovery responses for worker bootstrap
_cached_service_discoveries: Dict[str, Dict] = {}


def register_remote_service(base_url: str, timeout: int = 5) -> List[str]:
    """Query a remote service's discovery endpoint and register its datasets.

    Returns list of registered dataset names, or empty list on failure.
    Graceful degradation: catches all errors so atlas starts regardless.
    """
    import httpx

    base_url = base_url.rstrip("/")
    try:
        resp = httpx.get(f"{base_url}/api/v1/datasets", timeout=timeout)
        resp.raise_for_status()
        service_info = resp.json()

        # Cache for worker bootstrap (no HTTP in child processes)
        _cached_service_discoveries[base_url] = service_info

        registered = []
        for name, meta in service_info["datasets"].items():
            descriptor = RemoteDatasetDescriptor(
                name=meta["name"],
                description=meta["description"],
                args_schema=meta["args_schema"],
                base_url=base_url,
                result_kind=meta["result_kind"],
                supported_formats=meta["supported_formats"],
                source_service=service_info["service"],
                timeout_hint=meta.get("timeout_hint"),
                data_category=meta.get("data_category"),
            )
            register(name, descriptor)
            registered.append(name)
        return registered
    except Exception as e:
        logger.warning(f"Failed to discover service at {base_url}: {e}")
        return []


def register_remote_descriptors_from_cache(cached_discoveries: Dict[str, Dict]):
    """Register remote datasets from pre-fetched discovery data.

    Called by worker processes to avoid HTTP calls in child processes.
    """
    for base_url, service_info in cached_discoveries.items():
        for name, meta in service_info["datasets"].items():
            descriptor = RemoteDatasetDescriptor(
                name=meta["name"],
                description=meta["description"],
                args_schema=meta["args_schema"],
                base_url=base_url,
                result_kind=meta["result_kind"],
                supported_formats=meta["supported_formats"],
                source_service=service_info["service"],
                timeout_hint=meta.get("timeout_hint"),
                data_category=meta.get("data_category"),
            )
            register(name, descriptor)


def get_cached_discoveries() -> Dict[str, Dict]:
    """Return cached discovery data for passing to worker processes."""
    return dict(_cached_service_discoveries)


def _remote_data_category(
    data_category: Optional[str],
    source_service: str,
) -> str:
    if data_category:
        return data_category
    if source_service.lower() == "dtcc-sim":
        return "simulation"
    return "remote"
