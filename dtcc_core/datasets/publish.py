import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class PublishedFile:
    path: str
    original_filename: str | None
    size: int
    sha256: str
    media_type: str
    sniffed_media_type: str | None
    raw: Mapping[str, Any]

    @classmethod
    def from_response(cls, payload: Mapping[str, Any]) -> "PublishedFile":
        original_filename = payload.get("original_filename")
        return cls(
            path=str(payload["path"]),
            original_filename=(
                str(original_filename) if original_filename is not None else None
            ),
            size=int(payload["size"]),
            sha256=str(payload["sha256"]),
            media_type=str(payload["media_type"]),
            sniffed_media_type=(
                str(payload["sniffed_media_type"])
                if payload.get("sniffed_media_type") is not None
                else None
            ),
            raw=payload,
        )


@dataclass(frozen=True)
class DatasetPublication:
    dataset_key: str
    version_id: str
    version_number: int
    status: str
    owner: str
    manifest_sha256: str
    file_set_sha256: str
    files: tuple[PublishedFile, ...]
    upload_url: str | None = None
    manifest_url: str | None = None
    raw: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_response(
        cls, payload: Mapping[str, Any], *, upload_url: str | None = None
    ) -> "DatasetPublication":
        manifest_url = payload.get("manifest_url")
        return cls(
            dataset_key=str(payload["dataset_key"]),
            version_id=str(payload["version_id"]),
            version_number=int(payload["version_number"]),
            status=str(payload["status"]),
            owner=str(payload["owner"]),
            manifest_sha256=str(payload["manifest_sha256"]),
            file_set_sha256=str(payload["file_set_sha256"]),
            files=tuple(PublishedFile.from_response(item) for item in payload["files"]),
            upload_url=upload_url,
            manifest_url=str(manifest_url) if manifest_url is not None else None,
            raw=payload,
        )


class DatasetPublishError(Exception):
    def __init__(
        self,
        message: str,
        *,
        status_code: int | None = None,
        detail: Any = None,
        failure_class: str | None = None,
        retry_after: float | None = None,
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.detail = detail
        self.failure_class = failure_class
        self.retry_after = retry_after

    @property
    def is_transient(self) -> bool:
        return self.failure_class in {
            "connection",
            "timeout",
            "http_5xx",
            "rate_limited",
        }


class DatasetPublishConfigurationError(DatasetPublishError):
    pass


class DatasetPackageError(DatasetPublishError):
    pass


class DatasetUploadError(DatasetPublishError):
    pass


class DatasetUploadConflictError(DatasetUploadError):
    pass


class DatasetUploadInProgressError(DatasetUploadError):
    pass


class DatasetUploadRateLimitError(DatasetUploadError):
    pass


class DatasetUploadClient:
    def __init__(
        self,
        base_url: str,
        token: str,
        *,
        timeout: float | None = None,
        session: Any = None,
    ) -> None:
        self.base_url = str(base_url).rstrip("/")
        self.datasets_url = _normalize_datasets_url(base_url)
        self.token = token
        self.timeout = timeout
        self.session = session

    @classmethod
    def from_config(
        cls,
        upload_url: str | None = None,
        token: str | None = None,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
        session: Any = None,
    ) -> "DatasetUploadClient":
        config_env = os.environ if env is None else env
        resolved_url = (
            upload_url if upload_url is not None else config_env.get("DTCC_UPLOAD_URL")
        )
        resolved_token = (
            token if token is not None else config_env.get("DTCC_UPLOAD_TOKEN")
        )

        if not resolved_url or not resolved_url.strip():
            raise _config_error(
                "Missing dataset upload URL. Pass upload_url or set DTCC_UPLOAD_URL."
            )
        if not resolved_token or not resolved_token.strip():
            raise _config_error(
                "Missing dataset upload token. Pass token or set DTCC_UPLOAD_TOKEN."
            )

        return cls(resolved_url, resolved_token, timeout=timeout, session=session)


def build_publish_idempotency_key(
    dataset_key: str,
    manifest_path: str | Path,
    files: Sequence[str | Path],
    manifest: Mapping[str, Any],
) -> str:
    manifest_file = Path(manifest_path)
    if not manifest_file.is_file():
        raise _package_error(f"Manifest file does not exist: {manifest_file}")

    try:
        manifest_payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise _package_error(f"Unable to read manifest JSON: {manifest_file}") from error

    if not isinstance(manifest_payload, Mapping):
        raise _package_error("Manifest file must contain a JSON object.")
    if manifest_payload != manifest:
        raise _package_error(
            "Manifest mapping does not match manifest file content."
        )

    logical_path = manifest.get("file")
    manifest_file_path = manifest_payload.get("file")
    if not isinstance(logical_path, str) or not logical_path:
        raise _package_error("Manifest is missing required 'file' field.")
    if not isinstance(manifest_file_path, str) or not manifest_file_path:
        raise _package_error("Manifest file is missing required 'file' field.")
    if manifest_file_path != logical_path:
        raise _package_error(
            "Manifest mapping 'file' field does not match manifest file content."
        )

    if len(files) != 1:
        raise _package_error("Package must include exactly one file.")

    manifest_sha256 = _sha256_file(manifest_file)
    records = []
    for file in files:
        file_path = Path(file)
        if not file_path.is_file():
            raise _package_error(f"Package file does not exist: {file_path}")
        if file_path.name != logical_path:
            raise _package_error(
                f"Package file name {file_path.name!r} does not match manifest file "
                f"{logical_path!r}."
            )
        records.append(
            {
                "path": logical_path,
                "size": file_path.stat().st_size,
                "sha256": _sha256_file(file_path),
            }
        )

    records.sort(key=lambda record: record["path"])
    file_set_payload = json.dumps(
        records, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    file_set_sha256 = hashlib.sha256(file_set_payload).hexdigest()
    key_payload = json.dumps(
        {
            "dataset_key": dataset_key,
            "manifest_sha256": manifest_sha256,
            "file_set_sha256": file_set_sha256,
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")

    return f"dtcc-publish-v1:{hashlib.sha256(key_payload).hexdigest()}"


def _normalize_datasets_url(upload_url: str) -> str:
    base = upload_url.rstrip("/")
    if base.endswith("/v1/datasets"):
        return base
    if base.endswith("/v1"):
        return f"{base}/datasets"
    return f"{base}/v1/datasets"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _config_error(message: str) -> DatasetPublishConfigurationError:
    return DatasetPublishConfigurationError(message, failure_class="configuration")


def _package_error(message: str) -> DatasetPackageError:
    return DatasetPackageError(message, failure_class="invalid_package")
