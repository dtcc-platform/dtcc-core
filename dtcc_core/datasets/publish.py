import hashlib
import json
import os
import re
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import requests


MANIFEST_V2_SCHEMA_VERSION = "dtcc-dataset-manifest-v2"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
WINDOWS_DRIVE_RE = re.compile(r"^[A-Za-z]:")
WINDOWS_RESERVED_BASENAMES = {
    "con",
    "prn",
    "aux",
    "nul",
    "com1",
    "com2",
    "com3",
    "com4",
    "com5",
    "com6",
    "com7",
    "com8",
    "com9",
    "lpt1",
    "lpt2",
    "lpt3",
    "lpt4",
    "lpt5",
    "lpt6",
    "lpt7",
    "lpt8",
    "lpt9",
}


@dataclass(frozen=True)
class PublishedFile:
    path: str
    original_filename: str | None
    size: int
    sha256: str
    media_type: str
    sniffed_media_type: str
    raw: Mapping[str, Any]

    @classmethod
    def from_response(cls, payload: Mapping[str, Any]) -> "PublishedFile":
        original_filename = payload.get("original_filename")
        sniffed_media_type = payload["sniffed_media_type"]
        if sniffed_media_type is None:
            raise ValueError("Published file response missing sniffed_media_type")
        return cls(
            path=str(payload["path"]),
            original_filename=(
                str(original_filename) if original_filename is not None else None
            ),
            size=int(payload["size"]),
            sha256=str(payload["sha256"]),
            media_type=str(payload["media_type"]),
            sniffed_media_type=str(sniffed_media_type),
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


class DatasetUploadInProgressError(DatasetUploadConflictError):
    pass


class DatasetUploadRateLimitError(DatasetUploadError):
    pass


@dataclass(frozen=True)
class _PackageUploadFile:
    logical_path: str
    path: Path
    media_type: str
    size: int
    sha256: str


@dataclass(frozen=True)
class _ValidatedUploadPackage:
    manifest: dict[str, Any]
    files: tuple[_PackageUploadFile, ...]
    idempotency_prefix: str


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
        self.session = requests.Session() if session is None else session

    @classmethod
    def from_config(
        cls,
        *,
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

    def upload_package(
        self,
        *,
        dataset_key: str,
        manifest_path: str | Path,
        files: Sequence[str | Path],
        manifest: Mapping[str, Any] | None = None,
        idempotency_key: str | None = None,
    ) -> DatasetPublication:
        package = _validate_upload_package(
            manifest_path,
            files,
            manifest=manifest,
        )
        key = idempotency_key or build_publish_idempotency_key(
            dataset_key=dataset_key,
            manifest_path=manifest_path,
            files=files,
            manifest=package.manifest,
        )
        headers = {
            "Authorization": f"Bearer {self.token}",
            "Idempotency-Key": key,
        }

        try:
            with ExitStack() as stack:
                manifest_handle = stack.enter_context(Path(manifest_path).open("rb"))
                file_handles = [
                    stack.enter_context(upload_file.path.open("rb"))
                    for upload_file in package.files
                ]
                multipart_files: Any
                if package.idempotency_prefix == "dtcc-publish-v1":
                    upload_file = package.files[0]
                    multipart_files = {
                        "manifest": (
                            "manifest.json",
                            manifest_handle,
                            "application/json",
                        ),
                        "files": (
                            upload_file.logical_path,
                            file_handles[0],
                            upload_file.media_type,
                        ),
                    }
                else:
                    multipart_files = [
                        (
                            "manifest",
                            ("manifest.json", manifest_handle, "application/json"),
                        ),
                        *[
                            (
                                "files",
                                (
                                    upload_file.logical_path,
                                    file_handle,
                                    upload_file.media_type,
                                ),
                            )
                            for upload_file, file_handle in zip(
                                package.files, file_handles
                            )
                        ],
                    ]
                response = self.session.post(
                    self.datasets_url,
                    data={"dataset_key": dataset_key},
                    files=multipart_files,
                    headers=headers,
                    timeout=self.timeout,
                )
        except requests.Timeout as error:
            detail = _sanitize_error_text(str(error), token=self.token)
            raise DatasetUploadError(
                f"Dataset upload timed out: {detail}",
                detail=detail,
                failure_class="timeout",
            ) from None
        except requests.ConnectionError as error:
            detail = _sanitize_error_text(str(error), token=self.token)
            raise DatasetUploadError(
                f"Dataset upload connection failed: {detail}",
                detail=detail,
                failure_class="connection",
            ) from None
        except requests.RequestException as error:
            detail = _sanitize_error_text(str(error), token=self.token)
            raise DatasetUploadError(
                f"Dataset upload request failed: {detail}",
                detail=detail,
                failure_class="request",
            ) from None

        if response.status_code < 200 or response.status_code >= 300:
            raise _upload_error_for_response(response, token=self.token)

        try:
            payload = response.json()
        except ValueError as error:
            detail = _response_detail(response, token=self.token)
            raise _malformed_success_error(response, detail=detail) from None
        if not isinstance(payload, Mapping):
            detail = _response_detail(response, token=self.token)
            raise _malformed_success_error(response, detail=detail)

        try:
            return DatasetPublication.from_response(payload, upload_url=self.base_url)
        except (AttributeError, KeyError, TypeError, ValueError) as error:
            detail = _sanitize_error_text(str(error), token=self.token)
            raise _malformed_success_error(response, detail=detail) from None


def build_publish_idempotency_key(
    *,
    dataset_key: str,
    manifest_path: str | Path,
    files: Sequence[str | Path],
    manifest: Mapping[str, Any],
) -> str:
    package = _validate_upload_package(
        manifest_path,
        files,
        manifest=manifest,
    )
    manifest_sha256 = _sha256_file(Path(manifest_path))
    file_set_sha256 = _file_set_sha256(package.files)
    key_payload = json.dumps(
        {
            "dataset_key": dataset_key,
            "manifest_sha256": manifest_sha256,
            "file_set_sha256": file_set_sha256,
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")

    return f"{package.idempotency_prefix}:{hashlib.sha256(key_payload).hexdigest()}"


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


def _read_manifest(manifest_path: Path) -> dict[str, Any]:
    if not manifest_path.is_file():
        raise _package_error(f"Manifest file does not exist: {manifest_path}")
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise _package_error(
            f"Unable to read manifest JSON: {manifest_path}"
        ) from error
    if not isinstance(payload, dict):
        raise _package_error("Manifest file must contain a JSON object.")
    return payload


def _validate_upload_package(
    manifest_path: str | Path,
    files: Sequence[str | Path],
    *,
    manifest: Mapping[str, Any] | None = None,
) -> _ValidatedUploadPackage:
    manifest_file = Path(manifest_path)
    manifest_payload = _read_manifest(manifest_file)
    if manifest is not None and dict(manifest) != manifest_payload:
        raise _package_error("Manifest mapping does not match manifest file content.")

    schema_version = manifest_payload.get("schema_version")
    if schema_version == MANIFEST_V2_SCHEMA_VERSION:
        return _validate_v2_package(manifest_file, files, manifest_payload)
    if schema_version is not None:
        raise _package_error(
            f"Unsupported manifest schema_version: {schema_version!r}; "
            f"expected {MANIFEST_V2_SCHEMA_VERSION!r}"
        )
    if "artifacts" in manifest_payload or "identity" in manifest_payload:
        raise _package_error(
            f"Dataset Manifest v2 packages must declare schema_version={MANIFEST_V2_SCHEMA_VERSION!r}"
        )
    manifest_payload, file_path, logical_name = _validate_single_file_package(
        manifest_file,
        files,
        manifest=manifest_payload,
    )
    media_type = str(manifest_payload.get("media_type") or "application/octet-stream")
    return _ValidatedUploadPackage(
        manifest=manifest_payload,
        files=(
            _PackageUploadFile(
                logical_path=logical_name,
                path=file_path,
                media_type=media_type,
                size=file_path.stat().st_size,
                sha256=_sha256_file(file_path),
            ),
        ),
        idempotency_prefix="dtcc-publish-v1",
    )


def _validate_logical_filename(name: str) -> str:
    if (
        not name
        or name.startswith(".")
        or name.startswith("/")
        or "/" in name
        or "\\" in name
    ):
        raise _package_error(f"invalid manifest file name: {name!r}")
    return name


def _validate_single_file_package(
    manifest_path: str | Path,
    files: Sequence[str | Path],
    *,
    manifest: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], Path, str]:
    manifest_payload = _read_manifest(Path(manifest_path))
    if manifest is not None and dict(manifest) != manifest_payload:
        raise _package_error("Manifest mapping does not match manifest file content.")
    if len(files) != 1:
        raise _package_error(
            "V1 dataset publishing supports single-file packages only."
        )

    logical_file = manifest_payload.get("file")
    if not isinstance(logical_file, str):
        raise _package_error("Manifest is missing required 'file' field.")
    logical_name = _validate_logical_filename(logical_file)

    manifest_files = manifest_payload.get("files")
    if manifest_files is not None and manifest_files != [logical_name]:
        raise _package_error("manifest.files must match the manifest file field.")

    file_path = Path(files[0])
    if not file_path.is_file():
        raise _package_error(f"Package file does not exist: {file_path}")
    if file_path.name != logical_name:
        raise _package_error(
            f"Package file name {file_path.name!r} does not match manifest file "
            f"{logical_name!r}."
        )
    return manifest_payload, file_path, logical_name


def _validate_v2_package(
    manifest_path: Path,
    files: Sequence[str | Path],
    manifest: Mapping[str, Any],
) -> _ValidatedUploadPackage:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise _package_error("Dataset Manifest v2 must contain a non-empty artifacts list.")

    expected_artifacts: dict[str, Mapping[str, Any]] = {}
    for index, artifact in enumerate(artifacts):
        if not isinstance(artifact, Mapping):
            raise _package_error(f"Manifest artifact at index {index} must be an object.")
        artifact_path = _required_artifact_string(artifact, "path", index)
        logical_path = _validate_package_path(artifact_path)
        for field_name in ("role", "format", "media_type", "data_kind"):
            _required_artifact_string(artifact, field_name, index)
        if logical_path in expected_artifacts:
            raise _package_error(f"Duplicate artifact path in manifest: {logical_path}")
        expected_artifacts[logical_path] = artifact

    package_root = manifest_path.parent.resolve()
    provided: dict[str, Path] = {}
    for file in files:
        file_path = Path(file)
        if not file_path.is_file():
            raise _package_error(f"Package file does not exist: {file_path}")
        logical_path = _relative_package_path(package_root, file_path)
        if logical_path not in expected_artifacts:
            raise _package_error(f"Unexpected package artifact file: {logical_path}")
        if logical_path in provided:
            raise _package_error(f"Duplicate package artifact file: {logical_path}")
        provided[logical_path] = file_path

    missing_paths = [path for path in expected_artifacts if path not in provided]
    if missing_paths:
        raise _package_error(f"Missing package artifact file: {missing_paths[0]}")
    if len(provided) != len(expected_artifacts):
        raise _package_error("Package artifact files do not match manifest artifacts.")

    upload_files: list[_PackageUploadFile] = []
    for logical_path, artifact in expected_artifacts.items():
        file_path = provided[logical_path]
        size = file_path.stat().st_size
        declared_size = artifact.get("size")
        if declared_size is not None:
            if not isinstance(declared_size, int) or declared_size < 0:
                raise _package_error(
                    f"Manifest artifact {logical_path!r} has invalid size."
                )
            if size != declared_size:
                raise _package_error(
                    f"Package artifact size does not match manifest for {logical_path}."
                )
        sha256 = _sha256_file(file_path)
        declared_sha256 = artifact.get("sha256")
        if declared_sha256 is not None:
            if not isinstance(declared_sha256, str) or not SHA256_RE.fullmatch(
                declared_sha256
            ):
                raise _package_error(
                    f"Manifest artifact {logical_path!r} has invalid sha256."
                )
            if sha256 != declared_sha256:
                raise _package_error(
                    f"Package artifact sha256 does not match manifest for {logical_path}."
                )
        media_type = str(artifact["media_type"]).strip()
        if not media_type:
            raise _package_error(
                f"Manifest artifact {logical_path!r} has invalid media_type."
            )
        upload_files.append(
            _PackageUploadFile(
                logical_path=logical_path,
                path=file_path,
                media_type=media_type,
                size=size,
                sha256=sha256,
            )
        )

    return _ValidatedUploadPackage(
        manifest=dict(manifest),
        files=tuple(upload_files),
        idempotency_prefix="dtcc-publish-v2",
    )


def _required_artifact_string(
    artifact: Mapping[str, Any], field_name: str, index: int
) -> str:
    value = artifact.get(field_name)
    if not isinstance(value, str) or not value.strip():
        raise _package_error(
            f"Manifest artifact at index {index} is missing required {field_name!r}."
        )
    return value.strip()


def _validate_package_path(value: str) -> str:
    if not isinstance(value, str):
        raise _package_error("Invalid artifact path")
    normalized = value
    if not normalized or normalized in {".", ".."}:
        raise _package_error("Invalid artifact path")
    if normalized.startswith("/") or normalized.startswith("."):
        raise _package_error("Invalid artifact path")
    if WINDOWS_DRIVE_RE.match(normalized):
        raise _package_error("Invalid artifact path")
    if "\\" in normalized or "\x00" in normalized or "//" in normalized:
        raise _package_error("Invalid artifact path")
    if any(ord(char) < 32 or ord(char) == 127 for char in normalized):
        raise _package_error("Invalid artifact path")

    parts = normalized.split("/")
    for part in parts:
        if not part or part in {".", ".."} or part.startswith("."):
            raise _package_error("Invalid artifact path")
        basename = part.rsplit(".", 1)[0].lower()
        if basename in WINDOWS_RESERVED_BASENAMES:
            raise _package_error("Invalid artifact path")
    return normalized


def _relative_package_path(package_root: Path, file_path: Path) -> str:
    try:
        relative = file_path.resolve().relative_to(package_root).as_posix()
    except ValueError as error:
        raise _package_error(
            "V2 package artifact files must be under the manifest directory."
        ) from error
    return _validate_package_path(relative)


def _file_set_sha256(files: Sequence[_PackageUploadFile]) -> str:
    records = sorted(
        (
            {
                "path": upload_file.logical_path,
                "size": upload_file.size,
                "sha256": upload_file.sha256,
            }
            for upload_file in files
        ),
        key=lambda record: record["path"],
    )
    payload = json.dumps(records, separators=(",", ":"), sort_keys=True).encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def _retry_after(response: Any) -> float | None:
    value = response.headers.get("Retry-After")
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _sanitize_error_text(text: Any, *, token: str | None) -> str:
    rendered = str(text)
    if token:
        rendered = rendered.replace(token, "<redacted>")
    return rendered


def _response_detail(response: Any, *, token: str | None) -> str:
    try:
        payload = response.json()
    except ValueError:
        return _sanitize_error_text(response.text, token=token)

    if isinstance(payload, Mapping) and "detail" in payload:
        detail = payload["detail"]
        if isinstance(detail, str):
            return _sanitize_error_text(detail, token=token)
        detail_text = json.dumps(detail, sort_keys=True, default=str)
        return _sanitize_error_text(detail_text[:200], token=token)

    try:
        detail_text = json.dumps(payload, sort_keys=True, default=str)
    except (TypeError, ValueError):
        detail_text = response.text
    return _sanitize_error_text(detail_text[:200], token=token)


def _upload_error_for_response(
    response: Any, *, token: str | None
) -> DatasetUploadError:
    status_code = int(response.status_code)
    detail = _response_detail(response, token=token)
    message = f"Dataset upload failed with HTTP {status_code}: {detail}"

    if status_code == 409 and detail == "Idempotency-Key is already in progress":
        return DatasetUploadInProgressError(
            message,
            status_code=status_code,
            detail=detail,
            failure_class="conflict",
        )
    if status_code == 409:
        return DatasetUploadConflictError(
            message,
            status_code=status_code,
            detail=detail,
            failure_class="conflict",
        )
    if status_code == 429:
        return DatasetUploadRateLimitError(
            message,
            status_code=status_code,
            detail=detail,
            failure_class="rate_limited",
            retry_after=_retry_after(response),
        )
    if 400 <= status_code < 500:
        return DatasetUploadError(
            message,
            status_code=status_code,
            detail=detail,
            failure_class="http_4xx",
        )

    # upload_package always sends an idempotency key, so 5xx responses can be retried.
    if status_code >= 500:
        return DatasetUploadError(
            message,
            status_code=status_code,
            detail=detail,
            failure_class="http_5xx",
        )

    return DatasetUploadError(
        message,
        status_code=status_code,
        detail=detail,
        failure_class="request",
    )


def _malformed_success_error(response: Any, *, detail: str) -> DatasetUploadError:
    return DatasetUploadError(
        f"Dataset upload returned malformed success response: {detail}",
        status_code=int(response.status_code),
        detail=detail,
        failure_class="request",
    )
