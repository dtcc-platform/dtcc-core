"""Legacy v2 and canonical v3 Dataset Package I/O."""

from __future__ import annotations

import hashlib
from inspect import signature
import json
import os
import shutil
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path

from dtcc_core.datasets.dataset import DatasetDescriptor
from dtcc_core.datasets.schema import DatasetArtifact, DatasetManifest, DatasetContext


CANONICAL_MANIFEST_VERSION = "dtcc-dataset-manifest-v3"
CANONICAL_ARTIFACT = "artifacts/model.dtcc"


@dataclass(frozen=True)
class DatasetPackage:
    """Object-first Dataset export package."""

    path: Path
    manifest_path: Path
    manifest: DatasetManifest
    artifacts: tuple[DatasetArtifact, ...]
    files: tuple[Path, ...]
    package_format: str

    def publish(
        self,
        *,
        dataset_key: str,
        uploader=None,
        upload_url: str | None = None,
        token: str | None = None,
        idempotency_key: str | None = None,
    ):
        """Publish an artifact or canonical model package through dtcc-upload."""
        from dtcc_core.datasets.publish import DatasetUploadClient

        resolved_uploader = uploader or DatasetUploadClient.from_config(
            upload_url=upload_url,
            token=token,
        )
        manifest_payload = self.manifest.model_dump(mode="json")

        if self.package_format == "directory":
            artifact_files = self._directory_artifact_files()
            return resolved_uploader.upload_package(
                dataset_key=dataset_key,
                manifest_path=self.manifest_path,
                files=artifact_files,
                manifest=manifest_payload,
                idempotency_key=idempotency_key,
            )

        if self.package_format == "dtccpkg":
            with tempfile.TemporaryDirectory() as tmpdir:
                package_dir = Path(tmpdir) / _safe_stem(self.path.stem)
                manifest_path, artifact_files = self._extract_archive_package(
                    package_dir
                )
                return resolved_uploader.upload_package(
                    dataset_key=dataset_key,
                    manifest_path=manifest_path,
                    files=artifact_files,
                    manifest=manifest_payload,
                    idempotency_key=idempotency_key,
                )

        raise ValueError(f"Unsupported Dataset package format: {self.package_format}")

    def _directory_artifact_files(self) -> tuple[Path, ...]:
        package_dir = self.manifest_path.parent
        artifact_files = tuple(package_dir / artifact.path for artifact in self.artifacts)
        missing = [path for path in artifact_files if not path.is_file()]
        if missing:
            raise ValueError(f"Dataset package artifact file does not exist: {missing[0]}")
        return artifact_files

    def _extract_archive_package(self, package_dir: Path) -> tuple[Path, tuple[Path, ...]]:
        if not self.path.is_file():
            raise ValueError(f"Dataset package archive does not exist: {self.path}")
        package_dir.mkdir(parents=True, exist_ok=False)
        with zipfile.ZipFile(self.path) as archive:
            if self.manifest.schema_version == CANONICAL_MANIFEST_VERSION:
                from .publish import _validate_package_path
                entries = archive.infolist()
                expected = {"manifest.json", *(_validate_package_path(a.path) for a in self.artifacts)}
                names = [entry.filename for entry in entries]
                if len(names) != len(set(names)) or set(names) != expected:
                    raise ValueError("Canonical archive members must match the manifest exactly")
                if archive.getinfo('manifest.json').file_size > 4 * 1024 * 1024:
                    raise ValueError("Package manifest exceeds 4 MiB limit")
                for artifact in self.artifacts:
                    if archive.getinfo(artifact.path).file_size != artifact.size:
                        raise ValueError("Artifact size does not match manifest")
            manifest_path = package_dir / "manifest.json"
            _extract_zip_member(archive, "manifest.json", manifest_path)
            artifact_files = []
            for artifact in self.artifacts:
                artifact_path = package_dir / artifact.path
                _extract_zip_member(archive, artifact.path, artifact_path)
                artifact_files.append(artifact_path)
        return manifest_path, tuple(artifact_files)


def export_model_package(
    obj,
    path: str | Path,
    *,
    format: str | None = None,
    canonical: bool = False,
    validate_schema: bool = True,
) -> DatasetPackage:
    """Export a dataset realization; canonical=True opts into the v3 contract.

    In canonical mode, format selects an optional supplemental artifact. With no
    format, only the canonical model is written. Legacy v2 export remains the
    default until its producers and consumers have migrated.
    Canonical exports validate the standard semantic schema by default. Set
    validate_schema=False to bypass semantics while retaining native admission.
    """
    context = getattr(obj, "dataset_context", None)
    if context is None:
        raise ValueError(
            "Cannot export Dataset package: this object has no DatasetContext."
        )
    if canonical:
        return _export_canonical_package(obj, Path(path), _normalize_format(format), validate_schema)

    if validate_schema is not True:
        raise ValueError('validate_schema applies to canonical packages; pass canonical=True')

    target_path = Path(path)
    package_format = "dtccpkg" if target_path.suffix.lower() == ".dtccpkg" else "directory"
    artifact_format = _normalize_format(format) or _default_format(obj)
    if artifact_format is None:
        raise ValueError(
            f"Cannot infer a safe export format for {type(obj).__name__}. "
            "Pass format= explicitly when a serializer is available."
        )

    if package_format == "directory":
        package_dir = target_path
        if package_dir.exists() and package_dir.is_file():
            raise ValueError(f"Dataset package path is a file: {package_dir}")
        if package_dir.exists() and any(package_dir.iterdir()):
            raise ValueError(
                f"Dataset package directory already exists and is not empty: "
                f"{package_dir}"
            )
        return _write_directory_package(
            obj,
            package_dir,
            artifact_format=artifact_format,
            package_path=package_dir,
            package_format=package_format,
        )

    target_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        package_dir = Path(tmpdir) / _safe_stem(target_path.stem)
        package = _write_directory_package(
            obj,
            package_dir,
            artifact_format=artifact_format,
            package_path=target_path,
            package_format=package_format,
        )
        with zipfile.ZipFile(target_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for file_path in package.files:
                zf.write(file_path, file_path.relative_to(package_dir).as_posix())
        return DatasetPackage(
            path=target_path,
            manifest_path=Path("manifest.json"),
            manifest=package.manifest,
            artifacts=package.artifacts,
            files=(target_path,),
            package_format=package_format,
        )


def _write_directory_package(
    obj,
    package_dir: Path,
    *,
    artifact_format: str | None,
    package_path: Path,
    package_format: str,
    canonical_data: bytes | None = None,
) -> DatasetPackage:
    artifact_dir = package_dir / "artifacts"
    artifact_dir.mkdir(parents=True, exist_ok=True)

    context = obj.dataset_context
    canonical_path = package_dir / CANONICAL_ARTIFACT
    if canonical_data is not None:
        canonical_path.write_bytes(canonical_data)
    artifact_path = canonical_path
    if artifact_format is not None:
        extension = DatasetDescriptor.format_extension(artifact_format)
        artifact_stem = _artifact_stem(obj, context)
        artifact_path = artifact_dir / f"{artifact_stem}.{extension}"
        if canonical_data is not None and artifact_path == canonical_path:
            raise ValueError("Supplemental artifact collides with the canonical model path")
        _write_artifact(obj, artifact_path, artifact_format, crs=_crs_value(obj))

    artifact_files = _artifact_files(artifact_dir, artifact_path)
    artifacts = tuple(
        _artifact_metadata(
            file_path,
            package_dir=package_dir,
            requested_format=artifact_format if file_path == artifact_path else None,
            obj=obj,
        )
        for file_path in artifact_files
    )
    manifest = context.manifest(artifacts=list(artifacts))
    if canonical_data is not None:
        from ..model.exchange import VERSION
        manifest.schema_version = CANONICAL_MANIFEST_VERSION
        manifest.health = context.health
        manifest.warnings = list(context.warnings)
        for artifact in artifacts:
            if artifact.path == CANONICAL_ARTIFACT:
                artifact.role = "canonical_model"
                artifact.format = "dtcc"
                artifact.media_type = "application/vnd.dtcc.model+protobuf"
                artifact.data_kind = "model"
                artifact.model_type = type(obj).__name__
                artifact.model_schema_version = VERSION
                artifact.crs = getattr(getattr(obj, 'transform', None), 'srs', '') or None
                # Object/geometry transforms may differ; do not invent aggregate
                # global bounds from the request or from local cached bounds.
                artifact.bounds = None
            else:
                artifact.role = "derived"
                artifact.derived_from = CANONICAL_ARTIFACT
                artifact.bounds = None
                artifact.crs = None
    manifest_path = package_dir / "manifest.json"
    manifest_path.write_text(manifest.model_dump_json(indent=2) + "\n", encoding="utf-8")

    return DatasetPackage(
        path=package_path,
        manifest_path=manifest_path,
        manifest=manifest,
        artifacts=artifacts,
        files=(manifest_path, *artifact_files),
        package_format=package_format,
    )


def _write_artifact(
    obj,
    path: Path,
    artifact_format: str,
    *,
    crs: str | None = None,
) -> None:
    fmt = _normalize_format(artifact_format)
    write_artifact = getattr(obj, "write_artifact", None)
    if callable(write_artifact):
        write_artifact(path, format=fmt)
        return

    if fmt == "geojson" and hasattr(obj, "to_geojson"):
        to_geojson = obj.to_geojson
        parameters = signature(to_geojson).parameters
        if "crs" in parameters:
            payload = to_geojson(crs=crs)
        else:
            payload = to_geojson()
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return

    try:
        import dtcc_core.io  # noqa: F401  Registers save methods lazily.
    except Exception as exc:  # pragma: no cover - defensive import guard
        raise RuntimeError(f"Could not import dtcc_core.io serializers: {exc}") from exc

    save = getattr(obj, "save", None)
    if save is None:
        raise ValueError(
            f"Object type {type(obj).__name__} does not support Dataset v2 "
            f"artifact export in format {fmt!r}."
        )
    try:
        save(path)
    except Exception as exc:
        raise ValueError(
            f"Could not export {type(obj).__name__} as {fmt!r}: {exc}"
        ) from exc


def _artifact_metadata(
    path: Path,
    *,
    package_dir: Path,
    requested_format: str | None,
    obj,
) -> DatasetArtifact:
    fmt = requested_format or _format_from_path(path)
    relative_path = path.relative_to(package_dir).as_posix()
    return DatasetArtifact(
        path=relative_path,
        role="primary",
        format=fmt,
        media_type=DatasetDescriptor.format_media_type(fmt),
        data_kind=DatasetDescriptor.format_kind(fmt),
        crs=_crs_value(obj),
        bounds=_bounds_value(obj),
        size=path.stat().st_size,
        sha256=_sha256_file(path),
    )


def _artifact_files(artifact_dir: Path, primary_path: Path) -> tuple[Path, ...]:
    files = sorted(path for path in artifact_dir.iterdir() if path.is_file())
    if primary_path not in files:
        files.insert(0, primary_path)
    return tuple(files)


def _format_from_path(path: Path) -> str:
    suffixes = [suffix[1:].lower() for suffix in path.suffixes if suffix]
    if len(suffixes) >= 2 and suffixes[-2:] == ["json", "zip"]:
        return "json.zip"
    if not suffixes:
        return "unknown"
    return suffixes[-1]


def _default_format(obj) -> str | None:
    default_artifact_format = getattr(obj, "default_artifact_format", None)
    if callable(default_artifact_format):
        return _normalize_format(default_artifact_format())
    if isinstance(default_artifact_format, str):
        return _normalize_format(default_artifact_format)

    type_name = type(obj).__name__
    defaults = {
        "City": "json",
        "Mesh": "vtu",
        "VolumeMesh": "vtu",
        "PointCloud": "dtcc",
        "Raster": "tif",
        "FootprintCollection": "geojson",
        "CalibrationGrid": "geojson",
    }
    return defaults.get(type_name)


def _artifact_stem(obj, context) -> str:
    """Return a safe package artifact stem for a dataset-produced object.

    Object-first export prefers product/model-specific names when available:
    ``dataset_artifact_stem()``, then ``artifact_stem``, then ``name``. The
    dataset identity remains authoritative in the manifest and is used as the
    fallback stem. All candidates are sanitized through ``_safe_stem`` before
    being used under ``artifacts/``.
    """
    candidate = None
    dataset_artifact_stem = getattr(obj, "dataset_artifact_stem", None)
    if callable(dataset_artifact_stem):
        candidate = dataset_artifact_stem()
    elif isinstance(getattr(obj, "artifact_stem", None), str):
        candidate = obj.artifact_stem
    elif isinstance(getattr(obj, "name", None), str):
        candidate = obj.name
    return _safe_stem(candidate or context.identity.name)


def _bounds_value(obj) -> list[float] | None:
    context = getattr(obj, "dataset_context", None)
    if context is not None and context.request.bounds is not None:
        return list(context.request.bounds)

    bounds = getattr(obj, "bounds", None)
    if bounds is None:
        return None
    try:
        return [
            float(bounds.xmin),
            float(bounds.ymin),
            float(bounds.xmax),
            float(bounds.ymax),
        ]
    except AttributeError:
        return None


def _crs_value(obj) -> str | None:
    context = getattr(obj, "dataset_context", None)
    if context is not None and context.metadata.crs:
        return str(context.metadata.crs[0])

    transform = getattr(obj, "transform", None)
    srs = getattr(transform, "srs", None)
    return str(srs) if srs else None


def _safe_stem(value: str) -> str:
    stem = str(value).strip().replace(" ", "_").replace("-", "_").lower()
    sanitized = "".join(char if char.isalnum() or char == "_" else "_" for char in stem)
    while "__" in sanitized:
        sanitized = sanitized.replace("__", "_")
    return sanitized.strip("_") or "dataset"


def _normalize_format(format: str | None) -> str | None:
    if format is None:
        return None
    return str(format).lower().lstrip(".")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _extract_zip_member(archive: zipfile.ZipFile, member: str, target: Path) -> None:
    try:
        info = archive.getinfo(member)
    except KeyError as error:
        raise ValueError(f"Dataset package archive is missing {member!r}.") from error
    if info.is_dir():
        raise ValueError(f"Dataset package archive member is a directory: {member!r}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with archive.open(info) as source, target.open("wb") as destination:
        shutil.copyfileobj(source, destination)


def _export_canonical_package(obj, path, supplemental_format, validate_schema):
    """Stage the existing package writer and publish the local snapshot atomically."""
    from ..model import exchange
    if not isinstance(obj.dataset_context, DatasetContext):
        raise ValueError("Canonical package export requires a DatasetContext")
    from ..model.object.object import _validate_attributes
    _validate_attributes(obj.dataset_context.model_dump(mode='python'), max_depth=exchange.MAX_DEPTH)
    data = exchange.dumps(obj, validate_schema=validate_schema)  # Fail before touching the destination.
    if supplemental_format == 'dtcc':
        supplemental_format = None
    archive = path.suffix.lower() == '.dtccpkg'
    if not archive and path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError("Canonical package directory must be absent or empty")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent, prefix='.dtcc-package-') as temporary:
        stage = Path(temporary) / 'package'
        package = _write_directory_package(
            obj, stage, artifact_format=supplemental_format, package_path=stage,
            package_format='directory', canonical_data=data,
        )
        if (len(package.artifacts) > 9999
                or package.manifest_path.stat().st_size > 4 * 1024 * 1024):
            raise ValueError("Canonical package exceeds artifact or manifest limits")
        if archive:
            staged_archive = Path(temporary) / 'package.dtccpkg'
            with zipfile.ZipFile(staged_archive, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
                for file in package.files:
                    zf.write(file, file.relative_to(stage).as_posix())
            os.replace(staged_archive, path)
            return DatasetPackage(path, Path('manifest.json'), package.manifest,
                                  package.artifacts, (path,), 'dtccpkg')
        os.replace(stage, path)
    manifest_path = path / 'manifest.json'
    files = (manifest_path, *(path / artifact.path for artifact in package.artifacts))
    return DatasetPackage(path, manifest_path, package.manifest, package.artifacts, files, 'directory')


def load_model_package(path, *, validate_schema=True, max_bytes=None):
    """Read a canonical v3 directory/archive and restore its model and Context.

    All declared artifacts are integrity-checked. No archive extraction or network
    access occurs. Legacy v2 packages require their existing consumers.
    validate_schema=False bypasses semantic evaluation, not integrity/admission.
    max_bytes optionally bounds the total uncompressed artifact bytes before
    reading artifacts. The default has no package byte cap; the native model
    must still fit in one Protobuf message smaller than 2 GiB. Applications
    accepting untrusted packages should supply a budget appropriate to them.
    """
    from contextlib import ExitStack
    from ..model import exchange
    from .publish import _validate_package_path

    if max_bytes is not None and (type(max_bytes) is not int or max_bytes < 1):
        raise ValueError('max_bytes must be a positive integer or None')
    path = Path(path)
    with ExitStack() as stack:
        if path.is_dir():
            root = path.resolve()

            def open_member(name):
                member = (root / _validate_package_path(name)).resolve()
                if not member.is_relative_to(root):
                    raise ValueError("Package artifact resolves outside its directory")
                return member.open('rb')

            archive_names = None
        else:
            archive = stack.enter_context(zipfile.ZipFile(path))
            infos = archive.infolist()
            archive_names = [entry.filename for entry in infos]
            if len(infos) > 10000 or len(set(archive_names)) != len(archive_names):
                raise ValueError("Package archive has too many or duplicate members")
            for entry in infos:
                _validate_package_path(entry.filename)
                if entry.is_dir():
                    raise ValueError("Package archive contains a directory member")
            open_member = archive.open
        with open_member('manifest.json') as stream:
            manifest_bytes = stream.read(4 * 1024 * 1024 + 1)
        if len(manifest_bytes) > 4 * 1024 * 1024:
            raise ValueError("Package manifest exceeds 4 MiB limit")
        def unique_pairs(pairs):
            value = {}
            for key, item in pairs:
                if key in value:
                    raise ValueError(f"Duplicate manifest key {key!r}")
                value[key] = item
            return value

        def invalid_constant(value):
            raise ValueError(f"Invalid manifest number {value}")

        payload = json.loads(manifest_bytes, object_pairs_hook=unique_pairs,
                             parse_constant=invalid_constant)
        from ..model.object.object import _validate_attributes
        _validate_attributes(payload, max_depth=exchange.MAX_DEPTH)
        manifest = DatasetManifest.model_validate(payload, strict=True)
        if manifest.schema_version != CANONICAL_MANIFEST_VERSION:
            raise ValueError("Unsupported package version; expected canonical manifest v3")
        if not manifest.artifacts or len(manifest.artifacts) > 9999:
            raise ValueError("Canonical package requires a bounded, nonempty artifact list")
        paths = [_validate_package_path(artifact.path) for artifact in manifest.artifacts]
        if len(set(paths)) != len(paths) or 'manifest.json' in paths:
            raise ValueError("Duplicate or reserved package artifact path")
        if archive_names is not None and set(archive_names) != {'manifest.json', *paths}:
            raise ValueError("Archive members do not match the manifest")
        canonical = [x for x in manifest.artifacts if x.role == 'canonical_model']
        if len(canonical) != 1:
            raise ValueError("Package must contain exactly one canonical model artifact")
        canonical = canonical[0]
        total = 0
        for artifact in manifest.artifacts:
            if artifact.size is None or artifact.size < 0 or artifact.sha256 is None:
                raise ValueError("Canonical packages require artifact size and sha256")
            total += artifact.size
            if archive_names is not None and archive.getinfo(artifact.path).file_size != artifact.size:
                raise ValueError("Artifact size does not match manifest")
        if max_bytes is not None and total > max_bytes:
            raise ValueError(f"Package artifacts exceed max_bytes={max_bytes}")
        if canonical.size > exchange.MAX_PROTOBUF_BYTES:
            raise ValueError("Canonical Protobuf message must be smaller than 2 GiB")
        data = None
        for artifact in manifest.artifacts:
            if artifact is not canonical and (
                artifact.role != 'derived' or artifact.derived_from != canonical.path
                or artifact.bounds is not None or artifact.crs is not None
                or artifact.model_type is not None or artifact.model_schema_version is not None
            ):
                raise ValueError("Supplemental artifact must identify its canonical source")
            digest, size, chunks = hashlib.sha256(), 0, []
            with open_member(artifact.path) as stream:
                while chunk := stream.read(1024 * 1024):
                    size += len(chunk)
                    if size > artifact.size:
                        raise ValueError("Artifact size does not match manifest")
                    digest.update(chunk)
                    if artifact is canonical:
                        chunks.append(chunk)
            if size != artifact.size or digest.hexdigest() != artifact.sha256:
                raise ValueError("Artifact size/sha256 does not match manifest")
            if artifact is canonical:
                data = b''.join(chunks)
        model, model_version = exchange._decode_model(data, validate_schema=validate_schema)
        if (canonical.format != 'dtcc' or canonical.model_schema_version != model_version
                or canonical.model_type != type(model).__name__
                or canonical.media_type != 'application/vnd.dtcc.model+protobuf'
                or canonical.data_kind != 'model' or canonical.derived_from is not None):
            raise ValueError("Canonical artifact metadata does not match payload")
        crs = getattr(getattr(model, 'transform', None), 'srs', '') or None
        if canonical.crs != crs or canonical.bounds is not None:
            raise ValueError("Canonical spatial metadata does not match payload contract")
        model.dataset_context = DatasetContext(
            identity=manifest.identity, metadata=manifest.metadata, provenance=manifest.provenance,
            presentation=manifest.presentation, request=manifest.request,
            health=manifest.health, warnings=manifest.warnings or [],
        )
        return model
