"""Dataset Manifest v2 package export helpers."""

from __future__ import annotations

import hashlib
import json
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path

from dtcc_core.datasets.dataset import DatasetDescriptor
from dtcc_core.datasets.schema import DatasetArtifact, DatasetManifest


@dataclass(frozen=True)
class DatasetPackage:
    """Object-first Dataset v2 export package."""

    path: Path
    manifest_path: Path
    manifest: DatasetManifest
    artifacts: tuple[DatasetArtifact, ...]
    files: tuple[Path, ...]
    package_format: str


def export_model_package(
    obj,
    path: str | Path,
    *,
    format: str | None = None,
) -> DatasetPackage:
    """Export a dataset-produced model object as a Dataset Manifest v2 package."""
    context = getattr(obj, "dataset_context", None)
    if context is None:
        raise ValueError(
            "Cannot export Dataset v2 package: this object has no DatasetContext."
        )

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
    artifact_format: str,
    package_path: Path,
    package_format: str,
) -> DatasetPackage:
    artifact_dir = package_dir / "artifacts"
    artifact_dir.mkdir(parents=True, exist_ok=True)

    context = obj.dataset_context
    extension = DatasetDescriptor.format_extension(artifact_format)
    artifact_stem = _artifact_stem(obj, context)
    artifact_path = artifact_dir / f"{artifact_stem}.{extension}"
    _write_artifact(obj, artifact_path, artifact_format)

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


def _write_artifact(obj, path: Path, artifact_format: str) -> None:
    fmt = _normalize_format(artifact_format)
    write_artifact = getattr(obj, "write_artifact", None)
    if callable(write_artifact):
        write_artifact(path, format=fmt)
        return

    if fmt == "geojson" and hasattr(obj, "to_geojson"):
        path.write_text(
            json.dumps(obj.to_geojson(), indent=2, sort_keys=True) + "\n",
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
        "PointCloud": "pb",
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
