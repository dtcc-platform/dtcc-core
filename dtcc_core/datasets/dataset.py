from abc import ABC, abstractmethod
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import tempfile
from typing import Any, Optional, Sequence, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator

import dtcc_core
from dtcc_core.model import Bounds
from dtcc_core.model import Object as DTCCObject
from dtcc_core.model import Geometry as DTCCGeometry

from .context import attach_dataset_context
from .schema import (
    DatasetContext,
    DatasetIdentity,
    DatasetMetadata,
    DatasetPresentation,
    DatasetProvenance,
    DatasetRequest,
)


_FORMAT_KIND_MAP = {
    "tif": "raster",
    "tiff": "raster",
    "asc": "raster",
    "png": "raster",
    "jpg": "raster",
    "jpeg": "raster",
    "mp4": "video",
    "geojson": "vector",
    "gpkg": "vector",
    "shp.zip": "vector",
    "obj": "mesh",
    "stl": "mesh",
    "ply": "mesh",
    "vtk": "mesh",
    "vtu": "mesh",
    "xdmf": "mesh",
    "inp": "mesh",
    "bdf": "mesh",
    "las": "point_cloud",
    "laz": "point_cloud",
    "copc": "point_cloud",
    "cityjson": "city_model",
    "city.json": "city_model",
    "json.zip": "city_model",
    "pb": "protobuf",
}


_FORMAT_MEDIA_TYPE_MAP = {
    "png": "image/png",
    "jpg": "image/jpeg",
    "jpeg": "image/jpeg",
    "mp4": "video/mp4",
    "tif": "image/tiff",
    "tiff": "image/tiff",
    "geojson": "application/geo+json",
    "gpkg": "application/geopackage+sqlite3",
    "shp.zip": "application/zip",
    "obj": "model/obj",
    "stl": "model/stl",
    "vtk": "application/vnd.vtk",
    "vtu": "application/vnd.vtk.vtu+xml",
    "xdmf": "application/x-xdmf",
    "cityjson": "application/json",
    "city.json": "application/json",
    "json": "application/json",
    "json.zip": "application/zip",
    "tar.gz": "application/gzip",
    "pb": "application/x-protobuf",
}


_FORMAT_EXTENSION_MAP = {
    "cityjson": "city.json",
}


class DatasetBaseArgs(BaseModel):
    model_config = ConfigDict(extra="forbid")

    bounds: Sequence[float] = Field(
        ...,
        description="Bounding box [minx, miny, maxx, maxy] or [minx, miny, minz, maxx, maxy, maxz]",
    )
    strict_live: bool = Field(
        False,
        description=(
            "Raise a typed upstream error instead of silently degrading to "
            "partial or empty results when live network fetches fail."
        ),
    )

    @field_validator("bounds")
    @classmethod
    def validate_bounds(cls, v: Sequence[float]):
        """Validate bounds length and order."""
        if len(v) not in (4, 6):
            raise ValueError("Bounds must be 4 or 6 floats")

        # Validate order
        if len(v) == 4:
            if v[0] >= v[2] or v[1] >= v[3]:
                raise ValueError("Invalid bounds: xmin < xmax, ymin < ymax")
        elif len(v) == 6:
            if v[0] >= v[3] or v[1] >= v[4] or v[2] >= v[5]:
                raise ValueError("Invalid bounds: min < max for all dimensions")
        return v


class DatasetUpstreamError(RuntimeError):
    """Typed error for live upstream failures.

    Live test suites can use ``failure_class`` to distinguish transient
    upstream outages (connection, timeout, HTTP 5xx) from hard failures such
    as HTTP 4xx or invalid upstream payloads.
    """

    def __init__(
        self,
        *,
        dataset: str,
        operation: str,
        target: str,
        failure_class: str,
        status_code: Optional[int] = None,
        message: str,
    ):
        super().__init__(message)
        self.dataset = dataset
        self.operation = operation
        self.target = target
        self.failure_class = failure_class
        self.status_code = status_code
        self.message = message

    @property
    def is_transient(self) -> bool:
        """Whether the error represents a likely transient upstream issue."""
        return self.failure_class in {"connection", "timeout", "http_5xx"}


@dataclass(frozen=True)
class DatasetExportResult:
    """Result returned by :meth:`DatasetDescriptor.export`."""

    path: Path
    manifest_path: Optional[Path]
    manifest: Optional[dict[str, Any]]
    files: tuple[Path, ...] = ()
    format: str = ""

    def __post_init__(self):
        if not self.files:
            object.__setattr__(self, "files", (Path(self.path),))
        if not self.format:
            object.__setattr__(
                self, "format", self._infer_format_from_export_path(self.path)
            )

    @staticmethod
    def _infer_format_from_export_path(path: Path) -> str:
        # Direct DatasetExportResult construction has no descriptor context.
        # Use only global extension metadata as a compatibility fallback.
        format_extensions = {
            format_name: format_name
            for format_name in {*_FORMAT_KIND_MAP, *_FORMAT_MEDIA_TYPE_MAP}
        }
        for format_name, extension in _FORMAT_EXTENSION_MAP.items():
            format_extensions[extension] = format_name

        path_name = Path(path).name.lower()
        for extension, format_name in sorted(
            format_extensions.items(), key=lambda item: len(item[0]), reverse=True
        ):
            if path_name.endswith(f".{extension.lower()}"):
                return format_name
        return ""

    def publish(
        self,
        *,
        dataset_key: str,
        uploader=None,
        upload_url: Optional[str] = None,
        token: Optional[str] = None,
        idempotency_key: Optional[str] = None,
    ):
        from dtcc_core.datasets.publish import DatasetPackageError, DatasetUploadClient

        if self.manifest_path is None or self.manifest is None:
            raise DatasetPackageError(
                "Cannot publish an export result without a manifest.",
                failure_class="invalid_package",
            )

        resolved_uploader = uploader or DatasetUploadClient.from_config(
            upload_url=upload_url, token=token
        )
        return resolved_uploader.upload_package(
            dataset_key=dataset_key,
            manifest_path=self.manifest_path,
            files=self.files,
            manifest=self.manifest,
            idempotency_key=idempotency_key,
        )


class DatasetDescriptor(ABC):
    """Callable, self-describing dataset."""

    name: str
    title: Optional[str] = None
    description: str = ""
    ArgsModel: BaseModel
    data_category: str = "unknown"
    result_kind: str = "unknown"
    python_return_type: str = "object"
    timeout_hint: Optional[int] = None
    multi_file_formats: Sequence[str] = ()

    def __init_subclass__(cls, register=True, **kwargs):
        """
        Auto-register dataset subclasses when they're defined.

        Args:
            register: Whether to auto-register this dataset (default: True).
                     Set to False for abstract base classes.
            **kwargs: Additional keyword arguments passed to super().__init_subclass__
        """
        super().__init_subclass__(**kwargs)

        # Only register if:
        # - registration is enabled (register=True)
        # - class has a name attribute
        # - name is not empty
        if register and hasattr(cls, "name") and cls.name:
            from dtcc_core.datasets.registry import _register_dataset_class

            _register_dataset_class(cls.name, cls)

    def __call__(self, **kwargs):
        args = self.validate(kwargs)
        result = self.build(args)
        result = self.prepare_result(result, args)
        context = self.create_context(args)
        return attach_dataset_context(result, context)

    def validate(self, kwargs):
        if isinstance(kwargs.get("bounds"), Bounds):
            kwargs["bounds"] = kwargs["bounds"].tuple
        args = self.ArgsModel(**kwargs)
        return args

    @abstractmethod
    def build(self, validated_args):
        """Resolve the dataset and return the result."""
        raise NotImplementedError

    def prepare_result(self, result, validated_args):
        """Adapt a public dataset-call result before context is attached.

        Subclasses can override this to migrate selected bare list/dict
        returns to DatasetCollection or DatasetValue without changing internal
        build/export code paths.
        """
        return result

    # TODO(Dataset v2): keep DatasetDescriptor during migration; Dataset is a
    # public alias below so new code can use the Dataset name without breaking
    # existing registrations.
    def create_context(self, args) -> DatasetContext:
        """Build initial Dataset v2 context from descriptor metadata."""
        parameters = args.model_dump(mode="json")
        descriptor = self.describe()
        title = descriptor.get("title") or self._title_from_identifier(self.name)
        bounds = parameters.get("bounds")
        crs_values = self._context_crs_values(parameters)
        data_category = descriptor.get("data_category")
        result_kind = descriptor.get("result_kind")

        return DatasetContext(
            identity=DatasetIdentity(
                name=self.name,
                title=title,
                version=getattr(self, "version", None),
            ),
            metadata=DatasetMetadata(
                description=descriptor.get("description") or "",
                provider=self._context_list_attr("provider", "providers"),
                source=self._context_list_attr("source", "sources", "source_service"),
                crs=crs_values,
                lod=self._context_optional_string(parameters.get("lod")),
                data_types=self._dedupe_strings(
                    item for item in (result_kind,) if item
                ),
                formats=list(descriptor.get("supported_formats") or ()),
                data_category=str(data_category) if data_category else None,
                result_kind=str(result_kind) if result_kind else None,
                python_return_type=descriptor.get("python_return_type"),
            ),
            provenance=DatasetProvenance(
                sources=self._context_list_attr("source", "sources", "source_service"),
                processing_steps=[f"Build dataset '{self.name}'"],
                generated_by={
                    "package": "dtcc-core",
                    "version": self._package_version(),
                },
            ),
            presentation=DatasetPresentation(
                headline=title,
                summary=descriptor.get("description") or None,
            ),
            request=DatasetRequest(
                dataset_name=self.name,
                parameters=parameters,
                bounds=bounds,
            ),
        )

    @staticmethod
    def _dedupe_strings(values) -> list[str]:
        strings: list[str] = []
        for value in values:
            if value is None:
                continue
            text = str(value)
            if text and text not in strings:
                strings.append(text)
        return strings

    @staticmethod
    def _context_optional_string(value: Any) -> str | None:
        if value is None:
            return None
        if hasattr(value, "name"):
            return str(value.name)
        return str(value)

    @classmethod
    def _context_crs_values(cls, parameters: dict[str, Any]) -> list[str]:
        crs = parameters.get("crs")
        if crs is None:
            return []
        if isinstance(crs, (list, tuple, set)):
            return cls._dedupe_strings(crs)
        return cls._dedupe_strings((crs,))

    def _context_list_attr(self, *names: str) -> list[Any]:
        values: list[Any] = []
        for name in names:
            value = getattr(self, name, None)
            if value is None:
                continue
            if isinstance(value, (list, tuple, set)):
                values.extend(value)
            else:
                values.append(value)

        deduped: list[Any] = []
        seen: set[str] = set()
        for value in values:
            marker = (
                json.dumps(value, sort_keys=True)
                if isinstance(value, dict)
                else str(value)
            )
            if marker not in seen:
                deduped.append(value)
                seen.add(marker)
        return deduped

    def show_options(self):
        return self.ArgsModel.model_json_schema()

    @staticmethod
    def _normalize_format_values(value: Any) -> list[str]:
        if value is None:
            return []
        if not isinstance(value, list):
            value = [value]
        formats = []
        for item in value:
            if item is None:
                continue
            fmt = str(item).strip().lower()
            if fmt and fmt not in {"none", "null"} and fmt not in formats:
                formats.append(fmt)
        return formats

    @classmethod
    def extract_supported_formats_from_schema(cls, schema: dict[str, Any]) -> list[str]:
        """Extract supported ``format`` values from a Pydantic JSON schema."""
        properties = schema.get("properties", {}) if isinstance(schema, dict) else {}
        format_prop = properties.get("format", {}) if isinstance(properties, dict) else {}
        if not isinstance(format_prop, dict):
            return []

        formats = []
        formats.extend(cls._normalize_format_values(format_prop.get("enum")))
        if "const" in format_prop:
            formats.extend(cls._normalize_format_values(format_prop.get("const")))

        any_of = format_prop.get("anyOf")
        if isinstance(any_of, list):
            for variant in any_of:
                if not isinstance(variant, dict):
                    continue
                formats.extend(cls._normalize_format_values(variant.get("enum")))
                if "const" in variant:
                    formats.extend(cls._normalize_format_values(variant.get("const")))

        if not formats:
            formats.extend(cls._normalize_format_values(format_prop.get("default")))

        deduped = []
        for fmt in formats:
            if fmt not in deduped:
                deduped.append(fmt)
        return deduped

    def list_supported_formats(self) -> list[str]:
        """Return serialized output formats supported by this dataset."""
        explicit_formats = getattr(self, "supported_formats", None)
        if explicit_formats is not None and not callable(explicit_formats):
            return self._normalize_format_values(explicit_formats)
        return self.extract_supported_formats_from_schema(self.show_options())

    @staticmethod
    def format_kind(format: str) -> str:
        """Return a coarse data kind for a serialized output format."""
        return _FORMAT_KIND_MAP.get(str(format).lower(), "unknown")

    @staticmethod
    def format_media_type(format: str) -> str:
        """Return the default HTTP media type for a serialized output format."""
        return _FORMAT_MEDIA_TYPE_MAP.get(
            str(format).lower(), "application/octet-stream"
        )

    @staticmethod
    def format_extension(format: str) -> str:
        """Return the recommended filename extension for a format value."""
        fmt = str(format).lower()
        return _FORMAT_EXTENSION_MAP.get(fmt, fmt)

    def format_metadata(self) -> list[dict[str, Any]]:
        """Return JSON-safe metadata for each supported serialized format."""
        multi_file_formats = {
            str(fmt).lower() for fmt in getattr(self, "multi_file_formats", ())
        }
        return [
            {
                "format": fmt,
                "extension": self.format_extension(fmt),
                "media_type": self.format_media_type(fmt),
                "data_kind": self.format_kind(fmt),
                "multi_file": fmt in multi_file_formats,
            }
            for fmt in self.list_supported_formats()
        ]

    @staticmethod
    def _title_from_identifier(identifier: str) -> str:
        return str(identifier).replace("_", " ").replace("-", " ").title()

    @staticmethod
    def _package_version() -> str:
        try:
            return version("dtcc-core")
        except PackageNotFoundError:
            return "0.9.8dev"

    def describe(self) -> dict[str, Any]:
        """Return the dataset contract used by Python clients and web services."""
        args_schema = self.show_options()
        formats = self.list_supported_formats()
        schema_properties = (
            args_schema.get("properties", {}) if isinstance(args_schema, dict) else {}
        )
        return {
            "name": self.name,
            "title": getattr(self, "title", None)
            or self._title_from_identifier(self.name),
            "description": self.description,
            "data_category": getattr(self, "data_category", "unknown"),
            "result_kind": getattr(self, "result_kind", "unknown"),
            "python_return_type": getattr(self, "python_return_type", "object"),
            "args_schema": args_schema,
            "supported_formats": formats,
            "formats": self.format_metadata(),
            "multi_file_formats": list(getattr(self, "multi_file_formats", ())),
            "timeout_hint": getattr(self, "timeout_hint", None),
            "serialization": {
                "python_object_when_format_omitted": True,
                "bytes_when_format_is_set": bool(formats),
                "format_parameter": "format" in schema_properties,
            },
        }

    def _infer_format_from_path(self, path: Path) -> str:
        """Infer a dataset format from a path suffix."""
        candidates = sorted(
            (
                (item["extension"], item["format"])
                for item in self.format_metadata()
                if item.get("extension") and item.get("format")
            ),
            key=lambda item: len(item[0]),
            reverse=True,
        )
        path_name = path.name.lower()
        for extension, format_name in candidates:
            if path_name.endswith(f".{extension.lower()}"):
                return format_name

        supported_formats = ", ".join(self.list_supported_formats()) or "none"
        raise ValueError(
            f"Could not infer export format from '{path}'. "
            f"Pass format explicitly. Supported formats: {supported_formats}"
        )

    def _build_manifest(
        self,
        args,
        path: Path,
        *,
        manifest_id: Optional[str] = None,
        title: Optional[str] = None,
        description: Optional[str] = None,
    ) -> dict[str, Any]:
        """Build a manifest for a concrete exported dataset request."""
        manifest = self.describe()
        manifest["manifest_schema_version"] = "dtcc-dataset-manifest-v1"
        manifest["created_by"] = {
            "package": "dtcc-core",
            "version": self._package_version(),
        }
        parameters = args.model_dump(mode="json")
        product = getattr(args, "product", None)

        if manifest_id is not None:
            manifest["id"] = manifest_id
        if title is not None:
            manifest["title"] = title
        if description is not None:
            manifest["description"] = description

        manifest["file"] = path.name
        manifest["format"] = getattr(args, "format", None)
        if manifest["format"] is not None:
            manifest["media_type"] = self.format_media_type(manifest["format"])
            manifest["data_kind"] = self.format_kind(manifest["format"])
        if product is not None:
            manifest["product"] = product
        manifest["bounds"] = list(parameters["bounds"])
        manifest["parameters"] = parameters
        manifest.update(self.export_manifest_metadata(args, path))
        return manifest

    def export_manifest_metadata(self, args, path: Path) -> dict[str, Any]:
        """Return dataset-specific metadata for an exported artifact."""
        return {}

    @staticmethod
    def _sanitize_publish_filename(value: str) -> str:
        filename = str(value).strip().replace(" ", "_").replace("-", "_").lower()
        sanitized = "".join(
            char if char.isalnum() or char in {"_", "."} else "_"
            for char in filename
        )
        while "__" in sanitized:
            sanitized = sanitized.replace("__", "_")
        return sanitized.strip("_")

    def _default_publish_filename(
        self, *, format: str | None, kwargs: dict[str, Any]
    ) -> str:
        if format is None:
            raise ValueError("publish() requires format when filename is omitted")

        extension = self.format_extension(format)
        product = kwargs.get("product")
        stem = self.name if product in {None, "", "field"} else f"{self.name}_{product}"
        filename = f"{self._sanitize_publish_filename(stem)}.{extension}"
        if (
            filename.startswith(".")
            or filename.startswith("/")
            or filename.startswith("\\")
            or ".." in filename
        ):
            raise ValueError(f"Unsafe generated publish filename: {filename!r}")
        return filename

    @staticmethod
    def _validate_publish_filename(value: Union[str, Path]) -> Path:
        filename = str(value)
        path = Path(filename)
        if (
            not filename
            or path.is_absolute()
            or path.name != filename
            or filename in {".", ".."}
            or filename.startswith(".")
            or "/" in filename
            or "\\" in filename
            or ".." in filename
        ):
            raise ValueError(f"Unsafe publish filename: {filename!r}")
        return path

    def _reject_multi_file_publish_format(self, format: str) -> None:
        multi_file_formats = {
            str(fmt).lower() for fmt in getattr(self, "multi_file_formats", ())
        }
        if str(format).lower() in multi_file_formats:
            from dtcc_core.datasets.publish import DatasetPackageError

            raise DatasetPackageError(
                f"publish() does not support multi-file format {format!r}.",
                failure_class="invalid_package",
            )

    def export(
        self,
        path: Union[str, Path],
        *,
        format: Optional[str] = None,
        manifest: bool = True,
        manifest_path: Optional[Union[str, Path]] = None,
        manifest_id: Optional[str] = None,
        title: Optional[str] = None,
        description: Optional[str] = None,
        **kwargs,
    ) -> DatasetExportResult:
        """Export a serialized dataset artifact and optional manifest to disk.

        Args:
            path: Output path for the serialized dataset file.
            format: Dataset format. If omitted, inferred from ``path``.
            manifest: Whether to write an Atlas-style manifest sidecar.
            manifest_path: Optional manifest output path. Defaults to
                ``path`` with ``.manifest.json`` as suffix.
            manifest_id: Optional manifest ``id`` value to include.
            title: Optional manifest title override.
            description: Optional manifest description override.
            **kwargs: Dataset arguments passed to the dataset build call.

        Returns:
            Paths and manifest data for the exported artifact.
        """
        output_path = Path(path)
        output_format = format or self._infer_format_from_path(output_path)

        request_kwargs = dict(kwargs)
        request_kwargs["format"] = output_format
        args = self.validate(request_kwargs)
        payload = self.build(args)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(payload, str):
            output_path.write_text(payload, encoding="utf-8")
        elif isinstance(payload, (bytes, bytearray)):
            output_path.write_bytes(payload)
        else:
            raise TypeError(
                f"Dataset export requires a serialized bytes or string payload, "
                f"got {type(payload).__name__}. Pass a supported format."
            )

        exported_manifest = None
        exported_manifest_path = None
        if manifest:
            exported_manifest_path = (
                Path(manifest_path)
                if manifest_path is not None
                else output_path.with_suffix(".manifest.json")
            )
            exported_manifest = self._build_manifest(
                args,
                output_path,
                manifest_id=manifest_id,
                title=title,
                description=description,
            )
            exported_manifest_path.parent.mkdir(parents=True, exist_ok=True)
            exported_manifest_path.write_text(
                json.dumps(exported_manifest, indent=2) + "\n",
                encoding="utf-8",
            )

        return DatasetExportResult(
            path=output_path,
            manifest_path=exported_manifest_path,
            manifest=exported_manifest,
            files=(output_path,),
            format=output_format,
        )

    def publish(
        self,
        *,
        dataset_key: str,
        format: Optional[str] = None,
        filename: Optional[Union[str, Path]] = None,
        output_dir: Optional[Union[str, Path]] = None,
        keep_export: bool = False,
        manifest_id: Optional[str] = None,
        title: Optional[str] = None,
        description: Optional[str] = None,
        upload_url: Optional[str] = None,
        token: Optional[str] = None,
        idempotency_key: Optional[str] = None,
        uploader=None,
        **kwargs,
    ):
        if filename is None:
            output_filename = Path(
                self._default_publish_filename(format=format, kwargs=kwargs)
            )
        else:
            output_filename = self._validate_publish_filename(filename)

        output_format = format or self._infer_format_from_path(output_filename)
        self._reject_multi_file_publish_format(output_format)

        if keep_export:
            package_dir = Path(output_dir) if output_dir is not None else Path(".")
            package_dir.mkdir(parents=True, exist_ok=True)
            package = self.export(
                package_dir / output_filename,
                format=output_format,
                manifest_id=manifest_id,
                title=title,
                description=description,
                **kwargs,
            )
            return package.publish(
                dataset_key=dataset_key,
                uploader=uploader,
                upload_url=upload_url,
                token=token,
                idempotency_key=idempotency_key,
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            package = self.export(
                Path(tmpdir) / output_filename,
                format=output_format,
                manifest_id=manifest_id,
                title=title,
                description=description,
                **kwargs,
            )
            return package.publish(
                dataset_key=dataset_key,
                uploader=uploader,
                upload_url=upload_url,
                token=token,
                idempotency_key=idempotency_key,
            )

    def __str__(self):
        """Return a nicely formatted summary of the dataset."""
        lines = []
        lines.append("=" * 70)
        lines.append(f"Dataset: {self.name}")
        lines.append("=" * 70)

        if self.description:
            lines.append(f"\nDescription:")
            # Wrap long descriptions nicely
            desc_lines = self.description.split("\n")
            for desc_line in desc_lines:
                lines.append(f"  {desc_line}")

        lines.append(f"\nAvailable Parameters:")
        lines.append("-" * 70)

        # Get schema information from ArgsModel
        schema = self.ArgsModel.model_json_schema()
        properties = schema.get("properties", {})
        required_fields = schema.get("required", [])

        if properties:
            for param_name, param_info in properties.items():
                param_type = param_info.get("type", "any")
                param_desc = param_info.get("description", "")
                default_val = param_info.get("default")
                is_required = param_name in required_fields

                # Format parameter type
                if "anyOf" in param_info:
                    # Handle union types
                    types = [t.get("type", str(t)) for t in param_info["anyOf"]]
                    param_type = " | ".join(str(t) for t in types)
                elif "items" in param_info:
                    # Handle array types
                    item_type = param_info["items"].get("type", "any")
                    param_type = f"array of {item_type}"

                # Format the line
                required_marker = "*" if is_required else " "
                param_line = f"  {required_marker} {param_name} ({param_type})"

                # Add default value if present
                if default_val is not None and not is_required:
                    param_line += f" = {default_val}"

                lines.append(param_line)
                if param_desc:
                    lines.append(f"      {param_desc}")
        else:
            lines.append("  No parameters defined")

        lines.append("\n" + "=" * 70)
        lines.append("* = required parameter")

        return "\n".join(lines)

    @staticmethod
    def parse_bounds(bounds: Sequence[float]) -> Bounds:
        """Convert bounds list to a Bounds object.

        Args:
            bounds: [minx, miny, maxx, maxy] or
                   [minx, miny, minz, maxx, maxy, maxz]

        Returns:
            Bounds object
        """
        if len(bounds) == 4:
            return Bounds(
                xmin=bounds[0], ymin=bounds[1], xmax=bounds[2], ymax=bounds[3]
            )
        elif len(bounds) == 6:
            return Bounds(
                xmin=bounds[0],
                ymin=bounds[1],
                zmin=bounds[2],
                xmax=bounds[3],
                ymax=bounds[4],
                zmax=bounds[5],
            )
        else:
            raise ValueError(f"Bounds must be 4 or 6 floats, got {len(bounds)}")

    @staticmethod
    def point_within_bounds(
        x: float, y: float, bounds: Bounds, tol: float = 1e-6
    ) -> bool:
        """Check whether a 2D point lies within bounds in the same CRS."""
        return (
            bounds.xmin - tol <= x <= bounds.xmax + tol
            and bounds.ymin - tol <= y <= bounds.ymax + tol
        )

    @staticmethod
    def build_upstream_error(
        dataset: str, operation: str, target: str, exc: Exception
    ) -> DatasetUpstreamError:
        """Classify request-layer failures into a typed upstream error."""
        try:
            import requests
        except ImportError:
            requests = None

        failure_class = "request"
        status_code = None

        if requests is not None:
            if isinstance(exc, requests.Timeout):
                failure_class = "timeout"
            elif isinstance(exc, requests.ConnectionError):
                failure_class = "connection"
            elif isinstance(exc, requests.HTTPError):
                response = getattr(exc, "response", None)
                status_code = getattr(response, "status_code", None)
                if status_code is not None and 500 <= status_code <= 599:
                    failure_class = "http_5xx"
                elif status_code is not None and 400 <= status_code <= 499:
                    failure_class = "http_4xx"
                else:
                    failure_class = "http_error"
            elif isinstance(exc, requests.RequestException):
                failure_class = "request"

        message = f"{dataset} {operation} failed for {target}: {exc}"
        return DatasetUpstreamError(
            dataset=dataset,
            operation=operation,
            target=target,
            failure_class=failure_class,
            status_code=status_code,
            message=message,
        )

    @staticmethod
    def serialize_upstream_error(exc: DatasetUpstreamError) -> dict[str, Any]:
        """Convert a typed upstream error into JSON-safe metadata."""
        return {
            "dataset": exc.dataset,
            "operation": exc.operation,
            "target": exc.target,
            "failure_class": exc.failure_class,
            "status_code": exc.status_code,
            "message": exc.message,
            "is_transient": exc.is_transient,
        }

    @classmethod
    def apply_result_health_metadata(
        cls,
        attributes: dict[str, Any],
        *,
        upstream_errors: Sequence[DatasetUpstreamError],
        stations_skipped_upstream: int = 0,
        requested_parameters: Optional[Sequence[Any]] = None,
        fetched_parameters: Optional[Sequence[Any]] = None,
    ) -> None:
        """Attach a uniform graceful-degradation contract to result metadata."""
        serialized_errors = [cls.serialize_upstream_error(exc) for exc in upstream_errors]
        attributes["partial_result"] = bool(serialized_errors)
        attributes["upstream_error_count"] = len(serialized_errors)
        attributes["upstream_errors"] = serialized_errors
        attributes["stations_skipped_upstream"] = stations_skipped_upstream

        if requested_parameters is not None:
            attributes["requested_parameters"] = list(requested_parameters)
        if fetched_parameters is not None:
            attributes["fetched_parameters"] = list(fetched_parameters)

    @staticmethod
    def export_to_bytes(
        obj: Union[DTCCObject, DTCCGeometry, list[DTCCObject], list[DTCCGeometry]],
        format: str,
        save_callable=None,
        **save_kwargs,
    ) -> bytes:
        """Export object to bytes.

        Args:
            obj: Object with .save() method
            format: File format extension
            save_callable: Custom save function
            **save_kwargs: Passed to obj.save()

        Returns:
            File contents as bytes
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpfile = Path(tmpdir) / f"data.{format}"
            if save_callable is not None:
                save_callable(obj, tmpfile, **save_kwargs)
            else:
                obj.save(tmpfile, **save_kwargs)
            return Path(tmpfile).read_bytes()


Dataset = DatasetDescriptor
