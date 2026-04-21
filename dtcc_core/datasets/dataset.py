import dtcc_core
from dtcc_core.model import Bounds
from dtcc_core.model import Object as DTCCObject
from dtcc_core.model import Geometry as DTCCGeometry

from abc import ABC, abstractmethod
from typing import Any, Optional, Sequence, Union
from pydantic import BaseModel, Field, field_validator
from pathlib import Path
import tempfile


class DatasetBaseArgs(BaseModel):
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


class DatasetDescriptor(ABC):
    """Callable, self-describing dataset."""

    name: str
    description: str = ""
    ArgsModel: BaseModel

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
        return self.build(args)

    def validate(self, kwargs):
        if isinstance(kwargs.get("bounds"), Bounds):
            kwargs["bounds"] = kwargs["bounds"].tuple
        args = self.ArgsModel(**kwargs)
        return args

    @abstractmethod
    def build(self, validated_args):
        """Resolve the dataset and return the result."""
        raise NotImplementedError

    def show_options(self):
        return self.ArgsModel.model_json_schema()

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
        with tempfile.TemporaryDirectory(delete=True) as tmpdir:
            tmpfile = Path(tmpdir) / f"data.{format}"
            if save_callable is not None:
                save_callable(obj, tmpfile, **save_kwargs)
            else:
                obj.save(tmpfile, **save_kwargs)
            return Path(tmpfile).read_bytes()
