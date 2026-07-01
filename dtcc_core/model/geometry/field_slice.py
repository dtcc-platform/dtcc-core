"""Native simulation field slices."""

from __future__ import annotations

import json
from dataclasses import dataclass, field as dataclass_field, replace
from pathlib import Path
from typing import Any

import numpy as np

from ..values import Field
from .bounds import Bounds
from .pointcloud import PointCloud


@dataclass
class FieldSlice(PointCloud):
    """Planar point samples with simulation fields attached to the points."""

    slice_axis: str = "z"
    slice_position: float = 0.0
    resolution: int = 0
    axes: tuple[int, int] = (0, 1)
    name: str = "field_slice"
    time: float = 0.0
    period: float = 0.0
    crs: str | None = None
    include_z: bool = True
    domain_bounds: Bounds | None = None
    metadata_payload: dict[str, Any] = dataclass_field(default_factory=dict)
    default_artifact_format: str = "png"

    def field(self, name: str) -> Field | None:
        """Return a field by name."""
        for field in self.fields:
            if field.name == name:
                return field
        return None

    @property
    def field_names(self) -> list[str]:
        """Return attached field names in storage order."""
        return [field.name for field in self.fields]

    def to_geojson(
        self,
        include_z: bool | None = None,
        crs: str | None = None,
    ) -> dict[str, Any]:
        """Return a GeoJSON adapter representation for debug/vector workflows."""
        include_z = self.include_z if include_z is None else include_z
        crs = self.crs if crs is None else crs
        velocity = self.field("velocity")
        speed = self.field("speed")
        pressure = self.field("pressure")
        velocity_values = _field_values(velocity)
        speed_values = _field_values(speed)
        pressure_values = _field_values(pressure)

        features: list[dict[str, Any]] = []
        for index, point in enumerate(np.asarray(self.points, dtype=float)):
            properties: dict[str, Any] = {"sample_index": index}
            if velocity_values is not None and index < len(velocity_values):
                vector = np.asarray(velocity_values[index], dtype=float).reshape(-1)
                if len(vector) >= 3:
                    properties.update(
                        {
                            "u": float(vector[0]),
                            "v": float(vector[1]),
                            "w": float(vector[2]),
                        }
                    )
            if speed_values is not None and index < len(speed_values):
                properties["speed"] = float(
                    np.asarray(speed_values[index], dtype=float).reshape(-1)[0]
                )
            if pressure_values is not None and index < len(pressure_values):
                properties["pressure"] = float(
                    np.asarray(pressure_values[index], dtype=float).reshape(-1)[0]
                )

            features.append(
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": _coordinates(point, include_z),
                    },
                    "properties": properties,
                }
            )

        payload: dict[str, Any] = {
            "type": "FeatureCollection",
            "name": self.name,
            "features": features,
            "metadata": self._metadata(len(features), include_z=include_z),
        }
        if crs is not None:
            payload["crs"] = {
                "type": "name",
                "properties": {"name": crs},
            }
            payload["metadata"]["crs"] = crs
        return payload

    def to_python(self) -> dict[str, Any]:
        """Return a transition-friendly Python representation."""
        return self.to_geojson()

    def to_plot_product(self, field_name: str = "speed"):
        """Return the plotting product used by DTCC smoke renderers."""
        from dtcc_core.plotting.products import SliceProduct

        scalar = self.field(field_name)
        if scalar is None:
            raise ValueError(f"FieldSlice has no field named {field_name!r}.")
        values = np.asarray(scalar.values, dtype=float).reshape((len(self.points), -1))
        if values.shape[1] != 1:
            raise ValueError(
                f"FieldSlice plot field {field_name!r} must be scalar, "
                f"got dim={values.shape[1]}."
            )
        velocity = self.field("velocity")
        vector_values = None
        if velocity is not None:
            vector_values = np.asarray(velocity.values, dtype=float).reshape(
                (len(self.points), -1)
            )

        metadata = {
            "dataset": "smoke",
            "time": self.time,
            "time_period": self.period,
            "fields": self.field_names,
            **dict(self.metadata_payload),
        }
        return SliceProduct(
            name=self.name,
            bounds=self.domain_bounds or self.bounds,
            axis=self.slice_axis,
            position=self.slice_position,
            coordinates=np.asarray(self.points, dtype=float),
            values=values.reshape(-1),
            resolution=self.resolution,
            field_name=field_name,
            field_unit=scalar.unit,
            vector_values=vector_values,
            axes=self.axes,
            metadata=metadata,
        )

    def plot(
        self,
        ax=None,
        show: bool = True,
        *,
        presentation: bool = True,
        field_name: str = "speed",
        **kwargs,
    ):
        """Plot the slice, optionally with a Dataset v2 presentation panel."""
        from dtcc_core.datasets.presentation import plot_product_with_presentation
        from dtcc_core.plotting.renderers import plot_product

        product = self.to_plot_product(field_name=field_name)
        options = _raster_options_with_overrides(self, kwargs)
        if presentation:
            return plot_product_with_presentation(
                product,
                options,
                context=self.dataset_context,
                obj=self,
                ax=ax,
                show=show,
            )
        return plot_product(product, options, ax=ax, show=show)

    def write_artifact(self, path: str | Path, format: str | None = None) -> None:
        """Write a smoke visualization artifact for object-first package export."""
        fmt = _normalize_format(format or self.default_artifact_format)
        path = Path(path)
        if fmt == "png":
            from dtcc_core.plotting.renderers import render_product_png

            path.write_bytes(
                render_product_png(
                    self.to_plot_product(),
                    _raster_options_from_context(self),
                )
            )
            return
        if fmt == "geojson":
            path.write_text(
                json.dumps(self.to_geojson(), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            return
        if fmt == "mp4":
            raise ValueError(
                "MP4 object-first export for FieldSlice is planned. "
                "Use datasets.smoke.export(..., product='slice', format='mp4') "
                "for the existing dataset-level MP4 path."
            )
        raise ValueError(
            f"FieldSlice object-first export does not support format {fmt!r}."
        )

    def to_proto(self):
        raise NotImplementedError(
            "FieldSlice protobuf serialization is not implemented; use "
            "PointCloud or explicit GeoJSON/PNG export when needed."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "FieldSlice protobuf deserialization is not implemented."
        )

    def _metadata(self, sample_count: int, *, include_z: bool) -> dict[str, Any]:
        metadata = {
            "dataset": "smoke",
            "product": "slice",
            "geometry": "Point",
            "sample_count": sample_count,
            "bounds": _bounds_values(self.domain_bounds or self.bounds),
            "fields": self.field_names,
            "include_z": include_z,
            "slice_axis": self.slice_axis,
            "slice_position": self.slice_position,
            "resolution": self.resolution,
            "visual_axes": list(_axis_names(self.axes)),
            "time": self.time,
            "time_period": self.period,
            **dict(self.metadata_payload),
        }
        return metadata


def _field_values(field: Field | None) -> np.ndarray | None:
    if field is None:
        return None
    return np.asarray(field.values, dtype=float)


def _coordinates(point: np.ndarray, include_z: bool) -> list[float]:
    values = np.asarray(point, dtype=float).reshape(-1)
    if not include_z:
        values = values[:2]
    return values.tolist()


def _bounds_values(bounds: Bounds) -> list[float]:
    return [
        float(bounds.xmin),
        float(bounds.ymin),
        float(bounds.zmin),
        float(bounds.xmax),
        float(bounds.ymax),
        float(bounds.zmax),
    ]


def _axis_names(axes: tuple[int, int]) -> tuple[str, str]:
    names = ("x", "y", "z")
    return (names[axes[0]], names[axes[1]])


def _normalize_format(format: str | None) -> str:
    if format is None:
        return "png"
    return str(format).lower().lstrip(".")


def _raster_options_from_context(obj):
    from dtcc_core.plotting.options import RasterRenderOptions

    context = getattr(obj, "dataset_context", None)
    parameters = {} if context is None else dict(context.request.parameters)
    allowed = RasterRenderOptions.__dataclass_fields__
    option_values = {
        key: parameters[key]
        for key in allowed
        if key in parameters and parameters[key] is not None
    }
    return RasterRenderOptions(**option_values)


def _raster_options_with_overrides(obj, overrides: dict[str, Any]):
    from dtcc_core.plotting.options import RasterRenderOptions

    allowed = set(RasterRenderOptions.__dataclass_fields__)
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise TypeError(f"Unsupported plot option: {unknown[0]}")
    return replace(_raster_options_from_context(obj), **overrides)
