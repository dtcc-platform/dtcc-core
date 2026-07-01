"""Native streamline collections with per-vertex fields."""

from __future__ import annotations

import json
from collections.abc import Iterator
from dataclasses import dataclass, field as dataclass_field
from pathlib import Path
from typing import Any

import numpy as np

from ..model import Model
from ..values import Field
from .bounds import Bounds
from .field_slice import (
    _axis_names,
    _bounds_values,
    _normalize_format,
    _raster_options_from_context,
    _raster_options_with_overrides,
)
from .linestring import LineString


@dataclass
class StreamlineCollection(Model):
    """Streamlines represented as line geometry with fields on each vertex."""

    lines: list[LineString] = dataclass_field(default_factory=list)
    seed_axis: str = "z"
    seed_position: float = 0.0
    requested_line_count: int = 0
    streamline_steps: int = 0
    streamline_step_size: float = 0.0
    time: float = 0.0
    period: float = 0.0
    crs: str | None = None
    include_z: bool = True
    axes: tuple[int, int] = (0, 1)
    name: str = "streamlines"
    domain_bounds: Bounds | None = None
    metadata_payload: dict[str, Any] = dataclass_field(default_factory=dict)
    default_artifact_format: str = "png"

    def __len__(self) -> int:
        return len(self.lines)

    def __iter__(self) -> Iterator[LineString]:
        return iter(self.lines)

    def __getitem__(self, index):
        return self.lines[index]

    @property
    def bounds(self) -> Bounds:
        """Return bounds spanning all line geometry."""
        non_empty = [line.calculate_bounds() for line in self.lines if len(line.vertices)]
        if not non_empty:
            return Bounds()
        return Bounds(
            xmin=float(min(bounds.xmin for bounds in non_empty)),
            ymin=float(min(bounds.ymin for bounds in non_empty)),
            zmin=float(min(bounds.zmin for bounds in non_empty)),
            xmax=float(max(bounds.xmax for bounds in non_empty)),
            ymax=float(max(bounds.ymax for bounds in non_empty)),
            zmax=float(max(bounds.zmax for bounds in non_empty)),
        )

    @property
    def field_names(self) -> list[str]:
        """Return field names attached to streamline vertices."""
        names: list[str] = []
        for line in self.lines:
            for field in line.fields:
                if field.name not in names:
                    names.append(field.name)
        return names

    def to_geojson(
        self,
        include_z: bool = True,
        crs: str | None = None,
    ) -> dict[str, Any]:
        """Return a GeoJSON adapter representation for debug/vector workflows."""
        crs = self.crs if crs is None else crs
        features: list[dict[str, Any]] = []
        for seed_index, line in enumerate(self.lines):
            vertices = np.asarray(line.vertices, dtype=float)
            if len(vertices) < 2:
                continue
            features.append(
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": _coordinates(vertices, include_z),
                    },
                    "properties": {
                        "seed_index": seed_index,
                        "num_points": int(len(vertices)),
                        **_line_field_properties(line),
                    },
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
        return self.to_geojson(include_z=self.include_z)

    def to_plot_product(self, field_name: str = "speed"):
        """Return the plotting product used by DTCC smoke renderers."""
        from dtcc_core.plotting.products import StreamlineProduct

        lines = tuple(np.asarray(line.vertices, dtype=float) for line in self.lines)
        values: list[np.ndarray] = []
        for line in self.lines:
            field = _line_field(line, field_name)
            if field is None:
                values.append(np.empty(0, dtype=float))
                continue
            values.append(np.asarray(field.values, dtype=float).reshape(-1))

        scalar_unit = ""
        for line in self.lines:
            field = _line_field(line, field_name)
            if field is not None:
                scalar_unit = field.unit
                break

        metadata = {
            "dataset": "smoke",
            "requested_line_count": self.requested_line_count,
            "streamline_steps": self.streamline_steps,
            "streamline_step_size": self.streamline_step_size,
            "time": self.time,
            "time_period": self.period,
            "fields": self.field_names,
            **dict(self.metadata_payload),
        }
        return StreamlineProduct(
            name=self.name,
            bounds=self.domain_bounds or self.bounds,
            lines=lines,
            line_values=tuple(values),
            value_name=field_name,
            value_unit=scalar_unit,
            vector_field_name="velocity",
            seed_axis=self.seed_axis,
            seed_position=self.seed_position,
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
        """Plot the streamlines, optionally with a Dataset v2 presentation panel."""
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
                json.dumps(self.to_geojson(include_z=self.include_z), indent=2, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
            return
        if fmt == "mp4":
            raise ValueError(
                "MP4 object-first export for StreamlineCollection is planned. "
                "Use datasets.smoke.export(..., product='streamlines', format='mp4') "
                "for the existing dataset-level MP4 path."
            )
        raise ValueError(
            "StreamlineCollection object-first export does not support "
            f"format {fmt!r}."
        )

    def to_proto(self):
        raise NotImplementedError(
            "StreamlineCollection protobuf serialization is not implemented."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "StreamlineCollection protobuf deserialization is not implemented."
        )

    def _metadata(self, line_count: int, *, include_z: bool) -> dict[str, Any]:
        return {
            "dataset": "smoke",
            "product": "streamlines",
            "geometry": "LineString",
            "sample_count": line_count,
            "line_count": line_count,
            "point_count": int(sum(len(line.vertices) for line in self.lines)),
            "bounds": _bounds_values(self.domain_bounds or self.bounds),
            "fields": self.field_names,
            "include_z": include_z,
            "slice_axis": self.seed_axis,
            "slice_position": self.seed_position,
            "seed_axis": self.seed_axis,
            "seed_position": self.seed_position,
            "requested_line_count": self.requested_line_count,
            "streamline_steps": self.streamline_steps,
            "streamline_step_size": self.streamline_step_size,
            "visual_axes": list(_axis_names(self.axes)),
            "time": self.time,
            "time_period": self.period,
            **dict(self.metadata_payload),
        }


def _line_field(line: LineString, name: str) -> Field | None:
    for field in line.fields:
        if field.name == name:
            return field
    return None


def _line_field_properties(line: LineString) -> dict[str, Any]:
    properties: dict[str, Any] = {}
    for field in line.fields:
        values = np.asarray(field.values, dtype=float)
        if field.dim == 1:
            properties[field.name] = values.reshape(-1).tolist()
        else:
            properties[field.name] = values.reshape((-1, field.dim)).tolist()
    return properties


def _coordinates(points: np.ndarray, include_z: bool) -> list:
    values = np.asarray(points, dtype=float)
    if not include_z:
        values = values[:, :2]
    return values.tolist()
