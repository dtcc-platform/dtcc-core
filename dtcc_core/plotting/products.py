"""Visualization product containers shared by plot and export paths."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from dtcc_core.model import Bounds


AXIS_NAMES = ("x", "y", "z")


@dataclass(frozen=True)
class SliceProduct:
    """A scalar field sampled on a planar slice."""

    name: str
    bounds: Bounds
    axis: str
    position: float
    coordinates: np.ndarray
    values: np.ndarray
    resolution: int
    field_name: str = "value"
    field_unit: str = ""
    vector_values: np.ndarray | None = None
    axes: tuple[int, int] = (0, 1)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def axes_names(self) -> tuple[str, str]:
        return (AXIS_NAMES[self.axes[0]], AXIS_NAMES[self.axes[1]])

    @property
    def extent(self) -> tuple[float, float, float, float]:
        xy = self.coordinates[:, self.axes]
        return (
            float(np.min(xy[:, 0])),
            float(np.max(xy[:, 0])),
            float(np.min(xy[:, 1])),
            float(np.max(xy[:, 1])),
        )

    @property
    def image(self) -> np.ndarray:
        return np.asarray(self.values, dtype=float).reshape(
            (self.resolution, self.resolution)
        ).T

    @property
    def fields(self) -> list[str]:
        fields = []
        if self.vector_values is not None:
            fields.append("velocity")
        fields.append(self.field_name)
        return fields

    def manifest_dict(self) -> dict[str, Any]:
        return {
            "product_kind": "slice",
            "field": self.field_name,
            "field_unit": self.field_unit,
            "slice_axis": self.axis,
            "slice_position": self.position,
            "resolution": self.resolution,
            "visual_axes": list(self.axes_names),
            "extent": _bbox_extent(self.extent),
            "fields": self.fields,
            **dict(self.metadata),
        }


@dataclass(frozen=True)
class StreamlineProduct:
    """Streamlines derived from a vector field."""

    name: str
    bounds: Bounds
    lines: tuple[np.ndarray, ...]
    line_values: tuple[np.ndarray, ...] | None = None
    value_name: str = "speed"
    value_unit: str = ""
    vector_field_name: str = "velocity"
    seed_axis: str = "z"
    seed_position: float = 0.5
    axes: tuple[int, int] = (0, 1)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def axes_names(self) -> tuple[str, str]:
        return (AXIS_NAMES[self.axes[0]], AXIS_NAMES[self.axes[1]])

    @property
    def extent(self) -> tuple[float, float, float, float]:
        mins = np.array([self.bounds.xmin, self.bounds.ymin, self.bounds.zmin])
        maxs = np.array([self.bounds.xmax, self.bounds.ymax, self.bounds.zmax])
        return (
            float(mins[self.axes[0]]),
            float(maxs[self.axes[0]]),
            float(mins[self.axes[1]]),
            float(maxs[self.axes[1]]),
        )

    @property
    def fields(self) -> list[str]:
        return [self.vector_field_name, self.value_name]

    def manifest_dict(self) -> dict[str, Any]:
        point_count = int(sum(len(line) for line in self.lines))
        return {
            "product_kind": "streamlines",
            "vector_field": self.vector_field_name,
            "value_field": self.value_name,
            "value_unit": self.value_unit,
            "line_count": len(self.lines),
            "point_count": point_count,
            "seed_axis": self.seed_axis,
            "seed_position": self.seed_position,
            "visual_axes": list(self.axes_names),
            "extent": _bbox_extent(self.extent),
            "fields": self.fields,
            **dict(self.metadata),
        }


def _bbox_extent(extent: tuple[float, float, float, float]) -> list[float]:
    xmin, xmax, ymin, ymax = extent
    return [xmin, ymin, xmax, ymax]
