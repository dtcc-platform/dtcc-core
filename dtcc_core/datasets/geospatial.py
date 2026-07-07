"""Shared geospatial helpers for dataset providers."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..reproject.reproject import reproject_array


WGS84_CRS_ALIASES = {"CRS84", "EPSG:4326", "WGS84"}


def bounds_to_wgs84(
    bounds: Sequence[float],
    source_crs: str,
) -> tuple[float, float, float, float]:
    """Return ``bounds`` as a WGS84 lon/lat bounding box.

    The transformation uses all four 2D bbox corners, then returns the min/max
    extent of the transformed corner set. Using only diagonal corners can clip
    transformed extents for projected or rotated coordinate systems.
    """

    source_crs = _validate_source_crs(source_crs)
    xmin, ymin, xmax, ymax = validate_2d_bounds(bounds)
    if is_wgs84_crs(source_crs):
        return xmin, ymin, xmax, ymax

    corners = np.array(
        (
            (xmin, ymin, 0.0),
            (xmin, ymax, 0.0),
            (xmax, ymin, 0.0),
            (xmax, ymax, 0.0),
        ),
        dtype=float,
    )
    try:
        transformed = reproject_array(corners, source_crs, "EPSG:4326")
    except Exception as exc:
        raise ValueError(
            f"Could not transform bounds from {source_crs!r} to EPSG:4326: {exc}"
        ) from exc

    if transformed.shape[0] != 4 or transformed.shape[1] < 2:
        raise ValueError(
            "Bounds transformation returned an invalid coordinate array with "
            f"shape {transformed.shape}."
        )

    xs = transformed[:, 0]
    ys = transformed[:, 1]
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def is_wgs84_crs(crs: str) -> bool:
    """Whether ``crs`` is one of the provider-facing WGS84 aliases."""

    return _validate_source_crs(crs).upper() in WGS84_CRS_ALIASES


def validate_2d_bounds(bounds: Sequence[float]) -> tuple[float, float, float, float]:
    """Validate and normalize a four-number 2D bounding box."""

    if isinstance(bounds, (str, bytes)) or not isinstance(bounds, Sequence):
        raise ValueError("bounds must be a sequence of four numeric values.")
    if len(bounds) != 4:
        raise ValueError("bounds must contain exactly four values: xmin, ymin, xmax, ymax.")

    try:
        xmin, ymin, xmax, ymax = (float(value) for value in bounds)
    except (TypeError, ValueError) as exc:
        raise ValueError("bounds must contain only numeric values.") from exc

    values = (xmin, ymin, xmax, ymax)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("bounds must contain only finite numeric values.")
    if xmin >= xmax or ymin >= ymax:
        raise ValueError("bounds must satisfy xmin < xmax and ymin < ymax.")
    return values


def _validate_source_crs(source_crs: str) -> str:
    if not isinstance(source_crs, str) or not source_crs.strip():
        raise ValueError("source_crs must be a non-empty CRS string.")
    return source_crs.strip()


__all__ = [
    "WGS84_CRS_ALIASES",
    "bounds_to_wgs84",
    "is_wgs84_crs",
    "validate_2d_bounds",
]
