"""Building adapters for the footprint conditioning API."""

from __future__ import annotations

from typing import Sequence

from shapely.geometry.base import BaseGeometry

from dtcc_core.model import Building, GeometryType

from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    condition_polygon_coverage,
)


def condition_building_footprints(
    buildings: Sequence[Building],
    *,
    lod: GeometryType = GeometryType.LOD0,
    options: ConditioningOptions,
) -> ConditioningResult:
    """Condition building footprints as a regularized polygon coverage.

    The adapter extracts one footprint geometry per building with
    ``building.flatten_geometry(lod).to_polygon()``, ignores Z during
    conditioning, preserves original building indices in ``source_map``, and
    returns the same deterministic guarantees as
    :func:`condition_polygon_coverage`.
    """

    footprints: list[BaseGeometry] = []
    source_map: list[list[int]] = []

    for index, building in enumerate(buildings):
        geometry = building.flatten_geometry(lod)
        if geometry is None:
            continue
        polygon = geometry.to_polygon(simplify=0.0)
        footprints.append(polygon)
        source_map.append([index])

    return condition_polygon_coverage(
        footprints,
        source_map=source_map,
        options=options,
    )
