"""Deterministic Shapely-based footprint conditioning helpers."""

from .buildings import condition_building_footprints
from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    condition_polygon_coverage,
)

__all__ = [
    "ConditioningOptions",
    "ConditioningResult",
    "condition_polygon_coverage",
    "condition_building_footprints",
]
