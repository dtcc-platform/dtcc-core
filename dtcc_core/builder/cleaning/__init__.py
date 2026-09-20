"""Deterministic Shapely-based footprint conditioning helpers."""

from .buildings import condition_building_footprints
from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    UnresolvedFootprintCleaningError,
    condition_polygon_coverage,
)
from .plotting import plot_footprint_cleaning_comparison
from .selection import select_footprints

__all__ = [
    "ConditioningOptions",
    "ConditioningResult",
    "UnresolvedFootprintCleaningError",
    "condition_polygon_coverage",
    "condition_building_footprints",
    "select_footprints",
    "plot_footprint_cleaning_comparison",
]
