"""Deterministic Shapely-based footprint conditioning helpers."""

from .buildings import condition_building_footprints
from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    condition_polygon_coverage,
)
from .plotting import plot_footprint_cleaning_comparison

__all__ = [
    "ConditioningOptions",
    "ConditioningResult",
    "condition_polygon_coverage",
    "condition_building_footprints",
    "plot_footprint_cleaning_comparison",
]
