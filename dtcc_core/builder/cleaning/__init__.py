"""Footprint cleaning.

:func:`clean_footprints` is the single entry point and :func:`building_cells`
reads its result; the other modules are internal building blocks.
"""

from .clean import FootprintCell, building_cells, clean_footprints
from .footprints import UnresolvedFootprintCleaningError
from .plotting import plot_footprint_cleaning_comparison

__all__ = [
    "clean_footprints",
    "building_cells",
    "FootprintCell",
    "UnresolvedFootprintCleaningError",
    "plot_footprint_cleaning_comparison",
]
