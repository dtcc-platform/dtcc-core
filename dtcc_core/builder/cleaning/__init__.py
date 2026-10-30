"""Footprint cleaning.

:func:`clean_footprints` is the single entry point; the other modules are
internal building blocks.
"""

from .clean import clean_footprints
from .footprints import UnresolvedFootprintCleaningError
from .plotting import plot_footprint_cleaning_comparison

__all__ = [
    "clean_footprints",
    "UnresolvedFootprintCleaningError",
    "plot_footprint_cleaning_comparison",
]
