"""Explicit area selection after geometric cleaning."""

from __future__ import annotations

import math

from .footprints import ConditioningResult


def select_footprints(
    result: ConditioningResult, *, min_area: float
) -> ConditioningResult:
    """Keep cleaned regions of at least min_area, without changing their geometry.

    Selection is not geometric repair and does not alter the cleaning reference.
    Return a new result, recording exclusions separately from cleaning defects.
    Source indices remain in the original input's index space.
    """
    if not math.isfinite(min_area) or min_area < 0:
        raise ValueError("min_area must be finite and nonnegative")
    kept, excluded = [], []
    for i, polygon in enumerate(result.polygons):
        (kept if polygon.area >= min_area else excluded).append(i)
    represented = {j for i in kept for j in result.source_map[i]}
    removed_sources = {j for i in excluded for j in result.source_map[i]}
    policy_sources = {
        source
        for exclusion in result.diagnostics.get("policy_exclusions", [])
        for source in exclusion.get("source_indices", [])
    }
    all_sources = {source for group in result.source_map for source in group} | policy_sources
    selection = {
        "min_area": min_area,
        "input_count": len(result.polygons),
        "output_count": len(kept),
        "removed_count": len(excluded),
        "removed_area": sum(result.polygons[i].area for i in excluded),
        "removed_source_map": [list(result.source_map[i]) for i in excluded],
        "excluded_regions": [
            {
                "region_index": i,
                "source_indices": list(result.source_map[i]),
                "area": float(result.polygons[i].area),
                "reason": "minimum_area",
            }
            for i in excluded
        ],
        "represented_source_indices": sorted(represented),
        "unrepresented_source_indices": sorted(all_sources - represented),
    }
    return ConditioningResult(
        polygons=[result.polygons[i] for i in kept],
        source_map=[list(result.source_map[i]) for i in kept],
        diagnostics={**result.diagnostics, "selection": selection},
    )
