"""Roof levels of merged blocks at the vertical resolution delta.

Two parts of a block whose roofs differ by less than delta would get a step
wall lower than the resolution, and such differences are below the accuracy of
the height estimates. Before meshing, the roof levels of each block are
therefore grouped (:func:`group_levels`):

* sort the levels and group them greedily from the lowest, each group spanning
  at most delta;
* give the groups, in increasing order, the values

      eta_c = max(clamp(w_c, [max_c - delta / 2, min_c + delta / 2]), eta_{c-1} + delta),

  where w_c is the area-weighted mean of group c.

Every level moves by at most delta / 2, their order is kept, and distinct
levels end up at least delta apart; levels at least delta apart are never
unified. Walls between cells of the same level are then dissolved
(:func:`level_regions`), which never adds a feature pair or a small sector.

The union of the extruded regions is a manifold solid only if every level set
Q_t, the union of the regions with roofs at or above t, is a 2-manifold with
boundary (the level-union profile). :func:`pinched_levels` reports the level
sets where regions touch at a point only, which gives a vertical line contact
between buildings.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import shapely
from shapely.geometry import Polygon

from .construction import polygon_parts
from .contract import essential_boundary_graph


def group_levels(
    levels: Sequence[float], areas: Sequence[float], delta: float
) -> list[float]:
    """Group roof levels closer than delta (see the module documentation).

    Parameters
    ----------
    levels : sequence of float
        Roof levels of the parts of one block.
    areas : sequence of float
        Their footprint areas, the weights of the group means.
    delta : float
        Vertical resolution; zero groups equal levels only.

    Returns
    -------
    list[float]
        The grouped level of each input level.
    """
    levels = np.asarray(levels, dtype=float)
    weights = np.asarray(areas, dtype=float)
    weights = np.where(weights > 0, weights, 1.0)
    order = np.argsort(levels, kind="stable")
    sorted_levels = levels[order]
    grouped = np.empty(len(levels))
    previous = -np.inf
    start = 0
    for end in range(1, len(levels) + 1):
        if end < len(levels) and sorted_levels[end] <= sorted_levels[start] + delta:
            continue
        members = order[start:end]
        low, high = sorted_levels[start], sorted_levels[end - 1]
        if low == high:
            mean = low
        else:
            mean = float(np.average(levels[members], weights=weights[members]))
        value = max(min(max(mean, high - delta / 2), low + delta / 2), previous + delta)
        grouped[members] = value
        previous = value
        start = end
    return grouped.tolist()


def level_regions(
    cells: Sequence[tuple[Polygon, int]],
    levels: Mapping[int, float],
    delta: float,
) -> list[tuple[Polygon, float, list[int]]]:
    """Dissolve the cells of one block into one region per grouped level.

    Parameters
    ----------
    cells : sequence of (Polygon, label)
        Cells that tile the block; cells of one label may have several pieces.
    levels : mapping of label to float
        Roof level of each label.
    delta : float
        Vertical resolution (see :func:`group_levels`).

    Returns
    -------
    list of (Polygon, float, list[int])
        Each region (a connected component of the cells at one grouped
        level), its level and its labels, ordered by level and geometry.
    """
    labels = sorted({label for _, label in cells})
    areas = {label: 0.0 for label in labels}
    for polygon, label in cells:
        areas[label] += polygon.area
    grouped = dict(
        zip(
            labels,
            group_levels([levels[label] for label in labels], [areas[label] for label in labels], delta),
        )
    )
    by_level: dict[float, list[tuple[Polygon, int]]] = {}
    for polygon, label in cells:
        by_level.setdefault(grouped[label], []).append((polygon, label))
    regions = []
    for value, members in by_level.items():
        polygons = sorted((polygon.normalize() for polygon, _ in members), key=lambda p: p.wkb)
        union = polygons[0] if len(polygons) == 1 else shapely.coverage_union_all(polygons)
        if not union.is_valid:
            union = shapely.unary_union(polygons)
        for part in polygon_parts(union):
            part = part.normalize()
            inside = sorted(
                {label for polygon, label in members if polygon.intersection(part).area > 0}
            )
            regions.append((part, float(value), inside))
    regions.sort(key=lambda item: (item[1], item[0].wkb))
    return regions


def pinched_levels(regions: Sequence[tuple[Polygon, float]]) -> list[dict]:
    """Level sets of one block whose regions touch at a point only.

    Parameters
    ----------
    regions : sequence of (Polygon, float)
        The regions of one block and their roof levels.

    Returns
    -------
    list[dict]
        For each pinched level set above the lowest level, the level and the
        pinch points.
    """
    values = sorted({value for _, value in regions})
    found = []
    for value in values[1:]:
        polygons = [polygon for polygon, level in regions if level >= value]
        union = shapely.unary_union(polygons)
        _, edges = essential_boundary_graph([union.boundary])
        degree: dict[tuple[float, float], int] = {}
        for a, b in edges:
            degree[a] = degree.get(a, 0) + 1
            degree[b] = degree.get(b, 0) + 1
        points = sorted(vertex for vertex, count in degree.items() if count != 2)
        if points:
            found.append({"level": value, "points": [list(point) for point in points]})
    return found
