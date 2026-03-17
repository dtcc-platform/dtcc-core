"""Fixed-precision polygon coverage conditioning for meshing.

The public contract of :func:`condition_polygon_coverage` is intentionally
small:

- output polygons are valid single ``Polygon`` objects
- output polygons are pairwise interior-disjoint
- output order is deterministic
- output ``source_map`` entries are deterministic sorted unique source indices
- geometry is only dropped because of explicit scale rules in
  :class:`ConditioningOptions`

Data-induced geometry failures are converted into diagnostics whenever
possible; only malformed user arguments raise hard errors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, Sequence

import numpy as np
import shapely
from shapely import BufferJoinStyle, GeometryCollection
from shapely.errors import GEOSException
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient
from shapely.ops import polygonize_full, unary_union
from shapely.strtree import STRtree
from shapely.validation import make_valid


@dataclass(slots=True)
class ConditioningOptions:
    precision_grid: float | None = None
    min_feature_size: float = 0.5
    merge_distance: float = 0.5
    min_area: float = 15.0
    min_hole_area: float = 0.25


@dataclass(slots=True)
class ConditioningResult:
    polygons: list[Polygon]
    source_map: list[list[int]]
    diagnostics: dict[str, Any]


def _validate_options(options: ConditioningOptions) -> None:
    for name in (
        "min_feature_size",
        "merge_distance",
        "min_area",
        "min_hole_area",
    ):
        value = getattr(options, name)
        if value < 0:
            raise ValueError(f"{name} must be non-negative, got {value}.")
    if options.precision_grid is not None and options.precision_grid <= 0:
        raise ValueError(
            f"precision_grid must be positive when provided, got {options.precision_grid}."
        )


def _sanitize_source_indices(indices: Sequence[int]) -> list[int]:
    if not indices:
        raise ValueError("source_map entries must not be empty.")
    cleaned: list[int] = []
    for value in indices:
        if not isinstance(value, int):
            raise TypeError("source_map entries must contain integers.")
        cleaned.append(value)
    return sorted(set(cleaned))


def _empty_diagnostics(input_count: int) -> dict[str, Any]:
    return {
        "input_count": input_count,
        "atomic_input_count": 0,
        "output_count": 0,
        "merged_group_count": 0,
        "repaired_invalid_count": 0,
        "collapsed_count": 0,
        "dropped_small_count": 0,
        "overlap_area_before": 0.0,
        "overlap_area_after": 0.0,
        "min_clearance_before": None,
        "min_clearance_after": None,
        "polygonize_cut_edges": 0,
        "polygonize_dangles": 0,
        "polygonize_invalid_rings": 0,
        "geos_exception_count": 0,
        "geos_exception_messages": [],
    }


def _record_geos_exception(
    diagnostics: dict[str, Any],
    step: str,
    exc: BaseException,
) -> None:
    diagnostics["geos_exception_count"] += 1
    message = f"{step}: {exc.__class__.__name__}: {str(exc).strip()}"
    messages = diagnostics["geos_exception_messages"]
    if message not in messages:
        messages.append(message)


def _extract_polygon_parts(geom: BaseGeometry) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []

    polygons: list[Polygon] = []
    stack = [geom]
    while stack:
        current = stack.pop()
        if current.is_empty:
            continue
        if isinstance(current, Polygon):
            polygons.append(current)
            continue
        if hasattr(current, "geoms"):
            stack.extend(reversed(list(current.geoms)))
    return [poly for poly in polygons if not poly.is_empty and poly.area > 0]


def _extract_polygon_parts_with_counts(
    geom: BaseGeometry,
) -> tuple[list[Polygon], int]:
    if geom is None or geom.is_empty:
        return [], 0

    polygons: list[Polygon] = []
    discarded = 0
    stack = [geom]
    while stack:
        current = stack.pop()
        if current.is_empty:
            continue
        if isinstance(current, Polygon):
            polygons.append(current)
            continue
        if hasattr(current, "geoms"):
            stack.extend(reversed(list(current.geoms)))
            continue
        discarded += 1
    return [poly for poly in polygons if not poly.is_empty and poly.area > 0], discarded


def _as_geometry(polygons: Sequence[Polygon]) -> BaseGeometry:
    if not polygons:
        return GeometryCollection()
    if len(polygons) == 1:
        return polygons[0]
    return GeometryCollection(list(polygons))


def _derive_precision_grid(options: ConditioningOptions) -> float:
    if options.precision_grid is not None:
        return options.precision_grid
    positives = [
        value for value in (options.min_feature_size, options.merge_distance) if value > 0
    ]
    if positives:
        return min(positives) / 16.0
    return 0.01


def _make_valid(geom: BaseGeometry, diagnostics: dict[str, Any]) -> BaseGeometry:
    if geom.is_empty:
        return geom
    try:
        if not geom.is_valid:
            diagnostics["repaired_invalid_count"] += 1
        return make_valid(geom)
    except GEOSException as exc:
        _record_geos_exception(diagnostics, "make_valid", exc)
        return geom


def _canonicalize(
    geom: BaseGeometry,
    grid: float,
    diagnostics: dict[str, Any],
) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []

    valid = _make_valid(geom, diagnostics)
    polygon_parts, discarded = _extract_polygon_parts_with_counts(valid)
    diagnostics["discarded_non_polygon_parts"] = diagnostics.get(
        "discarded_non_polygon_parts", 0
    ) + discarded
    if not polygon_parts:
        if not valid.is_empty:
            diagnostics["collapsed_count"] += 1
        return []

    try:
        snapped = shapely.set_precision(
            _as_geometry(polygon_parts),
            grid,
            mode="valid_output",
        )
    except GEOSException as exc:
        _record_geos_exception(diagnostics, "set_precision", exc)
        snapped = _as_geometry(polygon_parts)

    snapped_parts, discarded_after = _extract_polygon_parts_with_counts(snapped)
    diagnostics["discarded_non_polygon_parts"] += discarded_after

    canonical: list[Polygon] = []
    for polygon in snapped_parts:
        if polygon.is_empty or polygon.area <= 0:
            continue
        canonical.append(orient(polygon, sign=1.0))

    if not canonical and not geom.is_empty:
        diagnostics["collapsed_count"] += 1
    return canonical


def _remove_small_holes(
    poly: Polygon,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> Polygon:
    if min_hole_area <= 0:
        return orient(poly, sign=1.0)

    kept_holes: list[list[tuple[float, float]]] = []
    removed = 0
    for ring in poly.interiors:
        hole = Polygon(ring)
        if hole.area < min_hole_area:
            removed += 1
            continue
        kept_holes.append(list(ring.coords))
    if removed:
        diagnostics["dropped_small_count"] += removed
    return orient(Polygon(poly.exterior.coords, kept_holes), sign=1.0)


def _run_constructive_step(
    step: str,
    geom: BaseGeometry,
    operation: Callable[[BaseGeometry], BaseGeometry],
    grid: float,
    diagnostics: dict[str, Any],
) -> BaseGeometry:
    canonical_input = _as_geometry(_canonicalize(geom, grid, diagnostics))
    if canonical_input.is_empty:
        return canonical_input

    for attempt in range(2):
        try:
            return operation(canonical_input)
        except GEOSException as exc:
            _record_geos_exception(diagnostics, step, exc)
            canonical_input = _as_geometry(_canonicalize(canonical_input, grid, diagnostics))
            if canonical_input.is_empty:
                return canonical_input
    return canonical_input


def _apply_opening(
    poly: Polygon,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> list[Polygon]:
    if radius <= 0:
        return _canonicalize(poly, grid, diagnostics)

    opened = _run_constructive_step(
        "opening",
        poly,
        lambda geom: geom.buffer(
            -radius,
            join_style=BufferJoinStyle.mitre,
        ).buffer(
            radius,
            join_style=BufferJoinStyle.mitre,
        ),
        grid,
        diagnostics,
    )
    return _canonicalize(opened, grid, diagnostics)


def _apply_closing(
    geom: BaseGeometry,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> list[Polygon]:
    if radius <= 0:
        return _canonicalize(geom, grid, diagnostics)

    closed = _run_constructive_step(
        "closing",
        geom,
        lambda value: value.buffer(
            radius,
            join_style=BufferJoinStyle.mitre,
        ).buffer(
            -radius,
            join_style=BufferJoinStyle.mitre,
        ),
        grid,
        diagnostics,
    )
    return _canonicalize(closed, grid, diagnostics)


def _build_merge_groups(
    polygons: list[Polygon],
    merge_distance: float,
) -> list[list[int]]:
    if not polygons:
        return []

    parent = list(range(len(polygons)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left == root_right:
            return
        if root_left < root_right:
            parent[root_right] = root_left
        else:
            parent[root_left] = root_right

    tree = STRtree(polygons)
    for index, polygon in enumerate(polygons):
        query_geom = polygon if merge_distance == 0 else polygon.buffer(merge_distance)
        candidate_indices = tree.query(query_geom)
        for candidate in candidate_indices:
            other = int(candidate)
            if other <= index:
                continue
            intersection_area = polygon.intersection(polygons[other]).area
            if intersection_area > 0:
                union(index, other)
                continue
            if merge_distance > 0 and polygon.distance(polygons[other]) <= merge_distance:
                union(index, other)

    groups: dict[int, list[int]] = {}
    for index in range(len(polygons)):
        root = find(index)
        groups.setdefault(root, []).append(index)
    return [sorted(group) for _, group in sorted(groups.items())]


def _count_parts(geom: BaseGeometry) -> int:
    return len(list(shapely.get_parts(geom)))


def _group_tiebreak(source_indices: Sequence[int]) -> int:
    return min(source_indices) if source_indices else 2**31 - 1


def _reconstruct_global_coverage(
    group_geometries: list[BaseGeometry],
    group_sources: list[list[int]],
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[BaseGeometry], list[list[int]]]:
    if not group_geometries:
        return [], []

    linework = [geom.boundary for geom in group_geometries if not geom.is_empty]
    if not linework:
        return [], []

    noded = _run_constructive_step(
        "node",
        _as_geometry(_extract_polygon_parts(unary_union(group_geometries))),
        lambda _: shapely.node(unary_union(linework)),
        grid,
        diagnostics,
    )
    faces, cut_edges, dangles, invalid_rings = polygonize_full(noded)

    diagnostics["polygonize_cut_edges"] = _count_parts(cut_edges)
    diagnostics["polygonize_dangles"] = _count_parts(dangles)
    diagnostics["polygonize_invalid_rings"] = _count_parts(invalid_rings)

    face_polygons = _canonicalize(faces, grid, diagnostics)
    if not face_polygons:
        return group_geometries, group_sources

    tree = STRtree(group_geometries)
    assigned_faces: dict[int, list[Polygon]] = {index: [] for index in range(len(group_geometries))}

    for face in face_polygons:
        representative = face.representative_point()
        candidate_indices = [int(i) for i in tree.query(representative)]
        matches: list[tuple[float, int]] = []

        for candidate in candidate_indices:
            group_geom = group_geometries[candidate]
            if group_geom.covers(representative):
                overlap_area = group_geom.intersection(face).area
                matches.append((overlap_area, candidate))
                continue
            overlap_area = group_geom.intersection(face).area
            if overlap_area > 0:
                matches.append((overlap_area, candidate))

        if not matches:
            candidate_indices = [int(i) for i in tree.query(face)]
            for candidate in candidate_indices:
                overlap_area = group_geometries[candidate].intersection(face).area
                if overlap_area > 0:
                    matches.append((overlap_area, candidate))
        if not matches:
            continue

        best = max(
            matches,
            key=lambda item: (item[0], -_group_tiebreak(group_sources[item[1]])),
        )
        assigned_faces[best[1]].append(face)

    rebuilt_geometries: list[BaseGeometry] = []
    rebuilt_sources: list[list[int]] = []

    for index, faces_for_group in assigned_faces.items():
        if not faces_for_group:
            continue
        try:
            rebuilt = shapely.coverage_union_all(faces_for_group)
        except (AttributeError, GEOSException) as exc:
            if isinstance(exc, GEOSException):
                _record_geos_exception(diagnostics, "coverage_union_all", exc)
            rebuilt = unary_union(faces_for_group)
            diagnostics.setdefault("coverage_union_fallback_groups", []).append(index)

        canonical = _canonicalize(rebuilt, grid, diagnostics)
        if not canonical:
            continue
        rebuilt_geometries.append(_as_geometry(canonical))
        rebuilt_sources.append(group_sources[index])

    return rebuilt_geometries, rebuilt_sources


def _stable_sort(
    polygons: list[Polygon],
    source_map: list[list[int]],
) -> tuple[list[Polygon], list[list[int]]]:
    paired = list(zip(polygons, source_map))
    paired.sort(
        key=lambda item: (
            item[0].bounds[0],
            item[0].bounds[1],
            -item[0].area,
            tuple(item[1]),
        )
    )
    if not paired:
        return [], []
    sorted_polygons, sorted_sources = zip(*paired)
    return list(sorted_polygons), list(sorted_sources)


def _coverage_overlap_area(polygons: Sequence[Polygon]) -> float:
    if not polygons:
        return 0.0
    total_area = sum(polygon.area for polygon in polygons)
    union_area = unary_union(polygons).area
    return float(max(total_area - union_area, 0.0))


def _minimum_clearance(polygons: Sequence[Polygon]) -> float | None:
    values: list[float] = []
    for polygon in polygons:
        try:
            clearance = shapely.minimum_clearance(polygon)
        except GEOSException:
            continue
        if np.isfinite(clearance):
            values.append(float(clearance))
    if not values:
        return None
    return min(values)


def condition_polygon_coverage(
    polygons: Sequence[BaseGeometry],
    *,
    source_map: Sequence[Sequence[int]] | None = None,
    options: ConditioningOptions,
) -> ConditioningResult:
    """Condition a noisy polygon coverage into a deterministic meshing input.

    The conditioner uses only Shapely / GEOS operations, regularizes geometry at
    an explicit spatial scale, reconstructs the final global coverage from
    noded boundary linework, and returns valid, non-overlapping single polygons.
    Only malformed user arguments raise hard errors; data-induced failures are
    converted into diagnostics and the pipeline continues with the best valid
    geometry available at each stage.
    """

    _validate_options(options)
    diagnostics = _empty_diagnostics(len(polygons))

    if source_map is not None and len(source_map) != len(polygons):
        raise ValueError(
            "source_map length must match the number of input geometries."
        )

    initial_sources = (
        [_sanitize_source_indices(indices) for indices in source_map]
        if source_map is not None
        else [[index] for index in range(len(polygons))]
    )

    grid = _derive_precision_grid(options)

    atomic_polygons: list[Polygon] = []
    atomic_sources: list[list[int]] = []

    for index, geometry in enumerate(polygons):
        if geometry is None:
            continue
        if not isinstance(geometry, BaseGeometry):
            raise TypeError("polygons must contain Shapely geometries.")

        valid = _make_valid(geometry, diagnostics)
        parts, discarded = _extract_polygon_parts_with_counts(valid)
        diagnostics["discarded_non_polygon_parts"] = diagnostics.get(
            "discarded_non_polygon_parts", 0
        ) + discarded
        if not parts and not valid.is_empty:
            diagnostics["collapsed_count"] += 1
        for part in parts:
            atomic_polygons.append(part)
            atomic_sources.append(initial_sources[index])

    diagnostics["atomic_input_count"] = len(atomic_polygons)
    diagnostics["overlap_area_before"] = _coverage_overlap_area(atomic_polygons)
    diagnostics["min_clearance_before"] = _minimum_clearance(atomic_polygons)

    opened_polygons: list[Polygon] = []
    opened_sources: list[list[int]] = []
    r_open = options.min_feature_size / 2.0

    for polygon, indices in zip(atomic_polygons, atomic_sources):
        parts = _apply_opening(polygon, r_open, grid, diagnostics)
        filtered_parts: list[Polygon] = []
        for part in parts:
            without_small_holes = _remove_small_holes(
                part,
                options.min_hole_area,
                diagnostics,
            )
            filtered_parts.extend(_canonicalize(without_small_holes, grid, diagnostics))
        for part in filtered_parts:
            opened_polygons.append(part)
            opened_sources.append(indices)

    merge_groups = _build_merge_groups(opened_polygons, options.merge_distance)
    diagnostics["merged_group_count"] = len(merge_groups)

    regularized_group_geometries: list[BaseGeometry] = []
    regularized_group_sources: list[list[int]] = []

    r_close = options.merge_distance / 2.0
    for group in merge_groups:
        group_polygons = [opened_polygons[index] for index in group]
        group_sources = sorted(
            set(source for index in group for source in opened_sources[index])
        )
        unioned = _run_constructive_step(
            "group_union",
            _as_geometry(group_polygons),
            lambda geom: unary_union(_extract_polygon_parts(geom)),
            grid,
            diagnostics,
        )
        closed_parts = _apply_closing(unioned, r_close, grid, diagnostics)

        filtered_parts: list[Polygon] = []
        for part in closed_parts:
            without_small_holes = _remove_small_holes(
                part,
                options.min_hole_area,
                diagnostics,
            )
            filtered_parts.extend(_canonicalize(without_small_holes, grid, diagnostics))

        if not filtered_parts:
            continue
        regularized_group_geometries.append(_as_geometry(filtered_parts))
        regularized_group_sources.append(group_sources)

    rebuilt_geometries, rebuilt_sources = _reconstruct_global_coverage(
        regularized_group_geometries,
        regularized_group_sources,
        grid,
        diagnostics,
    )

    final_polygons: list[Polygon] = []
    final_sources: list[list[int]] = []
    for geometry, indices in zip(rebuilt_geometries, rebuilt_sources):
        for polygon in _canonicalize(geometry, grid, diagnostics):
            if polygon.area < options.min_area:
                diagnostics["dropped_small_count"] += 1
                continue
            final_polygons.append(orient(polygon, sign=1.0))
            final_sources.append(sorted(set(indices)))

    final_polygons, final_sources = _stable_sort(final_polygons, final_sources)
    diagnostics["output_count"] = len(final_polygons)
    diagnostics["overlap_area_after"] = _coverage_overlap_area(final_polygons)
    diagnostics["min_clearance_after"] = _minimum_clearance(final_polygons)

    return ConditioningResult(
        polygons=final_polygons,
        source_map=final_sources,
        diagnostics=diagnostics,
    )
