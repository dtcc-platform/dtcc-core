"""Fixed-precision polygon coverage conditioning for meshing.

The cleaner treats the input as a noisy polygon coverage and applies one
deterministic scale contract:

- ``min_feature_size`` removes positive features below the declared scale via a
  regularized opening.
- ``merge_distance`` removes negative gaps below the declared scale via a
  regularized closing.
- ``precision_grid`` defines the fixed-precision model used after every
  constructive step.
- a local defect-repair stage handles courtyard passages and narrow
  self-clearance defects before shared-boundary simplification.
- a shared-boundary coverage-simplification stage removes meshing-hostile
  short edges coverage-wide at the declared scale.
- a final boundary-regularization stage only touches residual short edges that
  survive the shared-boundary simplification, and only if every accepted
  change stays local to the short edges it is meant to remove.
- a source-coordinate recovery stage restores original supported vertex
  positions whenever doing so preserves the cleaned topology and declared-scale
  quality contract.
- any output polygon that still violates the declared minimum-clearance scale is
  regularized once more with the same morphological radius before meshing.

The public contract of :func:`condition_polygon_coverage` is intentionally
small:

- output polygons are valid single ``Polygon`` objects
- output polygons are pairwise interior-disjoint
- output order is deterministic
- output ``source_map`` entries are deterministic sorted unique source indices
- geometry is only dropped because of explicit scale rules in
  :class:`ConditioningOptions`
- local defect repair is monotone: accepted local candidates must strictly
  improve the polygon defect signature at the declared scale
- simplification is locality-preserving: accepted simplification changes must
  stay inside a ``2 * min_feature_size`` neighborhood of edges shorter than
  ``min_feature_size``
- simplification is area-balanced: accepted simplification changes must keep
  signed area drift bounded relative to the total local boundary change

Data-induced geometry failures are converted into diagnostics whenever
possible; only malformed user arguments raise hard errors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Literal, Sequence

import numpy as np
import shapely
from shapely import BufferJoinStyle, GeometryCollection
from shapely.errors import GEOSException
from shapely.geometry import LineString, Point, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient
from shapely.ops import polygonize_full, unary_union
from shapely.strtree import STRtree
from shapely.validation import make_valid

try:
    from .. import _dtcc_builder
except Exception:  # pragma: no cover - import fallback for partial builds
    _dtcc_builder = None

from ..logging import info


@dataclass(slots=True)
class ConditioningOptions:
    precision_grid: float | None = None
    min_feature_size: float = 0.5
    merge_distance: float = 0.5
    min_area: float = 15.0
    min_hole_area: float = 0.25
    collect_stage_metrics: bool = True
    enable_logging: bool = True


@dataclass(slots=True)
class ConditioningResult:
    polygons: list[Polygon]
    source_map: list[list[int]]
    diagnostics: dict[str, Any]


@dataclass(slots=True)
class _PolygonDefectSignature:
    clearance: float | None
    clearance_deficit: float
    short_edge_count: int
    min_edge_length: float | None
    vertex_count: int


@dataclass(slots=True)
class _CoverageDefectSignature:
    min_clearance: float | None
    pair_issue_count: int
    point_touch_count: int
    close_pair_count: int
    min_pair_clearance: float | None
    short_edge_count: int
    min_edge_length: float | None
    vertex_count: int
    ring_contact_count: int = 0


@dataclass(slots=True)
class _RepairCandidate:
    polygon: Polygon
    edit_zone: BaseGeometry
    operator: str


@dataclass(slots=True)
class _CoverageSimplifyCandidate:
    label: str
    polygons: list[Polygon]
    source_map: list[list[int]]
    signature: _CoverageDefectSignature
    difference_metrics: dict[str, float]
    change_outside_edit_zone: float
    edit_zone_area: float
    area_balance_budget: float
    patch_count: int
    patch_applied_count: int
    operator_attempts: dict[str, int]
    operator_applied: dict[str, int]


@dataclass(slots=True)
class _PostCoverageBranch:
    label: str
    coverage_candidate: _CoverageSimplifyCandidate
    coverage_polygons: list[Polygon]
    coverage_source_map: list[list[int]]
    source_reclaimed_polygons: list[Polygon]
    source_reclaimed_source_map: list[list[int]]
    small_component_absorbed_polygons: list[Polygon]
    small_component_absorbed_source_map: list[list[int]]
    boundary_regularized_polygons: list[Polygon]
    boundary_regularized_source_map: list[list[int]]
    clearance_regularized_polygons: list[Polygon]
    clearance_regularized_source_map: list[list[int]]
    source_coordinate_recovered_polygons: list[Polygon]
    source_coordinate_recovered_source_map: list[list[int]]
    coverage_meshing_regularized_polygons: list[Polygon]
    coverage_meshing_regularized_source_map: list[list[int]]
    final_output_polygons: list[Polygon]
    final_output_source_map: list[list[int]]
    final_signature: _CoverageDefectSignature
    final_difference_metrics: dict[str, float]
    diagnostics: dict[str, Any]


@dataclass(slots=True)
class _CoverageEvalCache:
    union_cache: dict[tuple[int, ...], BaseGeometry] = field(default_factory=dict)
    overlap_area_cache: dict[tuple[int, ...], float] = field(default_factory=dict)
    edit_zone_cache: dict[tuple[tuple[int, ...], float, float], BaseGeometry] = (
        field(default_factory=dict)
    )
    signature_cache: dict[tuple[tuple[int, ...], float], _CoverageDefectSignature] = (
        field(default_factory=dict)
    )
    pair_issue_candidates_cache: dict[
        tuple[tuple[int, ...], float],
        list[tuple[int, int, float, Literal["point", "close"]]],
    ] = field(default_factory=dict)
    defect_cluster_cache: dict[
        tuple[tuple[int, ...], float, float],
        list[list[int]],
    ] = field(default_factory=dict)
    boundary_payload_cache: dict[
        tuple[int, ...],
        list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]],
    ] = field(default_factory=dict)
    boundary_descriptor_cache: dict[
        tuple[tuple[int, ...], float, float],
        list[dict[str, Any]],
    ] = field(default_factory=dict)
    difference_metrics_cache: dict[
        tuple[tuple[int, ...], tuple[int, ...]], dict[str, float]
    ] = field(default_factory=dict)


_AREA_BALANCE_RELATIVE_TOLERANCE = 0.35
_AREA_BALANCE_ABSOLUTE_GRID_MULTIPLIER = 4.0
_AREA_BALANCE_SHORT_EDGE_MULTIPLIER = 2.0
_COURTYARD_AREA_THRESHOLD = 10.0
_RECOVERY_TREE_REBUILD_THRESHOLD = 16


def _polygon_sequence_key(polygons: Sequence[Polygon]) -> tuple[int, ...]:
    return tuple(id(polygon) for polygon in polygons)


def _cached_union(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
) -> BaseGeometry:
    if cache is None:
        return unary_union(polygons)

    key = _polygon_sequence_key(polygons)
    cached = cache.union_cache.get(key)
    if cached is not None:
        return cached

    value = unary_union(polygons)
    cache.union_cache[key] = value
    return value


def _cached_overlap_area(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
) -> float:
    if cache is None:
        return _coverage_overlap_area(list(polygons))

    key = _polygon_sequence_key(polygons)
    cached = cache.overlap_area_cache.get(key)
    if cached is not None:
        return cached

    total_area = float(sum(polygon.area for polygon in polygons))
    overlap_area = float(max(total_area - _cached_union(cache, polygons).area, 0.0))
    cache.overlap_area_cache[key] = overlap_area
    return overlap_area


def _cached_coverage_defect_signature(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
) -> _CoverageDefectSignature:
    if cache is None:
        return _coverage_defect_signature(polygons, target_scale=target_scale)

    key = (_polygon_sequence_key(polygons), float(target_scale))
    cached = cache.signature_cache.get(key)
    if cached is not None:
        return cached

    value = _coverage_defect_signature(polygons, target_scale=target_scale)
    cache.signature_cache[key] = value
    return value


def _cached_pair_issue_candidates(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
) -> list[tuple[int, int, float, Literal["point", "close"]]]:
    if cache is None:
        return _pair_issue_candidates(polygons, target_scale=target_scale)

    key = (_polygon_sequence_key(polygons), float(target_scale))
    cached = cache.pair_issue_candidates_cache.get(key)
    if cached is not None:
        return cached

    value = _pair_issue_candidates(polygons, target_scale=target_scale)
    cache.pair_issue_candidates_cache[key] = value
    return value


def _cached_coverage_defect_clusters(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    cluster_radius: float,
) -> list[list[int]]:
    if cache is None:
        return _coverage_defect_clusters(
            polygons,
            target_scale=target_scale,
            cluster_radius=cluster_radius,
        )

    key = (_polygon_sequence_key(polygons), float(target_scale), float(cluster_radius))
    cached = cache.defect_cluster_cache.get(key)
    if cached is not None:
        return cached

    value = _coverage_defect_clusters(
        polygons,
        target_scale=target_scale,
        cluster_radius=cluster_radius,
        cache=cache,
    )
    cache.defect_cluster_cache[key] = value
    return value


def _cached_coverage_defect_cluster_descriptors(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    cluster_radius: float,
) -> list[dict[str, Any]]:
    # Large coverage snapshots are where Python-side defect clustering starts to
    # become noticeable. Use the native helper only there and keep the existing
    # Python path for smaller inputs, where it is simpler and just as fast.
    if len(polygons) >= 64:
        builder_clusters = _builder_boundary_defect_clusters(
            polygons,
            target_scale=target_scale,
            pair_tolerance=cluster_radius,
            cache=cache,
        )
        if builder_clusters is not None:
            return builder_clusters

    clusters = _cached_coverage_defect_clusters(
        cache,
        polygons,
        target_scale=target_scale,
        cluster_radius=cluster_radius,
    )
    descriptors: list[dict[str, Any]] = []
    for indices in clusters:
        subset_polygons = [polygons[index] for index in indices]
        signature = _cached_coverage_defect_signature(
            cache,
            subset_polygons,
            target_scale=target_scale,
        )
        kind = "short_edge_only"
        if signature.pair_issue_count > 0 and signature.short_edge_count > 0:
            kind = "mixed_pair_short_edge"
        elif signature.point_touch_count > 0 and signature.close_pair_count == 0:
            kind = "point_touch_pair"
        elif signature.pair_issue_count > 0:
            kind = "close_pair"
        descriptors.append(
            {
                "indices": list(indices),
                "kind": kind,
                "short_edge_count": signature.short_edge_count,
                "pair_issue_count": signature.pair_issue_count,
            }
        )
    return descriptors


def _zero_difference_metrics() -> dict[str, float]:
    return {
        "reference_minus_candidate_area": 0.0,
        "candidate_minus_reference_area": 0.0,
        "symmetric_difference_area": 0.0,
        "union_area_delta": 0.0,
    }


def _cached_difference_area_metrics(
    cache: _CoverageEvalCache | None,
    reference_polygons: Sequence[Polygon],
    candidate_polygons: Sequence[Polygon],
) -> dict[str, float]:
    if cache is None:
        return _difference_area_metrics(
            unary_union(reference_polygons),
            unary_union(candidate_polygons),
        )

    reference_key = _polygon_sequence_key(reference_polygons)
    candidate_key = _polygon_sequence_key(candidate_polygons)
    if reference_key == candidate_key:
        return _zero_difference_metrics()

    key = (reference_key, candidate_key)
    cached = cache.difference_metrics_cache.get(key)
    if cached is not None:
        return cached

    value = _difference_area_metrics(
        _cached_union(cache, reference_polygons),
        _cached_union(cache, candidate_polygons),
    )
    cache.difference_metrics_cache[key] = value
    return value


def _cached_coverage_edit_zone(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    radius: float,
) -> BaseGeometry:
    if cache is None:
        return _coverage_edit_zone(
            polygons,
            target_scale=target_scale,
            radius=radius,
        )

    key = (_polygon_sequence_key(polygons), float(target_scale), float(radius))
    cached = cache.edit_zone_cache.get(key)
    if cached is not None:
        return cached

    value = _coverage_edit_zone(
        polygons,
        target_scale=target_scale,
        radius=radius,
    )
    cache.edit_zone_cache[key] = value
    return value


def _polygon_boundary_payload(
    polygon: Polygon,
) -> tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]:
    return (
        [(float(x), float(y)) for x, y, *_ in polygon.exterior.coords],
        [
            [(float(x), float(y)) for x, y, *_ in ring.coords]
            for ring in polygon.interiors
        ],
    )


def _cached_polygon_boundary_payload(
    cache: _CoverageEvalCache | None,
    polygons: Sequence[Polygon],
) -> list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]]:
    if cache is None:
        return [_polygon_boundary_payload(polygon) for polygon in polygons]

    key = _polygon_sequence_key(polygons)
    cached = cache.boundary_payload_cache.get(key)
    if cached is not None:
        return cached

    payload = [_polygon_boundary_payload(polygon) for polygon in polygons]
    cache.boundary_payload_cache[key] = payload
    return payload


def _builder_boundary_defect_clusters(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    pair_tolerance: float,
    cache: _CoverageEvalCache | None = None,
) -> list[dict[str, Any]] | None:
    if _dtcc_builder is None:
        return None
    if cache is not None:
        key = (_polygon_sequence_key(polygons), float(target_scale), float(pair_tolerance))
        cached = cache.boundary_descriptor_cache.get(key)
        if cached is not None:
            return cached
    try:
        clusters = list(
            _dtcc_builder.boundary_defect_clusters(
                _cached_polygon_boundary_payload(cache, polygons),
                float(target_scale),
                float(pair_tolerance),
            )
        )
    except Exception:
        return None

    normalized_clusters: list[dict[str, Any]] = []
    for cluster in clusters:
        normalized_clusters.append(
            {
                "indices": [int(index) for index in cluster.get("indices", [])],
                "kind": str(cluster.get("kind", "short_edge_only")),
                "short_edge_count": int(cluster.get("short_edge_count", 0)),
                "pair_issue_count": int(cluster.get("pair_issue_count", 0)),
            }
        )
    if cache is not None:
        cache.boundary_descriptor_cache[key] = normalized_clusters
    return normalized_clusters


def _builder_rewrite_defect_cluster(
    polygons: Sequence[Polygon],
    cluster: dict[str, Any],
    *,
    target_scale: float,
    grid: float,
    cache: _CoverageEvalCache | None = None,
) -> tuple[list[int], list[Polygon]] | None:
    if _dtcc_builder is None:
        return None
    try:
        payload = _dtcc_builder.rewrite_defect_cluster(
            _cached_polygon_boundary_payload(cache, polygons),
            cluster,
            float(target_scale),
            float(grid),
        )
    except Exception:
        return None
    if not payload or not bool(payload.get("supported", False)):
        return None

    changed_indices = [int(index) for index in payload.get("changed_indices", [])]
    polygon_payloads = list(payload.get("polygons", []))
    if not changed_indices or len(changed_indices) != len(polygon_payloads):
        return None

    rewritten_polygons: list[Polygon] = []
    for shell, holes in polygon_payloads:
        rewritten_polygons.append(Polygon(shell, holes))
    return changed_indices, rewritten_polygons


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
        "discarded_non_polygon_parts": 0,
        "enable_logging": True,
        "precision_grid": None,
        "output_grid": None,
        "opening_radius": 0.0,
        "closing_radius": 0.0,
        "output_lattice_applied": False,
        "output_canonicalization_candidate_count": 0,
        "output_canonicalization_output_grid_applied_count": 0,
        "output_canonicalization_fallback_count": 0,
        "output_canonicalization_reference_minus_candidate_area": 0.0,
        "output_canonicalization_candidate_minus_reference_area": 0.0,
        "output_canonicalization_signed_area_delta": 0.0,
        "output_canonicalization_area_balance_budget": 0.0,
        "final_min_area_filter_removed_count": 0,
        "final_min_area_filter_removed_area": 0.0,
        "global_reconstruction_applied": False,
        "global_reconstruction_reason": "not_evaluated",
        "global_reconstruction_overlap_threshold": None,
        "multi_component_group_count": 0,
        "component_source_reassignment_count": 0,
        "clearance_regularization_threshold": None,
        "clearance_regularization_applied": False,
        "clearance_regularization_candidate_count": 0,
        "clearance_regularization_improved_count": 0,
        "clearance_regularization_failed_count": 0,
        "clearance_regularization_overlap_area": 0.0,
        "source_coordinate_recovery_applied": False,
        "source_coordinate_recovery_candidate_count": 0,
        "source_coordinate_recovery_applied_count": 0,
        "source_coordinate_recovery_exact_count": 0,
        "source_coordinate_recovery_vertex_count": 0,
        "source_coordinate_recovery_short_edge_count_before": 0,
        "source_coordinate_recovery_short_edge_count_after": 0,
        "source_coordinate_recovery_reference_minus_candidate_area": 0.0,
        "source_coordinate_recovery_candidate_minus_reference_area": 0.0,
        "source_coordinate_recovery_signed_area_delta": 0.0,
        "source_coordinate_recovery_support_reference_minus_candidate_area": 0.0,
        "source_coordinate_recovery_support_candidate_minus_reference_area": 0.0,
        "source_coordinate_recovery_support_signed_area_delta": 0.0,
        "source_coordinate_recovery_rejected_overlap_count": 0,
        "source_coordinate_recovery_rejected_non_improving_count": 0,
        "source_coordinate_recovery_operator_attempts": {},
        "source_coordinate_recovery_operator_applied": {},
        "coverage_meshing_regularization_applied": False,
        "coverage_meshing_regularization_selected_branch": "identity",
        "coverage_meshing_regularization_short_edge_count_before": 0,
        "coverage_meshing_regularization_short_edge_count_after": 0,
        "coverage_meshing_regularization_pair_issue_count_before": 0,
        "coverage_meshing_regularization_pair_issue_count_after": 0,
        "coverage_meshing_regularization_ring_contact_count_before": 0,
        "coverage_meshing_regularization_ring_contact_count_after": 0,
        "coverage_meshing_regularization_ring_contact_polygon_count": 0,
        "coverage_meshing_regularization_ring_contact_component_count": 0,
        "coverage_meshing_regularization_ring_contact_failed_count": 0,
        "coverage_meshing_regularization_ring_contact_area_delta": 0.0,
        "coverage_meshing_regularization_reference_minus_candidate_area": 0.0,
        "coverage_meshing_regularization_candidate_minus_reference_area": 0.0,
        "coverage_meshing_regularization_signed_area_delta": 0.0,
        "coverage_meshing_regularization_operator_attempts": {},
        "coverage_meshing_regularization_operator_applied": {},
        "small_component_absorb_applied": False,
        "small_component_absorb_candidate_count": 0,
        "small_component_absorb_applied_count": 0,
        "small_component_absorb_component_count": 0,
        "small_component_absorb_edit_zone_area": 0.0,
        "small_component_absorb_change_outside_edit_zone": 0.0,
        "small_component_absorb_reference_minus_candidate_area": 0.0,
        "small_component_absorb_candidate_minus_reference_area": 0.0,
        "small_component_absorb_signed_area_delta": 0.0,
        "small_component_absorb_rejected_nonlocal_count": 0,
        "small_component_absorb_rejected_non_improving_count": 0,
        "small_component_absorb_operator_attempts": {},
        "small_component_absorb_operator_applied": {},
        "coverage_simplify_short_edge_count_before": 0,
        "coverage_simplify_short_edge_count_after": 0,
        "coverage_simplify_edit_zone_area": 0.0,
        "coverage_simplify_change_outside_edit_zone": 0.0,
        "coverage_simplify_reference_minus_candidate_area": 0.0,
        "coverage_simplify_candidate_minus_reference_area": 0.0,
        "coverage_simplify_signed_area_delta": 0.0,
        "coverage_simplify_area_balance_budget": 0.0,
        "coverage_simplify_rejected_nonlocal": False,
        "coverage_simplify_rejected_area_imbalance": False,
        "coverage_simplify_patch_count": 0,
        "coverage_simplify_patch_applied_count": 0,
        "coverage_simplify_operator_attempts": {},
        "coverage_simplify_operator_applied": {},
        "coverage_simplify_selected_branch": "identity",
        "coverage_simplify_branch_scores": {},
        "local_defect_repair_candidate_count": 0,
        "local_defect_repair_applied_count": 0,
        "local_defect_repair_short_edge_count_before": 0,
        "local_defect_repair_short_edge_count_after": 0,
        "local_defect_repair_edit_zone_area": 0.0,
        "local_defect_repair_change_outside_edit_zone": 0.0,
        "local_defect_repair_reference_minus_candidate_area": 0.0,
        "local_defect_repair_candidate_minus_reference_area": 0.0,
        "local_defect_repair_signed_area_delta": 0.0,
        "local_defect_repair_area_balance_budget": 0.0,
        "local_defect_repair_rejected_nonlocal_count": 0,
        "local_defect_repair_rejected_area_imbalance_count": 0,
        "local_defect_repair_rejected_non_improving_count": 0,
        "local_defect_repair_operator_attempts": {},
        "local_defect_repair_operator_applied": {},
        "source_reclaim_candidate_count": 0,
        "source_reclaim_applied_count": 0,
        "source_reclaim_component_count": 0,
        "source_reclaim_short_edge_count_before": 0,
        "source_reclaim_short_edge_count_after": 0,
        "source_reclaim_edit_zone_area": 0.0,
        "source_reclaim_change_outside_edit_zone": 0.0,
        "source_reclaim_reference_minus_candidate_area": 0.0,
        "source_reclaim_candidate_minus_reference_area": 0.0,
        "source_reclaim_signed_area_delta": 0.0,
        "source_reclaim_rejected_nonlocal_count": 0,
        "source_reclaim_rejected_non_improving_count": 0,
        "source_reclaim_operator_attempts": {},
        "source_reclaim_operator_applied": {},
        "polygon_simplify_candidate_count": 0,
        "polygon_simplify_applied_count": 0,
        "polygon_simplify_short_edge_count_before": 0,
        "polygon_simplify_short_edge_count_after": 0,
        "polygon_simplify_edit_zone_area": 0.0,
        "polygon_simplify_change_outside_edit_zone": 0.0,
        "polygon_simplify_reference_minus_candidate_area": 0.0,
        "polygon_simplify_candidate_minus_reference_area": 0.0,
        "polygon_simplify_signed_area_delta": 0.0,
        "polygon_simplify_area_balance_budget": 0.0,
        "polygon_simplify_rejected_nonlocal_count": 0,
        "polygon_simplify_rejected_area_imbalance_count": 0,
        "polygon_simplify_rejected_non_improving_count": 0,
        "polygon_simplify_operator_attempts": {},
        "polygon_simplify_operator_applied": {},
        "collect_stage_metrics": True,
        "stage_metrics": {},
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


def _derive_output_grid(options: ConditioningOptions, base_grid: float) -> float:
    return base_grid


def _derive_coverage_simplify_tolerance(options: ConditioningOptions) -> float:
    if options.min_feature_size <= 0:
        return 0.0
    return 2.0 * options.min_feature_size


def _derive_local_edit_radius(scale: float) -> float:
    if scale <= 0:
        return 0.0
    return 2.0 * scale


def _derive_patch_neighborhood_radius(scale: float) -> float:
    if scale <= 0:
        return 0.0
    return scale


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


def _clean_ring_coords(
    coords: Sequence[Sequence[float]],
    grid: float,
) -> list[tuple[float, float]] | None:
    ring = [(float(x), float(y)) for x, y, *_ in coords]
    if ring:
        if ring[0] != ring[-1]:
            ring.append(ring[0])
    if len(ring) < 4:
        return None

    unique_count = len(ring) - 1
    if unique_count < 3:
        return None

    closure_tolerance = max(4.0 * grid, 1e-9)
    tolerance = max(grid * grid * 0.25, 1e-12)
    needs_cleanup = False
    for index in range(unique_count):
        point = ring[index]
        next_point = ring[index + 1]
        if point == next_point and index + 1 < len(ring) - 1:
            needs_cleanup = True
            break
        prev_point = ring[index - 1] if index > 0 else ring[-2]
        next_unique = ring[index + 1] if index + 1 < unique_count else ring[0]
        cross = (
            (point[0] - prev_point[0]) * (next_unique[1] - point[1])
            - (point[1] - prev_point[1]) * (next_unique[0] - point[0])
        )
        if abs(cross) <= tolerance:
            needs_cleanup = True
            break
    if not needs_cleanup:
        closing_length = float(
            np.hypot(
                ring[0][0] - ring[-2][0],
                ring[0][1] - ring[-2][1],
            )
        )
        if closing_length > closure_tolerance:
            return ring

    unique: list[tuple[float, float]] = []
    for point in ring:
        if unique and point == unique[-1]:
            continue
        unique.append(point)

    if len(unique) >= 2 and unique[0] == unique[-1]:
        unique.pop()
    if len(unique) < 3:
        return None

    closure_tolerance = max(4.0 * grid, 1e-9)
    while len(unique) >= 4:
        closing_length = float(
            np.hypot(
                unique[0][0] - unique[-1][0],
                unique[0][1] - unique[-1][1],
            )
        )
        if closing_length > closure_tolerance:
            break

        drop_first_cost = _point_to_segment_distance(
            unique[0],
            unique[-1],
            unique[1],
        )
        drop_last_cost = _point_to_segment_distance(
            unique[-1],
            unique[-2],
            unique[0],
        )
        if drop_first_cost <= drop_last_cost:
            unique.pop(0)
        else:
            unique.pop()
        if len(unique) < 3:
            return None

    tolerance = max(grid * grid * 0.25, 1e-12)
    changed = True
    while changed and len(unique) >= 3:
        changed = False
        cleaned: list[tuple[float, float]] = []
        count = len(unique)
        for index in range(count):
            prev_point = unique[index - 1]
            point = unique[index]
            next_point = unique[(index + 1) % count]
            cross = (
                (point[0] - prev_point[0]) * (next_point[1] - point[1])
                - (point[1] - prev_point[1]) * (next_point[0] - point[0])
            )
            if abs(cross) <= tolerance:
                changed = True
                continue
            cleaned.append(point)
        unique = cleaned
        if len(unique) < 3:
            return None

    return unique + [unique[0]]


def _clean_polygon_vertices(poly: Polygon, grid: float) -> Polygon | None:
    shell = _clean_ring_coords(poly.exterior.coords, grid)
    if shell is None:
        return None

    holes: list[list[tuple[float, float]]] = []
    for ring in poly.interiors:
        cleaned = _clean_ring_coords(ring.coords, grid)
        if cleaned is None:
            continue
        holes.append(cleaned)

    cleaned_polygon = Polygon(shell, holes)
    if cleaned_polygon.is_empty or cleaned_polygon.area <= 0:
        return None
    return cleaned_polygon


def _point_to_segment_distance(
    point: tuple[float, float],
    start: tuple[float, float],
    end: tuple[float, float],
) -> float:
    dx = float(end[0] - start[0])
    dy = float(end[1] - start[1])
    base = float(np.hypot(dx, dy))
    if base <= 1e-12:
        return float(np.hypot(point[0] - start[0], point[1] - start[1]))
    return abs((point[0] - start[0]) * dy - (point[1] - start[1]) * dx) / base


def _find_short_collinear_vertex(
    points: Sequence[tuple[float, float]],
    *,
    target_scale: float,
    grid: float,
) -> int | None:
    count = len(points)
    if count < 3:
        return None

    line_tolerance = max(2.0 * grid, 0.05 * target_scale, 1e-9)
    for index in range(count):
        prev_point = points[index - 1]
        point = points[index]
        next_point = points[(index + 1) % count]
        prev_length = float(
            np.hypot(point[0] - prev_point[0], point[1] - prev_point[1])
        )
        next_length = float(
            np.hypot(next_point[0] - point[0], next_point[1] - point[1])
        )
        if min(prev_length, next_length) + 1e-12 >= target_scale:
            continue
        offset = _point_to_segment_distance(point, prev_point, next_point)
        if offset <= line_tolerance:
            return index
    return None


def _find_short_step_pair(
    points: Sequence[tuple[float, float]],
    *,
    target_scale: float,
    grid: float,
) -> tuple[int, int] | None:
    count = len(points)
    if count < 4:
        return None

    step_width_tolerance = max(4.0 * grid, 0.5 * target_scale, 1e-9)
    parallel_tolerance = 0.1
    for index in range(count):
        a = points[index - 1]
        b = points[index]
        c = points[(index + 1) % count]
        d = points[(index + 2) % count]

        ab = (float(b[0] - a[0]), float(b[1] - a[1]))
        bc = (float(c[0] - b[0]), float(c[1] - b[1]))
        cd = (float(d[0] - c[0]), float(d[1] - c[1]))
        ad = (float(d[0] - a[0]), float(d[1] - a[1]))

        len_ab = float(np.hypot(ab[0], ab[1]))
        len_bc = float(np.hypot(bc[0], bc[1]))
        len_cd = float(np.hypot(cd[0], cd[1]))
        len_ad = float(np.hypot(ad[0], ad[1]))
        if len_ab + 1e-12 >= target_scale or len_cd + 1e-12 >= target_scale:
            continue
        if len_bc <= max(grid, 1e-9) or len_ad <= max(grid, 1e-9):
            continue

        ab_cd_cross = abs(ab[0] * cd[1] - ab[1] * cd[0]) / max(len_ab * len_cd, 1e-12)
        ab_cd_dot = (ab[0] * cd[0] + ab[1] * cd[1]) / max(len_ab * len_cd, 1e-12)
        ad_bc_cross = abs(ad[0] * bc[1] - ad[1] * bc[0]) / max(len_ad * len_bc, 1e-12)
        if ab_cd_cross > parallel_tolerance or ab_cd_dot > -0.5:
            continue
        if ad_bc_cross > parallel_tolerance:
            continue

        step_width = max(
            _point_to_segment_distance(b, a, d),
            _point_to_segment_distance(c, a, d),
        )
        if step_width > step_width_tolerance:
            continue
        return index, (index + 1) % count
    return None


def _clean_ring_short_edge_chains(
    coords: Sequence[Sequence[float]],
    *,
    target_scale: float,
    grid: float,
) -> tuple[list[tuple[float, float]] | None, int]:
    unique: list[tuple[float, float]] = []
    for x, y, *_ in coords:
        point = (float(x), float(y))
        if unique and point == unique[-1]:
            continue
        unique.append(point)

    if len(unique) >= 2 and unique[0] == unique[-1]:
        unique.pop()
    if len(unique) < 3:
        return None, 0

    removed_count = 0
    while len(unique) >= 3:
        pair = _find_short_step_pair(
            unique,
            target_scale=target_scale,
            grid=grid,
        )
        if pair is not None:
            remove = {pair[0] % len(unique), pair[1] % len(unique)}
            unique = [point for index, point in enumerate(unique) if index not in remove]
            removed_count += 2
            continue

        vertex = _find_short_collinear_vertex(
            unique,
            target_scale=target_scale,
            grid=grid,
        )
        if vertex is None:
            break
        unique.pop(vertex)
        removed_count += 1

    if len(unique) < 3:
        return None, removed_count
    return unique + [unique[0]], removed_count


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
        cleaned = _clean_polygon_vertices(polygon, grid)
        if cleaned is None:
            diagnostics["collapsed_count"] += 1
            continue
        canonical.append(orient(cleaned, sign=1.0))

    if not canonical and not geom.is_empty:
        diagnostics["collapsed_count"] += 1
    return canonical


def _canonicalize_for_output(
    geom: BaseGeometry,
    *,
    output_grid: float,
    min_area: float,
    diagnostics: dict[str, Any],
) -> list[Polygon]:
    reference_parts: list[Polygon] = []
    reference_dropped = 0
    for polygon in _extract_polygon_parts(geom):
        if polygon.area + 1e-12 < min_area:
            reference_dropped += 1
            continue
        reference_parts.append(orient(polygon, sign=1.0))
    diagnostics["output_canonicalization_candidate_count"] += 1

    candidate_parts: list[Polygon] = []
    candidate_dropped = 0
    for polygon in _canonicalize(geom, output_grid, diagnostics):
        if polygon.area + 1e-12 < min_area:
            candidate_dropped += 1
            continue
        candidate_parts.append(orient(polygon, sign=1.0))
    if not candidate_parts:
        diagnostics["output_canonicalization_fallback_count"] += 1
        diagnostics["dropped_small_count"] += reference_dropped
        return reference_parts

    difference_metrics = _difference_area_metrics(
        _as_geometry(reference_parts),
        _as_geometry(candidate_parts),
    )
    balance_budget = _area_balance_budget(
        difference_metrics["symmetric_difference_area"],
        grid=output_grid,
    )
    diagnostics["output_canonicalization_reference_minus_candidate_area"] += (
        difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics["output_canonicalization_candidate_minus_reference_area"] += (
        difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics["output_canonicalization_signed_area_delta"] += difference_metrics[
        "union_area_delta"
    ]
    diagnostics["output_canonicalization_area_balance_budget"] += balance_budget

    if abs(difference_metrics["union_area_delta"]) > balance_budget:
        diagnostics["output_canonicalization_fallback_count"] += 1
        diagnostics["dropped_small_count"] += reference_dropped
        return reference_parts

    diagnostics["output_canonicalization_output_grid_applied_count"] += 1
    diagnostics["dropped_small_count"] += candidate_dropped
    return candidate_parts


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


def _filter_small_output_polygons(
    polygons: Sequence[Polygon],
    source_map: Sequence[Sequence[int]],
    *,
    min_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    if min_area <= 0:
        return list(polygons), [list(indices) for indices in source_map]

    kept_polygons: list[Polygon] = []
    kept_sources: list[list[int]] = []
    removed_count = 0
    removed_area = 0.0

    for polygon, indices in zip(polygons, source_map):
        if polygon.area + 1e-12 >= min_area:
            kept_polygons.append(polygon)
            kept_sources.append(list(indices))
            continue
        removed_count += 1
        removed_area += float(polygon.area)

    diagnostics["final_min_area_filter_removed_count"] += removed_count
    diagnostics["final_min_area_filter_removed_area"] += removed_area
    diagnostics["dropped_small_count"] += removed_count
    return kept_polygons, kept_sources


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
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
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
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            -radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
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


def _assign_component_sources(
    components: Sequence[Polygon],
    member_polygons: Sequence[Polygon],
    member_sources: Sequence[Sequence[int]],
    *,
    grid: float,
    diagnostics: dict[str, Any],
) -> list[list[int]]:
    if not components:
        return []

    group_sources = sorted(
        {
            source
            for indices in member_sources
            for source in indices
        }
    )
    if len(components) == 1:
        return [group_sources]

    diagnostics["multi_component_group_count"] += 1
    area_tolerance = max(grid * grid, 1e-9)
    distance_tolerance = grid
    component_sources: list[list[int]] = []

    for component in components:
        assigned: list[int] = []
        for member_polygon, indices in zip(member_polygons, member_sources):
            try:
                overlap_area = component.intersection(member_polygon).area
            except GEOSException as exc:
                _record_geos_exception(diagnostics, "component_source_overlap", exc)
                continue
            if overlap_area > area_tolerance:
                assigned.extend(indices)

        if not assigned:
            representative = component.representative_point()
            for member_polygon, indices in zip(member_polygons, member_sources):
                try:
                    if member_polygon.covers(representative):
                        assigned.extend(indices)
                        continue
                    if component.distance(member_polygon) <= distance_tolerance:
                        assigned.extend(indices)
                except GEOSException as exc:
                    _record_geos_exception(diagnostics, "component_source_distance", exc)

        cleaned = sorted(set(assigned)) if assigned else group_sources
        diagnostics["component_source_reassignment_count"] += 1
        component_sources.append(cleaned)

    return component_sources


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


def _should_reconstruct_global_coverage(
    polygons: Sequence[Polygon],
    *,
    grid: float,
    diagnostics: dict[str, Any],
) -> bool:
    overlap_tolerance = max(grid * grid, 1e-9)
    diagnostics["global_reconstruction_overlap_threshold"] = overlap_tolerance

    overlap_area = _coverage_overlap_area(polygons)
    if overlap_area > overlap_tolerance:
        diagnostics["global_reconstruction_applied"] = True
        diagnostics["global_reconstruction_reason"] = "overlap_detected"
        return True

    diagnostics["global_reconstruction_applied"] = False
    diagnostics["global_reconstruction_reason"] = "coverage_already_disjoint"
    return False


def _expanded_query_geometry(
    geometry: BaseGeometry,
    *,
    radius: float,
) -> BaseGeometry:
    if radius <= 0:
        return geometry.envelope
    minx, miny, maxx, maxy = geometry.bounds
    return shapely.box(minx - radius, miny - radius, maxx + radius, maxy + radius)


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


def _count_ring_boundary_contacts(
    intersection: BaseGeometry,
    *,
    tolerance: float,
) -> int:
    if intersection.is_empty:
        return 0

    try:
        if intersection.length > tolerance:
            return 1
    except (AttributeError, TypeError):
        return 1

    if intersection.geom_type == "Point":
        return 1
    if intersection.geom_type == "MultiPoint":
        return sum(1 for _ in shapely.get_parts(intersection))
    if intersection.geom_type == "GeometryCollection":
        return sum(
            _count_ring_boundary_contacts(part, tolerance=tolerance)
            for part in shapely.get_parts(intersection)
        )
    return 1


def _polygon_ring_boundary_contact_count(
    polygon: Polygon,
    *,
    tolerance: float = 1e-12,
) -> int:
    rings = [polygon.exterior, *polygon.interiors]
    if len(rings) < 2:
        return 0

    count = 0
    for ring_index, ring in enumerate(rings):
        for other in rings[ring_index + 1 :]:
            try:
                boundary_intersection = ring.intersection(other)
            except GEOSException:
                return count + 1
            count += _count_ring_boundary_contacts(
                boundary_intersection,
                tolerance=tolerance,
            )
    return count


def _polygon_has_ring_boundary_contacts(
    polygon: Polygon,
    *,
    tolerance: float = 1e-12,
) -> bool:
    return _polygon_ring_boundary_contact_count(
        polygon,
        tolerance=tolerance,
    ) > 0


def _normalize_mesher_ready_polygon(
    polygon: Polygon,
    *,
    declared_scale: float,
    diagnostics: dict[str, Any] | None = None,
) -> list[Polygon]:
    polygon = orient(polygon, sign=1.0)
    if polygon.is_empty or not polygon.interiors:
        return [polygon]

    if not _polygon_has_ring_boundary_contacts(polygon):
        return [polygon]

    shell = Polygon(np.asarray(polygon.exterior.coords, dtype=np.float64))
    hole_polygons = [
        Polygon(np.asarray(ring.coords, dtype=np.float64))
        for ring in polygon.interiors
    ]
    max_extent = max(
        polygon.bounds[2] - polygon.bounds[0],
        polygon.bounds[3] - polygon.bounds[1],
        1.0,
    )
    distance = max(
        float(declared_scale) * 5e-5,
        float(max_extent) * 5e-7,
        1e-9,
    )

    for _attempt in range(8):
        expanded_holes = unary_union(
            [
                hole.buffer(
                    distance,
                    quad_segs=1,
                    join_style=BufferJoinStyle.mitre,
                )
                for hole in hole_polygons
            ]
        )
        candidate_geometry = shell.difference(expanded_holes)
        candidate_polygons = [
            orient(candidate, sign=1.0)
            for candidate in _extract_polygon_parts(candidate_geometry)
        ]
        if candidate_polygons and all(
            not _polygon_has_ring_boundary_contacts(
                candidate,
                tolerance=max(distance * 0.25, 1e-12),
            )
            for candidate in candidate_polygons
        ):
            if diagnostics is not None:
                diagnostics["mesher_regularized_polygon_count"] = (
                    diagnostics.get("mesher_regularized_polygon_count", 0) + 1
                )
                diagnostics["mesher_regularized_component_count"] = (
                    diagnostics.get("mesher_regularized_component_count", 0)
                    + len(candidate_polygons)
                )
                diagnostics["mesher_regularization_area_delta_total"] = (
                    diagnostics.get("mesher_regularization_area_delta_total", 0.0)
                    + float(sum(part.area for part in candidate_polygons) - polygon.area)
                )
                diagnostics["mesher_regularization_max_distance"] = max(
                    float(diagnostics.get("mesher_regularization_max_distance", 0.0)),
                    float(distance),
                )
            return candidate_polygons
        distance *= 2.0

    if diagnostics is not None:
        diagnostics["mesher_regularization_failed_count"] = (
            diagnostics.get("mesher_regularization_failed_count", 0) + 1
        )
    return [polygon]


def _regularize_ring_contact_polygon(
    polygon: Polygon,
    *,
    min_distance: float,
    target_clearance: float,
    grid: float,
) -> list[Polygon]:
    polygon = orient(polygon, sign=1.0)
    if polygon.is_empty or not polygon.interiors:
        return [polygon]

    if not _polygon_has_ring_boundary_contacts(polygon):
        return [polygon]

    shell = Polygon(np.asarray(polygon.exterior.coords, dtype=np.float64))
    hole_polygons = [
        Polygon(np.asarray(ring.coords, dtype=np.float64))
        for ring in polygon.interiors
    ]
    distance = max(float(min_distance), 1e-9)

    for _attempt in range(8):
        expanded_holes = unary_union(
            [
                hole.buffer(
                    distance,
                    quad_segs=1,
                    join_style=BufferJoinStyle.mitre,
                )
                for hole in hole_polygons
            ]
        )
        candidate_geometry = shell.difference(expanded_holes)
        candidate_polygons = [
            orient(candidate, sign=1.0)
            for candidate in _extract_polygon_parts(candidate_geometry)
        ]
        if not candidate_polygons:
            distance *= 2.0
            continue

        if any(_polygon_has_ring_boundary_contacts(candidate) for candidate in candidate_polygons):
            distance *= 2.0
            continue

        if all(
            _signature_satisfies_scale_contract(
                _polygon_defect_signature(
                    candidate,
                    target_scale=target_clearance,
                ),
                target_scale=target_clearance,
                grid=grid,
            )
            for candidate in candidate_polygons
        ):
            return candidate_polygons

        distance *= 2.0

    return [polygon]


def _coverage_pairwise_metrics(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
) -> dict[str, float | int | None]:
    if target_scale <= 0 or len(polygons) < 2:
        return {
            "pair_issue_count": 0,
            "point_touch_count": 0,
            "close_pair_count": 0,
            "min_pair_clearance": None,
        }

    point_tolerance = max(target_scale * 1e-6, 1e-9)
    line_tolerance = point_tolerance
    tree = STRtree(polygons)
    point_touch_count = 0
    close_pair_count = 0
    min_pair_clearance: float | None = None

    for index, polygon in enumerate(polygons):
        for candidate in tree.query(
            _expanded_query_geometry(polygon, radius=target_scale)
        ):
            other_index = int(candidate)
            if other_index <= index:
                continue
            other = polygons[other_index]
            try:
                distance = float(polygon.distance(other))
            except GEOSException:
                continue
            if distance > target_scale + point_tolerance:
                continue

            if distance <= point_tolerance:
                try:
                    boundary_intersection = polygon.boundary.intersection(other.boundary)
                except GEOSException:
                    boundary_intersection = GeometryCollection()
                if boundary_intersection.length > line_tolerance:
                    continue
                point_touch_count += 1
                min_pair_clearance = 0.0
                continue

            close_pair_count += 1
            if min_pair_clearance is None or distance < min_pair_clearance:
                min_pair_clearance = distance

    return {
        "pair_issue_count": point_touch_count + close_pair_count,
        "point_touch_count": point_touch_count,
        "close_pair_count": close_pair_count,
        "min_pair_clearance": min_pair_clearance,
    }


def _segment_length_stats(
    polygons: Sequence[Polygon],
    *,
    short_edge_threshold: float = 0.0,
) -> dict[str, float | int | None]:
    edge_count = 0
    vertex_count = 0
    short_edge_count = 0
    min_edge_length: float | None = None
    total_edge_length = 0.0

    for polygon in polygons:
        rings = [polygon.exterior, *polygon.interiors]
        for ring in rings:
            coords = np.asarray(ring.coords, dtype=float)
            if len(coords) < 2:
                continue
            vertex_count += max(len(coords) - 1, 0)
            deltas = np.diff(coords[:, :2], axis=0)
            lengths = np.hypot(deltas[:, 0], deltas[:, 1])
            if lengths.size == 0:
                continue
            edge_count += int(lengths.size)
            total_edge_length += float(lengths.sum())
            ring_min = float(lengths.min())
            min_edge_length = (
                ring_min if min_edge_length is None else min(min_edge_length, ring_min)
            )
            if short_edge_threshold > 0:
                short_edge_count += int(
                    np.count_nonzero(lengths + 1e-12 < short_edge_threshold)
                )

    return {
        "edge_count": edge_count,
        "vertex_count": vertex_count,
        "short_edge_count": short_edge_count,
        "min_edge_length": min_edge_length,
        "mean_edge_length": (total_edge_length / edge_count) if edge_count else None,
    }


def _short_edge_segments(
    polygons: Sequence[Polygon],
    *,
    short_edge_threshold: float,
) -> list[LineString]:
    if short_edge_threshold <= 0:
        return []

    segments: list[LineString] = []
    for polygon in polygons:
        rings = [polygon.exterior, *polygon.interiors]
        for ring in rings:
            coords = list(ring.coords)
            for start, end in zip(coords, coords[1:]):
                dx = float(end[0] - start[0])
                dy = float(end[1] - start[1])
                length = float(np.hypot(dx, dy))
                if length + 1e-12 >= short_edge_threshold:
                    continue
                segments.append(LineString([start, end]))
    return segments


def _short_edge_edit_zone(
    polygons: Sequence[Polygon],
    *,
    short_edge_threshold: float,
    radius: float,
) -> BaseGeometry:
    if short_edge_threshold <= 0 or radius <= 0:
        return GeometryCollection()

    segments = _short_edge_segments(
        polygons,
        short_edge_threshold=short_edge_threshold,
    )
    if not segments:
        return GeometryCollection()

    return unary_union(segments).buffer(
        radius,
        quad_segs=1,
        join_style=BufferJoinStyle.mitre,
        mitre_limit=1000.0,
    )


def _coverage_pair_issue_edit_zone(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    radius: float,
) -> BaseGeometry:
    if target_scale <= 0 or radius <= 0 or len(polygons) < 2:
        return GeometryCollection()

    point_tolerance = max(target_scale * 1e-6, 1e-9)
    line_tolerance = point_tolerance
    tree = STRtree(polygons)
    issue_geometries: list[BaseGeometry] = []

    for index, polygon in enumerate(polygons):
        for candidate in tree.query(
            _expanded_query_geometry(polygon, radius=target_scale)
        ):
            other_index = int(candidate)
            if other_index <= index:
                continue
            other = polygons[other_index]
            try:
                distance = float(polygon.distance(other))
            except GEOSException:
                continue
            if distance > target_scale + point_tolerance:
                continue

            if distance <= point_tolerance:
                try:
                    boundary_intersection = polygon.boundary.intersection(other.boundary)
                except GEOSException:
                    boundary_intersection = GeometryCollection()
                if boundary_intersection.length > line_tolerance:
                    continue
                issue_geometries.append(
                    boundary_intersection
                    if not boundary_intersection.is_empty
                    else polygon.representative_point()
                )
                continue

            try:
                issue_geometries.append(shapely.shortest_line(polygon, other))
            except GEOSException:
                continue

    if not issue_geometries:
        return GeometryCollection()
    return unary_union(issue_geometries).buffer(
        radius,
        quad_segs=1,
        join_style=BufferJoinStyle.mitre,
        mitre_limit=1000.0,
    )


def _coverage_edit_zone(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    radius: float,
) -> BaseGeometry:
    zones: list[BaseGeometry] = []
    short_edge_zone = _short_edge_edit_zone(
        polygons,
        short_edge_threshold=target_scale,
        radius=radius,
    )
    if not short_edge_zone.is_empty:
        zones.append(short_edge_zone)
    pair_issue_zone = _coverage_pair_issue_edit_zone(
        polygons,
        target_scale=target_scale,
        radius=radius,
    )
    if not pair_issue_zone.is_empty:
        zones.append(pair_issue_zone)
    if not zones:
        return GeometryCollection()
    return unary_union(zones)


def _change_outside_edit_zone(
    reference: BaseGeometry,
    candidate: BaseGeometry,
    *,
    edit_zone: BaseGeometry,
) -> float:
    change = reference.symmetric_difference(candidate)
    if edit_zone.is_empty:
        return float(change.area)
    return float(change.difference(edit_zone).area)


def _difference_area_metrics(
    reference: BaseGeometry,
    candidate: BaseGeometry,
) -> dict[str, float]:
    symmetric_difference_area = float(reference.symmetric_difference(candidate).area)
    union_area_delta = float(candidate.area - reference.area)
    missing_area = max(
        0.0,
        0.5 * (symmetric_difference_area - union_area_delta),
    )
    extra_area = max(
        0.0,
        0.5 * (symmetric_difference_area + union_area_delta),
    )
    return {
        "reference_minus_candidate_area": missing_area,
        "candidate_minus_reference_area": extra_area,
        "symmetric_difference_area": symmetric_difference_area,
        "union_area_delta": union_area_delta,
    }


def _area_balance_budget(
    symmetric_difference_area: float,
    *,
    grid: float,
    short_edge_count: int = 0,
    short_edge_threshold: float = 0.0,
) -> float:
    return max(
        _AREA_BALANCE_RELATIVE_TOLERANCE * symmetric_difference_area,
        (
            _AREA_BALANCE_SHORT_EDGE_MULTIPLIER
            * short_edge_count
            * short_edge_threshold
            * grid
        ),
        _AREA_BALANCE_ABSOLUTE_GRID_MULTIPLIER * grid * grid,
        1e-9,
    )


def _operator_count_increment(
    diagnostics: dict[str, Any],
    key: str,
    operator: str,
) -> None:
    counts = diagnostics.setdefault(key, {})
    counts[operator] = counts.get(operator, 0) + 1


def _safe_minimum_clearance(geometry: BaseGeometry) -> float | None:
    try:
        clearance = float(shapely.minimum_clearance(geometry))
    except GEOSException:
        return None
    if not np.isfinite(clearance):
        return None
    return clearance


def _polygon_defect_signature(
    polygon: Polygon,
    *,
    target_scale: float,
) -> _PolygonDefectSignature:
    segment_stats = _segment_length_stats(
        [polygon],
        short_edge_threshold=target_scale,
    )
    clearance = _safe_minimum_clearance(polygon)
    clearance_deficit = 0.0
    if target_scale > 0 and clearance is not None:
        clearance_deficit = max(target_scale - clearance, 0.0)
    return _PolygonDefectSignature(
        clearance=clearance,
        clearance_deficit=clearance_deficit,
        short_edge_count=int(segment_stats["short_edge_count"]),
        min_edge_length=segment_stats["min_edge_length"],
        vertex_count=int(segment_stats["vertex_count"]),
    )


def _coverage_defect_signature(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
) -> _CoverageDefectSignature:
    segment_stats = _segment_length_stats(
        polygons,
        short_edge_threshold=target_scale,
    )
    pair_metrics = _coverage_pairwise_metrics(
        polygons,
        target_scale=target_scale,
    )
    polygon_clearance = _minimum_clearance(polygons)
    pair_clearance = pair_metrics["min_pair_clearance"]
    clearances = [
        value
        for value in (polygon_clearance, pair_clearance)
        if value is not None
    ]
    return _CoverageDefectSignature(
        min_clearance=min(clearances) if clearances else None,
        ring_contact_count=sum(
            _polygon_ring_boundary_contact_count(polygon)
            for polygon in polygons
        ),
        pair_issue_count=int(pair_metrics["pair_issue_count"]),
        point_touch_count=int(pair_metrics["point_touch_count"]),
        close_pair_count=int(pair_metrics["close_pair_count"]),
        min_pair_clearance=(
            float(pair_clearance) if pair_clearance is not None else None
        ),
        short_edge_count=int(segment_stats["short_edge_count"]),
        min_edge_length=segment_stats["min_edge_length"],
        vertex_count=int(segment_stats["vertex_count"]),
    )


def _signature_improves(
    reference: _PolygonDefectSignature,
    candidate: _PolygonDefectSignature,
    *,
    grid: float,
) -> bool:
    tolerance = max(grid, 1e-9)
    if (
        reference.clearance is not None
        and candidate.clearance is not None
        and candidate.clearance + tolerance < reference.clearance
    ):
        return False

    if candidate.clearance_deficit + tolerance < reference.clearance_deficit:
        return True
    if candidate.clearance_deficit > reference.clearance_deficit + tolerance:
        return False

    if candidate.short_edge_count < reference.short_edge_count:
        return True
    if candidate.short_edge_count > reference.short_edge_count:
        return False

    reference_min_edge = reference.min_edge_length or 0.0
    candidate_min_edge = candidate.min_edge_length or 0.0
    if candidate_min_edge > reference_min_edge + tolerance:
        return True

    if candidate.vertex_count < reference.vertex_count:
        if (
            reference.clearance is None
            or candidate.clearance is None
            or candidate.clearance + tolerance >= reference.clearance
        ):
            return True

    return False


def _coverage_signature_improves(
    reference: _CoverageDefectSignature,
    candidate: _CoverageDefectSignature,
    *,
    grid: float,
    target_scale: float,
) -> bool:
    tolerance = max(grid, 1e-9)
    reference_clearance_deficit = max(
        target_scale - (reference.min_clearance or 0.0),
        0.0,
    )
    candidate_clearance_deficit = max(
        target_scale - (candidate.min_clearance or 0.0),
        0.0,
    )
    if candidate_clearance_deficit + tolerance < reference_clearance_deficit:
        return True
    if candidate_clearance_deficit > reference_clearance_deficit + tolerance:
        return False
    if candidate.ring_contact_count < reference.ring_contact_count:
        return True
    if candidate.ring_contact_count > reference.ring_contact_count:
        return False
    if candidate.pair_issue_count < reference.pair_issue_count:
        return True
    if candidate.pair_issue_count > reference.pair_issue_count:
        return False
    if (
        reference.min_clearance is not None
        and candidate.min_clearance is not None
        and candidate.min_clearance + tolerance < reference.min_clearance
    ):
        return False

    if candidate.short_edge_count < reference.short_edge_count:
        return True
    if candidate.short_edge_count > reference.short_edge_count:
        return False

    reference_min_edge = reference.min_edge_length or 0.0
    candidate_min_edge = candidate.min_edge_length or 0.0
    if candidate_min_edge > reference_min_edge + tolerance:
        return True

    if candidate.vertex_count < reference.vertex_count:
        if (
            reference.min_clearance is None
            or candidate.min_clearance is None
            or candidate.min_clearance + tolerance >= reference.min_clearance
        ):
            return True

    return False


def _signature_not_worse(
    reference: _PolygonDefectSignature,
    candidate: _PolygonDefectSignature,
    *,
    grid: float,
) -> bool:
    tolerance = max(grid, 1e-9)
    if (
        reference.clearance is not None
        and candidate.clearance is not None
        and candidate.clearance + tolerance < reference.clearance
    ):
        return False
    if candidate.clearance_deficit > reference.clearance_deficit + tolerance:
        return False
    if candidate.short_edge_count > reference.short_edge_count:
        return False
    reference_min_edge = reference.min_edge_length or 0.0
    candidate_min_edge = candidate.min_edge_length or 0.0
    if candidate_min_edge + tolerance < reference_min_edge:
        return False
    return True


def _signature_satisfies_scale_contract(
    signature: _PolygonDefectSignature,
    *,
    target_scale: float,
    grid: float,
) -> bool:
    if target_scale <= 0:
        return True

    tolerance = max(grid, 1e-9)
    if signature.short_edge_count > 0:
        return False
    if signature.clearance_deficit > tolerance:
        return False

    min_edge_length = signature.min_edge_length or 0.0
    if min_edge_length + tolerance < target_scale:
        return False
    return True


def _coverage_signature_not_worse(
    reference: _CoverageDefectSignature,
    candidate: _CoverageDefectSignature,
    *,
    grid: float,
    target_scale: float,
) -> bool:
    tolerance = max(grid, 1e-9)
    reference_clearance_deficit = max(
        target_scale - (reference.min_clearance or 0.0),
        0.0,
    )
    candidate_clearance_deficit = max(
        target_scale - (candidate.min_clearance or 0.0),
        0.0,
    )
    if candidate_clearance_deficit > reference_clearance_deficit + tolerance:
        return False
    if candidate.ring_contact_count > reference.ring_contact_count:
        return False
    if candidate.pair_issue_count > reference.pair_issue_count:
        return False
    if (
        reference.min_clearance is not None
        and candidate.min_clearance is not None
        and candidate.min_clearance + tolerance < reference.min_clearance
    ):
        return False
    if candidate.short_edge_count > reference.short_edge_count:
        return False
    reference_min_edge = reference.min_edge_length or 0.0
    candidate_min_edge = candidate.min_edge_length or 0.0
    if candidate_min_edge + tolerance < reference_min_edge:
        return False
    return True


def _candidate_metric_key(prefix: str, suffix: str) -> str:
    return f"{prefix}_{suffix}"


def _coverage_signature_score(
    signature: _CoverageDefectSignature,
    *,
    target_scale: float,
) -> tuple[float, int, int, int, float, int]:
    return (
        max(target_scale - (signature.min_clearance or 0.0), 0.0),
        signature.ring_contact_count,
        signature.pair_issue_count,
        signature.short_edge_count,
        -(signature.min_edge_length or 0.0),
        signature.vertex_count,
    )


def _coverage_signature_satisfies_scale_contract(
    signature: _CoverageDefectSignature,
    *,
    target_scale: float,
    grid: float,
) -> bool:
    if target_scale <= 0:
        return True

    tolerance = max(grid, 1e-9)
    if signature.ring_contact_count > 0:
        return False
    if signature.pair_issue_count > 0:
        return False
    if signature.short_edge_count > 0:
        return False
    if max(target_scale - (signature.min_clearance or 0.0), 0.0) > tolerance:
        return False

    min_edge_length = signature.min_edge_length or 0.0
    if min_edge_length + tolerance < target_scale:
        return False
    return True


def _coverage_candidate_small_output_metrics(
    candidate: _CoverageSimplifyCandidate,
    *,
    output_min_area: float,
) -> tuple[int, float]:
    if output_min_area <= 0:
        return 0, 0.0

    small_count = 0
    small_area = 0.0
    for polygon in candidate.polygons:
        if polygon.area + 1e-9 < output_min_area:
            small_count += 1
            small_area += float(polygon.area)
    return small_count, small_area


def _coverage_candidate_prescore(
    candidate: _CoverageSimplifyCandidate,
    *,
    target_scale: float,
    output_min_area: float,
) -> tuple[float, ...]:
    small_output_count, small_output_area = _coverage_candidate_small_output_metrics(
        candidate,
        output_min_area=output_min_area,
    )
    clearance_deficit = max(target_scale - (candidate.signature.min_clearance or 0.0), 0.0)
    return (
        float(candidate.signature.ring_contact_count),
        float(candidate.signature.pair_issue_count),
        float(candidate.signature.short_edge_count),
        clearance_deficit,
        -(candidate.signature.min_edge_length or 0.0),
        float(candidate.signature.vertex_count),
        float(small_output_count),
        small_output_area,
        candidate.difference_metrics["reference_minus_candidate_area"],
        candidate.difference_metrics["candidate_minus_reference_area"],
        candidate.difference_metrics["symmetric_difference_area"],
        abs(candidate.difference_metrics["union_area_delta"]),
    )


def _coverage_candidate_has_small_fidelity_drift(
    candidate: _CoverageSimplifyCandidate,
    *,
    target_scale: float,
    grid: float,
) -> bool:
    return _difference_metrics_have_small_fidelity_drift(
        candidate.difference_metrics,
        edit_zone_area=candidate.edit_zone_area,
        target_scale=target_scale,
        grid=grid,
    )


def _difference_metrics_have_small_fidelity_drift(
    difference_metrics: dict[str, float],
    *,
    edit_zone_area: float,
    target_scale: float,
    grid: float,
) -> bool:
    area_unit = max(target_scale * target_scale, 16.0 * grid * grid, 1e-9)
    reference_area = max(edit_zone_area, area_unit)
    return (
        difference_metrics["reference_minus_candidate_area"]
        <= max(reference_area * 0.025, area_unit)
        and difference_metrics["candidate_minus_reference_area"]
        <= max(reference_area * 0.025, area_unit)
        and difference_metrics["symmetric_difference_area"]
        <= max(reference_area * 0.05, 2.0 * area_unit)
        and abs(difference_metrics["union_area_delta"])
        <= max(reference_area * 0.0125, 0.5 * area_unit)
    )


def _should_attempt_local_coverage_candidate(
    reference_signature: _CoverageDefectSignature,
    global_candidate: _CoverageSimplifyCandidate | None,
    *,
    target_scale: float,
    grid: float,
    purpose: Literal["coverage", "meshing"],
) -> bool:
    if target_scale <= 0:
        return False
    if global_candidate is None:
        return (
            reference_signature.pair_issue_count > 0
            or reference_signature.short_edge_count > 0
        )

    residual_short_edge_budget = max(12, reference_signature.short_edge_count // 12)
    if global_candidate.signature.pair_issue_count > reference_signature.pair_issue_count:
        return True
    if not _coverage_signature_satisfies_scale_contract(
        global_candidate.signature,
        target_scale=target_scale,
        grid=grid,
    ):
        if (
            purpose == "coverage"
            and global_candidate.signature.short_edge_count <= residual_short_edge_budget
        ):
            return False
        if purpose == "meshing":
            # The final meshing output must satisfy the declared scale contract.
            # Even a small residual short-edge cluster is worth a local cleanup
            # attempt here because those few defects can still create very hard
            # 2D/3D meshing failures downstream.
            return True
        return True
    if purpose == "meshing":
        return False
    return not _coverage_candidate_has_small_fidelity_drift(
        global_candidate,
        target_scale=target_scale,
        grid=grid,
    )


def _should_attempt_meshing_local_candidate(
    reference_signature: _CoverageDefectSignature,
    global_candidate: _CoverageSimplifyCandidate | None,
    *,
    target_scale: float,
    grid: float,
) -> bool:
    if not _should_attempt_local_coverage_candidate(
        reference_signature,
        global_candidate,
        target_scale=target_scale,
        grid=grid,
        purpose="meshing",
    ):
        return False

    # This used to skip expensive local search for mixed residual defects that
    # often collapsed back to identity. For the meshing stage we still need to
    # try, because a few surviving residual defects can violate the declared
    # min_feature_size contract and later break 2D/3D meshing.
    if (
        reference_signature.point_touch_count > 0
        and reference_signature.close_pair_count > 0
        and reference_signature.short_edge_count > 0
        and reference_signature.short_edge_count <= 8
    ):
        if global_candidate is not None and _coverage_signature_satisfies_scale_contract(
            global_candidate.signature,
            target_scale=target_scale,
            grid=grid,
        ):
            if (
                global_candidate.signature.pair_issue_count
                >= reference_signature.pair_issue_count
                and global_candidate.signature.point_touch_count
                >= reference_signature.point_touch_count
            ):
                return False

    return True


def _should_evaluate_additional_post_coverage_candidate(
    primary: _CoverageSimplifyCandidate,
    secondary: _CoverageSimplifyCandidate,
    *,
    target_scale: float,
    grid: float,
    output_min_area: float,
) -> bool:
    area_tolerance = max(target_scale * target_scale, 16.0 * grid * grid, 1e-9)
    primary_small_count, primary_small_area = _coverage_candidate_small_output_metrics(
        primary,
        output_min_area=output_min_area,
    )
    secondary_small_count, secondary_small_area = (
        _coverage_candidate_small_output_metrics(
            secondary,
            output_min_area=output_min_area,
        )
    )
    if secondary.signature.pair_issue_count < primary.signature.pair_issue_count:
        return True
    if secondary.signature.short_edge_count + 1 < primary.signature.short_edge_count:
        return True
    if (
        secondary.signature.pair_issue_count == primary.signature.pair_issue_count
        and secondary.signature.short_edge_count < primary.signature.short_edge_count
        and secondary.difference_metrics["reference_minus_candidate_area"]
        <= primary.difference_metrics["reference_minus_candidate_area"] + area_tolerance
    ):
        return True
    if (
        secondary_small_count < primary_small_count
        and secondary_small_area + area_tolerance < primary_small_area
    ):
        return True
    if (
        secondary.difference_metrics["reference_minus_candidate_area"] + area_tolerance
        < primary.difference_metrics["reference_minus_candidate_area"]
        and secondary.signature.pair_issue_count <= primary.signature.pair_issue_count
        and secondary.signature.short_edge_count <= primary.signature.short_edge_count
    ):
        return True
    if (
        secondary.difference_metrics["symmetric_difference_area"]
        + 2.0 * area_tolerance
        < primary.difference_metrics["symmetric_difference_area"]
        and secondary.signature.pair_issue_count <= primary.signature.pair_issue_count
        and secondary.signature.short_edge_count <= primary.signature.short_edge_count
    ):
        return True
    return False


def _should_finalize_identity_post_coverage_candidate(
    identity_candidate: _CoverageSimplifyCandidate,
    best_candidate: _CoverageSimplifyCandidate | None,
    *,
    target_scale: float,
    grid: float,
    output_min_area: float,
) -> bool:
    if best_candidate is None:
        return True

    identity_prescore = _coverage_candidate_prescore(
        identity_candidate,
        target_scale=target_scale,
        output_min_area=output_min_area,
    )
    best_prescore = _coverage_candidate_prescore(
        best_candidate,
        target_scale=target_scale,
        output_min_area=output_min_area,
    )
    if best_prescore >= identity_prescore:
        return True

    area_tolerance = max(target_scale * target_scale, 16.0 * grid * grid, 1e-9)
    identity_small_count, identity_small_area = _coverage_candidate_small_output_metrics(
        identity_candidate,
        output_min_area=output_min_area,
    )
    best_small_count, best_small_area = _coverage_candidate_small_output_metrics(
        best_candidate,
        output_min_area=output_min_area,
    )
    if best_small_count > identity_small_count:
        return True
    if best_small_area > identity_small_area + area_tolerance:
        return True

    if _coverage_candidate_has_small_fidelity_drift(
        best_candidate,
        target_scale=target_scale,
        grid=grid,
    ):
        strong_pair_win = (
            identity_candidate.signature.pair_issue_count > 0
            and best_candidate.signature.pair_issue_count == 0
        )
        strong_short_edge_win = (
            identity_candidate.signature.short_edge_count >= 8
            and best_candidate.signature.short_edge_count
            <= max(0, identity_candidate.signature.short_edge_count // 4)
        )
        return not (strong_pair_win or strong_short_edge_win)

    identity_clearance_deficit = max(
        target_scale - (identity_candidate.signature.min_clearance or 0.0),
        0.0,
    )
    best_clearance_deficit = max(
        target_scale - (best_candidate.signature.min_clearance or 0.0),
        0.0,
    )
    strong_pair_win = (
        best_candidate.signature.pair_issue_count
        < identity_candidate.signature.pair_issue_count
    )
    strong_short_edge_win = (
        identity_candidate.signature.short_edge_count >= 24
        and best_candidate.signature.short_edge_count
        <= max(8, identity_candidate.signature.short_edge_count // 6)
    )
    if (
        best_candidate.signature.pair_issue_count
        <= identity_candidate.signature.pair_issue_count
        and strong_short_edge_win
        and best_clearance_deficit <= identity_clearance_deficit + area_tolerance
    ):
        return False
    return not strong_pair_win


def _select_post_coverage_candidates_for_evaluation(
    identity_candidate: _CoverageSimplifyCandidate,
    candidates: Sequence[_CoverageSimplifyCandidate],
    *,
    target_scale: float,
    grid: float,
    output_min_area: float,
) -> list[_CoverageSimplifyCandidate]:
    selected: list[_CoverageSimplifyCandidate] = []
    if not candidates:
        return [identity_candidate]

    ranked_candidates = sorted(
        candidates,
        key=lambda candidate: _coverage_candidate_prescore(
            candidate,
            target_scale=target_scale,
            output_min_area=output_min_area,
        ),
    )
    identity_prescore = _coverage_candidate_prescore(
        identity_candidate,
        target_scale=target_scale,
        output_min_area=output_min_area,
    )
    best_candidate = ranked_candidates[0]
    best_prescore = _coverage_candidate_prescore(
        best_candidate,
        target_scale=target_scale,
        output_min_area=output_min_area,
    )
    if _should_finalize_identity_post_coverage_candidate(
        identity_candidate,
        best_candidate if best_prescore < identity_prescore else None,
        target_scale=target_scale,
        grid=grid,
        output_min_area=output_min_area,
    ):
        selected.append(identity_candidate)
    if best_prescore < identity_prescore:
        selected.append(best_candidate)
        if (
            len(ranked_candidates) > 1
            and ranked_candidates[1] is not best_candidate
            and (
                ranked_candidates[1].signature.pair_issue_count
                <= best_candidate.signature.pair_issue_count + 1
            )
            and (
                ranked_candidates[1].signature.short_edge_count
                <= best_candidate.signature.short_edge_count + 4
            )
            and ranked_candidates[1].difference_metrics["reference_minus_candidate_area"]
            <= best_candidate.difference_metrics["reference_minus_candidate_area"]
            + max(target_scale * target_scale, 16.0 * grid * grid, 1e-9)
            and _should_evaluate_additional_post_coverage_candidate(
                best_candidate,
                ranked_candidates[1],
                target_scale=target_scale,
                grid=grid,
                output_min_area=output_min_area,
            )
        ):
            selected.append(ranked_candidates[1])
    if not selected:
        selected.append(identity_candidate)
    return selected


def _combined_edit_zone(
    polygon: Polygon,
    *,
    short_edge_threshold: float,
    radius: float,
    diagnostics: dict[str, Any],
) -> BaseGeometry:
    zone = _short_edge_edit_zone(
        [polygon],
        short_edge_threshold=short_edge_threshold,
        radius=radius,
    )
    if radius <= 0:
        return zone
    try:
        clearance_line = shapely.minimum_clearance_line(polygon)
    except GEOSException as exc:
        _record_geos_exception(diagnostics, "minimum_clearance_line", exc)
        clearance_line = GeometryCollection()
    if clearance_line is None or clearance_line.is_empty:
        return zone
    line_zone = clearance_line.buffer(
        radius,
        quad_segs=1,
        join_style=BufferJoinStyle.mitre,
        mitre_limit=1000.0,
    )
    return line_zone if zone.is_empty else zone.union(line_zone)


def _ring_classification(
    polygon: Polygon,
    point: Point,
    *,
    tolerance: float,
) -> tuple[str | None, int | None]:
    if polygon.exterior.distance(point) <= tolerance:
        return "exterior", None
    for index, ring in enumerate(polygon.interiors):
        if LineString(ring.coords).distance(point) <= tolerance:
            return "hole", index
    return None, None


def _nearest_vertex_index(coords: np.ndarray, point: np.ndarray) -> int:
    distances = np.linalg.norm(coords - point, axis=1)
    return int(np.argmin(distances))


def _path_polygon(coords: np.ndarray, indices: Sequence[int]) -> Polygon | None:
    if len(indices) < 3:
        return None
    path = coords[list(indices)]
    closed = np.vstack([path, path[0]])
    polygon = Polygon(closed)
    if polygon.is_empty or polygon.area <= 0:
        return None
    return polygon


def _smaller_enclosed_ring_path(
    polygon: Polygon,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[int], Polygon] | None:
    try:
        clearance_line = shapely.minimum_clearance_line(polygon)
    except GEOSException as exc:
        _record_geos_exception(diagnostics, "minimum_clearance_line", exc)
        return None
    if clearance_line is None or clearance_line.is_empty:
        return None

    coords = np.asarray(polygon.exterior.coords[:-1], dtype=float)
    if len(coords) < 6:
        return None

    clearance_coords = np.asarray(clearance_line.coords, dtype=float)
    if len(clearance_coords) < 2:
        return None

    tolerance = max(grid, min_clearance * 0.25, 1e-6)
    start_type, _ = _ring_classification(
        polygon,
        Point(clearance_coords[0]),
        tolerance=tolerance,
    )
    end_type, _ = _ring_classification(
        polygon,
        Point(clearance_coords[1]),
        tolerance=tolerance,
    )
    if start_type != "exterior" or end_type != "exterior":
        return None

    raw_idx_a = _nearest_vertex_index(coords, clearance_coords[0])
    raw_idx_b = _nearest_vertex_index(coords, clearance_coords[1])
    if raw_idx_a == raw_idx_b:
        return None

    idx_a, idx_b = (
        (raw_idx_a, raw_idx_b)
        if raw_idx_a < raw_idx_b
        else (raw_idx_b, raw_idx_a)
    )
    forward = list(range(idx_a, idx_b + 1))
    backward = list(range(idx_b, len(coords))) + list(range(0, idx_a + 1))

    forward_polygon = _path_polygon(coords, forward)
    backward_polygon = _path_polygon(coords, backward)
    if forward_polygon is None or backward_polygon is None:
        return None

    courtyard_indices, courtyard_polygon = (
        (forward, forward_polygon)
        if forward_polygon.area <= backward_polygon.area
        else (backward, backward_polygon)
    )
    if courtyard_polygon.area < _COURTYARD_AREA_THRESHOLD * min_clearance * min_clearance:
        return None
    return courtyard_indices, courtyard_polygon


def _try_close_courtyard_passage(
    polygon: Polygon,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    enclosure = _smaller_enclosed_ring_path(
        polygon,
        min_clearance=min_clearance,
        grid=grid,
        diagnostics=diagnostics,
    )
    if enclosure is None:
        return None

    courtyard_indices, courtyard_polygon = enclosure
    coords = np.asarray(polygon.exterior.coords[:-1], dtype=float)
    idx_a = courtyard_indices[0]
    idx_b = courtyard_indices[-1]
    midpoint = (coords[idx_a] + coords[idx_b]) / 2.0

    snapped_coords = coords.copy()
    snapped_coords[idx_a] = midpoint
    snapped_coords[idx_b] = midpoint
    snapped_ring = np.vstack([snapped_coords, snapped_coords[0]])
    holes = [list(hole.coords) for hole in polygon.interiors]

    try:
        self_touching = Polygon(snapped_ring, holes) if holes else Polygon(snapped_ring)
        repaired = make_valid(self_touching)
    except (GEOSException, ValueError) as exc:
        _record_geos_exception(diagnostics, "courtyard_passage_close", exc)
        return None

    candidates = _extract_polygon_parts(repaired)
    if not candidates:
        return None
    best = max(candidates, key=lambda value: value.area)
    if len(best.interiors) <= len(polygon.interiors):
        return None

    passage_zone = LineString([tuple(coords[idx_a]), tuple(coords[idx_b])]).buffer(
        max(min_clearance, grid),
        quad_segs=1,
        join_style=BufferJoinStyle.mitre,
        mitre_limit=1000.0,
    )
    edit_zone = courtyard_polygon.union(passage_zone)
    return _RepairCandidate(
        polygon=orient(best, sign=1.0),
        edit_zone=edit_zone,
        operator="courtyard_passage",
    )


def _try_polygon_clearance_opening(
    polygon: Polygon,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    radius = min_clearance / 2.0
    if radius <= 0:
        return None
    opened_parts = _apply_opening(polygon, radius, grid, diagnostics)
    if len(opened_parts) != 1:
        return None
    closed_parts = _apply_closing(opened_parts[0], radius, grid, diagnostics)
    if len(closed_parts) != 1:
        return None
    return _RepairCandidate(
        polygon=orient(closed_parts[0], sign=1.0),
        edit_zone=_combined_edit_zone(
            polygon,
            short_edge_threshold=min_clearance,
            radius=_derive_local_edit_radius(min_clearance),
            diagnostics=diagnostics,
        ),
        operator="clearance_open_close",
    )


def _try_polygon_local_simplify(
    polygon: Polygon,
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    if tolerance <= 0:
        return None
    try:
        simplified = shapely.simplify(
            polygon,
            tolerance,
            preserve_topology=True,
        )
    except (GEOSException, ValueError, TypeError) as exc:
        _record_geos_exception(diagnostics, "local_simplify", exc)
        return None
    return _RepairCandidate(
        polygon=orient(simplified, sign=1.0),
        edit_zone=_combined_edit_zone(
            polygon,
            short_edge_threshold=tolerance,
            radius=_derive_local_edit_radius(tolerance),
            diagnostics=diagnostics,
        ),
        operator=f"local_simplify_{tolerance:.3f}",
    )


def _accept_local_candidate(
    reference_polygon: Polygon,
    candidate: _RepairCandidate,
    *,
    target_scale: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[Polygon | None, str | None, dict[str, float]]:
    canonical_parts: list[Polygon] = []
    for part in _canonicalize(candidate.polygon, grid, diagnostics):
        without_small_holes = _remove_small_holes(
            part,
            min_hole_area,
            diagnostics,
        )
        canonical_parts.extend(_canonicalize(without_small_holes, grid, diagnostics))
    filtered_parts = [
        orient(part, sign=1.0)
        for part in canonical_parts
        if part.area + 1e-12 >= min_area
    ]
    if len(filtered_parts) != 1:
        return None, "non_improving", {}

    accepted = filtered_parts[0]
    reference_signature = _polygon_defect_signature(
        reference_polygon,
        target_scale=target_scale,
    )
    candidate_signature = _polygon_defect_signature(
        accepted,
        target_scale=target_scale,
    )
    if not _signature_improves(
        reference_signature,
        candidate_signature,
        grid=grid,
    ):
        return None, "non_improving", {}

    difference_metrics = _difference_area_metrics(reference_polygon, accepted)
    area_balance_budget = _area_balance_budget(
        difference_metrics["symmetric_difference_area"],
        grid=grid,
        short_edge_count=reference_signature.short_edge_count,
        short_edge_threshold=target_scale,
    )
    polygon_change_outside_edit_zone = _change_outside_edit_zone(
        reference_polygon,
        accepted,
        edit_zone=candidate.edit_zone,
    )
    edit_zone_budget = max(float(candidate.edit_zone.area), grid * grid, 1e-9)
    if polygon_change_outside_edit_zone > edit_zone_budget:
        return None, "nonlocal", {
            **difference_metrics,
            "change_outside_edit_zone": polygon_change_outside_edit_zone,
            "area_balance_budget": area_balance_budget,
            "edit_zone_area": float(candidate.edit_zone.area),
        }
    if abs(difference_metrics["union_area_delta"]) > area_balance_budget:
        return None, "area_imbalance", {
            **difference_metrics,
            "change_outside_edit_zone": polygon_change_outside_edit_zone,
            "area_balance_budget": area_balance_budget,
            "edit_zone_area": float(candidate.edit_zone.area),
        }
    return accepted, None, {
        **difference_metrics,
        "change_outside_edit_zone": polygon_change_outside_edit_zone,
        "area_balance_budget": area_balance_budget,
        "edit_zone_area": float(candidate.edit_zone.area),
    }


def _build_source_geometry_lookup(
    polygons: Sequence[Polygon],
    source_map: Sequence[Sequence[int]],
) -> dict[int, BaseGeometry]:
    grouped: dict[int, list[Polygon]] = {}
    for polygon, indices in zip(polygons, source_map):
        for source_index in indices:
            grouped.setdefault(source_index, []).append(polygon)
    return {
        source_index: unary_union(parts)
        for source_index, parts in grouped.items()
    }


def _derive_source_recovery_distance(
    min_segment_length: float,
    grid: float,
) -> float:
    return max(2.0 * grid, 0.05 * max(min_segment_length, 0.0), 1e-6)


@dataclass(frozen=True)
class _SupportVertexIndex:
    cell_size: float
    buckets: dict[tuple[int, int], list[tuple[float, float]]]


def _build_support_vertex_index(
    vertices: Sequence[tuple[float, float]],
    *,
    distance_tolerance: float,
) -> _SupportVertexIndex:
    cell_size = max(distance_tolerance, 1e-9)
    buckets: dict[tuple[int, int], list[tuple[float, float]]] = {}
    inverse = 1.0 / cell_size
    for x, y in vertices:
        key = (int(np.floor(x * inverse)), int(np.floor(y * inverse)))
        buckets.setdefault(key, []).append((x, y))
    return _SupportVertexIndex(cell_size=cell_size, buckets=buckets)


def _boundary_vertices(geometry: BaseGeometry) -> list[tuple[float, float]]:
    vertices: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()
    for polygon in _extract_polygon_parts(geometry):
        for ring in [polygon.exterior, *polygon.interiors]:
            for x, y, *_ in list(ring.coords)[:-1]:
                vertex = (float(x), float(y))
                if vertex in seen:
                    continue
                seen.add(vertex)
                vertices.append(vertex)
    return vertices


def _recover_ring_vertices_from_support(
    coords: Sequence[Sequence[float]],
    *,
    support_vertices: Sequence[tuple[float, float]],
    support_vertex_index: _SupportVertexIndex | None = None,
    distance_tolerance: float,
) -> tuple[list[tuple[float, float]], int]:
    recovered: list[tuple[float, float]] = []
    recovered_count = 0
    ring = list(coords)
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    for x, y, *_ in ring:
        current = (float(x), float(y))
        best_vertex: tuple[float, float] | None = None
        best_distance: float | None = None
        if support_vertex_index is None:
            candidate_vertices = support_vertices
        else:
            inverse = 1.0 / support_vertex_index.cell_size
            cell_x = int(np.floor(current[0] * inverse))
            cell_y = int(np.floor(current[1] * inverse))
            bucketed: list[tuple[float, float]] = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    bucketed.extend(
                        support_vertex_index.buckets.get((cell_x + dx, cell_y + dy), ())
                    )
            candidate_vertices = bucketed

        for support_vertex in candidate_vertices:
            distance = float(
                np.hypot(
                    support_vertex[0] - current[0],
                    support_vertex[1] - current[1],
                )
            )
            if distance > distance_tolerance:
                continue
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_vertex = support_vertex
        if best_vertex is not None:
            recovered.append(best_vertex)
            if best_vertex != current:
                recovered_count += 1
        else:
            recovered.append(current)
    if recovered:
        recovered.append(recovered[0])
    return recovered, recovered_count


def _recover_polygon_vertices_from_support(
    polygon: Polygon,
    *,
    support_vertices: Sequence[tuple[float, float]],
    support_vertex_index: _SupportVertexIndex | None = None,
    grid: float,
) -> tuple[Polygon | None, int]:
    if not support_vertices:
        return None, 0

    shell, shell_count = _recover_ring_vertices_from_support(
        polygon.exterior.coords,
        support_vertices=support_vertices,
        support_vertex_index=support_vertex_index,
        distance_tolerance=_derive_source_recovery_distance(0.0, grid),
    )
    if shell_count == 0:
        return None, 0
    cleaned_shell = _clean_ring_coords(shell, grid)
    if cleaned_shell is None:
        return None, 0

    holes: list[list[tuple[float, float]]] = []
    recovered_count = shell_count
    for ring in polygon.interiors:
        recovered_hole, hole_count = _recover_ring_vertices_from_support(
            ring.coords,
            support_vertices=support_vertices,
            support_vertex_index=support_vertex_index,
            distance_tolerance=_derive_source_recovery_distance(0.0, grid),
        )
        recovered_count += hole_count
        cleaned_hole = _clean_ring_coords(recovered_hole, grid)
        if cleaned_hole is None:
            continue
        holes.append(cleaned_hole)

    candidate = Polygon(cleaned_shell, holes)
    if candidate.is_empty or candidate.area <= 0 or not candidate.is_valid:
        return None, 0
    return orient(candidate, sign=1.0), recovered_count


def _recover_polygon_source_coordinates(
    polygon: Polygon,
    *,
    support: BaseGeometry,
    support_vertices: Sequence[tuple[float, float]] | None = None,
    support_vertex_index: _SupportVertexIndex | None = None,
    support_source_count: int = 1,
    min_segment_length: float,
    grid: float,
) -> tuple[Polygon | None, str | None, dict[str, float] | None]:
    if support.is_empty:
        return None, None, None

    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=min_segment_length,
    )
    if not _signature_satisfies_scale_contract(
        reference_signature,
        target_scale=min_segment_length,
        grid=grid,
    ):
        return None, None, None

    distance_tolerance = _derive_source_recovery_distance(min_segment_length, grid)
    exact_tolerance = max(distance_tolerance, 2.0 * grid)
    exact_parts: list[Polygon] = []
    for part in _extract_polygon_parts(support):
        if part.distance(polygon) > exact_tolerance:
            continue
        cleaned_part = _clean_polygon_vertices(orient(part, sign=1.0), grid)
        if cleaned_part is None:
            continue
        exact_parts.append(cleaned_part)
    best_exact: tuple[tuple[float, float], Polygon, dict[str, float]] | None = None
    for part in exact_parts:
        candidate_signature = _polygon_defect_signature(
            part,
            target_scale=min_segment_length,
        )
        if not _signature_not_worse(reference_signature, candidate_signature, grid=grid):
            continue
        if not _signature_satisfies_scale_contract(
            candidate_signature,
            target_scale=min_segment_length,
            grid=grid,
        ):
            continue
        if candidate_signature.vertex_count > reference_signature.vertex_count:
            continue
        hausdorff = float(shapely.hausdorff_distance(polygon.boundary, part.boundary))
        if hausdorff > exact_tolerance:
            continue
        support_metrics = _difference_area_metrics(part, polygon)
        if support_metrics["symmetric_difference_area"] <= max(grid * grid, 1e-9):
            continue
        score = (hausdorff, support_metrics["symmetric_difference_area"])
        if best_exact is None or score < best_exact[0]:
            best_exact = (score, part, _difference_area_metrics(part, part))
    if best_exact is not None:
        return best_exact[1], "exact_source_polygon", best_exact[2]

    if support_source_count <= 1 and exact_parts:
        return None, None, None

    if support_vertices is None:
        support_vertices = _boundary_vertices(support)
    candidate, recovered_count = _recover_polygon_vertices_from_support(
        polygon,
        support_vertices=support_vertices,
        support_vertex_index=support_vertex_index,
        grid=grid,
    )
    if candidate is None or recovered_count == 0:
        return None, None, None

    candidate_signature = _polygon_defect_signature(
        candidate,
        target_scale=min_segment_length,
    )
    if not _signature_not_worse(reference_signature, candidate_signature, grid=grid):
        return None, None, None
    if not _signature_satisfies_scale_contract(
        candidate_signature,
        target_scale=min_segment_length,
        grid=grid,
    ):
        return None, None, None

    current_support_metrics = _difference_area_metrics(support, polygon)
    candidate_support_metrics = _difference_area_metrics(support, candidate)
    tolerance = max(grid * grid, 1e-9)
    if (
        candidate_support_metrics["reference_minus_candidate_area"]
        > current_support_metrics["reference_minus_candidate_area"] + tolerance
    ):
        return None, None, None
    if (
        candidate_support_metrics["candidate_minus_reference_area"]
        > current_support_metrics["candidate_minus_reference_area"] + tolerance
    ):
        return None, None, None
    if (
        candidate_support_metrics["symmetric_difference_area"]
        + tolerance
        >= current_support_metrics["symmetric_difference_area"]
    ):
        return None, None, None

    return candidate, "support_vertex_restore", candidate_support_metrics


def _normalize_single_polygon_candidate(
    geometry: BaseGeometry,
    *,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> Polygon | None:
    canonical_parts: list[Polygon] = []
    for part in _canonicalize(geometry, grid, diagnostics):
        without_small_holes = _remove_small_holes(
            part,
            min_hole_area,
            diagnostics,
        )
        canonical_parts.extend(_canonicalize(without_small_holes, grid, diagnostics))
    filtered_parts = [
        orient(part, sign=1.0)
        for part in canonical_parts
        if part.area + 1e-12 >= min_area
    ]
    if len(filtered_parts) != 1:
        return None
    return filtered_parts[0]


def _characteristic_width(polygon: Polygon) -> float:
    perimeter = float(polygon.length)
    if perimeter <= 1e-9:
        return 0.0
    return 2.0 * float(polygon.area) / perimeter


def _reclaim_source_supported_area(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    source_lookup: dict[int, BaseGeometry],
    min_segment_length: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    before_stats = _segment_length_stats(
        polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics["source_reclaim_short_edge_count_before"] = before_stats["short_edge_count"]
    if min_segment_length <= 0 or not polygons:
        diagnostics["source_reclaim_applied"] = False
        diagnostics["source_reclaim_short_edge_count_after"] = before_stats["short_edge_count"]
        return polygons, source_map

    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []
    applied_count = 0
    candidate_count = 0
    component_count = 0
    rejected_nonlocal_count = 0
    rejected_non_improving_count = 0
    edit_zone_area = 0.0
    change_outside_edit_zone = 0.0

    for polygon, indices in zip(polygons, source_map):
        support_parts = [source_lookup[index] for index in indices if index in source_lookup]
        if not support_parts:
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue
        support = unary_union(support_parts)
        missing = support.difference(polygon)
        missing_parts = [
            part
            for part in _extract_polygon_parts(missing)
            if part.area > max(0.5 * min_segment_length * min_segment_length, grid * grid, 1e-9)
            and _characteristic_width(part) + 1e-12 >= 0.5 * min_segment_length
        ]
        if not missing_parts:
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        candidate_count += 1
        component_count += len(missing_parts)
        current = polygon
        current_missing = _difference_area_metrics(support, current)[
            "reference_minus_candidate_area"
        ]

        for missing_part in sorted(missing_parts, key=lambda value: value.area, reverse=True):
            edit_zone = missing_part.buffer(
                max(min_segment_length, grid),
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
            unioned = current.union(missing_part)
            normalized = _normalize_single_polygon_candidate(
                unioned,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
            )
            if normalized is None:
                rejected_non_improving_count += 1
                continue

            current_signature = _polygon_defect_signature(
                current,
                target_scale=min_segment_length,
            )
            candidate_signature = _polygon_defect_signature(
                normalized,
                target_scale=min_segment_length,
            )
            if not _signature_not_worse(
                current_signature,
                candidate_signature,
                grid=grid,
            ):
                rejected_non_improving_count += 1
                continue

            candidate_missing = _difference_area_metrics(support, normalized)[
                "reference_minus_candidate_area"
            ]
            if candidate_missing + max(grid * grid, 1e-9) >= current_missing:
                rejected_non_improving_count += 1
                continue

            nonlocal_change = _change_outside_edit_zone(
                current,
                normalized,
                edit_zone=edit_zone,
            )
            if nonlocal_change > max(float(edit_zone.area), grid * grid, 1e-9):
                edit_zone_area += float(edit_zone.area)
                change_outside_edit_zone += nonlocal_change
                rejected_nonlocal_count += 1
                continue

            edit_zone_area += float(edit_zone.area)
            change_outside_edit_zone += nonlocal_change
            current = normalized
            current_missing = candidate_missing
            applied_count += 1
            _operator_count_increment(
                diagnostics,
                "source_reclaim_operator_applied",
                "raw_supported_union",
            )

        if current is not polygon:
            _operator_count_increment(
                diagnostics,
                "source_reclaim_operator_attempts",
                "raw_supported_union",
            )
        elif missing_parts:
            _operator_count_increment(
                diagnostics,
                "source_reclaim_operator_attempts",
                "raw_supported_union",
            )

        candidate_polygons.append(current)
        candidate_sources.append(indices)

    diagnostics["source_reclaim_candidate_count"] = candidate_count
    diagnostics["source_reclaim_applied_count"] = applied_count
    diagnostics["source_reclaim_component_count"] = component_count
    diagnostics["source_reclaim_edit_zone_area"] = edit_zone_area
    diagnostics["source_reclaim_change_outside_edit_zone"] = change_outside_edit_zone
    diagnostics["source_reclaim_rejected_nonlocal_count"] = rejected_nonlocal_count
    diagnostics["source_reclaim_rejected_non_improving_count"] = (
        rejected_non_improving_count
    )
    if applied_count == 0:
        diagnostics["source_reclaim_short_edge_count_after"] = before_stats[
            "short_edge_count"
        ]
        diagnostics["source_reclaim_applied"] = False
        diagnostics["source_reclaim_reference_minus_candidate_area"] = 0.0
        diagnostics["source_reclaim_candidate_minus_reference_area"] = 0.0
        diagnostics["source_reclaim_signed_area_delta"] = 0.0
        return polygons, source_map

    after_stats = _segment_length_stats(
        candidate_polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics["source_reclaim_short_edge_count_after"] = after_stats["short_edge_count"]
    diagnostics["source_reclaim_applied"] = True
    difference_metrics = _difference_area_metrics(
        unary_union(polygons),
        unary_union(candidate_polygons),
    )
    diagnostics["source_reclaim_reference_minus_candidate_area"] = (
        difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics["source_reclaim_candidate_minus_reference_area"] = (
        difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics["source_reclaim_signed_area_delta"] = difference_metrics["union_area_delta"]
    return _stable_sort(candidate_polygons, candidate_sources)


def _absorb_small_supported_components(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    source_lookup: dict[int, BaseGeometry],
    raw_support_union: BaseGeometry | None,
    min_area: float,
    min_segment_length: float,
    grid: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    if min_area <= 0 or not polygons:
        diagnostics["small_component_absorb_applied"] = False
        return polygons, source_map

    current_polygons = list(polygons)
    current_sources = [list(indices) for indices in source_map]
    applied_count = 0
    candidate_count = 0
    component_count = 0
    rejected_nonlocal_count = 0
    rejected_non_improving_count = 0
    edit_zone_area = 0.0
    change_outside_edit_zone = 0.0
    support_radius = max(min_segment_length, grid)
    bridge_search_radius = max(
        4.0 * support_radius,
        float(np.sqrt(max(min_area, 0.0))),
        support_radius + grid,
    )
    support_tolerance = max(grid * grid, 1e-9)
    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}

    while True:
        small_indices = [
            index
            for index, polygon in enumerate(current_polygons)
            if polygon.area + 1e-12 < min_area
        ]
        if not small_indices:
            break

        accepted = False
        for small_index in sorted(
            small_indices,
            key=lambda index: current_polygons[index].area,
        ):
            small_polygon = current_polygons[small_index]
            small_sources = current_sources[small_index]
            component_count += 1

            best_candidate: tuple[
                tuple[float, ...],
                int,
                Polygon,
                list[int],
                BaseGeometry,
                dict[str, float],
            ] | None = None

            for neighbor_index, neighbor_polygon in enumerate(current_polygons):
                if neighbor_index == small_index:
                    continue
                if neighbor_polygon.area + 1e-12 < min_area:
                    continue
                if small_polygon.distance(neighbor_polygon) > bridge_search_radius:
                    continue

                operator_attempts["source_supported_absorb"] = (
                    operator_attempts.get("source_supported_absorb", 0) + 1
                )
                candidate_count += 1

                local_sources = sorted(
                    set(small_sources).union(current_sources[neighbor_index])
                )
                source_support_parts = [
                    source_lookup[index]
                    for index in local_sources
                    if index in source_lookup
                ]

                edit_zone = unary_union([small_polygon, neighbor_polygon]).convex_hull.buffer(
                    support_radius,
                    quad_segs=1,
                    join_style=BufferJoinStyle.mitre,
                    mitre_limit=1000.0,
                )
                support_geometries: list[BaseGeometry] = []
                if source_support_parts:
                    support_geometries.append(unary_union(source_support_parts))
                if raw_support_union is not None and not raw_support_union.is_empty:
                    support_geometries.append(raw_support_union)
                if not support_geometries:
                    continue
                support_patch = unary_union(support_geometries).intersection(edit_zone)
                if support_patch.is_empty:
                    continue
                buffered_small = small_polygon.buffer(grid)
                buffered_neighbor = neighbor_polygon.buffer(grid)
                if not (
                    support_patch.intersects(buffered_small)
                    and support_patch.intersects(buffered_neighbor)
                ):
                    continue

                candidate_geometry = neighbor_polygon.union(small_polygon).union(
                    support_patch
                )
                normalized = _normalize_single_polygon_candidate(
                    candidate_geometry,
                    grid=grid,
                    min_area=0.0,
                    min_hole_area=min_hole_area,
                    diagnostics=diagnostics,
                )
                if normalized is None:
                    rejected_non_improving_count += 1
                    continue

                local_diagnostics = _empty_diagnostics(1)
                local_diagnostics["collect_stage_metrics"] = False
                regularized_polygons, regularized_sources = _simplify_polygons_for_meshing(
                    [normalized],
                    [local_sources],
                    min_segment_length=min_segment_length,
                    grid=grid,
                    min_area=0.0,
                    min_hole_area=min_hole_area,
                    diagnostics=local_diagnostics,
                )
                regularized_polygons, regularized_sources = _regularize_low_clearance_polygons(
                    regularized_polygons,
                    regularized_sources,
                    min_clearance=min_segment_length,
                    grid=grid,
                    min_area=0.0,
                    min_hole_area=min_hole_area,
                    diagnostics=local_diagnostics,
                )
                diagnostics["geos_exception_count"] += local_diagnostics["geos_exception_count"]
                diagnostics["geos_exception_messages"].extend(
                    local_diagnostics["geos_exception_messages"]
                )
                if len(regularized_polygons) != 1:
                    rejected_non_improving_count += 1
                    continue
                normalized = regularized_polygons[0]
                local_sources = regularized_sources[0]

                reference_signature = _coverage_defect_signature(
                    [neighbor_polygon, small_polygon],
                    target_scale=min_segment_length,
                )
                candidate_signature = _coverage_defect_signature(
                    [normalized],
                    target_scale=min_segment_length,
                )
                if not _coverage_signature_not_worse(
                    reference_signature,
                    candidate_signature,
                    grid=grid,
                    target_scale=min_segment_length,
                ):
                    rejected_non_improving_count += 1
                    continue

                unaffected_polygons = [
                    polygon
                    for index, polygon in enumerate(current_polygons)
                    if index not in (small_index, neighbor_index)
                ]
                unaffected_union = unary_union(unaffected_polygons)
                if normalized.intersection(unaffected_union).area > support_tolerance:
                    rejected_non_improving_count += 1
                    continue

                reference_local_union = unary_union([neighbor_polygon, small_polygon])
                difference_metrics = _difference_area_metrics(
                    reference_local_union,
                    normalized,
                )
                nonlocal_change = _change_outside_edit_zone(
                    reference_local_union,
                    normalized,
                    edit_zone=edit_zone,
                )
                if nonlocal_change > max(float(edit_zone.area), support_tolerance):
                    rejected_nonlocal_count += 1
                    continue
                supported_reference = reference_local_union.union(support_patch)
                unsupported_extra_area = float(
                    normalized.difference(supported_reference).area
                )
                if (
                    difference_metrics["reference_minus_candidate_area"]
                    > support_tolerance
                ):
                    rejected_non_improving_count += 1
                    continue
                if unsupported_extra_area > support_tolerance:
                    rejected_non_improving_count += 1
                    continue

                score = (
                    candidate_signature.short_edge_count,
                    -(candidate_signature.min_edge_length or 0.0),
                    unsupported_extra_area,
                    difference_metrics["candidate_minus_reference_area"],
                    difference_metrics["symmetric_difference_area"],
                    normalized.area,
                )
                if best_candidate is None or score < best_candidate[0]:
                    best_candidate = (
                        score,
                        neighbor_index,
                        normalized,
                        local_sources,
                        edit_zone,
                        difference_metrics,
                    )

            if best_candidate is None:
                continue

            (
                _,
                neighbor_index,
                normalized,
                local_sources,
                edit_zone,
                difference_metrics,
            ) = best_candidate
            retained_polygons: list[Polygon] = []
            retained_sources: list[list[int]] = []
            for index, polygon in enumerate(current_polygons):
                if index == small_index:
                    continue
                if index == neighbor_index:
                    retained_polygons.append(normalized)
                    retained_sources.append(local_sources)
                    continue
                retained_polygons.append(polygon)
                retained_sources.append(current_sources[index])

            current_polygons, current_sources = _stable_sort(
                retained_polygons,
                retained_sources,
            )
            operator_applied["source_supported_absorb"] = (
                operator_applied.get("source_supported_absorb", 0) + 1
            )
            applied_count += 1
            edit_zone_area += float(edit_zone.area)
            change_outside_edit_zone += 0.0
            diagnostics["small_component_absorb_reference_minus_candidate_area"] += (
                difference_metrics["reference_minus_candidate_area"]
            )
            diagnostics["small_component_absorb_candidate_minus_reference_area"] += (
                difference_metrics["candidate_minus_reference_area"]
            )
            diagnostics["small_component_absorb_signed_area_delta"] += (
                difference_metrics["union_area_delta"]
            )
            accepted = True
            break

        if not accepted:
            break

    diagnostics["small_component_absorb_candidate_count"] = candidate_count
    diagnostics["small_component_absorb_applied_count"] = applied_count
    diagnostics["small_component_absorb_component_count"] = component_count
    diagnostics["small_component_absorb_edit_zone_area"] = edit_zone_area
    diagnostics["small_component_absorb_change_outside_edit_zone"] = (
        change_outside_edit_zone
    )
    diagnostics["small_component_absorb_rejected_nonlocal_count"] = (
        rejected_nonlocal_count
    )
    diagnostics["small_component_absorb_rejected_non_improving_count"] = (
        rejected_non_improving_count
    )
    diagnostics["small_component_absorb_operator_attempts"] = operator_attempts
    diagnostics["small_component_absorb_operator_applied"] = operator_applied
    diagnostics["small_component_absorb_applied"] = applied_count > 0
    return current_polygons, current_sources


def _recover_source_supported_coordinates(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    source_lookup: dict[int, BaseGeometry],
    min_segment_length: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    before_stats = _segment_length_stats(
        polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics["source_coordinate_recovery_short_edge_count_before"] = before_stats[
        "short_edge_count"
    ]
    if min_segment_length <= 0 or not polygons:
        diagnostics["source_coordinate_recovery_applied"] = False
        diagnostics["source_coordinate_recovery_short_edge_count_after"] = before_stats[
            "short_edge_count"
        ]
        return polygons, source_map

    current_polygons = list(polygons)
    current_sources = [list(indices) for indices in source_map]
    candidate_count = 0
    applied_count = 0
    exact_count = 0
    vertex_count = 0
    rejected_overlap_count = 0
    rejected_non_improving_count = 0
    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}
    overlap_tolerance = max(grid * grid, 1e-9)
    recovery_query_radius = max(min_segment_length, grid, 1e-9)
    support_cache: dict[
        tuple[int, ...],
        tuple[BaseGeometry, list[tuple[float, float]], _SupportVertexIndex],
    ] = {}
    proposal_map: dict[int, tuple[Polygon, str, dict[str, float]]] = {}
    recovery_signature_cache = _CoverageEvalCache()

    for polygon_index, (polygon, indices) in enumerate(zip(polygons, source_map)):
        support_key = tuple(indices)
        cached_support = support_cache.get(support_key)
        if cached_support is None:
            support_parts = [
                source_lookup[source_index]
                for source_index in indices
                if source_index in source_lookup
            ]
            if not support_parts:
                continue
            if len(support_parts) == 1:
                support = support_parts[0]
            else:
                support = unary_union(support_parts)
            support_vertices = _boundary_vertices(support)
            support_vertex_index = _build_support_vertex_index(
                support_vertices,
                distance_tolerance=_derive_source_recovery_distance(
                    min_segment_length,
                    grid,
                ),
            )
            cached_support = (support, support_vertices, support_vertex_index)
            support_cache[support_key] = cached_support
        support, support_vertices, support_vertex_index = cached_support
        if support.is_empty:
            continue
        candidate_count += 1

        candidate, operator, support_metrics = _recover_polygon_source_coordinates(
            polygon,
            support=support,
            support_vertices=support_vertices,
            support_vertex_index=support_vertex_index,
            support_source_count=len(indices),
            min_segment_length=min_segment_length,
            grid=grid,
        )
        if operator is None or candidate is None or support_metrics is None:
            rejected_non_improving_count += 1
            continue
        operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
        proposal_map[polygon_index] = (candidate, operator, support_metrics)

    def _apply_recovery_candidate(
        polygon_index: int,
        candidate: Polygon,
        operator: str,
        support_metrics: dict[str, float],
    ) -> None:
        nonlocal applied_count, exact_count, vertex_count
        current_polygons[polygon_index] = candidate
        current_sources[polygon_index] = list(source_map[polygon_index])
        operator_applied[operator] = operator_applied.get(operator, 0) + 1
        applied_count += 1
        if operator == "exact_source_polygon":
            exact_count += 1
        elif operator == "support_vertex_restore":
            vertex_count += 1
        diagnostics["source_coordinate_recovery_support_reference_minus_candidate_area"] += (
            support_metrics["reference_minus_candidate_area"]
        )
        diagnostics["source_coordinate_recovery_support_candidate_minus_reference_area"] += (
            support_metrics["candidate_minus_reference_area"]
        )
        diagnostics["source_coordinate_recovery_support_signed_area_delta"] += (
            support_metrics["union_area_delta"]
        )

    def _candidate_preserves_local_coverage(
        polygon_index: int,
        candidate: Polygon,
        neighbor_indices: set[int],
    ) -> bool:
        relevant_indices = [
            other_index
            for other_index in sorted(neighbor_indices)
            if other_index != polygon_index
        ]
        if not relevant_indices:
            return True

        reference_subset = [current_polygons[polygon_index]]
        candidate_subset = [candidate]
        for other_index in relevant_indices:
            other = current_polygons[other_index]
            if other.is_empty:
                continue
            reference_subset.append(other)
            candidate_subset.append(other)

        if len(reference_subset) <= 1:
            return True

        reference_signature = _cached_coverage_defect_signature(
            recovery_signature_cache,
            reference_subset,
            target_scale=min_segment_length,
        )
        candidate_signature = _cached_coverage_defect_signature(
            recovery_signature_cache,
            candidate_subset,
            target_scale=min_segment_length,
        )
        return _coverage_signature_not_worse(
            reference_signature,
            candidate_signature,
            grid=grid,
            target_scale=min_segment_length,
        )

    if proposal_map:
        recovery_tree: STRtree | None = None
        dirty_indices: set[int] = set()
        for polygon_index in sorted(proposal_map):
            candidate, operator, support_metrics = proposal_map[polygon_index]
            if (
                recovery_tree is None
                or len(dirty_indices) > _RECOVERY_TREE_REBUILD_THRESHOLD
            ):
                recovery_tree = STRtree(np.asarray(current_polygons, dtype=object))
                dirty_indices.clear()

            candidate_indices = {
                int(index)
                for index in recovery_tree.query(
                    _expanded_query_geometry(
                        current_polygons[polygon_index],
                        radius=recovery_query_radius,
                    )
                )
            }
            candidate_indices.update(
                int(index)
                for index in recovery_tree.query(
                    _expanded_query_geometry(
                        candidate,
                        radius=recovery_query_radius,
                    )
                )
            )
            candidate_indices.update(dirty_indices)
            candidate_indices.discard(polygon_index)

            overlap_area = 0.0
            for other_index in candidate_indices:
                other = current_polygons[other_index]
                if other.is_empty or not candidate.intersects(other):
                    continue
                overlap_area += candidate.intersection(other).area
                if overlap_area > overlap_tolerance:
                    break
            if overlap_area > overlap_tolerance:
                rejected_overlap_count += 1
                continue

            if not _candidate_preserves_local_coverage(
                polygon_index,
                candidate,
                candidate_indices,
            ):
                rejected_non_improving_count += 1
                continue

            tree_geometry_changed = not candidate.equals_exact(
                current_polygons[polygon_index],
                tolerance=0.0,
            )
            _apply_recovery_candidate(
                polygon_index,
                candidate,
                operator,
                support_metrics,
            )
            if tree_geometry_changed:
                dirty_indices.add(polygon_index)

    diagnostics["source_coordinate_recovery_candidate_count"] = candidate_count
    diagnostics["source_coordinate_recovery_applied_count"] = applied_count
    diagnostics["source_coordinate_recovery_exact_count"] = exact_count
    diagnostics["source_coordinate_recovery_vertex_count"] = vertex_count
    diagnostics["source_coordinate_recovery_rejected_overlap_count"] = (
        rejected_overlap_count
    )
    diagnostics["source_coordinate_recovery_rejected_non_improving_count"] = (
        rejected_non_improving_count
    )
    diagnostics["source_coordinate_recovery_operator_attempts"] = operator_attempts
    diagnostics["source_coordinate_recovery_operator_applied"] = operator_applied
    if applied_count == 0:
        diagnostics["source_coordinate_recovery_short_edge_count_after"] = before_stats[
            "short_edge_count"
        ]
        diagnostics["source_coordinate_recovery_applied"] = False
        diagnostics["source_coordinate_recovery_reference_minus_candidate_area"] = 0.0
        diagnostics["source_coordinate_recovery_candidate_minus_reference_area"] = 0.0
        diagnostics["source_coordinate_recovery_signed_area_delta"] = 0.0
        return polygons, source_map

    current_polygons, current_sources = _stable_sort(current_polygons, current_sources)
    after_stats = _segment_length_stats(
        current_polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics["source_coordinate_recovery_short_edge_count_after"] = after_stats[
        "short_edge_count"
    ]
    diagnostics["source_coordinate_recovery_applied"] = True
    difference_metrics = _difference_area_metrics(
        unary_union(polygons),
        unary_union(current_polygons),
    )
    diagnostics["source_coordinate_recovery_reference_minus_candidate_area"] = (
        difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics["source_coordinate_recovery_candidate_minus_reference_area"] = (
        difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics["source_coordinate_recovery_signed_area_delta"] = difference_metrics[
        "union_area_delta"
    ]
    return current_polygons, current_sources


def _regularize_coverage_for_meshing(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_segment_length: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> tuple[list[Polygon], list[list[int]]]:
    def _regularize_ring_contacts(
        input_polygons: list[Polygon],
        input_sources: list[list[int]],
    ) -> tuple[list[Polygon], list[list[int]]] | None:
        normalization_distance = max(
            grid,
            min_segment_length * 0.0625,
            1e-6,
        )
        candidate_polygons: list[Polygon] = []
        candidate_sources: list[list[int]] = []
        changed_count = 0
        component_count = 0
        failed_count = 0
        area_delta = 0.0

        for polygon, indices in zip(input_polygons, input_sources):
            ring_contact_count = _polygon_ring_boundary_contact_count(polygon)
            if ring_contact_count == 0:
                candidate_polygons.append(polygon)
                candidate_sources.append(indices)
                continue

            normalized_parts = _regularize_ring_contact_polygon(
                polygon,
                min_distance=normalization_distance,
                target_clearance=min_segment_length,
                grid=grid,
            )
            remaining_contacts = sum(
                _polygon_ring_boundary_contact_count(part)
                for part in normalized_parts
            )
            if remaining_contacts > 0:
                failed_count += 1
                candidate_polygons.append(polygon)
                candidate_sources.append(indices)
                continue

            changed_count += 1
            component_count += len(normalized_parts)
            area_delta += float(sum(part.area for part in normalized_parts) - polygon.area)
            for part in normalized_parts:
                candidate_polygons.append(part)
                candidate_sources.append(list(indices))

        diagnostics["coverage_meshing_regularization_ring_contact_polygon_count"] = (
            changed_count
        )
        diagnostics["coverage_meshing_regularization_ring_contact_component_count"] = (
            component_count
        )
        diagnostics["coverage_meshing_regularization_ring_contact_failed_count"] = (
            failed_count
        )
        diagnostics["coverage_meshing_regularization_ring_contact_area_delta"] = (
            area_delta
        )

        if changed_count == 0:
            return None

        return _stable_sort(candidate_polygons, candidate_sources)

    if cache is None:
        cache = _CoverageEvalCache()
    before_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=min_segment_length,
    )
    diagnostics["coverage_meshing_regularization_short_edge_count_before"] = (
        before_signature.short_edge_count
    )
    diagnostics["coverage_meshing_regularization_pair_issue_count_before"] = (
        before_signature.pair_issue_count
    )
    diagnostics["coverage_meshing_regularization_ring_contact_count_before"] = (
        before_signature.ring_contact_count
    )

    if min_segment_length <= 0 or not polygons:
        diagnostics["coverage_meshing_regularization_applied"] = False
        diagnostics["coverage_meshing_regularization_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        diagnostics["coverage_meshing_regularization_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["coverage_meshing_regularization_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        return polygons, source_map

    ring_contact_applied = False
    ring_contact_candidate = None
    if before_signature.ring_contact_count > 0:
        ring_contact_candidate = _regularize_ring_contacts(polygons, source_map)
    if ring_contact_candidate is not None:
        polygons, source_map = ring_contact_candidate
        before_signature = _cached_coverage_defect_signature(
            cache,
            polygons,
            target_scale=min_segment_length,
        )
        ring_contact_applied = True

    if (
        before_signature.short_edge_count == 0
        and before_signature.ring_contact_count == 0
        and before_signature.pair_issue_count == 0
        and max(min_segment_length - (before_signature.min_clearance or 0.0), 0.0)
        <= max(grid, 1e-9)
    ):
        diagnostics["coverage_meshing_regularization_applied"] = ring_contact_applied
        diagnostics["coverage_meshing_regularization_selected_branch"] = (
            "ring_contacts" if ring_contact_applied else "identity"
        )
        diagnostics["coverage_meshing_regularization_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        diagnostics["coverage_meshing_regularization_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["coverage_meshing_regularization_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        return polygons, source_map

    graphical_candidate = _simplify_coverage_short_edge_graphically(
        polygons,
        source_map,
        tolerance=min_segment_length,
        grid=grid,
        diagnostics=diagnostics,
        cache=cache,
    )
    if (
        before_signature.pair_issue_count == 0
        and graphical_candidate is not None
        and _coverage_signature_satisfies_scale_contract(
            graphical_candidate.signature,
            target_scale=min_segment_length,
            grid=grid,
        )
        and _coverage_candidate_has_small_fidelity_drift(
            graphical_candidate,
            target_scale=min_segment_length,
            grid=grid,
        )
    ):
        diagnostics["coverage_meshing_regularization_local_candidate_attempted"] = False
        diagnostics["coverage_meshing_regularization_applied"] = True
        diagnostics["coverage_meshing_regularization_selected_branch"] = (
            graphical_candidate.label
        )
        diagnostics["coverage_meshing_regularization_short_edge_count_after"] = (
            graphical_candidate.signature.short_edge_count
        )
        diagnostics["coverage_meshing_regularization_pair_issue_count_after"] = (
            graphical_candidate.signature.pair_issue_count
        )
        diagnostics["coverage_meshing_regularization_reference_minus_candidate_area"] = (
            graphical_candidate.difference_metrics["reference_minus_candidate_area"]
        )
        diagnostics["coverage_meshing_regularization_candidate_minus_reference_area"] = (
            graphical_candidate.difference_metrics["candidate_minus_reference_area"]
        )
        diagnostics["coverage_meshing_regularization_signed_area_delta"] = (
            graphical_candidate.difference_metrics["union_area_delta"]
        )
        diagnostics["coverage_meshing_regularization_operator_attempts"] = dict(
            graphical_candidate.operator_attempts
        )
        diagnostics["coverage_meshing_regularization_operator_applied"] = dict(
            graphical_candidate.operator_applied
        )
        return _stable_sort(
            list(graphical_candidate.polygons),
            [list(indices) for indices in graphical_candidate.source_map],
        )

    global_candidate = _simplify_coverage_globally(
        polygons,
        source_map,
        tolerance=min_segment_length,
        grid=grid,
        diagnostics=diagnostics,
        cache=cache,
    )
    attempt_local_candidate = _should_attempt_meshing_local_candidate(
        before_signature,
        global_candidate,
        target_scale=min_segment_length,
        grid=grid,
    )
    diagnostics["coverage_meshing_regularization_local_candidate_attempted"] = (
        attempt_local_candidate
    )
    local_candidate = (
        _simplify_coverage_locally(
            polygons,
            source_map,
            tolerance=min_segment_length,
            grid=grid,
            diagnostics=diagnostics,
            enable_patch_union_fallback=True,
            enable_pair_cluster_rescue=True,
            cache=cache,
        )
        if attempt_local_candidate
        else None
    )
    residual_candidate = _simplify_coverage_residual_scale_polygons(
        polygons,
        source_map,
        tolerance=min_segment_length,
        grid=grid,
        min_area=min_area,
        min_hole_area=min_hole_area,
        diagnostics=diagnostics,
        cache=cache,
    )

    if (
        graphical_candidate is None
        and global_candidate is None
        and local_candidate is None
        and residual_candidate is None
    ):
        diagnostics["coverage_meshing_regularization_applied"] = ring_contact_applied
        diagnostics["coverage_meshing_regularization_selected_branch"] = (
            "ring_contacts" if ring_contact_applied else "identity"
        )
        diagnostics["coverage_meshing_regularization_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        diagnostics["coverage_meshing_regularization_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["coverage_meshing_regularization_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        diagnostics["coverage_meshing_regularization_reference_minus_candidate_area"] = (
            0.0
        )
        diagnostics["coverage_meshing_regularization_candidate_minus_reference_area"] = (
            0.0
        )
        diagnostics["coverage_meshing_regularization_signed_area_delta"] = 0.0
        diagnostics["coverage_meshing_regularization_operator_attempts"] = {}
        diagnostics["coverage_meshing_regularization_operator_applied"] = {}
        return polygons, source_map

    reference_union = _cached_union(cache, polygons)

    best_polygons = list(polygons)
    best_sources = [list(indices) for indices in source_map]
    best_signature = before_signature
    best_difference_metrics = _zero_difference_metrics()
    best_label = "identity"
    best_score = (
        *_coverage_signature_score(
            best_signature,
            target_scale=min_segment_length,
        ),
        best_difference_metrics["reference_minus_candidate_area"],
        best_difference_metrics["candidate_minus_reference_area"],
        best_difference_metrics["symmetric_difference_area"],
        abs(best_difference_metrics["union_area_delta"]),
    )
    best_operator_attempts: dict[str, int] = {}
    best_operator_applied: dict[str, int] = {}

    for candidate in (
        graphical_candidate,
        global_candidate,
        local_candidate,
        residual_candidate,
    ):
        if candidate is None:
            continue

        if _coverage_signature_satisfies_scale_contract(
            candidate.signature,
            target_scale=min_segment_length,
            grid=grid,
        ):
            candidate_polygons = candidate.polygons
            candidate_sources = candidate.source_map
            candidate_signature = candidate.signature
            difference_metrics = candidate.difference_metrics
        else:
            local_diagnostics = _empty_diagnostics(len(candidate.polygons))
            local_diagnostics["collect_stage_metrics"] = False
            local_diagnostics["enable_logging"] = False

            candidate_polygons, candidate_sources = _simplify_polygons_for_meshing(
                candidate.polygons,
                candidate.source_map,
                min_segment_length=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=local_diagnostics,
            )
            candidate_polygons, candidate_sources = _regularize_low_clearance_polygons(
                candidate_polygons,
                candidate_sources,
                min_clearance=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=local_diagnostics,
            )
            diagnostics["geos_exception_count"] += local_diagnostics["geos_exception_count"]
            diagnostics["geos_exception_messages"].extend(
                local_diagnostics["geos_exception_messages"]
            )
            candidate_signature = _cached_coverage_defect_signature(
                cache,
                candidate_polygons,
                target_scale=min_segment_length,
            )
            difference_metrics = _cached_difference_area_metrics(
                cache,
                polygons,
                candidate_polygons,
            )
        score = (
            *_coverage_signature_score(
                candidate_signature,
                target_scale=min_segment_length,
            ),
            difference_metrics["reference_minus_candidate_area"],
            difference_metrics["candidate_minus_reference_area"],
            difference_metrics["symmetric_difference_area"],
            abs(difference_metrics["union_area_delta"]),
        )
        if score >= best_score:
            continue

        best_polygons = candidate_polygons
        best_sources = candidate_sources
        best_signature = candidate_signature
        best_difference_metrics = difference_metrics
        best_label = candidate.label
        best_score = score
        best_operator_attempts = dict(candidate.operator_attempts)
        best_operator_applied = dict(candidate.operator_applied)

    def _repair_nested_courtyard_pair_issues(
        input_polygons: list[Polygon],
        input_sources: list[list[int]],
    ) -> tuple[list[Polygon], list[list[int]]] | None:
        current_polygons = list(input_polygons)
        current_sources = [list(indices) for indices in input_sources]
        progress = False

        while True:
            pair_candidates = _cached_pair_issue_candidates(
                cache,
                current_polygons,
                target_scale=min_segment_length,
            )
            if not pair_candidates:
                break

            repaired_this_round = False
            for left_index, right_index, _, _ in pair_candidates:
                pair_polygons = [current_polygons[left_index], current_polygons[right_index]]
                pair_sources = [current_sources[left_index], current_sources[right_index]]
                if (
                    _courtyard_pair_shrink_budget_context(
                        pair_polygons,
                        tolerance=min_segment_length,
                        grid=grid,
                    )
                    is None
                ):
                    continue

                reference_signature = _cached_coverage_defect_signature(
                    cache,
                    pair_polygons,
                    target_scale=min_segment_length,
                )
                best_pair_candidate: tuple[
                    list[Polygon],
                    list[list[int]],
                    _CoverageDefectSignature,
                ] | None = None
                best_pair_area_loss: float | None = None
                for shrink_radius in _iter_pair_issue_shrink_radii(
                    tolerance=min_segment_length,
                    grid=grid,
                ):
                    candidate = _apply_local_smaller_polygon_shrink_operator(
                        pair_polygons,
                        pair_sources,
                        radius=shrink_radius,
                        grid=grid,
                        diagnostics=diagnostics,
                    )
                    if candidate is None:
                        continue
                    candidate_polygons, candidate_sources = candidate
                    candidate_signature = _cached_coverage_defect_signature(
                        cache,
                        candidate_polygons,
                        target_scale=min_segment_length,
                    )
                    if candidate_signature.pair_issue_count >= reference_signature.pair_issue_count:
                        continue
                    candidate_difference = _cached_difference_area_metrics(
                        cache,
                        pair_polygons,
                        candidate_polygons,
                    )
                    area_loss = candidate_difference["reference_minus_candidate_area"]
                    if (
                        best_pair_candidate is None
                        or candidate_signature.pair_issue_count
                        < best_pair_candidate[2].pair_issue_count
                        or (
                            candidate_signature.pair_issue_count
                            == best_pair_candidate[2].pair_issue_count
                            and (
                                best_pair_area_loss is None
                                or area_loss < best_pair_area_loss
                            )
                        )
                    ):
                        best_pair_candidate = (
                            candidate_polygons,
                            candidate_sources,
                            candidate_signature,
                        )
                        best_pair_area_loss = area_loss

                if best_pair_candidate is None:
                    continue

                replacement_polygons, replacement_sources, _ = best_pair_candidate
                unaffected_polygons = [
                    polygon
                    for index, polygon in enumerate(current_polygons)
                    if index not in (left_index, right_index)
                ]
                unaffected_sources = [
                    indices
                    for index, indices in enumerate(current_sources)
                    if index not in (left_index, right_index)
                ]
                current_polygons, current_sources = _stable_sort(
                    unaffected_polygons + replacement_polygons,
                    unaffected_sources + replacement_sources,
                )
                repaired_this_round = True
                progress = True
                break

            if not repaired_this_round:
                break

        if not progress:
            return None
        return current_polygons, current_sources

    if best_signature.pair_issue_count > 0:
        rescue_candidate = _repair_nested_courtyard_pair_issues(
            best_polygons,
            best_sources,
        )
        if rescue_candidate is not None:
            rescue_polygons, rescue_sources = rescue_candidate
            rescue_diagnostics = _empty_diagnostics(len(rescue_polygons))
            rescue_diagnostics["collect_stage_metrics"] = False
            rescue_diagnostics["enable_logging"] = False
            rescue_polygons, rescue_sources = _simplify_polygons_for_meshing(
                rescue_polygons,
                rescue_sources,
                min_segment_length=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=rescue_diagnostics,
            )
            rescue_polygons, rescue_sources = _regularize_low_clearance_polygons(
                rescue_polygons,
                rescue_sources,
                min_clearance=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=rescue_diagnostics,
            )
            rescue_signature = _cached_coverage_defect_signature(
                cache,
                rescue_polygons,
                target_scale=min_segment_length,
            )
            rescue_difference_metrics = _cached_difference_area_metrics(
                cache,
                polygons,
                rescue_polygons,
            )
            if _coverage_signature_improves(
                best_signature,
                rescue_signature,
                grid=grid,
                target_scale=min_segment_length,
            ):
                best_polygons = rescue_polygons
                best_sources = rescue_sources
                best_signature = rescue_signature
                best_difference_metrics = rescue_difference_metrics
                best_label = "nested_courtyard_pair_rescue"
                best_operator_attempts = dict(best_operator_attempts)
                best_operator_attempts["coverage_pair_issue_shrink_nested_courtyard"] = (
                    best_operator_attempts.get(
                        "coverage_pair_issue_shrink_nested_courtyard",
                        0,
                    )
                    + 1
                )
                best_operator_applied = dict(best_operator_applied)
                best_operator_applied["coverage_pair_issue_shrink_nested_courtyard"] = (
                    best_operator_applied.get(
                        "coverage_pair_issue_shrink_nested_courtyard",
                        0,
                    )
                    + 1
                )

    diagnostics["coverage_meshing_regularization_applied"] = (
        ring_contact_applied or best_label != "identity"
    )
    diagnostics["coverage_meshing_regularization_selected_branch"] = (
        "ring_contacts"
        if ring_contact_applied and best_label == "identity"
        else best_label
    )
    diagnostics["coverage_meshing_regularization_short_edge_count_after"] = (
        best_signature.short_edge_count
    )
    diagnostics["coverage_meshing_regularization_pair_issue_count_after"] = (
        best_signature.pair_issue_count
    )
    diagnostics["coverage_meshing_regularization_ring_contact_count_after"] = (
        best_signature.ring_contact_count
    )
    diagnostics["coverage_meshing_regularization_reference_minus_candidate_area"] = (
        best_difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics["coverage_meshing_regularization_candidate_minus_reference_area"] = (
        best_difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics["coverage_meshing_regularization_signed_area_delta"] = (
        best_difference_metrics["union_area_delta"]
    )
    diagnostics["coverage_meshing_regularization_operator_attempts"] = (
        best_operator_attempts
    )
    diagnostics["coverage_meshing_regularization_operator_applied"] = (
        best_operator_applied
    )
    return _stable_sort(best_polygons, best_sources)


def _vertex_count(polygons: Sequence[Polygon]) -> int:
    count = 0
    for polygon in polygons:
        rings = [polygon.exterior, *polygon.interiors]
        for ring in rings:
            count += max(len(ring.coords) - 1, 0)
    return count


def _coverage_quality_metrics(
    polygons: Sequence[Polygon],
    *,
    short_edge_threshold: float = 0.0,
) -> dict[str, float | int | None]:
    total_area = float(sum(polygon.area for polygon in polygons))
    union_area = float(unary_union(polygons).area) if polygons else 0.0
    segment_stats = _segment_length_stats(
        polygons,
        short_edge_threshold=short_edge_threshold,
    )
    pair_metrics = _coverage_pairwise_metrics(
        polygons,
        target_scale=short_edge_threshold,
    )
    polygon_clearance = _minimum_clearance(polygons)
    pair_clearance = pair_metrics["min_pair_clearance"]
    clearances = [
        value
        for value in (polygon_clearance, pair_clearance)
        if value is not None
    ]
    return {
        "polygon_count": len(polygons),
        "vertex_count": int(segment_stats["vertex_count"]),
        "total_area": total_area,
        "union_area": union_area,
        "overlap_area": float(max(total_area - union_area, 0.0)),
        "min_clearance": min(clearances) if clearances else None,
        "ring_contact_count": sum(
            _polygon_ring_boundary_contact_count(polygon)
            for polygon in polygons
        ),
        **pair_metrics,
        **segment_stats,
    }


def _record_stage_metrics(
    diagnostics: dict[str, Any],
    stage: str,
    polygons: Sequence[Polygon],
    *,
    short_edge_threshold: float = 0.0,
    extra: dict[str, Any] | None = None,
) -> dict[str, float | int | None] | None:
    if not diagnostics.get("collect_stage_metrics", True):
        return None
    metrics = _coverage_quality_metrics(
        polygons,
        short_edge_threshold=short_edge_threshold,
    )
    if extra:
        metrics = {**metrics, **extra}
    diagnostics["stage_metrics"][stage] = metrics
    return metrics


def _fmt_metric(value: float | int | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def _format_stage_metrics(
    stage: str,
    metrics: dict[str, Any],
    *,
    short_edge_threshold: float,
) -> str:
    line = (
        f"{stage}: polys={metrics['polygon_count']} verts={metrics['vertex_count']} "
        f"area={_fmt_metric(metrics['union_area'])} m2 "
        f"overlap={_fmt_metric(metrics['overlap_area'])} m2 "
        f"clear={_fmt_metric(metrics['min_clearance'])} m "
        f"edge_min={_fmt_metric(metrics['min_edge_length'])} m "
        f"edge_mean={_fmt_metric(metrics['mean_edge_length'])} m"
    )
    if short_edge_threshold > 0:
        line += (
            f" short<{short_edge_threshold:.2f}m={metrics['short_edge_count']}"
        )
        if "ring_contact_count" in metrics:
            line += f" ring_touch={metrics['ring_contact_count']}"
        if "pair_issue_count" in metrics:
            line += f" pair<{short_edge_threshold:.2f}m={metrics['pair_issue_count']}"
    if "group_count" in metrics:
        line += f" groups={metrics['group_count']}"
    return line


def _log_conditioning_start(
    input_count: int,
    options: ConditioningOptions,
    *,
    grid: float,
    output_grid: float,
    opening_radius: float,
    closing_radius: float,
) -> None:
    if not options.enable_logging:
        return
    info(
        "\n".join(
            [
                "Footprint cleaning started",
                (
                    f"  inputs={input_count} precision_grid={grid:.4f} m "
                    f"output_grid={output_grid:.4f} m"
                ),
                (
                    f"  opening_radius={opening_radius:.3f} m "
                    f"closing_radius={closing_radius:.3f} m"
                ),
                (
                    f"  min_feature_size={options.min_feature_size:.3f} m "
                    f"merge_distance={options.merge_distance:.3f} m "
                    f"min_area={options.min_area:.3f} m2 "
                    f"min_hole_area={options.min_hole_area:.3f} m2"
                ),
            ]
        )
    )


def _log_conditioning_stage(
    label: str,
    metrics: dict[str, Any],
    *,
    short_edge_threshold: float,
    enabled: bool,
) -> None:
    if not enabled:
        return
    info(
        "Footprint cleaning stage\n"
        f"  {_format_stage_metrics(label, metrics, short_edge_threshold=short_edge_threshold)}"
    )


def _log_conditioning_summary(
    diagnostics: dict[str, Any],
    *,
    short_edge_threshold: float,
) -> None:
    if not diagnostics.get("enable_logging", True):
        return
    stage_metrics = diagnostics["stage_metrics"]
    atomic = stage_metrics.get("atomic_input", {})
    final = stage_metrics.get("final_output", {})
    lines = ["Footprint cleaning summary"]
    if diagnostics.get("collect_stage_metrics", True):
        lines.extend(
            [
                (
                    f"  {_format_stage_metrics('atomic_input', atomic, short_edge_threshold=short_edge_threshold)}"
                    if atomic
                    else "  atomic_input: n/a"
                ),
                (
                    f"  {_format_stage_metrics('final_output', final, short_edge_threshold=short_edge_threshold)}"
                    if final
                    else "  final_output: n/a"
                ),
            ]
        )
    else:
        lines.append("  stage_metrics=disabled")

    lines.extend(
        [
        (
            f"  polygonize cut/dangle/invalid="
            f"{diagnostics['polygonize_cut_edges']}/"
            f"{diagnostics['polygonize_dangles']}/"
            f"{diagnostics['polygonize_invalid_rings']}"
        ),
        (
            f"  repaired_invalid={diagnostics['repaired_invalid_count']} "
            f"collapsed={diagnostics['collapsed_count']} "
            f"dropped_small={diagnostics['dropped_small_count']} "
            f"discarded_non_polygon={diagnostics['discarded_non_polygon_parts']}"
        ),
        (
            f"  multi_component_groups={diagnostics['multi_component_group_count']} "
            f"component_source_reassignments="
            f"{diagnostics['component_source_reassignment_count']}"
        ),
        (
            f"  overlap={_fmt_metric(diagnostics['overlap_area_before'])} -> "
            f"{_fmt_metric(diagnostics['overlap_area_after'])} m2 "
            f"clearance={_fmt_metric(diagnostics['min_clearance_before'])} -> "
            f"{_fmt_metric(diagnostics['min_clearance_after'])} m"
        ),
        (
            f"  output_canonicalization candidates="
            f"{diagnostics['output_canonicalization_candidate_count']} "
            f"output_grid_applied="
            f"{diagnostics['output_canonicalization_output_grid_applied_count']} "
            f"fallbacks={diagnostics['output_canonicalization_fallback_count']} "
            f"missing="
            f"{_fmt_metric(diagnostics['output_canonicalization_reference_minus_candidate_area'])} "
            f"extra="
            f"{_fmt_metric(diagnostics['output_canonicalization_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['output_canonicalization_signed_area_delta'])} "
            f"budget={_fmt_metric(diagnostics['output_canonicalization_area_balance_budget'])}"
        ),
        (
            f"  final_min_area_filter removed="
            f"{diagnostics['final_min_area_filter_removed_count']} "
            f"removed_area={_fmt_metric(diagnostics['final_min_area_filter_removed_area'])} m2"
        ),
        (
            f"  clearance_regularization candidates="
            f"{diagnostics['clearance_regularization_candidate_count']} "
            f"improved={diagnostics['clearance_regularization_improved_count']} "
            f"failed={diagnostics['clearance_regularization_failed_count']}"
        ),
        (
            f"  local_defect_repair candidates="
            f"{diagnostics['local_defect_repair_candidate_count']} "
            f"applied={diagnostics['local_defect_repair_applied_count']} "
            f"short_edges={diagnostics['local_defect_repair_short_edge_count_before']} -> "
            f"{diagnostics['local_defect_repair_short_edge_count_after']} "
            f"outside_zone={_fmt_metric(diagnostics['local_defect_repair_change_outside_edit_zone'])} m2 "
            f"missing={_fmt_metric(diagnostics['local_defect_repair_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['local_defect_repair_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['local_defect_repair_signed_area_delta'])} "
            f"budget={_fmt_metric(diagnostics['local_defect_repair_area_balance_budget'])} "
            f"rejected_nonlocal={diagnostics['local_defect_repair_rejected_nonlocal_count']} "
            f"rejected_area_bias={diagnostics['local_defect_repair_rejected_area_imbalance_count']}"
        ),
        (
            f"  local_defect_repair operators "
            f"attempted={diagnostics.get('local_defect_repair_operator_attempts', {})} "
            f"applied={diagnostics.get('local_defect_repair_operator_applied', {})}"
        ),
        (
            f"  coverage_simplify branch={diagnostics['coverage_simplify_selected_branch']} "
            f"patches={diagnostics['coverage_simplify_patch_count']} "
            f"applied={diagnostics['coverage_simplify_patch_applied_count']} short_edges="
            f"{diagnostics['coverage_simplify_short_edge_count_before']} -> "
            f"{diagnostics['coverage_simplify_short_edge_count_after']} "
            f"edit_zone={_fmt_metric(diagnostics['coverage_simplify_edit_zone_area'])} m2 "
            f"outside_zone={_fmt_metric(diagnostics['coverage_simplify_change_outside_edit_zone'])} m2 "
            f"missing={_fmt_metric(diagnostics['coverage_simplify_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['coverage_simplify_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['coverage_simplify_signed_area_delta'])} "
            f"budget={_fmt_metric(diagnostics['coverage_simplify_area_balance_budget'])} "
            f"rejected_nonlocal={diagnostics['coverage_simplify_rejected_nonlocal']} "
            f"rejected_area_bias={diagnostics['coverage_simplify_rejected_area_imbalance']}"
        ),
        (
            f"  coverage_simplify operators "
            f"attempted={diagnostics.get('coverage_simplify_operator_attempts', {})} "
            f"applied={diagnostics.get('coverage_simplify_operator_applied', {})}"
        ),
        (
            f"  coverage_simplify scores="
            f"{diagnostics.get('coverage_simplify_branch_scores', {})}"
        ),
        (
            f"  source_reclaim candidates={diagnostics['source_reclaim_candidate_count']} "
            f"components={diagnostics['source_reclaim_component_count']} "
            f"applied={diagnostics['source_reclaim_applied_count']} "
            f"short_edges={diagnostics['source_reclaim_short_edge_count_before']} -> "
            f"{diagnostics['source_reclaim_short_edge_count_after']} "
            f"outside_zone={_fmt_metric(diagnostics['source_reclaim_change_outside_edit_zone'])} m2 "
            f"missing={_fmt_metric(diagnostics['source_reclaim_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['source_reclaim_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['source_reclaim_signed_area_delta'])} "
            f"rejected_nonlocal={diagnostics['source_reclaim_rejected_nonlocal_count']} "
            f"rejected_non_improving={diagnostics['source_reclaim_rejected_non_improving_count']}"
        ),
        (
            f"  source_reclaim operators "
            f"attempted={diagnostics.get('source_reclaim_operator_attempts', {})} "
            f"applied={diagnostics.get('source_reclaim_operator_applied', {})}"
        ),
        (
            f"  small_component_absorb candidates="
            f"{diagnostics['small_component_absorb_candidate_count']} "
            f"components={diagnostics['small_component_absorb_component_count']} "
            f"applied={diagnostics['small_component_absorb_applied_count']} "
            f"outside_zone={_fmt_metric(diagnostics['small_component_absorb_change_outside_edit_zone'])} m2 "
            f"missing={_fmt_metric(diagnostics['small_component_absorb_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['small_component_absorb_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['small_component_absorb_signed_area_delta'])} "
            f"rejected_nonlocal={diagnostics['small_component_absorb_rejected_nonlocal_count']} "
            f"rejected_non_improving={diagnostics['small_component_absorb_rejected_non_improving_count']}"
        ),
        (
            f"  small_component_absorb operators "
            f"attempted={diagnostics.get('small_component_absorb_operator_attempts', {})} "
            f"applied={diagnostics.get('small_component_absorb_operator_applied', {})}"
        ),
        (
            f"  source_coordinate_recovery candidates="
            f"{diagnostics['source_coordinate_recovery_candidate_count']} "
            f"applied={diagnostics['source_coordinate_recovery_applied_count']} "
            f"exact={diagnostics['source_coordinate_recovery_exact_count']} "
            f"vertex={diagnostics['source_coordinate_recovery_vertex_count']} "
            f"short_edges={diagnostics['source_coordinate_recovery_short_edge_count_before']} -> "
            f"{diagnostics['source_coordinate_recovery_short_edge_count_after']} "
            f"missing={_fmt_metric(diagnostics['source_coordinate_recovery_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['source_coordinate_recovery_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['source_coordinate_recovery_signed_area_delta'])} "
            f"support_missing={_fmt_metric(diagnostics['source_coordinate_recovery_support_reference_minus_candidate_area'])} "
            f"support_extra={_fmt_metric(diagnostics['source_coordinate_recovery_support_candidate_minus_reference_area'])} "
            f"rejected_overlap={diagnostics['source_coordinate_recovery_rejected_overlap_count']} "
            f"rejected_non_improving={diagnostics['source_coordinate_recovery_rejected_non_improving_count']}"
        ),
        (
            f"  source_coordinate_recovery operators "
            f"attempted={diagnostics.get('source_coordinate_recovery_operator_attempts', {})} "
            f"applied={diagnostics.get('source_coordinate_recovery_operator_applied', {})}"
        ),
        (
            f"  coverage_meshing_regularization branch="
            f"{diagnostics['coverage_meshing_regularization_selected_branch']} "
            f"short_edges={diagnostics['coverage_meshing_regularization_short_edge_count_before']} -> "
            f"{diagnostics['coverage_meshing_regularization_short_edge_count_after']} "
            f"pair_issues={diagnostics['coverage_meshing_regularization_pair_issue_count_before']} -> "
            f"{diagnostics['coverage_meshing_regularization_pair_issue_count_after']} "
            f"missing={_fmt_metric(diagnostics['coverage_meshing_regularization_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['coverage_meshing_regularization_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['coverage_meshing_regularization_signed_area_delta'])}"
        ),
        (
            f"  coverage_meshing_regularization operators "
            f"attempted={diagnostics.get('coverage_meshing_regularization_operator_attempts', {})} "
            f"applied={diagnostics.get('coverage_meshing_regularization_operator_applied', {})}"
        ),
        (
            f"  boundary_regularization candidates="
            f"{diagnostics['polygon_simplify_candidate_count']} "
            f"applied={diagnostics['polygon_simplify_applied_count']} "
            f"short_edges={diagnostics['polygon_simplify_short_edge_count_before']} -> "
            f"{diagnostics['polygon_simplify_short_edge_count_after']} "
            f"outside_zone={_fmt_metric(diagnostics['polygon_simplify_change_outside_edit_zone'])} m2 "
            f"missing={_fmt_metric(diagnostics['polygon_simplify_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['polygon_simplify_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['polygon_simplify_signed_area_delta'])} "
            f"budget={_fmt_metric(diagnostics['polygon_simplify_area_balance_budget'])} "
            f"rejected_nonlocal={diagnostics['polygon_simplify_rejected_nonlocal_count']} "
            f"rejected_area_bias={diagnostics['polygon_simplify_rejected_area_imbalance_count']}"
        ),
        (
            f"  boundary_regularization operators "
            f"attempted={diagnostics.get('polygon_simplify_operator_attempts', {})} "
            f"applied={diagnostics.get('polygon_simplify_operator_applied', {})}"
        ),
        (
            f"  geos_exceptions={diagnostics['geos_exception_count']} "
            f"output_count={diagnostics['output_count']}"
        ),
        ]
    )
    if diagnostics["geos_exception_messages"]:
        lines.append(
            "  messages=" + " | ".join(diagnostics["geos_exception_messages"][:3])
        )
    info("\n".join(lines))


def _apply_coverage_candidate_diagnostics(
    diagnostics: dict[str, Any],
    candidate: _CoverageSimplifyCandidate | None,
) -> None:
    if candidate is None:
        diagnostics["coverage_simplify_applied"] = False
        diagnostics["coverage_simplify_patch_count"] = 0
        diagnostics["coverage_simplify_patch_applied_count"] = 0
        diagnostics["coverage_simplify_operator_attempts"] = {}
        diagnostics["coverage_simplify_operator_applied"] = {}
        return

    diagnostics["coverage_simplify_applied"] = candidate.label != "identity"
    diagnostics["coverage_simplify_short_edge_count_after"] = (
        candidate.signature.short_edge_count
    )
    diagnostics["coverage_simplify_reference_minus_candidate_area"] = (
        candidate.difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics["coverage_simplify_candidate_minus_reference_area"] = (
        candidate.difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics["coverage_simplify_signed_area_delta"] = candidate.difference_metrics[
        "union_area_delta"
    ]
    diagnostics["coverage_simplify_area_balance_budget"] = (
        candidate.area_balance_budget
    )
    diagnostics["coverage_simplify_change_outside_edit_zone"] = (
        candidate.change_outside_edit_zone
    )
    diagnostics["coverage_simplify_edit_zone_area"] = candidate.edit_zone_area
    diagnostics["coverage_simplify_patch_count"] = candidate.patch_count
    diagnostics["coverage_simplify_patch_applied_count"] = (
        candidate.patch_applied_count
    )
    diagnostics["coverage_simplify_operator_attempts"] = dict(
        candidate.operator_attempts
    )
    diagnostics["coverage_simplify_operator_applied"] = dict(
        candidate.operator_applied
    )


def _simplify_coverage_globally(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> _CoverageSimplifyCandidate | None:
    reference_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=tolerance,
    )
    reference_union = _cached_union(cache, polygons)
    edit_zone = _cached_coverage_edit_zone(
        cache,
        polygons,
        target_scale=tolerance,
        radius=_derive_local_edit_radius(tolerance),
    )
    operator_attempts = {"coverage_simplify_global": 1}

    try:
        simplified = shapely.coverage_simplify(
            np.asarray(polygons, dtype=object),
            tolerance,
            simplify_boundary=True,
        )
    except (AttributeError, GEOSException, TypeError, ValueError) as exc:
        diagnostics["coverage_simplify_failed"] = True
        _record_geos_exception(diagnostics, "coverage_simplify", exc)
        return None

    simplified_polygons: list[Polygon] = []
    simplified_sources: list[list[int]] = []
    for geometry, indices in zip(simplified, source_map):
        for polygon in _canonicalize(geometry, grid, diagnostics):
            simplified_polygons.append(orient(polygon, sign=1.0))
            simplified_sources.append(indices)
    if not simplified_polygons:
        return None
    if _cached_overlap_area(cache, simplified_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_signature = _cached_coverage_defect_signature(
        cache,
        simplified_polygons,
        target_scale=tolerance,
    )
    if not _coverage_signature_improves(
        reference_signature,
        candidate_signature,
        grid=grid,
        target_scale=tolerance,
    ):
        return None

    candidate_union = _cached_union(cache, simplified_polygons)
    difference_metrics = _cached_difference_area_metrics(
        cache,
        polygons,
        simplified_polygons,
    )
    area_balance_budget = _area_balance_budget(
        difference_metrics["symmetric_difference_area"],
        grid=grid,
        short_edge_count=reference_signature.short_edge_count,
        short_edge_threshold=tolerance,
    )
    change_outside_edit_zone = _change_outside_edit_zone(
        reference_union,
        candidate_union,
        edit_zone=edit_zone,
    )
    if change_outside_edit_zone > max(float(edit_zone.area), grid * grid, 1e-9):
        return None
    if abs(difference_metrics["union_area_delta"]) > area_balance_budget:
        return None

    return _CoverageSimplifyCandidate(
        label="global",
        polygons=simplified_polygons,
        source_map=simplified_sources,
        signature=candidate_signature,
        difference_metrics=difference_metrics,
        change_outside_edit_zone=change_outside_edit_zone,
        edit_zone_area=float(edit_zone.area),
        area_balance_budget=area_balance_budget,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts=operator_attempts,
        operator_applied={"coverage_simplify_global": 1},
    )


def _apply_open_close(
    geom: BaseGeometry,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> list[Polygon]:
    if radius <= 0:
        return _canonicalize(geom, grid, diagnostics)

    regularized = _run_constructive_step(
        "patch_open_close",
        geom,
        lambda value: value.buffer(
            -radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            -radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ),
        grid,
        diagnostics,
    )
    return _canonicalize(regularized, grid, diagnostics)


def _apply_close_open(
    geom: BaseGeometry,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> list[Polygon]:
    if radius <= 0:
        return _canonicalize(geom, grid, diagnostics)

    regularized = _run_constructive_step(
        "patch_close_open",
        geom,
        lambda value: value.buffer(
            radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            -radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            -radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ).buffer(
            radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        ),
        grid,
        diagnostics,
    )
    return _canonicalize(regularized, grid, diagnostics)


def _apply_local_patch_union_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    operator: str,
    grid: float,
    diagnostics: dict[str, Any],
    merge_components: bool = False,
    patch_geometry: BaseGeometry | None = None,
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or not subset_polygons:
        return None

    if patch_geometry is None:
        patch_geometry = unary_union(subset_polygons)
    if patch_geometry.is_empty:
        return None

    candidate_polygons: list[Polygon] = []
    patch_components = (
        [patch_geometry] if merge_components else _extract_polygon_parts(patch_geometry)
    )
    if not patch_components:
        return None

    for component in patch_components:
        if operator == "open_close":
            parts = _apply_open_close(component, radius, grid, diagnostics)
        elif operator == "close_open":
            parts = _apply_close_open(component, radius, grid, diagnostics)
        else:
            raise ValueError(f"Unsupported local patch operator: {operator}")
        for polygon in parts:
            candidate_polygons.append(orient(polygon, sign=1.0))

    if not candidate_polygons:
        return None
    if not merge_components and len(candidate_polygons) < len(subset_polygons):
        return None
    if _coverage_overlap_area(candidate_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_sources = _assign_component_sources(
        candidate_polygons,
        subset_polygons,
        subset_sources,
        grid=grid,
        diagnostics=diagnostics,
    )
    return _stable_sort(candidate_polygons, candidate_sources)


def _apply_local_pair_merge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
    patch_geometry: BaseGeometry | None = None,
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or not subset_polygons:
        return None

    if patch_geometry is None:
        patch_geometry = unary_union(subset_polygons)
    if patch_geometry.is_empty:
        return None

    candidate_polygons = [
        orient(polygon, sign=1.0)
        for polygon in _apply_closing(
            patch_geometry,
            radius,
            grid,
            diagnostics,
        )
    ]
    if not candidate_polygons:
        return None
    if len(candidate_polygons) > len(subset_polygons):
        return None
    if _coverage_overlap_area(candidate_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_sources = _assign_component_sources(
        candidate_polygons,
        subset_polygons,
        subset_sources,
        grid=grid,
        diagnostics=diagnostics,
    )
    return _stable_sort(candidate_polygons, candidate_sources)


def _apply_local_point_touch_bridge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or len(subset_polygons) < 2:
        return None

    touch_points: list[BaseGeometry] = []
    point_tolerance = max(grid, radius) * 1e-3
    tree = STRtree(subset_polygons)
    for index, polygon in enumerate(subset_polygons):
        for candidate in tree.query(polygon):
            other_index = int(candidate)
            if other_index <= index:
                continue
            other = subset_polygons[other_index]
            try:
                if polygon.distance(other) > point_tolerance:
                    continue
                boundary_intersection = polygon.boundary.intersection(other.boundary)
            except GEOSException:
                continue
            if boundary_intersection.is_empty or boundary_intersection.length > point_tolerance:
                continue
            touch_points.append(boundary_intersection)

    if not touch_points:
        return None

    bridge = unary_union(touch_points).buffer(max(radius, grid), quad_segs=1)
    patch_geometry = unary_union([*subset_polygons, bridge])
    candidate_polygons = [
        orient(polygon, sign=1.0)
        for polygon in _canonicalize(patch_geometry, grid, diagnostics)
    ]
    if not candidate_polygons:
        return None
    if len(candidate_polygons) > len(subset_polygons):
        return None
    if _coverage_overlap_area(candidate_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_sources = _assign_component_sources(
        candidate_polygons,
        subset_polygons,
        subset_sources,
        grid=grid,
        diagnostics=diagnostics,
    )
    return _stable_sort(candidate_polygons, candidate_sources)


def _apply_local_close_pair_bridge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or len(subset_polygons) < 2:
        return None

    connectors: list[BaseGeometry] = []
    tree = STRtree(subset_polygons)
    point_tolerance = max(grid, radius) * 1e-3
    for index, polygon in enumerate(subset_polygons):
        for candidate in tree.query(_expanded_query_geometry(polygon, radius=radius)):
            other_index = int(candidate)
            if other_index <= index:
                continue
            other = subset_polygons[other_index]
            try:
                distance = float(polygon.distance(other))
            except GEOSException:
                continue
            if distance <= point_tolerance or distance > radius + point_tolerance:
                continue
            try:
                connectors.append(shapely.shortest_line(polygon, other))
            except GEOSException:
                continue

    if not connectors:
        return None

    bridge = unary_union(connectors).buffer(max(radius, grid), quad_segs=1)
    patch_geometry = unary_union([*subset_polygons, bridge])
    candidate_polygons = [
        orient(polygon, sign=1.0)
        for polygon in _canonicalize(patch_geometry, grid, diagnostics)
    ]
    if not candidate_polygons:
        return None
    if len(candidate_polygons) > len(subset_polygons):
        return None
    if _coverage_overlap_area(candidate_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_sources = _assign_component_sources(
        candidate_polygons,
        subset_polygons,
        subset_sources,
        grid=grid,
        diagnostics=diagnostics,
    )
    return _stable_sort(candidate_polygons, candidate_sources)


def _iter_pair_issue_shrink_radii(
    *,
    tolerance: float,
    grid: float,
) -> tuple[float, ...]:
    if tolerance <= 0:
        return ()

    radii: list[float] = []
    for factor in (0.5, 0.625, 0.75):
        radius = max(grid, tolerance * factor)
        if any(abs(radius - existing) <= 1e-12 for existing in radii):
            continue
        radii.append(radius)
    return tuple(radii)


def _apply_local_smaller_polygon_shrink_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or len(subset_polygons) < 2:
        return None

    areas = [float(polygon.area) for polygon in subset_polygons]
    if not areas:
        return None
    small_index = min(
        range(len(subset_polygons)),
        key=lambda index: (areas[index], index),
    )

    try:
        shrunken = subset_polygons[small_index].buffer(
            -radius,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        )
    except GEOSException:
        return None

    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []
    for index, (polygon, indices) in enumerate(zip(subset_polygons, subset_sources)):
        geometry = shrunken if index == small_index else polygon
        parts = [
            orient(part, sign=1.0)
            for part in _canonicalize(geometry, grid, diagnostics)
        ]
        if len(parts) != 1:
            return None
        candidate_polygons.append(parts[0])
        candidate_sources.append(list(indices))

    if _coverage_overlap_area(candidate_polygons) > max(grid * grid, 1e-9):
        return None

    return _stable_sort(candidate_polygons, candidate_sources)


def _courtyard_pair_shrink_budget_context(
    polygons: Sequence[Polygon],
    *,
    tolerance: float,
    grid: float,
) -> tuple[BaseGeometry, float] | None:
    if len(polygons) != 2 or tolerance <= 0:
        return None

    larger_index = max(range(len(polygons)), key=lambda index: float(polygons[index].area))
    smaller_index = 1 - larger_index
    larger = polygons[larger_index]
    smaller = polygons[smaller_index]
    if not larger.interiors:
        return None

    representative = smaller.representative_point()
    if not any(Polygon(ring).covers(representative) for ring in larger.interiors):
        return None

    edit_zone = smaller.buffer(
        max(tolerance, grid),
        quad_segs=1,
        join_style=BufferJoinStyle.mitre,
        mitre_limit=1000.0,
    )
    shrink_area_budget = max(
        0.2 * float(smaller.area),
        tolerance * tolerance,
        16.0 * grid * grid,
    )
    return edit_zone, shrink_area_budget


def _direct_pair_issue_cluster_candidates(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> list[tuple[str, list[Polygon], list[list[int]]]]:
    candidates: list[tuple[str, list[Polygon], list[list[int]]]] = []
    if tolerance <= 0 or len(subset_polygons) < 2:
        return candidates

    seen_signatures: set[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]] = set()
    for left_index, right_index, distance, issue_kind in _cached_pair_issue_candidates(
        cache,
        subset_polygons,
        target_scale=tolerance,
    ):
        pair_polygons = [subset_polygons[left_index], subset_polygons[right_index]]
        pair_sources = [subset_sources[left_index], subset_sources[right_index]]
        radius = max(grid, min(tolerance * 0.5, max(distance, grid) * 2.0))
        operator_candidates: list[
            tuple[str, tuple[list[Polygon], list[list[int]]] | None]
        ] = []
        if issue_kind == "point":
            operator_candidates.append(
                (
                    f"coverage_pair_issue_point_{radius:.3f}",
                    _apply_local_point_touch_bridge_operator(
                        pair_polygons,
                        pair_sources,
                        radius=radius,
                        grid=grid,
                        diagnostics=diagnostics,
                    ),
                )
            )
        else:
            operator_candidates.append(
                (
                    f"coverage_pair_issue_bridge_{radius:.3f}",
                    _apply_local_close_pair_bridge_operator(
                        pair_polygons,
                        pair_sources,
                        radius=radius,
                        grid=grid,
                        diagnostics=diagnostics,
                    ),
                )
            )

        for shrink_radius in _iter_pair_issue_shrink_radii(
            tolerance=tolerance,
            grid=grid,
        ):
            operator_candidates.append(
                (
                    f"coverage_pair_issue_shrink_{shrink_radius:.3f}",
                    _apply_local_smaller_polygon_shrink_operator(
                        pair_polygons,
                        pair_sources,
                        radius=shrink_radius,
                        grid=grid,
                        diagnostics=diagnostics,
                    ),
                )
            )

        for operator, candidate in operator_candidates:
            if candidate is None:
                continue

            repaired_pair_polygons, repaired_pair_sources = candidate
            candidate_polygons = [
                polygon
                for index, polygon in enumerate(subset_polygons)
                if index not in (left_index, right_index)
            ] + repaired_pair_polygons
            candidate_sources = [
                indices
                for index, indices in enumerate(subset_sources)
                if index not in (left_index, right_index)
            ] + repaired_pair_sources
            candidate_polygons, candidate_sources = _stable_sort(
                candidate_polygons,
                candidate_sources,
            )
            signature = (
                _polygon_sequence_key(candidate_polygons),
                tuple(tuple(indices) for indices in candidate_sources),
            )
            if signature in seen_signatures:
                continue
            seen_signatures.add(signature)
            candidates.append((operator, candidate_polygons, candidate_sources))

    return candidates


def _iter_nonpair_patch_fallback_configs(
    signature: _CoverageDefectSignature,
    *,
    subset_size: int,
    patch_radii: Sequence[float],
) -> tuple[tuple[str, float], ...]:
    if not patch_radii or signature.pair_issue_count > 0:
        return ()

    # Light residual short-edge clusters rarely benefit from the full
    # morphology ladder. Keep one deterministic rescue for the small
    # two-polygon case and otherwise skip straight to identity.
    if signature.short_edge_count <= 1:
        return ()
    if signature.short_edge_count <= 2:
        if subset_size <= 2:
            return (("open_close", patch_radii[0]),)
        return ()
    if signature.short_edge_count <= 4 and subset_size <= 3:
        return (
            ("open_close", patch_radii[0]),
            ("close_open", patch_radii[0]),
        )

    configs: list[tuple[str, float]] = []
    for radius in patch_radii:
        for operator_name in ("open_close", "close_open"):
            configs.append((operator_name, radius))
    return tuple(configs)


def _iter_local_coverage_simplify_tolerances(
    signature: _CoverageDefectSignature,
    *,
    patch_tolerances: Sequence[float],
) -> tuple[float, ...]:
    if not patch_tolerances:
        return ()

    # Single-short-edge residuals without pair defects consistently favor the
    # coarser simplify pass. Skip the finer trial and avoid paying for a
    # second canonicalize/score cycle on the same cluster.
    if (
        signature.pair_issue_count == 0
        and signature.point_touch_count == 0
        and signature.close_pair_count == 0
        and signature.short_edge_count <= 1
    ):
        return (patch_tolerances[-1],)

    return tuple(patch_tolerances)


def _simplify_polygon_short_edge_graphically(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[Polygon | None, int]:
    shell, shell_removed = _clean_ring_short_edge_chains(
        polygon.exterior.coords,
        target_scale=target_scale,
        grid=grid,
    )
    if shell is None:
        return None, 0

    holes: list[list[tuple[float, float]]] = []
    removed_count = shell_removed
    for ring in polygon.interiors:
        cleaned_hole, hole_removed = _clean_ring_short_edge_chains(
            ring.coords,
            target_scale=target_scale,
            grid=grid,
        )
        removed_count += hole_removed
        if cleaned_hole is None:
            continue
        holes.append(cleaned_hole)

    if removed_count == 0:
        return None, 0

    candidate = _normalize_single_polygon_candidate(
        Polygon(shell, holes),
        grid=grid,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )
    if candidate is None:
        return None, 0
    return candidate, removed_count


def _simplify_coverage_short_edge_graphically(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> _CoverageSimplifyCandidate | None:
    if tolerance <= 0 or not polygons:
        return None

    reference_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=tolerance,
    )
    if reference_signature.short_edge_count == 0:
        return None

    reference_union = _cached_union(cache, polygons)
    edit_zone = _cached_coverage_edit_zone(
        cache,
        polygons,
        target_scale=tolerance,
        radius=_derive_local_edit_radius(tolerance),
    )
    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []
    changed_count = 0

    for polygon, indices in zip(polygons, source_map):
        candidate, removed_count = _simplify_polygon_short_edge_graphically(
            polygon,
            target_scale=tolerance,
            grid=grid,
            diagnostics=diagnostics,
        )
        if candidate is None or removed_count == 0:
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue
        candidate_polygons.append(candidate)
        candidate_sources.append(indices)
        changed_count += 1

    if changed_count == 0:
        return None
    if _cached_overlap_area(cache, candidate_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_signature = _cached_coverage_defect_signature(
        cache,
        candidate_polygons,
        target_scale=tolerance,
    )
    if not _coverage_signature_improves(
        reference_signature,
        candidate_signature,
        grid=grid,
        target_scale=tolerance,
    ):
        return None

    candidate_union = _cached_union(cache, candidate_polygons)
    difference_metrics = _cached_difference_area_metrics(
        cache,
        polygons,
        candidate_polygons,
    )
    area_balance_budget = _area_balance_budget(
        difference_metrics["symmetric_difference_area"],
        grid=grid,
        short_edge_count=reference_signature.short_edge_count,
        short_edge_threshold=tolerance,
    )
    change_outside_edit_zone = _change_outside_edit_zone(
        reference_union,
        candidate_union,
        edit_zone=edit_zone,
    )
    if change_outside_edit_zone > max(float(edit_zone.area), grid * grid, 1e-9):
        return None
    if abs(difference_metrics["union_area_delta"]) > area_balance_budget:
        return None

    operator = "coverage_graph_short_edges"
    return _CoverageSimplifyCandidate(
        label="local",
        polygons=candidate_polygons,
        source_map=candidate_sources,
        signature=candidate_signature,
        difference_metrics=difference_metrics,
        change_outside_edit_zone=change_outside_edit_zone,
        edit_zone_area=float(edit_zone.area),
        area_balance_budget=area_balance_budget,
        patch_count=changed_count,
        patch_applied_count=changed_count,
        operator_attempts={operator: changed_count},
        operator_applied={operator: changed_count},
    )


def _simplify_coverage_residual_scale_polygons(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    tolerance: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> _CoverageSimplifyCandidate | None:
    if tolerance <= 0 or not polygons:
        return None

    reference_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=tolerance,
    )
    if _coverage_signature_satisfies_scale_contract(
        reference_signature,
        target_scale=tolerance,
        grid=grid,
    ):
        return None

    candidate_polygons = list(polygons)
    candidate_sources = [list(indices) for indices in source_map]
    changed_count = 0
    changed_area = 0.0

    for index, polygon in enumerate(polygons):
        reference_polygon_signature = _polygon_defect_signature(
            polygon,
            target_scale=tolerance,
        )
        if _signature_satisfies_scale_contract(
            reference_polygon_signature,
            target_scale=tolerance,
            grid=grid,
        ):
            continue

        best_candidate: Polygon | None = None
        best_metrics: dict[str, float] | None = None
        best_score: tuple[float, float, float, float] | None = None
        for factor in (0.75, 1.0, 1.25):
            simplify_tolerance = max(factor * tolerance, grid, 1e-9)
            try:
                simplified = shapely.simplify(
                    polygon,
                    simplify_tolerance,
                    preserve_topology=True,
                )
            except (GEOSException, ValueError, TypeError) as exc:
                _record_geos_exception(diagnostics, "residual_polygon_simplify", exc)
                continue

            normalized = _normalize_single_polygon_candidate(
                simplified,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
            )
            if normalized is None:
                continue

            candidate_signature = _polygon_defect_signature(
                normalized,
                target_scale=tolerance,
            )
            if not _signature_satisfies_scale_contract(
                candidate_signature,
                target_scale=tolerance,
                grid=grid,
            ):
                continue

            difference_metrics = _difference_area_metrics(polygon, normalized)
            if not _difference_metrics_have_small_fidelity_drift(
                difference_metrics,
                edit_zone_area=float(polygon.area),
                target_scale=tolerance,
                grid=grid,
            ):
                continue

            score = (
                difference_metrics["symmetric_difference_area"],
                abs(difference_metrics["union_area_delta"]),
                difference_metrics["reference_minus_candidate_area"],
                difference_metrics["candidate_minus_reference_area"],
            )
            if best_score is None or score < best_score:
                best_candidate = normalized
                best_metrics = difference_metrics
                best_score = score

        if best_candidate is None or best_metrics is None:
            continue

        candidate_polygons[index] = best_candidate
        changed_count += 1
        changed_area += float(polygon.area)

    if changed_count == 0:
        return None
    if _cached_overlap_area(cache, candidate_polygons) > max(grid * grid, 1e-9):
        return None

    candidate_signature = _cached_coverage_defect_signature(
        cache,
        candidate_polygons,
        target_scale=tolerance,
    )
    if not _coverage_signature_improves(
        reference_signature,
        candidate_signature,
        grid=grid,
        target_scale=tolerance,
    ):
        return None
    if not _coverage_signature_satisfies_scale_contract(
        candidate_signature,
        target_scale=tolerance,
        grid=grid,
    ):
        return None

    difference_metrics = _cached_difference_area_metrics(
        cache,
        polygons,
        candidate_polygons,
    )
    if not _difference_metrics_have_small_fidelity_drift(
        difference_metrics,
        edit_zone_area=changed_area,
        target_scale=tolerance,
        grid=grid,
    ):
        return None

    operator = "coverage_residual_polygon_simplify"
    return _CoverageSimplifyCandidate(
        label="residual_polygon_simplify",
        polygons=candidate_polygons,
        source_map=candidate_sources,
        signature=candidate_signature,
        difference_metrics=difference_metrics,
        change_outside_edit_zone=0.0,
        edit_zone_area=changed_area,
        area_balance_budget=max(changed_area, grid * grid, 1e-9),
        patch_count=changed_count,
        patch_applied_count=changed_count,
        operator_attempts={operator: changed_count},
        operator_applied={operator: changed_count},
    )


def _rewrite_coverage_cluster_graphically(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    cluster: dict[str, Any],
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> tuple[str, list[Polygon], list[list[int]]] | None:
    rewrite = _builder_rewrite_defect_cluster(
        subset_polygons,
        cluster,
        target_scale=tolerance,
        grid=grid,
        cache=cache,
    )
    if rewrite is None:
        return None

    changed_indices, rewritten_polygons = rewrite
    if not changed_indices or len(changed_indices) != len(rewritten_polygons):
        return None

    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []
    changed_index_set = set(changed_indices)
    rewritten_iter = iter(rewritten_polygons)
    for index, (polygon, indices) in enumerate(zip(subset_polygons, subset_sources)):
        if index in changed_index_set:
            rewritten = next(rewritten_iter)
            normalized = _normalize_single_polygon_candidate(
                rewritten,
                grid=grid,
                min_area=0.0,
                min_hole_area=0.0,
                diagnostics=diagnostics,
            )
            if normalized is None:
                return None
            candidate_polygons.append(normalized)
            candidate_sources.append(list(indices))
        else:
            candidate_polygons.append(polygon)
            candidate_sources.append(list(indices))

    return "coverage_cpp_cluster_rewrite", candidate_polygons, candidate_sources


def _pair_issue_candidates(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
) -> list[tuple[int, int, float, Literal["point", "close"]]]:
    if target_scale <= 0 or len(polygons) < 2:
        return []

    point_tolerance = max(target_scale * 1e-6, 1e-9)
    line_tolerance = point_tolerance
    tree = STRtree(polygons)
    candidates: list[tuple[int, int, float, Literal["point", "close"]]] = []
    for index, polygon in enumerate(polygons):
        for candidate in tree.query(
            _expanded_query_geometry(polygon, radius=target_scale)
        ):
            other_index = int(candidate)
            if other_index <= index:
                continue
            other = polygons[other_index]
            try:
                distance = float(polygon.distance(other))
            except GEOSException:
                continue
            if distance > target_scale + point_tolerance:
                continue
            if distance <= point_tolerance:
                try:
                    boundary_intersection = polygon.boundary.intersection(other.boundary)
                except GEOSException:
                    boundary_intersection = GeometryCollection()
                if boundary_intersection.length > line_tolerance:
                    continue
                candidates.append((index, other_index, distance, "point"))
                continue
            candidates.append((index, other_index, distance, "close"))
    candidates.sort(key=lambda item: (item[2], item[0], item[1]))
    return candidates


def _polygon_has_short_edges(
    polygon: Polygon,
    *,
    short_edge_threshold: float,
) -> bool:
    if short_edge_threshold <= 0:
        return False

    for ring in [polygon.exterior, *polygon.interiors]:
        coords = np.asarray(ring.coords, dtype=float)
        if len(coords) < 2:
            continue
        deltas = np.diff(coords[:, :2], axis=0)
        lengths = np.hypot(deltas[:, 0], deltas[:, 1])
        if lengths.size == 0:
            continue
        if np.any(lengths + 1e-12 < short_edge_threshold):
            return True
    return False


def _merge_index_clusters(clusters: list[set[int]]) -> list[list[int]]:
    if not clusters:
        return []

    merged: list[set[int]] = []
    for cluster in clusters:
        if not cluster:
            continue
        overlaps: list[int] = []
        for idx, existing in enumerate(merged):
            if existing & cluster:
                overlaps.append(idx)
        if not overlaps:
            merged.append(set(cluster))
            continue
        combined = set(cluster)
        for idx in reversed(overlaps):
            combined.update(merged.pop(idx))
        merged.append(combined)

    merged.sort(key=lambda values: (min(values), len(values), tuple(sorted(values))))
    return [sorted(values) for values in merged]


def _coverage_defect_clusters(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    cluster_radius: float,
    cache: _CoverageEvalCache | None = None,
) -> list[list[int]]:
    if target_scale <= 0 or cluster_radius <= 0 or not polygons:
        return []

    defect_indices: set[int] = {
        index
        for index, polygon in enumerate(polygons)
        if _polygon_has_short_edges(
            polygon,
            short_edge_threshold=target_scale,
        )
    }

    adjacency: dict[int, set[int]] = {index: set() for index in defect_indices}
    for left_index, right_index, _, _ in _cached_pair_issue_candidates(
        cache,
        polygons,
        target_scale=target_scale,
    ):
        defect_indices.add(left_index)
        defect_indices.add(right_index)
        adjacency.setdefault(left_index, set()).add(right_index)
        adjacency.setdefault(right_index, set()).add(left_index)

    if not defect_indices:
        return []

    search_radius = max(cluster_radius, target_scale, 1e-9)
    tree = STRtree(polygons)
    tolerance = max(target_scale * 1e-6, 1e-9)
    for index in list(defect_indices):
        polygon = polygons[index]
        query_geometry = _expanded_query_geometry(polygon, radius=search_radius)
        neighbors = adjacency.setdefault(index, set())
        for candidate in tree.query(query_geometry):
            other_index = int(candidate)
            if other_index == index or other_index not in defect_indices:
                continue
            other = polygons[other_index]
            try:
                if polygon.distance(other) > search_radius + tolerance:
                    continue
            except GEOSException:
                continue
            neighbors.add(other_index)
            adjacency.setdefault(other_index, set()).add(index)

    connected_components: list[set[int]] = []
    remaining = set(defect_indices)
    while remaining:
        seed = min(remaining)
        stack = [seed]
        component: set[int] = set()
        while stack:
            index = stack.pop()
            if index in component:
                continue
            component.add(index)
            remaining.discard(index)
            stack.extend(adjacency.get(index, ()))
        connected_components.append(component)

    expanded_components: list[set[int]] = []
    expansion_radius = max(cluster_radius, 1e-9)
    for component in connected_components:
        expanded = set(component)
        for index in component:
            query_geometry = _expanded_query_geometry(
                polygons[index],
                radius=expansion_radius,
            )
            for candidate in tree.query(query_geometry):
                expanded.add(int(candidate))
        expanded_components.append(expanded)

    return _merge_index_clusters(expanded_components)


def _simplify_coverage_pair_issues_graphically(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> _CoverageSimplifyCandidate | None:
    if tolerance <= 0 or len(polygons) < 2:
        return None

    current_polygons = list(polygons)
    current_sources = [list(indices) for indices in source_map]
    reference_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=tolerance,
    )
    if reference_signature.pair_issue_count == 0:
        return None

    reference_union = _cached_union(cache, polygons)
    local_radius = _derive_local_edit_radius(tolerance)
    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}
    patch_count = 0
    patch_applied_count = 0

    while True:
        accepted = False
        for (
            left_index,
            right_index,
            distance,
            issue_kind,
        ) in _cached_pair_issue_candidates(
            cache,
            current_polygons,
            target_scale=tolerance,
        ):
            subset_polygons = [current_polygons[left_index], current_polygons[right_index]]
            subset_sources = [current_sources[left_index], current_sources[right_index]]
            reference_subset_signature = _cached_coverage_defect_signature(
                cache,
                subset_polygons,
                target_scale=tolerance,
            )
            if reference_subset_signature.pair_issue_count == 0:
                continue

            patch_count += 1
            radius = max(grid, min(tolerance * 0.5, max(distance, grid) * 2.0))
            shrink_budget_context = _courtyard_pair_shrink_budget_context(
                subset_polygons,
                tolerance=tolerance,
                grid=grid,
            )
            operator_candidates: list[
                tuple[str, tuple[list[Polygon], list[list[int]]] | None]
            ] = []
            if issue_kind == "point":
                operator_candidates.append(
                    (
                        f"coverage_pair_issue_point_{radius:.3f}",
                        _apply_local_point_touch_bridge_operator(
                            subset_polygons,
                            subset_sources,
                            radius=radius,
                            grid=grid,
                            diagnostics=diagnostics,
                        ),
                    )
                )
            else:
                operator_candidates.append(
                    (
                        f"coverage_pair_issue_bridge_{radius:.3f}",
                        _apply_local_close_pair_bridge_operator(
                            subset_polygons,
                            subset_sources,
                            radius=radius,
                            grid=grid,
                            diagnostics=diagnostics,
                        ),
                    )
                )

            for shrink_radius in _iter_pair_issue_shrink_radii(
                tolerance=tolerance,
                grid=grid,
            ):
                operator_candidates.append(
                    (
                        f"coverage_pair_issue_shrink_{shrink_radius:.3f}",
                        _apply_local_smaller_polygon_shrink_operator(
                            subset_polygons,
                            subset_sources,
                            radius=shrink_radius,
                            grid=grid,
                            diagnostics=diagnostics,
                        ),
                    )
                )

            for operator, candidate in operator_candidates:
                operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                if candidate is None:
                    continue

                replacement_polygons, replacement_sources = candidate
                candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    replacement_polygons,
                    target_scale=tolerance,
                )
                if not _coverage_signature_improves(
                    reference_subset_signature,
                    candidate_signature,
                    grid=grid,
                    target_scale=tolerance,
                ):
                    continue

                reference_subset_union = _cached_union(cache, subset_polygons)
                candidate_union = _cached_union(cache, replacement_polygons)
                edit_zone = _cached_coverage_edit_zone(
                    cache,
                    subset_polygons,
                    target_scale=tolerance,
                    radius=local_radius,
                )
                difference_metrics = _cached_difference_area_metrics(
                    cache,
                    subset_polygons,
                    replacement_polygons,
                )
                area_balance_budget = _area_balance_budget(
                    difference_metrics["symmetric_difference_area"],
                    grid=grid,
                    short_edge_count=reference_subset_signature.short_edge_count,
                    short_edge_threshold=tolerance,
                )
                if (
                    operator.startswith("coverage_pair_issue_shrink_")
                    and shrink_budget_context is not None
                ):
                    edit_zone, shrink_area_budget = shrink_budget_context
                    area_balance_budget = max(area_balance_budget, shrink_area_budget)
                else:
                    edit_zone = _cached_coverage_edit_zone(
                        cache,
                        subset_polygons,
                        target_scale=tolerance,
                        radius=local_radius,
                    )
                change_outside_edit_zone = _change_outside_edit_zone(
                    reference_subset_union,
                    candidate_union,
                    edit_zone=edit_zone,
                )
                if change_outside_edit_zone > max(float(edit_zone.area), grid * grid, 1e-9):
                    continue
                extra_area_budget = max(
                    area_balance_budget,
                    tolerance * tolerance,
                    0.1 * float(edit_zone.area),
                    16.0 * grid * grid,
                )
                if (
                    difference_metrics["reference_minus_candidate_area"]
                    > area_balance_budget
                ):
                    continue
                if (
                    difference_metrics["candidate_minus_reference_area"]
                    > extra_area_budget
                ):
                    continue

                unaffected_polygons = [
                    polygon
                    for index, polygon in enumerate(current_polygons)
                    if index not in (left_index, right_index)
                ]
                unaffected_sources = [
                    indices
                    for index, indices in enumerate(current_sources)
                    if index not in (left_index, right_index)
                ]
                current_polygons, current_sources = _stable_sort(
                    unaffected_polygons + replacement_polygons,
                    unaffected_sources + replacement_sources,
                )
                operator_applied[operator] = operator_applied.get(operator, 0) + 1
                patch_applied_count += 1
                accepted = True
                break
            if accepted:
                break

        if not accepted:
            break

    if patch_applied_count == 0:
        return None

    candidate_signature = _cached_coverage_defect_signature(
        cache,
        current_polygons,
        target_scale=tolerance,
    )
    if not _coverage_signature_improves(
        reference_signature,
        candidate_signature,
        grid=grid,
        target_scale=tolerance,
    ):
        return None

    edit_zone = _cached_coverage_edit_zone(
        cache,
        polygons,
        target_scale=tolerance,
        radius=local_radius,
    )
    candidate_union = _cached_union(cache, current_polygons)
    difference_metrics = _cached_difference_area_metrics(
        cache,
        polygons,
        current_polygons,
    )
    area_balance_budget = _area_balance_budget(
        difference_metrics["symmetric_difference_area"],
        grid=grid,
        short_edge_count=reference_signature.short_edge_count,
        short_edge_threshold=tolerance,
    )
    change_outside_edit_zone = _change_outside_edit_zone(
        reference_union,
        candidate_union,
        edit_zone=edit_zone,
    )
    if change_outside_edit_zone > max(float(edit_zone.area), grid * grid, 1e-9):
        return None
    extra_area_budget = max(
        area_balance_budget,
        tolerance * tolerance,
        0.1 * float(edit_zone.area),
        16.0 * grid * grid,
    )
    if difference_metrics["reference_minus_candidate_area"] > area_balance_budget:
        return None
    if difference_metrics["candidate_minus_reference_area"] > extra_area_budget:
        return None

    return _CoverageSimplifyCandidate(
        label="local",
        polygons=current_polygons,
        source_map=current_sources,
        signature=candidate_signature,
        difference_metrics=difference_metrics,
        change_outside_edit_zone=change_outside_edit_zone,
        edit_zone_area=float(edit_zone.area),
        area_balance_budget=area_balance_budget,
        patch_count=patch_count,
        patch_applied_count=patch_applied_count,
        operator_attempts=operator_attempts,
        operator_applied=operator_applied,
    )


def _simplify_coverage_locally(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
    enable_patch_union_fallback: bool = True,
    enable_pair_cluster_rescue: bool = True,
    cache: _CoverageEvalCache | None = None,
) -> _CoverageSimplifyCandidate | None:
    if cache is None:
        cache = _CoverageEvalCache()
    current_polygons = list(polygons)
    current_sources = [list(indices) for indices in source_map]
    reference_union = _cached_union(cache, polygons)
    reference_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=tolerance,
    )
    local_radius = _derive_local_edit_radius(tolerance)
    patch_radius = _derive_patch_neighborhood_radius(tolerance)
    ordered_patch_tolerances: list[float] = []
    for factor in (1.0, 1.25):
        candidate_tolerance = max(tolerance * factor, grid)
        if candidate_tolerance <= 0:
            continue
        if candidate_tolerance in ordered_patch_tolerances:
            continue
        ordered_patch_tolerances.append(candidate_tolerance)
    patch_tolerances = tuple(ordered_patch_tolerances)
    patch_radii = tuple(
        sorted(
            {
                max(tolerance * factor, grid)
                for factor in (0.25, 0.375, 0.5)
                if tolerance * factor > 0
            }
        )
    )
    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}
    patch_count = 0
    patch_applied_count = 0

    pair_candidate = _simplify_coverage_pair_issues_graphically(
        polygons,
        source_map,
        tolerance=tolerance,
        grid=grid,
        diagnostics=diagnostics,
        cache=cache,
    )
    if pair_candidate is not None:
        current_polygons = list(pair_candidate.polygons)
        current_sources = [list(indices) for indices in pair_candidate.source_map]
        operator_attempts.update(pair_candidate.operator_attempts)
        operator_applied.update(pair_candidate.operator_applied)
        patch_count += pair_candidate.patch_count
        patch_applied_count += pair_candidate.patch_applied_count

    fast_candidate = _simplify_coverage_short_edge_graphically(
        current_polygons,
        current_sources,
        tolerance=tolerance,
        grid=grid,
        diagnostics=diagnostics,
        cache=cache,
    )
    if fast_candidate is not None:
        current_polygons = list(fast_candidate.polygons)
        current_sources = [list(indices) for indices in fast_candidate.source_map]
        operator_attempts.update(fast_candidate.operator_attempts)
        operator_applied.update(fast_candidate.operator_applied)
        patch_count += fast_candidate.patch_count
        patch_applied_count += fast_candidate.patch_applied_count

        fast_residual_short_edge_budget = max(
            48,
            reference_signature.short_edge_count // 3,
        )
        if (
            fast_candidate.signature.pair_issue_count == 0
            and fast_candidate.signature.short_edge_count
            <= fast_residual_short_edge_budget
        ):
            return fast_candidate

    while True:
        defect_clusters = _cached_coverage_defect_cluster_descriptors(
            cache,
            current_polygons,
            target_scale=tolerance,
            cluster_radius=max(local_radius + patch_radius, grid),
        )
        diagnostics["coverage_simplify_defect_cluster_count"] = len(defect_clusters)
        if not defect_clusters:
            break

        accepted_patch = False
        for cluster in defect_clusters:
            affected_indices = cluster["indices"]
            subset_polygons = [current_polygons[index] for index in affected_indices]
            subset_sources = [current_sources[index] for index in affected_indices]
            reference_subset_signature = _cached_coverage_defect_signature(
                cache,
                subset_polygons,
                target_scale=tolerance,
            )
            if (
                reference_subset_signature.short_edge_count == 0
                and reference_subset_signature.pair_issue_count == 0
                and max(tolerance - (reference_subset_signature.min_clearance or 0.0), 0.0)
                <= max(grid, 1e-9)
            ):
                continue

            patch_count += 1
            reference_subset_union = _cached_union(cache, subset_polygons)
            patch_neighborhood = reference_subset_union.buffer(
                max(grid, patch_radius, 1e-9),
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
            shrink_budget_context = _courtyard_pair_shrink_budget_context(
                subset_polygons,
                tolerance=tolerance,
                grid=grid,
            )
            cluster_kind = str(cluster.get("kind", "short_edge_only"))
            simple_nonpair_cluster = (
                reference_subset_signature.pair_issue_count == 0
                and len(subset_polygons) == 1
                and reference_subset_signature.short_edge_count <= 4
            )
            simple_pair_cluster = (
                reference_subset_signature.pair_issue_count > 0
                and len(subset_polygons) <= 2
                and reference_subset_signature.short_edge_count <= 4
            )
            light_close_pair_cluster = (
                reference_subset_signature.point_touch_count == 0
                and reference_subset_signature.close_pair_count <= 1
                and reference_subset_signature.short_edge_count <= 4
            )
            severe_nonpair_cluster = (
                len(subset_polygons) > 1
                or reference_subset_signature.short_edge_count > 8
                or max(
                    tolerance - (reference_subset_signature.min_clearance or 0.0),
                    0.0,
                )
                > max(tolerance * 0.5, grid)
            )
            best_candidate: tuple[
                tuple[float, ...],
                list[Polygon],
                list[list[int]],
                str,
                _CoverageDefectSignature,
                dict[str, float],
                dict[str, int],
            ] | None = None

            def best_candidate_has_operator_prefix(
                expected_prefixes: tuple[str, ...],
            ) -> bool:
                if best_candidate is None:
                    return False
                operator = best_candidate[3]
                return operator.startswith(expected_prefixes)

            def best_candidate_satisfies_local_contract(
                expected_prefixes: tuple[str, ...],
            ) -> bool:
                if best_candidate is None:
                    return False
                operator = best_candidate[3]
                signature = best_candidate[4]
                difference_metrics = best_candidate[5]
                if not operator.startswith(expected_prefixes):
                    return False
                if not _coverage_signature_satisfies_scale_contract(
                    signature,
                    target_scale=tolerance,
                    grid=grid,
                ):
                    return False
                return _difference_metrics_have_small_fidelity_drift(
                    difference_metrics,
                    edit_zone_area=float(patch_neighborhood.area),
                    target_scale=tolerance,
                    grid=grid,
                )

            direct_pair_cluster = (
                enable_pair_cluster_rescue
                and reference_subset_signature.pair_issue_count > 0
                and len(subset_polygons) <= 3
                and reference_subset_signature.short_edge_count <= 4
                and not light_close_pair_cluster
            )

            def consider_candidate(
                candidate_polygons: list[Polygon],
                candidate_sources: list[list[int]],
                *,
                operator: str,
                applied_counts: dict[str, int] | None = None,
            ) -> None:
                nonlocal best_candidate
                if not candidate_polygons:
                    return

                candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    candidate_polygons,
                    target_scale=tolerance,
                )
                if not _coverage_signature_improves(
                    reference_subset_signature,
                    candidate_signature,
                    grid=grid,
                    target_scale=tolerance,
                ):
                    return

                candidate_signature_score = _coverage_signature_score(
                    candidate_signature,
                    target_scale=tolerance,
                )
                if (
                    best_candidate is not None
                    and candidate_signature_score > best_candidate[0][:6]
                ):
                    return

                if _cached_overlap_area(cache, candidate_polygons) > max(grid * grid, 1e-9):
                    return

                candidate_union = _cached_union(cache, candidate_polygons)
                difference_metrics = _cached_difference_area_metrics(
                    cache,
                    subset_polygons,
                    candidate_polygons,
                )
                area_balance_budget = _area_balance_budget(
                    difference_metrics["symmetric_difference_area"],
                    grid=grid,
                    short_edge_count=reference_subset_signature.short_edge_count,
                    short_edge_threshold=tolerance,
                )
                edit_zone = patch_neighborhood
                if (
                    operator.startswith("coverage_pair_issue_shrink_")
                    and shrink_budget_context is not None
                ):
                    edit_zone, shrink_area_budget = shrink_budget_context
                    area_balance_budget = max(area_balance_budget, shrink_area_budget)
                change_outside_patch = _change_outside_edit_zone(
                    reference_subset_union,
                    candidate_union,
                    edit_zone=edit_zone,
                )
                if change_outside_patch > max(
                    float(edit_zone.area),
                    grid * grid,
                    1e-9,
                ):
                    return
                if abs(difference_metrics["union_area_delta"]) > area_balance_budget:
                    return

                score = (
                    *candidate_signature_score,
                    difference_metrics["reference_minus_candidate_area"],
                    difference_metrics["candidate_minus_reference_area"],
                    difference_metrics["symmetric_difference_area"],
                    abs(difference_metrics["union_area_delta"]),
                )
                if best_candidate is None or score < best_candidate[0]:
                    best_candidate = (
                        score,
                        candidate_polygons,
                        candidate_sources,
                        operator,
                        candidate_signature,
                        difference_metrics,
                        dict(applied_counts or {operator: 1}),
                    )

            if cluster_kind == "short_edge_only":
                direct_graph_rewrite = _rewrite_coverage_cluster_graphically(
                    subset_polygons,
                    subset_sources,
                    cluster=cluster,
                    tolerance=tolerance,
                    grid=grid,
                    diagnostics=diagnostics,
                    cache=cache,
                )
                if direct_graph_rewrite is not None:
                    operator, candidate_polygons, candidate_sources = direct_graph_rewrite
                    operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                    consider_candidate(
                        candidate_polygons,
                        candidate_sources,
                        operator=operator,
                    )

            if direct_pair_cluster:
                for (
                    operator,
                    candidate_polygons,
                    candidate_sources,
                ) in _direct_pair_issue_cluster_candidates(
                    subset_polygons,
                    subset_sources,
                    tolerance=tolerance,
                    grid=grid,
                    diagnostics=diagnostics,
                    cache=cache,
                ):
                    operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                    consider_candidate(
                        candidate_polygons,
                        candidate_sources,
                        operator=operator,
                    )

            if not best_candidate_satisfies_local_contract(("coverage_pair_issue_",)):
                for candidate_tolerance in _iter_local_coverage_simplify_tolerances(
                    reference_subset_signature,
                    patch_tolerances=patch_tolerances,
                ):
                    operator = f"coverage_simplify_{candidate_tolerance:.3f}"
                    operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                    try:
                        simplified = shapely.coverage_simplify(
                            np.asarray(subset_polygons, dtype=object),
                            candidate_tolerance,
                            simplify_boundary=True,
                        )
                    except (AttributeError, GEOSException, TypeError, ValueError) as exc:
                        diagnostics["coverage_simplify_failed"] = True
                        _record_geos_exception(diagnostics, "coverage_simplify", exc)
                        continue

                    candidate_polygons: list[Polygon] = []
                    candidate_sources: list[list[int]] = []
                    for geometry, indices in zip(simplified, subset_sources):
                        for polygon in _canonicalize(geometry, grid, diagnostics):
                            candidate_polygons.append(orient(polygon, sign=1.0))
                            candidate_sources.append(indices)
                    consider_candidate(
                        candidate_polygons,
                        candidate_sources,
                        operator=operator,
                    )
                    if best_candidate_satisfies_local_contract(
                        ("coverage_simplify_", "coverage_pair_issue_"),
                    ):
                        break

            simplify_candidate_satisfies_local_contract = (
                best_candidate_satisfies_local_contract(
                    ("coverage_simplify_", "coverage_pair_issue_")
                )
            )

            if (
                best_candidate is None
                and reference_subset_signature.pair_issue_count == 0
                and not simple_nonpair_cluster
                and severe_nonpair_cluster
                and enable_patch_union_fallback
            ):
                for operator_name, radius in _iter_nonpair_patch_fallback_configs(
                    reference_subset_signature,
                    subset_size=len(subset_polygons),
                    patch_radii=patch_radii,
                ):
                    operator = f"coverage_patch_{operator_name}_{radius:.3f}"
                    operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                    candidate = _apply_local_patch_union_operator(
                        subset_polygons,
                        subset_sources,
                        radius=radius,
                        operator=operator_name,
                        grid=grid,
                        diagnostics=diagnostics,
                        patch_geometry=reference_subset_union,
                    )
                    if candidate is None:
                        continue
                    candidate_polygons, candidate_sources = candidate
                    consider_candidate(
                        candidate_polygons,
                        candidate_sources,
                        operator=operator,
                    )

            if (
                reference_subset_signature.pair_issue_count > 0
                and not simplify_candidate_satisfies_local_contract
                and not light_close_pair_cluster
            ):
                if enable_pair_cluster_rescue and not direct_pair_cluster:
                    for (
                        operator,
                        candidate_polygons,
                        candidate_sources,
                    ) in _direct_pair_issue_cluster_candidates(
                        subset_polygons,
                        subset_sources,
                        tolerance=tolerance,
                        grid=grid,
                        diagnostics=diagnostics,
                        cache=cache,
                    ):
                        operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                        consider_candidate(
                            candidate_polygons,
                            candidate_sources,
                            operator=operator,
                        )

                if (
                    not simple_pair_cluster
                    and len(subset_polygons) > 2
                    and reference_subset_signature.short_edge_count > 8
                    and enable_patch_union_fallback
                    and enable_pair_cluster_rescue
                    and not best_candidate_has_operator_prefix(
                        ("coverage_simplify_", "coverage_point_bridge_", "coverage_pair_merge_"),
                    )
                ):
                    for radius in patch_radii:
                        for operator_name in ("open_close", "close_open"):
                            operator = f"coverage_pair_{operator_name}_{radius:.3f}"
                            operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                            candidate = _apply_local_patch_union_operator(
                                subset_polygons,
                                subset_sources,
                                radius=radius,
                                operator=operator_name,
                                grid=grid,
                                diagnostics=diagnostics,
                                merge_components=True,
                                patch_geometry=reference_subset_union,
                            )
                            if candidate is None:
                                continue
                            candidate_polygons, candidate_sources = candidate
                            consider_candidate(
                                candidate_polygons,
                                candidate_sources,
                                operator=operator,
                            )

            if best_candidate is None:
                continue

            (
                _,
                replacement_polygons,
                replacement_sources,
                operator,
                _,
                _,
                applied_counts,
            ) = best_candidate
            unaffected_polygons = [
                polygon
                for index, polygon in enumerate(current_polygons)
                if index not in affected_indices
            ]
            unaffected_sources = [
                indices
                for index, indices in enumerate(current_sources)
                if index not in affected_indices
            ]
            current_polygons, current_sources = _stable_sort(
                unaffected_polygons + replacement_polygons,
                unaffected_sources + replacement_sources,
            )
            for applied_operator, count in applied_counts.items():
                operator_applied[applied_operator] = (
                    operator_applied.get(applied_operator, 0) + count
                )
            patch_applied_count += 1
            accepted_patch = True
            break

        if not accepted_patch:
            break

    if patch_applied_count == 0:
        return None

    candidate_signature = _cached_coverage_defect_signature(
        cache,
        current_polygons,
        target_scale=tolerance,
    )
    if not _coverage_signature_improves(
        reference_signature,
        candidate_signature,
        grid=grid,
        target_scale=tolerance,
    ):
        return None

    edit_zone = _cached_coverage_edit_zone(
        cache,
        polygons,
        target_scale=tolerance,
        radius=local_radius,
    )
    candidate_union = _cached_union(cache, current_polygons)
    difference_metrics = _cached_difference_area_metrics(
        cache,
        polygons,
        current_polygons,
    )
    area_balance_budget = _area_balance_budget(
        difference_metrics["symmetric_difference_area"],
        grid=grid,
        short_edge_count=reference_signature.short_edge_count,
        short_edge_threshold=tolerance,
    )
    change_outside_edit_zone = _change_outside_edit_zone(
        reference_union,
        candidate_union,
        edit_zone=edit_zone,
    )
    if change_outside_edit_zone > max(float(edit_zone.area), grid * grid, 1e-9):
        return None
    if abs(difference_metrics["union_area_delta"]) > area_balance_budget:
        return None

    return _CoverageSimplifyCandidate(
        label="local",
        polygons=current_polygons,
        source_map=current_sources,
        signature=candidate_signature,
        difference_metrics=difference_metrics,
        change_outside_edit_zone=change_outside_edit_zone,
        edit_zone_area=float(edit_zone.area),
        area_balance_budget=area_balance_budget,
        patch_count=patch_count,
        patch_applied_count=patch_applied_count,
        operator_attempts=operator_attempts,
        operator_applied=operator_applied,
    )


def _simplify_coverage(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    tolerance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    if tolerance <= 0 or not polygons:
        diagnostics["coverage_simplify_applied"] = False
        diagnostics["coverage_simplify_tolerance"] = 0.0
        return polygons, source_map

    cache = _CoverageEvalCache()
    before_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=tolerance,
    )
    diagnostics["coverage_simplify_short_edge_count_before"] = (
        before_signature.short_edge_count
    )
    if (
        before_signature.short_edge_count == 0
        and before_signature.pair_issue_count == 0
        and max(tolerance - (before_signature.min_clearance or 0.0), 0.0)
        <= max(grid, 1e-9)
    ):
        diagnostics["coverage_simplify_applied"] = False
        diagnostics["coverage_simplify_tolerance"] = 0.0
        diagnostics["coverage_simplify_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        return polygons, source_map

    diagnostics["coverage_simplify_tolerance"] = tolerance
    global_candidate = _simplify_coverage_globally(
        polygons,
        source_map,
        tolerance=tolerance,
        grid=grid,
        diagnostics=diagnostics,
        cache=cache,
    )
    local_attempted = _should_attempt_local_coverage_candidate(
        before_signature,
        global_candidate,
        target_scale=tolerance,
        grid=grid,
        purpose="coverage",
    )
    diagnostics["coverage_simplify_local_candidate_attempted"] = local_attempted
    local_candidate = (
        _simplify_coverage_locally(
            polygons,
            source_map,
            tolerance=tolerance,
            grid=grid,
            diagnostics=diagnostics,
            enable_patch_union_fallback=False,
            enable_pair_cluster_rescue=False,
            cache=cache,
        )
        if local_attempted
        else None
    )

    chosen_candidate = global_candidate
    if global_candidate is None:
        chosen_candidate = local_candidate
    elif local_candidate is not None:
        local_score = (
            *_coverage_signature_score(
                local_candidate.signature,
                target_scale=tolerance,
            ),
            local_candidate.difference_metrics["reference_minus_candidate_area"],
            local_candidate.difference_metrics["candidate_minus_reference_area"],
            local_candidate.difference_metrics["symmetric_difference_area"],
            abs(local_candidate.difference_metrics["union_area_delta"]),
        )
        global_score = (
            *_coverage_signature_score(
                global_candidate.signature,
                target_scale=tolerance,
            ),
            global_candidate.difference_metrics["reference_minus_candidate_area"],
            global_candidate.difference_metrics["candidate_minus_reference_area"],
            global_candidate.difference_metrics["symmetric_difference_area"],
            abs(global_candidate.difference_metrics["union_area_delta"]),
        )
        if local_score < global_score:
            chosen_candidate = local_candidate

    _apply_coverage_candidate_diagnostics(diagnostics, chosen_candidate)
    if chosen_candidate is None:
        diagnostics["coverage_simplify_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        return polygons, source_map

    return _stable_sort(chosen_candidate.polygons, chosen_candidate.source_map)


def _identity_coverage_candidate(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    target_scale: float,
    signature: _CoverageDefectSignature | None = None,
) -> _CoverageSimplifyCandidate:
    return _CoverageSimplifyCandidate(
        label="identity",
        polygons=list(polygons),
        source_map=[list(indices) for indices in source_map],
        signature=signature
        if signature is not None
        else _coverage_defect_signature(
            polygons,
            target_scale=target_scale,
        ),
        difference_metrics=_zero_difference_metrics(),
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )


def _evaluate_post_coverage_branch(
    coverage_candidate: _CoverageSimplifyCandidate,
    *,
    reference_union: BaseGeometry,
    source_lookup: dict[int, BaseGeometry],
    raw_support_union: BaseGeometry,
    min_feature_size: float,
    source_recovery_scale: float,
    grid: float,
    min_area: float,
    output_min_area: float,
    min_hole_area: float,
    apply_coverage_meshing_regularization: bool = True,
    cache: _CoverageEvalCache | None = None,
) -> _PostCoverageBranch:
    diagnostics = _empty_diagnostics(len(coverage_candidate.polygons))
    diagnostics["collect_stage_metrics"] = False
    reclaim_area_threshold = max(grid * grid, 0.25 * min_feature_size * min_feature_size)
    should_reclaim = (
        coverage_candidate.difference_metrics["reference_minus_candidate_area"]
        > reclaim_area_threshold
    )
    diagnostics["source_reclaim_skipped"] = not should_reclaim
    if should_reclaim:
        source_reclaimed_polygons, source_reclaimed_sources = _reclaim_source_supported_area(
            coverage_candidate.polygons,
            coverage_candidate.source_map,
            source_lookup=source_lookup,
            min_segment_length=min_feature_size,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=diagnostics,
        )
    else:
        source_reclaimed_polygons = coverage_candidate.polygons
        source_reclaimed_sources = coverage_candidate.source_map

    should_absorb_small_components = output_min_area > 0 and any(
        polygon.area + 1e-12 < output_min_area for polygon in source_reclaimed_polygons
    )
    diagnostics["small_component_absorb_skipped"] = not should_absorb_small_components
    if should_absorb_small_components:
        small_component_absorbed_polygons, small_component_absorbed_sources = (
            _absorb_small_supported_components(
                source_reclaimed_polygons,
                source_reclaimed_sources,
                source_lookup=source_lookup,
                raw_support_union=raw_support_union,
                min_area=output_min_area,
                min_segment_length=min_feature_size,
                grid=grid,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
            )
        )
    else:
        small_component_absorbed_polygons = source_reclaimed_polygons
        small_component_absorbed_sources = source_reclaimed_sources

    post_absorb_signature = _cached_coverage_defect_signature(
        cache,
        small_component_absorbed_polygons,
        target_scale=min_feature_size,
    )
    should_simplify_for_meshing = post_absorb_signature.short_edge_count > 0
    diagnostics["polygon_simplify_skipped"] = not should_simplify_for_meshing
    if should_simplify_for_meshing:
        boundary_regularized_polygons, boundary_regularized_sources = (
            _simplify_polygons_for_meshing(
                small_component_absorbed_polygons,
                small_component_absorbed_sources,
                min_segment_length=min_feature_size,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
            )
        )
    else:
        boundary_regularized_polygons = small_component_absorbed_polygons
        boundary_regularized_sources = small_component_absorbed_sources

    post_boundary_signature = _cached_coverage_defect_signature(
        cache,
        boundary_regularized_polygons,
        target_scale=min_feature_size,
    )
    clearance_deficit = max(
        min_feature_size - (post_boundary_signature.min_clearance or 0.0),
        0.0,
    )
    should_regularize_clearance = clearance_deficit > max(grid, 1e-9)
    diagnostics["clearance_regularization_skipped"] = not should_regularize_clearance
    if should_regularize_clearance:
        clearance_regularized_polygons, clearance_regularized_sources = (
            _regularize_low_clearance_polygons(
                boundary_regularized_polygons,
                boundary_regularized_sources,
                min_clearance=min_feature_size,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
            )
        )
    else:
        clearance_regularized_polygons = boundary_regularized_polygons
        clearance_regularized_sources = boundary_regularized_sources

    source_coordinate_recovered_polygons, source_coordinate_recovered_sources = (
        _recover_source_supported_coordinates(
            clearance_regularized_polygons,
            clearance_regularized_sources,
            source_lookup=source_lookup,
            min_segment_length=source_recovery_scale,
            grid=grid,
            diagnostics=diagnostics,
        )
    )
    if apply_coverage_meshing_regularization:
        coverage_meshing_regularized_polygons, coverage_meshing_regularized_sources = (
            _regularize_coverage_for_meshing(
                source_coordinate_recovered_polygons,
                source_coordinate_recovered_sources,
                min_segment_length=min_feature_size,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
                cache=cache,
            )
        )
    else:
        coverage_meshing_regularized_polygons = source_coordinate_recovered_polygons
        coverage_meshing_regularized_sources = source_coordinate_recovered_sources

    final_output_polygons, final_output_sources = _filter_small_output_polygons(
        coverage_meshing_regularized_polygons,
        coverage_meshing_regularized_sources,
        min_area=output_min_area,
        diagnostics=diagnostics,
    )
    final_polygons, final_sources = _stable_sort(
        final_output_polygons,
        final_output_sources,
    )
    final_signature = _cached_coverage_defect_signature(
        cache,
        final_polygons,
        target_scale=min_feature_size,
    )
    final_difference_metrics = _difference_area_metrics(
        reference_union,
        _cached_union(cache, final_polygons),
    )
    return _PostCoverageBranch(
        label=coverage_candidate.label,
        coverage_candidate=coverage_candidate,
        coverage_polygons=coverage_candidate.polygons,
        coverage_source_map=coverage_candidate.source_map,
        source_reclaimed_polygons=source_reclaimed_polygons,
        source_reclaimed_source_map=source_reclaimed_sources,
        small_component_absorbed_polygons=small_component_absorbed_polygons,
        small_component_absorbed_source_map=small_component_absorbed_sources,
        boundary_regularized_polygons=boundary_regularized_polygons,
        boundary_regularized_source_map=boundary_regularized_sources,
        clearance_regularized_polygons=clearance_regularized_polygons,
        clearance_regularized_source_map=clearance_regularized_sources,
        source_coordinate_recovered_polygons=source_coordinate_recovered_polygons,
        source_coordinate_recovered_source_map=source_coordinate_recovered_sources,
        coverage_meshing_regularized_polygons=coverage_meshing_regularized_polygons,
        coverage_meshing_regularized_source_map=coverage_meshing_regularized_sources,
        final_output_polygons=final_polygons,
        final_output_source_map=final_sources,
        final_signature=final_signature,
        final_difference_metrics=final_difference_metrics,
        diagnostics=diagnostics,
    )


def _post_coverage_branch_score(
    branch: _PostCoverageBranch,
    *,
    target_scale: float,
) -> tuple[float, ...]:
    return (
        *_coverage_signature_score(
            branch.final_signature,
            target_scale=target_scale,
        ),
        branch.final_difference_metrics["reference_minus_candidate_area"],
        branch.final_difference_metrics["candidate_minus_reference_area"],
        branch.final_difference_metrics["symmetric_difference_area"],
        abs(branch.final_difference_metrics["union_area_delta"]),
    )


def _copy_post_coverage_diagnostics(
    diagnostics: dict[str, Any],
    branch_diagnostics: dict[str, Any],
) -> None:
    for key, value in branch_diagnostics.items():
        if key.startswith("source_reclaim_") or key.startswith("polygon_simplify_"):
            diagnostics[key] = value
        elif key.startswith("small_component_absorb_"):
            diagnostics[key] = value
        elif key.startswith("clearance_regularization_"):
            diagnostics[key] = value
        elif key.startswith("source_coordinate_recovery_"):
            diagnostics[key] = value
        elif key.startswith("coverage_meshing_regularization_"):
            diagnostics[key] = value
        elif key.startswith("final_min_area_filter_"):
            diagnostics[key] = value
    diagnostics["geos_exception_count"] += branch_diagnostics["geos_exception_count"]
    diagnostics["geos_exception_messages"].extend(
        branch_diagnostics["geos_exception_messages"]
    )


def _choose_post_coverage_branch(
    identity_branch: _PostCoverageBranch | None,
    global_branch: _PostCoverageBranch | None,
    local_branch: _PostCoverageBranch | None,
    *,
    target_scale: float,
    grid: float,
) -> tuple[_PostCoverageBranch, tuple[float, ...] | None, tuple[float, ...] | None]:
    chosen_branch = global_branch or local_branch or identity_branch
    if chosen_branch is None:
        raise ValueError("at least one post-coverage branch must be available")
    global_score = (
        _post_coverage_branch_score(
            global_branch,
            target_scale=target_scale,
        )
        if global_branch is not None
        else None
    )
    local_score = (
        _post_coverage_branch_score(
            local_branch,
            target_scale=target_scale,
        )
        if local_branch is not None
        else None
    )
    if (
        identity_branch is not None
        and chosen_branch is identity_branch
        and local_branch is not None
    ):
        chosen_branch = local_branch
    if global_branch is not None and local_branch is not None:
        fidelity_tolerance = max(grid * grid, 1e-9)
        if (
            local_score is not None
            and global_score is not None
            and local_score < global_score
            and local_branch.final_difference_metrics["reference_minus_candidate_area"]
            <= global_branch.final_difference_metrics["reference_minus_candidate_area"]
            * 1.1
            + fidelity_tolerance
            and local_branch.final_difference_metrics["symmetric_difference_area"]
            <= global_branch.final_difference_metrics["symmetric_difference_area"]
            * 1.1
            + fidelity_tolerance
            and abs(local_branch.final_difference_metrics["union_area_delta"])
            <= abs(global_branch.final_difference_metrics["union_area_delta"])
            + fidelity_tolerance
        ):
            chosen_branch = local_branch
    return chosen_branch, global_score, local_score


def _apply_local_polygon_repairs(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_segment_length: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
    stage_prefix: str,
    enable_defect_operators: bool,
    enable_simplify_operators: bool,
) -> tuple[list[Polygon], list[list[int]]]:
    before_stats = _segment_length_stats(
        polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics[_candidate_metric_key(stage_prefix, "segment_length_before")] = (
        before_stats
    )

    if min_segment_length <= 0 or not polygons:
        diagnostics[_candidate_metric_key(stage_prefix, "applied")] = False
        diagnostics[_candidate_metric_key(stage_prefix, "tolerance")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "segment_length_after")] = (
            before_stats
        )
        return polygons, source_map

    diagnostics[_candidate_metric_key(stage_prefix, "short_edge_count_before")] = (
        before_stats["short_edge_count"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "target_min_edge_length")] = (
        min_segment_length
    )
    diagnostics[_candidate_metric_key(stage_prefix, "tolerance")] = min_segment_length
    if (
        enable_simplify_operators
        and not enable_defect_operators
        and before_stats["short_edge_count"] == 0
    ):
        diagnostics[_candidate_metric_key(stage_prefix, "candidate_count")] = 0
        diagnostics[_candidate_metric_key(stage_prefix, "applied_count")] = 0
        diagnostics[_candidate_metric_key(stage_prefix, "edit_zone_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "change_outside_edit_zone")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "rejected_nonlocal_count")] = 0
        diagnostics[
            _candidate_metric_key(stage_prefix, "rejected_area_imbalance_count")
        ] = 0
        diagnostics[
            _candidate_metric_key(stage_prefix, "rejected_non_improving_count")
        ] = 0
        diagnostics[_candidate_metric_key(stage_prefix, "overlap_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "short_edge_count_after")] = 0
        diagnostics[_candidate_metric_key(stage_prefix, "applied")] = False
        diagnostics[_candidate_metric_key(stage_prefix, "segment_length_after")] = (
            before_stats
        )
        diagnostics[_candidate_metric_key(stage_prefix, "area_balance_budget")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "symmetric_difference_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "reference_minus_candidate_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "candidate_minus_reference_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "signed_area_delta")] = 0.0
        return polygons, source_map

    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []
    edit_zone_area = 0.0
    change_outside_edit_zone = 0.0
    candidate_count = 0
    applied_count = 0
    rejected_nonlocal_count = 0
    rejected_area_imbalance_count = 0
    rejected_non_improving_count = 0
    total_area_balance_budget = 0.0

    for polygon, indices in zip(polygons, source_map):
        polygon_signature = _polygon_defect_signature(
            polygon,
            target_scale=min_segment_length,
        )
        needs_clearance_repair = enable_defect_operators and (
            polygon_signature.clearance_deficit > max(grid, 1e-9)
        )
        needs_short_edge_repair = enable_simplify_operators and (
            polygon_signature.short_edge_count > 0
        )
        if not (needs_clearance_repair or needs_short_edge_repair):
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        candidate_count += 1
        accepted_candidate: Polygon | None = None
        local_candidates: list[_RepairCandidate] = []

        if needs_clearance_repair:
            courtyard_candidate = _try_close_courtyard_passage(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if courtyard_candidate is not None:
                local_candidates.append(courtyard_candidate)

            clearance_candidate = _try_polygon_clearance_opening(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if clearance_candidate is not None:
                local_candidates.append(clearance_candidate)

        if needs_short_edge_repair:
            simplify_tolerances = (
                min_segment_length * 0.5,
                min_segment_length * 0.75,
                min_segment_length,
            )
            for tolerance in simplify_tolerances:
                simplify_candidate = _try_polygon_local_simplify(
                    polygon,
                    tolerance=tolerance,
                    grid=grid,
                    diagnostics=diagnostics,
                )
                if simplify_candidate is not None:
                    local_candidates.append(simplify_candidate)

        for candidate in local_candidates:
            _operator_count_increment(
                diagnostics,
                _candidate_metric_key(stage_prefix, "operator_attempts"),
                candidate.operator,
            )
            accepted_candidate, reason, candidate_metrics = _accept_local_candidate(
                polygon,
                candidate,
                target_scale=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
            )
            if accepted_candidate is None:
                if candidate_metrics:
                    edit_zone_area += candidate_metrics.get("edit_zone_area", 0.0)
                    change_outside_edit_zone += candidate_metrics.get(
                        "change_outside_edit_zone",
                        0.0,
                    )
                if reason == "nonlocal":
                    rejected_nonlocal_count += 1
                elif reason == "area_imbalance":
                    rejected_area_imbalance_count += 1
                else:
                    rejected_non_improving_count += 1
                continue

            edit_zone_area += candidate_metrics["edit_zone_area"]
            change_outside_edit_zone += candidate_metrics["change_outside_edit_zone"]
            total_area_balance_budget += candidate_metrics["area_balance_budget"]
            _operator_count_increment(
                diagnostics,
                _candidate_metric_key(stage_prefix, "operator_applied"),
                candidate.operator,
            )
            break

        if accepted_candidate is not None:
            applied_count += 1
            candidate_polygons.append(accepted_candidate)
            candidate_sources.append(indices)
            continue

        candidate_polygons.append(polygon)
        candidate_sources.append(indices)

    diagnostics[_candidate_metric_key(stage_prefix, "candidate_count")] = candidate_count
    diagnostics[_candidate_metric_key(stage_prefix, "applied_count")] = applied_count
    diagnostics[_candidate_metric_key(stage_prefix, "edit_zone_area")] = edit_zone_area
    diagnostics[_candidate_metric_key(stage_prefix, "change_outside_edit_zone")] = (
        change_outside_edit_zone
    )
    diagnostics[_candidate_metric_key(stage_prefix, "rejected_nonlocal_count")] = (
        rejected_nonlocal_count
    )
    diagnostics[
        _candidate_metric_key(stage_prefix, "rejected_area_imbalance_count")
    ] = rejected_area_imbalance_count
    diagnostics[
        _candidate_metric_key(stage_prefix, "rejected_non_improving_count")
    ] = rejected_non_improving_count

    if applied_count == 0:
        diagnostics[_candidate_metric_key(stage_prefix, "overlap_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "short_edge_count_after")] = (
            before_stats["short_edge_count"]
        )
        diagnostics[_candidate_metric_key(stage_prefix, "applied")] = False
        diagnostics[_candidate_metric_key(stage_prefix, "segment_length_after")] = (
            before_stats
        )
        diagnostics[_candidate_metric_key(stage_prefix, "area_balance_budget")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "symmetric_difference_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "reference_minus_candidate_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "candidate_minus_reference_area")] = 0.0
        diagnostics[_candidate_metric_key(stage_prefix, "signed_area_delta")] = 0.0
        return polygons, source_map

    reference_union = unary_union(polygons)
    overlap_area = _coverage_overlap_area(candidate_polygons)
    diagnostics[_candidate_metric_key(stage_prefix, "overlap_area")] = overlap_area

    if overlap_area > max(grid * grid, 1e-9):
        diagnostics[_candidate_metric_key(stage_prefix, "failed")] = True
        diagnostics[_candidate_metric_key(stage_prefix, "short_edge_count_after")] = (
            before_stats["short_edge_count"]
        )
        diagnostics[_candidate_metric_key(stage_prefix, "segment_length_after")] = (
            before_stats
        )
        diagnostics[_candidate_metric_key(stage_prefix, "applied")] = False
        return polygons, source_map

    after_stats = _segment_length_stats(
        candidate_polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics[_candidate_metric_key(stage_prefix, "short_edge_count_after")] = (
        after_stats["short_edge_count"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "applied")] = applied_count > 0
    difference_metrics = _difference_area_metrics(
        reference_union,
        unary_union(candidate_polygons),
    )
    diagnostics[_candidate_metric_key(stage_prefix, "symmetric_difference_area")] = (
        difference_metrics["symmetric_difference_area"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "reference_minus_candidate_area")] = (
        difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "candidate_minus_reference_area")] = (
        difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "signed_area_delta")] = (
        difference_metrics["union_area_delta"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "area_balance_budget")] = (
        total_area_balance_budget
    )
    diagnostics[_candidate_metric_key(stage_prefix, "segment_length_after")] = (
        after_stats
    )
    return _stable_sort(candidate_polygons, candidate_sources)


def _repair_local_defects(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_clearance: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    return _apply_local_polygon_repairs(
        polygons,
        source_map,
        min_segment_length=min_clearance,
        grid=grid,
        min_area=min_area,
        min_hole_area=min_hole_area,
        diagnostics=diagnostics,
        stage_prefix="local_defect_repair",
        enable_defect_operators=True,
        enable_simplify_operators=False,
    )


def _simplify_polygons_for_meshing(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_segment_length: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    return _apply_local_polygon_repairs(
        polygons,
        source_map,
        min_segment_length=min_segment_length,
        grid=grid,
        min_area=min_area,
        min_hole_area=min_hole_area,
        diagnostics=diagnostics,
        stage_prefix="polygon_simplify",
        enable_defect_operators=False,
        enable_simplify_operators=True,
    )


def _regularize_low_clearance_polygons(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_clearance: float,
    grid: float,
    min_area: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    diagnostics["clearance_regularization_threshold"] = min_clearance
    if min_clearance <= 0 or not polygons:
        diagnostics["clearance_regularization_applied"] = False
        return polygons, source_map

    radius = min_clearance / 2.0
    overlap_tolerance = max(grid * grid, 1e-9)
    candidate_count = 0
    improved_count = 0
    failed_count = 0
    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []

    for polygon, indices in zip(polygons, source_map):
        try:
            clearance_before = float(shapely.minimum_clearance(polygon))
        except GEOSException as exc:
            _record_geos_exception(diagnostics, "minimum_clearance", exc)
            clearance_before = np.inf

        if not np.isfinite(clearance_before) or clearance_before + 1e-12 >= min_clearance:
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        candidate_count += 1
        opened_parts = _apply_opening(polygon, radius, grid, diagnostics)
        if len(opened_parts) != 1:
            failed_count += 1
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        closed_parts = _apply_closing(opened_parts[0], radius, grid, diagnostics)
        if len(closed_parts) != 1:
            failed_count += 1
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        cleaned_parts: list[Polygon] = []
        for part in closed_parts:
            without_small_holes = _remove_small_holes(
                part,
                min_hole_area,
                diagnostics,
            )
            cleaned_parts.extend(_canonicalize(without_small_holes, grid, diagnostics))

        if len(cleaned_parts) != 1 or cleaned_parts[0].area < min_area:
            failed_count += 1
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        candidate = orient(cleaned_parts[0], sign=1.0)
        try:
            clearance_after = float(shapely.minimum_clearance(candidate))
        except GEOSException as exc:
            _record_geos_exception(diagnostics, "minimum_clearance_candidate", exc)
            clearance_after = clearance_before

        if not np.isfinite(clearance_after) or clearance_after + 1e-12 < clearance_before:
            failed_count += 1
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        if clearance_after > clearance_before + 1e-12:
            improved_count += 1
        candidate_polygons.append(candidate)
        candidate_sources.append(indices)

    overlap_area = _coverage_overlap_area(candidate_polygons)
    diagnostics["clearance_regularization_overlap_area"] = overlap_area
    diagnostics["clearance_regularization_candidate_count"] = candidate_count
    diagnostics["clearance_regularization_improved_count"] = improved_count
    diagnostics["clearance_regularization_failed_count"] = failed_count
    if improved_count == 0:
        diagnostics["clearance_regularization_overlap_area"] = 0.0
        diagnostics["clearance_regularization_applied"] = False
        return polygons, source_map

    if overlap_area > overlap_tolerance:
        diagnostics["clearance_regularization_applied"] = False
        return polygons, source_map

    diagnostics["clearance_regularization_applied"] = improved_count > 0
    return _stable_sort(candidate_polygons, candidate_sources)


def condition_polygon_coverage(
    polygons: Sequence[BaseGeometry],
    *,
    source_map: Sequence[Sequence[int]] | None = None,
    options: ConditioningOptions,
) -> ConditioningResult:
    """Condition a noisy polygon coverage into a deterministic meshing input.

    The pipeline is:

    1. extract and canonicalize all polygonal parts on an explicit precision
       grid
    2. apply per-polygon opening at ``min_feature_size / 2``
    3. build merge groups and apply group-wise closing at
       ``merge_distance / 2``
    4. reconstruct the authoritative global coverage from noded linework
    5. repair local courtyard / clearance defects before any shared-boundary
       simplification
    6. regularize shared boundaries coverage-wide
    7. reclaim source-supported area locally where doing so does not worsen the
       declared-scale defect signature
    8. remove any remaining short meshing-hostile edges, but only if the
       candidate output still satisfies the coverage invariants
    9. regularize any residual low-clearance polygons with the same declared
       scale before meshing
    10. recover original supported coordinates wherever the cleaned topology
        already matches a valid source-supported boundary

    The conditioner uses only Shapely / GEOS operations and returns valid,
    non-overlapping single polygons. Only malformed user arguments raise hard
    errors; data-induced failures are converted into diagnostics and the
    pipeline continues with the best valid geometry available at each stage.
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
    output_grid = _derive_output_grid(options, grid)
    meshing_scale = max(options.min_feature_size, output_grid)
    r_open = options.min_feature_size / 2.0
    r_close = options.merge_distance / 2.0
    diagnostics["precision_grid"] = grid
    diagnostics["output_grid"] = output_grid
    diagnostics["opening_radius"] = r_open
    diagnostics["closing_radius"] = r_close
    diagnostics["collect_stage_metrics"] = options.collect_stage_metrics
    diagnostics["enable_logging"] = options.enable_logging

    _log_conditioning_start(
        len(polygons),
        options,
        grid=grid,
        output_grid=output_grid,
        opening_radius=r_open,
        closing_radius=r_close,
    )

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
    source_lookup = _build_source_geometry_lookup(atomic_polygons, atomic_sources)
    raw_support_union = unary_union(atomic_polygons)
    atomic_metrics = _record_stage_metrics(
        diagnostics,
        "atomic_input",
        atomic_polygons,
        short_edge_threshold=meshing_scale,
    )
    if atomic_metrics is not None:
        _log_conditioning_stage(
            "atomic_input",
            atomic_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    opened_polygons: list[Polygon] = []
    opened_sources: list[list[int]] = []

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

    opened_metrics = _record_stage_metrics(
        diagnostics,
        "opened",
        opened_polygons,
        short_edge_threshold=meshing_scale,
    )
    if opened_metrics is not None:
        _log_conditioning_stage(
            "opened",
            opened_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    merge_groups = _build_merge_groups(opened_polygons, options.merge_distance)
    diagnostics["merged_group_count"] = len(merge_groups)

    regularized_group_geometries: list[BaseGeometry] = []
    regularized_group_sources: list[list[int]] = []

    for group in merge_groups:
        group_polygons = [opened_polygons[index] for index in group]
        group_member_sources = [opened_sources[index] for index in group]
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

        component_sources = _assign_component_sources(
            filtered_parts,
            group_polygons,
            group_member_sources,
            grid=grid,
            diagnostics=diagnostics,
        )
        for component, indices in zip(filtered_parts, component_sources):
            regularized_group_geometries.append(component)
            regularized_group_sources.append(indices)

    regularized_polygons = [
        polygon
        for geometry in regularized_group_geometries
        for polygon in _extract_polygon_parts(geometry)
    ]
    regularized_metrics = _record_stage_metrics(
        diagnostics,
        "regularized_groups",
        regularized_polygons,
        short_edge_threshold=meshing_scale,
        extra={"group_count": len(regularized_group_geometries)},
    )
    if regularized_metrics is not None:
        _log_conditioning_stage(
            "regularized_groups",
            regularized_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    if _should_reconstruct_global_coverage(
        regularized_polygons,
        grid=grid,
        diagnostics=diagnostics,
    ):
        rebuilt_geometries, rebuilt_sources = _reconstruct_global_coverage(
            regularized_group_geometries,
            regularized_group_sources,
            grid,
            diagnostics,
        )
    else:
        rebuilt_geometries = regularized_group_geometries
        rebuilt_sources = regularized_group_sources

    rebuilt_polygons = [
        polygon
        for geometry in rebuilt_geometries
        for polygon in _extract_polygon_parts(geometry)
    ]
    reconstructed_metrics = _record_stage_metrics(
        diagnostics,
        "reconstructed",
        rebuilt_polygons,
        short_edge_threshold=meshing_scale,
        extra={"group_count": len(rebuilt_geometries)},
    )
    if reconstructed_metrics is not None:
        _log_conditioning_stage(
            "reconstructed",
            reconstructed_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons: list[Polygon] = []
    final_sources: list[list[int]] = []
    for geometry, indices in zip(rebuilt_geometries, rebuilt_sources):
        for polygon in _canonicalize_for_output(
            geometry,
            output_grid=output_grid,
            min_area=0.0,
            diagnostics=diagnostics,
        ):
            final_polygons.append(orient(polygon, sign=1.0))
            final_sources.append(sorted(set(indices)))

    presimplify_metrics = _record_stage_metrics(
        diagnostics,
        "presimplify",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if presimplify_metrics is not None:
        _log_conditioning_stage(
            "presimplify",
            presimplify_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons, final_sources = _repair_local_defects(
        final_polygons,
        final_sources,
        min_clearance=meshing_scale,
        grid=output_grid,
        min_area=0.0,
        min_hole_area=options.min_hole_area,
        diagnostics=diagnostics,
    )
    local_defect_metrics = _record_stage_metrics(
        diagnostics,
        "local_defect_repaired",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if local_defect_metrics is not None:
        _log_conditioning_stage(
            "local_defect_repaired",
            local_defect_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    coverage_eval_cache = _CoverageEvalCache()
    coverage_reference_union = _cached_union(coverage_eval_cache, final_polygons)
    coverage_simplify_tolerance = max(
        _derive_coverage_simplify_tolerance(options),
        output_grid,
    )
    identity_signature = _cached_coverage_defect_signature(
        coverage_eval_cache,
        final_polygons,
        target_scale=coverage_simplify_tolerance,
    )
    identity_candidate = _identity_coverage_candidate(
        final_polygons,
        final_sources,
        target_scale=coverage_simplify_tolerance,
        signature=identity_signature,
    )
    coverage_candidates = [identity_candidate]
    global_candidate = _simplify_coverage_globally(
        final_polygons,
        final_sources,
        tolerance=coverage_simplify_tolerance,
        grid=output_grid,
        diagnostics=diagnostics,
        cache=coverage_eval_cache,
    )
    if global_candidate is not None:
        coverage_candidates.append(global_candidate)
    local_candidate = None
    local_candidate_attempted = _should_attempt_local_coverage_candidate(
        identity_candidate.signature,
        global_candidate,
        target_scale=coverage_simplify_tolerance,
        grid=output_grid,
        purpose="coverage",
    )
    diagnostics["coverage_simplify_local_candidate_attempted"] = (
        local_candidate_attempted
    )
    if local_candidate_attempted:
        local_candidate = _simplify_coverage_locally(
            final_polygons,
            final_sources,
            tolerance=coverage_simplify_tolerance,
            grid=output_grid,
            diagnostics=diagnostics,
            enable_patch_union_fallback=False,
            enable_pair_cluster_rescue=False,
            cache=coverage_eval_cache,
        )
    if local_candidate is not None:
        coverage_candidates.append(local_candidate)

    diagnostics["coverage_simplify_branch_prescores"] = {
        candidate.label: list(
            _coverage_candidate_prescore(
                candidate,
                target_scale=meshing_scale,
                output_min_area=options.min_area,
            )
        )
        for candidate in coverage_candidates
    }
    candidates_to_evaluate = _select_post_coverage_candidates_for_evaluation(
        identity_candidate,
        [candidate for candidate in (global_candidate, local_candidate) if candidate is not None],
        target_scale=meshing_scale,
        grid=output_grid,
        output_min_area=options.min_area,
    )
    diagnostics["coverage_simplify_branches_evaluated"] = [
        candidate.label for candidate in candidates_to_evaluate
    ]
    finalized_branches = {
        candidate.label: _evaluate_post_coverage_branch(
            candidate,
            reference_union=coverage_reference_union,
            source_lookup=source_lookup,
            raw_support_union=raw_support_union,
            min_feature_size=meshing_scale,
            source_recovery_scale=options.min_feature_size,
            grid=output_grid,
            min_area=0.0,
            output_min_area=options.min_area,
            min_hole_area=options.min_hole_area,
            cache=coverage_eval_cache,
        )
        for candidate in candidates_to_evaluate
    }
    identity_branch = finalized_branches.get("identity")
    global_branch = finalized_branches.get("global")
    local_branch = finalized_branches.get("local")
    chosen_branch, global_score, local_score = _choose_post_coverage_branch(
        identity_branch,
        global_branch,
        local_branch,
        target_scale=meshing_scale,
        grid=output_grid,
    )
    finalized_labels = list(finalized_branches.keys())
    diagnostics["coverage_simplify_branches_finalized"] = finalized_labels
    diagnostics["coverage_simplify_branch_scores"] = {
        **(
            {
                "identity": list(
                    _post_coverage_branch_score(
                        identity_branch,
                        target_scale=meshing_scale,
                    )
                ),
            }
            if identity_branch is not None
            else {}
        ),
        **(
            {
                "global": list(global_score),
            }
            if global_score is not None
            else {}
        ),
        **(
            {
                "local": list(local_score),
            }
            if local_score is not None
            else {}
        ),
    }
    diagnostics["coverage_simplify_selected_branch"] = chosen_branch.label
    _apply_coverage_candidate_diagnostics(
        diagnostics,
        chosen_branch.coverage_candidate,
    )
    _copy_post_coverage_diagnostics(diagnostics, chosen_branch.diagnostics)

    final_polygons = chosen_branch.coverage_polygons
    final_sources = chosen_branch.coverage_source_map
    coverage_simplified_metrics = _record_stage_metrics(
        diagnostics,
        "coverage_simplified",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if coverage_simplified_metrics is not None:
        _log_conditioning_stage(
            "coverage_simplified",
            coverage_simplified_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.source_reclaimed_polygons
    final_sources = chosen_branch.source_reclaimed_source_map
    source_reclaimed_metrics = _record_stage_metrics(
        diagnostics,
        "source_reclaimed",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if source_reclaimed_metrics is not None:
        _log_conditioning_stage(
            "source_reclaimed",
            source_reclaimed_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.small_component_absorbed_polygons
    final_sources = chosen_branch.small_component_absorbed_source_map
    small_component_absorbed_metrics = _record_stage_metrics(
        diagnostics,
        "small_component_absorbed",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if small_component_absorbed_metrics is not None:
        _log_conditioning_stage(
            "small_component_absorbed",
            small_component_absorbed_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.boundary_regularized_polygons
    final_sources = chosen_branch.boundary_regularized_source_map
    boundary_regularized_metrics = _record_stage_metrics(
        diagnostics,
        "boundary_regularized",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if boundary_regularized_metrics is not None:
        _log_conditioning_stage(
            "boundary_regularized",
            boundary_regularized_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.clearance_regularized_polygons
    final_sources = chosen_branch.clearance_regularized_source_map
    clearance_regularized_metrics = _record_stage_metrics(
        diagnostics,
        "clearance_regularized",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if clearance_regularized_metrics is not None:
        _log_conditioning_stage(
            "clearance_regularized",
            clearance_regularized_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.source_coordinate_recovered_polygons
    final_sources = chosen_branch.source_coordinate_recovered_source_map
    source_coordinate_recovered_metrics = _record_stage_metrics(
        diagnostics,
        "source_coordinate_recovered",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if source_coordinate_recovered_metrics is not None:
        _log_conditioning_stage(
            "source_coordinate_recovered",
            source_coordinate_recovered_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.coverage_meshing_regularized_polygons
    final_sources = chosen_branch.coverage_meshing_regularized_source_map
    coverage_meshing_regularized_metrics = _record_stage_metrics(
        diagnostics,
        "coverage_meshing_regularized",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if coverage_meshing_regularized_metrics is not None:
        _log_conditioning_stage(
            "coverage_meshing_regularized",
            coverage_meshing_regularized_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )

    final_polygons = chosen_branch.final_output_polygons
    final_sources = chosen_branch.final_output_source_map
    final_polygons, final_sources = _stable_sort(final_polygons, final_sources)
    diagnostics["output_count"] = len(final_polygons)
    diagnostics["overlap_area_after"] = _coverage_overlap_area(final_polygons)
    diagnostics["min_clearance_after"] = _minimum_clearance(final_polygons)
    _record_stage_metrics(
        diagnostics,
        "final_output",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    _log_conditioning_summary(
        diagnostics,
        short_edge_threshold=meshing_scale,
    )

    return ConditioningResult(
        polygons=final_polygons,
        source_map=final_sources,
        diagnostics=diagnostics,
    )
