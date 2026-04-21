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
import math
import time
from typing import Any, Callable, Iterable, Literal, Sequence

import numpy as np
import shapely
from shapely import BufferCapStyle, BufferJoinStyle, GeometryCollection
from shapely.errors import GEOSException
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient
from shapely.ops import polygonize_full, unary_union
from shapely.strtree import STRtree
from shapely.validation import make_valid

try:
    from .. import _dtcc_builder
except Exception:  # pragma: no cover - import fallback for partial builds
    _dtcc_builder = None

from ..logging import debug, info


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
    ring_contact_count: int = 0


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
    area_balance_budget_override: float | None = None


def _record_stage_seconds(
    diagnostics: dict[str, Any],
    stage_name: str,
    started_at: float,
) -> float:
    elapsed = float(time.perf_counter() - started_at)
    diagnostics.setdefault("stage_seconds", {})[stage_name] = elapsed
    return elapsed


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
    post_recovery_regularized_polygons: list[Polygon]
    post_recovery_regularized_source_map: list[list[int]]
    coverage_contact_regularized_polygons: list[Polygon]
    coverage_contact_regularized_source_map: list[list[int]]
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
    # Keep the Python clustering path as the canonical implementation. The
    # native helper is still available for targeted experiments, but it can
    # change cluster membership across identical runs, which then changes the
    # accepted patch sequence upstream of 2D/3D meshing.

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
        indices = sorted(int(index) for index in cluster.get("indices", []))
        normalized_clusters.append(
            {
                "indices": indices,
                "kind": str(cluster.get("kind", "short_edge_only")),
                "short_edge_count": int(cluster.get("short_edge_count", 0)),
                "pair_issue_count": int(cluster.get("pair_issue_count", 0)),
            }
        )
    normalized_clusters.sort(
        key=lambda cluster: (
            tuple(cluster["indices"]),
            str(cluster["kind"]),
            int(cluster["pair_issue_count"]),
            int(cluster["short_edge_count"]),
        )
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
        "source_coordinate_recovery_contract_reject_count": 0,
        "source_coordinate_recovery_short_edge_count_before": 0,
        "source_coordinate_recovery_short_edge_count_after": 0,
        "source_coordinate_recovery_pair_issue_count_before": 0,
        "source_coordinate_recovery_pair_issue_count_after": 0,
        "source_coordinate_recovery_ring_contact_count_before": 0,
        "source_coordinate_recovery_ring_contact_count_after": 0,
        "source_coordinate_recovery_min_clearance_before": None,
        "source_coordinate_recovery_min_clearance_after": None,
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
        "coverage_void_regularization_applied": False,
        "coverage_void_regularization_threshold": 0.0,
        "coverage_void_regularization_hole_cleanup_count": 0,
        "coverage_void_regularization_gap_patch_count": 0,
        "coverage_void_regularization_gap_patch_applied_count": 0,
        "coverage_void_regularization_notch_simplify_count": 0,
        "coverage_void_regularization_operator_attempts": {},
        "coverage_void_regularization_operator_applied": {},
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
        "coverage_contact_regularization_applied": False,
        "coverage_contact_regularization_selected_branch": "identity",
        "coverage_contact_regularization_short_edge_count_before": 0,
        "coverage_contact_regularization_short_edge_count_after": 0,
        "coverage_contact_regularization_pair_issue_count_before": 0,
        "coverage_contact_regularization_pair_issue_count_after": 0,
        "coverage_contact_regularization_ring_contact_count_before": 0,
        "coverage_contact_regularization_ring_contact_count_after": 0,
        "coverage_contact_regularization_reference_minus_candidate_area": 0.0,
        "coverage_contact_regularization_candidate_minus_reference_area": 0.0,
        "coverage_contact_regularization_signed_area_delta": 0.0,
        "coverage_contact_regularization_operator_attempts": {},
        "coverage_contact_regularization_operator_applied": {},
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


def _remove_meshing_hostile_holes(
    poly: Polygon,
    *,
    min_hole_area: float,
    min_clearance: float,
    diagnostics: dict[str, Any],
    sliver_width_factor: float = 3.0,
    compactness_threshold: float = 0.1,
) -> Polygon:
    if poly.is_empty or not poly.interiors:
        return orient(poly, sign=1.0)

    kept_holes: list[list[tuple[float, float]]] = []
    removed_small = 0
    removed_sliver = 0
    width_limit = (
        sliver_width_factor * min_clearance if min_clearance > 0 else 0.0
    )

    for ring in poly.interiors:
        hole = Polygon(ring)
        if hole.is_empty:
            removed_small += 1
            continue
        if min_hole_area > 0 and hole.area < min_hole_area:
            removed_small += 1
            continue

        hole_clearance = _safe_minimum_clearance(hole)
        remove_hole = False
        if (
            min_clearance > 0
            and hole_clearance is not None
            and hole_clearance + 1.0e-12 < min_clearance
        ):
            remove_hole = True
        elif (
            min_clearance > 0
            and hole_clearance is not None
            and hole_clearance <= width_limit + 1.0e-12
        ):
            perimeter = float(hole.length)
            compactness = 1.0
            if perimeter > 0.0:
                compactness = (4.0 * math.pi * float(hole.area)) / (
                    perimeter * perimeter
                )
            if compactness < compactness_threshold:
                remove_hole = True

        if remove_hole:
            removed_sliver += 1
            continue

        kept_holes.append(list(ring.coords))

    if removed_small:
        diagnostics["dropped_small_count"] += removed_small
    if removed_sliver:
        diagnostics["dropped_meshing_hole_count"] = (
            diagnostics.get("dropped_meshing_hole_count", 0) + removed_sliver
        )
    return orient(Polygon(poly.exterior.coords, kept_holes), sign=1.0)


def _compactness_ratio(geometry: BaseGeometry) -> float:
    if geometry.is_empty:
        return 1.0
    perimeter = float(geometry.length)
    if perimeter <= 0.0:
        return 1.0
    return float((4.0 * math.pi * float(geometry.area)) / (perimeter * perimeter))


def _apply_meshing_hostile_hole_cleanup(
    polygons: Sequence[Polygon],
    source_map: Sequence[Sequence[int]],
    *,
    min_hole_area: float,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]], int]:
    cleaned_polygons: list[Polygon] = []
    cleaned_sources: list[list[int]] = []
    applied_count = 0

    for polygon, indices in zip(polygons, source_map):
        cleaned = _remove_meshing_hostile_holes(
            polygon,
            min_hole_area=min_hole_area,
            min_clearance=min_clearance,
            diagnostics=diagnostics,
        )
        if cleaned.equals_exact(polygon, tolerance=0.0):
            cleaned_polygons.append(polygon)
            cleaned_sources.append(list(indices))
            continue
        cleaned_parts = _canonicalize(cleaned, grid, diagnostics)
        if len(cleaned_parts) != 1:
            cleaned_polygons.append(polygon)
            cleaned_sources.append(list(indices))
            continue
        cleaned_polygon = orient(cleaned_parts[0], sign=1.0)
        if not cleaned_polygon.equals_exact(polygon, tolerance=0.0):
            applied_count += 1
        cleaned_polygons.append(cleaned_polygon)
        cleaned_sources.append(list(indices))

    return cleaned_polygons, cleaned_sources, applied_count


def _meshing_hostile_gap_candidates(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    grid: float,
    slit_width_factor: float = 1.5,
    compactness_threshold: float = 0.1,
) -> list[tuple[int, int, float, float]]:
    if target_scale <= 0 or len(polygons) < 2:
        return []

    width_limit = max(target_scale * slit_width_factor, grid)
    tree = STRtree(polygons)
    point_tolerance = max(target_scale * 1.0e-6, 1.0e-9)
    candidates: list[tuple[int, int, float, float]] = []

    for index, polygon in enumerate(polygons):
        for candidate in tree.query(_expanded_query_geometry(polygon, radius=width_limit)):
            other_index = int(candidate)
            if other_index <= index:
                continue
            other = polygons[other_index]
            try:
                distance = float(polygon.distance(other))
            except GEOSException:
                continue
            if distance <= target_scale + point_tolerance or distance > width_limit + point_tolerance:
                continue
            overlap_envelope = _local_offset_overlap_envelope(
                [polygon, other],
                radius=width_limit,
                grid=grid,
            )
            if overlap_envelope is None or overlap_envelope.is_empty:
                continue
            envelope_compactness = min(
                _compactness_ratio(component)
                for component in shapely.get_parts(overlap_envelope)
                if not component.is_empty and component.area > 0.0
            )
            if envelope_compactness >= compactness_threshold:
                continue
            candidates.append((index, other_index, distance, envelope_compactness))

    candidates.sort(key=lambda item: (item[3], item[2], item[0], item[1]))
    return candidates


def _regularize_meshing_hostile_voids(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_segment_length: float,
    grid: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> tuple[list[Polygon], list[list[int]]]:
    diagnostics["coverage_void_regularization_threshold"] = min_segment_length
    if min_segment_length <= 0 or not polygons:
        diagnostics["coverage_void_regularization_applied"] = False
        return polygons, source_map

    if cache is None:
        cache = _CoverageEvalCache()

    current_polygons, current_sources, hole_cleanup_count = _apply_meshing_hostile_hole_cleanup(
        polygons,
        source_map,
        min_hole_area=min_hole_area,
        min_clearance=min_segment_length,
        grid=grid,
        diagnostics=diagnostics,
    )

    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}
    gap_patch_count = 0
    gap_patch_applied_count = 0
    notch_simplify_count = 0

    while True:
        gap_candidates = _meshing_hostile_gap_candidates(
            current_polygons,
            target_scale=min_segment_length,
            grid=grid,
        )
        gap_patch_count += len(gap_candidates)
        if not gap_candidates:
            break

        current_signature = _cached_coverage_defect_signature(
            cache,
            current_polygons,
            target_scale=min_segment_length,
        )
        current_gap_count = len(gap_candidates)
        reference_union = _cached_union(cache, current_polygons)
        best_candidate: tuple[
            tuple[float, ...],
            list[Polygon],
            list[list[int]],
            str,
        ] | None = None

        for left_index, right_index, gap_distance, envelope_compactness in gap_candidates:
            operator = "coverage_void_gap_bridge"
            operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
            subset_polygons = [current_polygons[left_index], current_polygons[right_index]]
            subset_sources = [current_sources[left_index], current_sources[right_index]]
            candidate = _apply_local_close_pair_bridge_operator(
                subset_polygons,
                subset_sources,
                radius=min_segment_length,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                continue

            replacement_polygons, replacement_sources = candidate
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
            candidate_polygons, candidate_sources = _stable_sort(
                unaffected_polygons + replacement_polygons,
                unaffected_sources + replacement_sources,
            )
            candidate_polygons, candidate_sources, candidate_hole_cleanup_count = (
                _apply_meshing_hostile_hole_cleanup(
                    candidate_polygons,
                    candidate_sources,
                    min_hole_area=min_hole_area,
                    min_clearance=min_segment_length,
                    grid=grid,
                    diagnostics=diagnostics,
                )
            )
            candidate_signature = _cached_coverage_defect_signature(
                cache,
                candidate_polygons,
                target_scale=min_segment_length,
            )
            if not _coverage_signature_not_worse(
                current_signature,
                candidate_signature,
                grid=grid,
                target_scale=min_segment_length,
            ):
                continue

            candidate_gap_count = len(
                _meshing_hostile_gap_candidates(
                    candidate_polygons,
                    target_scale=min_segment_length,
                    grid=grid,
                )
            )
            if candidate_gap_count >= current_gap_count:
                continue

            difference_metrics = _cached_difference_area_metrics(
                cache,
                current_polygons,
                candidate_polygons,
            )
            overlap_envelope = _local_offset_overlap_envelope(
                subset_polygons,
                radius=max(1.5 * min_segment_length, grid),
                grid=grid,
            )
            if overlap_envelope is None or overlap_envelope.is_empty:
                continue
            edit_zone = overlap_envelope.buffer(
                max(min_segment_length, grid),
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
            change_outside_edit_zone = _change_outside_edit_zone(
                reference_union,
                _cached_union(cache, candidate_polygons),
                edit_zone=edit_zone,
            )
            if change_outside_edit_zone > max(float(edit_zone.area), grid * grid, 1.0e-9):
                continue
            if not _difference_metrics_have_small_fidelity_drift(
                difference_metrics,
                edit_zone_area=float(edit_zone.area),
                target_scale=min_segment_length,
                grid=grid,
            ):
                continue

            extra_area_budget = max(
                2.0 * float(overlap_envelope.area),
                min_segment_length * min_segment_length,
                16.0 * grid * grid,
            )
            if difference_metrics["candidate_minus_reference_area"] > extra_area_budget:
                continue

            score = (
                candidate_gap_count,
                -candidate_hole_cleanup_count,
                candidate_signature.short_edge_count,
                candidate_signature.vertex_count,
                difference_metrics["candidate_minus_reference_area"],
                difference_metrics["symmetric_difference_area"],
                gap_distance,
                envelope_compactness,
            )
            if best_candidate is None or score < best_candidate[0]:
                best_candidate = (
                    score,
                    candidate_polygons,
                    candidate_sources,
                    operator,
                )

        if best_candidate is None:
            break

        _, current_polygons, current_sources, applied_operator = best_candidate
        operator_applied[applied_operator] = operator_applied.get(applied_operator, 0) + 1
        gap_patch_applied_count += 1

    refined_polygons: list[Polygon] = []
    refined_sources: list[list[int]] = []
    for polygon, indices in zip(current_polygons, current_sources):
        signature = _polygon_defect_signature(
            polygon,
            target_scale=min_segment_length,
        )
        should_try_micro_detour_simplify = (
            signature.short_edge_count == 0
            and signature.ring_contact_count == 0
            and signature.clearance is not None
            and signature.min_edge_length is not None
        )
        should_try_notch_simplify = (
            should_try_micro_detour_simplify
            and signature.min_edge_length <= min_segment_length + 2.0 * grid + 1.0e-9
            and signature.clearance <= min_segment_length + 2.0 * grid + 1.0e-9
        )

        best_polygon = polygon
        best_signature = signature
        if should_try_micro_detour_simplify:
            micro_detour_operator = "coverage_void_micro_detour_chain"
            operator_attempts[micro_detour_operator] = (
                operator_attempts.get(micro_detour_operator, 0) + 1
            )
            micro_detour_candidate = _try_polygon_micro_detour_chain_simplify(
                polygon,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if micro_detour_candidate is not None:
                normalized = _remove_meshing_hostile_holes(
                    micro_detour_candidate.polygon,
                    min_hole_area=min_hole_area,
                    min_clearance=min_segment_length,
                    diagnostics=diagnostics,
                )
                parts = _canonicalize(normalized, grid, diagnostics)
                if len(parts) == 1:
                    candidate_polygon = orient(parts[0], sign=1.0)
                    candidate_signature = _polygon_defect_signature(
                        candidate_polygon,
                        target_scale=min_segment_length,
                    )
                    if _signature_not_worse(
                        signature,
                        candidate_signature,
                        grid=grid,
                    ) and (
                        candidate_signature.min_edge_length is not None
                        and best_signature.min_edge_length is not None
                        and candidate_signature.min_edge_length
                        > best_signature.min_edge_length + grid
                        and (
                            candidate_signature.clearance is None
                            or best_signature.clearance is None
                            or candidate_signature.clearance + grid
                            >= best_signature.clearance
                        )
                    ):
                        difference_metrics = _difference_area_metrics(
                            polygon,
                            candidate_polygon,
                        )
                        change_outside_edit_zone = _change_outside_edit_zone(
                            polygon,
                            candidate_polygon,
                            edit_zone=micro_detour_candidate.edit_zone,
                        )
                        if (
                            change_outside_edit_zone
                            > max(
                                float(micro_detour_candidate.edit_zone.area),
                                grid * grid,
                                1.0e-9,
                            )
                            or not _difference_metrics_have_small_fidelity_drift(
                                difference_metrics,
                                edit_zone_area=float(micro_detour_candidate.edit_zone.area),
                                target_scale=min_segment_length,
                                grid=grid,
                            )
                        ):
                            continue
                        best_polygon = candidate_polygon
                        best_signature = candidate_signature
                        operator_applied[micro_detour_operator] = (
                            operator_applied.get(micro_detour_operator, 0) + 1
                        )

        if not should_try_notch_simplify:
            if not best_polygon.equals_exact(polygon, tolerance=0.0):
                notch_simplify_count += 1
            refined_polygons.append(best_polygon)
            refined_sources.append(list(indices))
            continue

        for tolerance in (
            min_segment_length,
            min_segment_length * 1.25,
            min_segment_length * 1.5,
            min_segment_length * 2.0,
        ):
            operator = f"coverage_void_notch_simplify_{tolerance:.3f}"
            operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
            candidate = _try_polygon_local_simplify(
                polygon,
                tolerance=tolerance,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                continue
            normalized = _remove_meshing_hostile_holes(
                candidate.polygon,
                min_hole_area=min_hole_area,
                min_clearance=min_segment_length,
                diagnostics=diagnostics,
            )
            parts = _canonicalize(normalized, grid, diagnostics)
            if len(parts) != 1:
                continue
            candidate_polygon = orient(parts[0], sign=1.0)
            candidate_signature = _polygon_defect_signature(
                candidate_polygon,
                target_scale=min_segment_length,
            )
            if not _signature_not_worse(
                signature,
                candidate_signature,
                grid=grid,
            ):
                continue
            if (
                candidate_signature.min_edge_length is None
                or best_signature.min_edge_length is None
                or candidate_signature.min_edge_length
                <= best_signature.min_edge_length + grid
            ):
                continue
            if (
                candidate_signature.clearance is None
                or best_signature.clearance is None
                or candidate_signature.clearance + grid < best_signature.clearance
            ):
                continue
            difference_metrics = _difference_area_metrics(
                polygon,
                candidate_polygon,
            )
            change_outside_edit_zone = _change_outside_edit_zone(
                polygon,
                candidate_polygon,
                edit_zone=candidate.edit_zone,
            )
            if (
                change_outside_edit_zone
                > max(float(candidate.edit_zone.area), grid * grid, 1.0e-9)
                or not _difference_metrics_have_small_fidelity_drift(
                    difference_metrics,
                    edit_zone_area=float(candidate.edit_zone.area),
                    target_scale=min_segment_length,
                    grid=grid,
                )
            ):
                continue
            best_polygon = candidate_polygon
            best_signature = candidate_signature
            operator_applied[operator] = operator_applied.get(operator, 0) + 1

        if not best_polygon.equals_exact(polygon, tolerance=0.0):
            notch_simplify_count += 1
        refined_polygons.append(best_polygon)
        refined_sources.append(list(indices))

    current_polygons, current_sources = _stable_sort(refined_polygons, refined_sources)

    diagnostics["coverage_void_regularization_hole_cleanup_count"] = hole_cleanup_count
    diagnostics["coverage_void_regularization_gap_patch_count"] = gap_patch_count
    diagnostics["coverage_void_regularization_gap_patch_applied_count"] = (
        gap_patch_applied_count
    )
    diagnostics["coverage_void_regularization_notch_simplify_count"] = (
        notch_simplify_count
    )
    diagnostics["coverage_void_regularization_operator_attempts"] = operator_attempts
    diagnostics["coverage_void_regularization_operator_applied"] = operator_applied
    diagnostics["coverage_void_regularization_applied"] = (
        hole_cleanup_count > 0
        or gap_patch_applied_count > 0
        or notch_simplify_count > 0
    )
    return current_polygons, current_sources


def _regularize_final_polygon_shapes(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_segment_length: float,
    grid: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]]:
    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}
    refined_polygons: list[Polygon] = []
    refined_sources: list[list[int]] = []
    near_threshold_limit = min_segment_length + 2.0 * grid + 1.0e-9
    local_shape_edge_limit = 2.0 * min_segment_length + 2.0 * grid + 1.0e-9

    for polygon, indices in zip(polygons, source_map):
        signature = _polygon_defect_signature(
            polygon,
            target_scale=min_segment_length,
        )
        if (
            signature.short_edge_count != 0
            or signature.ring_contact_count != 0
            or signature.clearance is None
            or signature.min_edge_length is None
        ):
            refined_polygons.append(polygon)
            refined_sources.append(list(indices))
            continue
        if (
            signature.clearance > near_threshold_limit
            and signature.min_edge_length > local_shape_edge_limit
        ):
            refined_polygons.append(polygon)
            refined_sources.append(list(indices))
            continue

        def _normalize_final_shape_candidate(
            candidate_polygon: Polygon,
        ) -> tuple[Polygon, _PolygonDefectSignature] | None:
            normalized = _remove_meshing_hostile_holes(
                candidate_polygon,
                min_hole_area=min_hole_area,
                min_clearance=min_segment_length,
                diagnostics=diagnostics,
            )
            parts = _canonicalize(normalized, grid, diagnostics)
            if len(parts) != 1:
                return None
            finalized = orient(parts[0], sign=1.0)
            finalized_signature = _polygon_defect_signature(
                finalized,
                target_scale=min_segment_length,
            )
            return finalized, finalized_signature

        current_polygon = polygon
        current_signature = signature
        micro_detour_applied_count = 0
        micro_detour_operator = "final_shape_micro_detour_chain"
        while True:
            operator_attempts[micro_detour_operator] = (
                operator_attempts.get(micro_detour_operator, 0) + 1
            )
            candidate = _try_polygon_micro_detour_chain_simplify(
                current_polygon,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                break
            normalized_candidate = _normalize_final_shape_candidate(candidate.polygon)
            if normalized_candidate is None:
                break
            candidate_polygon, candidate_signature = normalized_candidate
            if not _signature_not_worse(
                current_signature,
                candidate_signature,
                grid=grid,
            ):
                break
            if (
                candidate_signature.clearance is None
                or current_signature.clearance is None
                or candidate_signature.clearance + grid < current_signature.clearance
            ):
                break
            if (
                candidate_signature.min_edge_length is None
                or current_signature.min_edge_length is None
                or candidate_signature.min_edge_length + grid
                < current_signature.min_edge_length
            ):
                break
            if candidate_signature.vertex_count >= current_signature.vertex_count:
                break

            difference_metrics = _difference_area_metrics(
                current_polygon,
                candidate_polygon,
            )
            relative_change_budget = max(
                0.0025 * float(current_polygon.area),
                32.0 * grid * grid,
                1.0e-9,
            )
            area_growth_budget = max(
                0.002 * float(current_polygon.area),
                16.0 * grid * grid,
                1.0e-9,
            )
            fill_loss_budget = max(
                16.0 * grid * grid,
                0.02 * difference_metrics["candidate_minus_reference_area"] + grid * grid,
                1.0e-9,
            )
            if difference_metrics["symmetric_difference_area"] > relative_change_budget:
                break
            if difference_metrics["union_area_delta"] > area_growth_budget:
                break
            if difference_metrics["reference_minus_candidate_area"] > fill_loss_budget:
                break

            current_polygon = candidate_polygon
            current_signature = candidate_signature
            micro_detour_applied_count += 1

        if micro_detour_applied_count > 0:
            operator_applied[micro_detour_operator] = (
                operator_applied.get(micro_detour_operator, 0)
                + micro_detour_applied_count
            )

        short_walk_applied_count = 0
        short_walk_added_area_total = 0.0
        short_walk_added_area_budget = max(
            0.005 * float(current_polygon.area),
            32.0 * grid * grid,
            1.0e-9,
        )
        short_walk_operator = "final_shape_same_turn_short_walk"
        while True:
            operator_attempts[short_walk_operator] = (
                operator_attempts.get(short_walk_operator, 0) + 1
            )
            candidate = _try_polygon_same_turn_short_walk_collapse(
                current_polygon,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                break
            normalized_candidate = _normalize_final_shape_candidate(candidate.polygon)
            if normalized_candidate is None:
                break
            candidate_polygon, candidate_signature = normalized_candidate
            if not _signature_not_worse(
                current_signature,
                candidate_signature,
                grid=grid,
            ):
                break
            if (
                candidate_signature.clearance is None
                or current_signature.clearance is None
                or candidate_signature.clearance + grid < current_signature.clearance
            ):
                break
            if (
                candidate_signature.min_edge_length is None
                or current_signature.min_edge_length is None
                or candidate_signature.min_edge_length + grid
                < current_signature.min_edge_length
            ):
                break
            if candidate_signature.vertex_count >= current_signature.vertex_count:
                break

            difference_metrics = _difference_area_metrics(
                current_polygon,
                candidate_polygon,
            )
            missing_tolerance = max(16.0 * grid * grid, 1.0e-9)
            area_growth_budget = max(
                0.003 * float(current_polygon.area),
                0.5 * float(candidate.edit_zone.area),
                16.0 * grid * grid,
                1.0e-9,
            )
            if difference_metrics["reference_minus_candidate_area"] > missing_tolerance:
                break
            if difference_metrics["candidate_minus_reference_area"] > area_growth_budget:
                break
            if (
                short_walk_added_area_total
                + difference_metrics["candidate_minus_reference_area"]
                > short_walk_added_area_budget
            ):
                break

            current_polygon = candidate_polygon
            current_signature = candidate_signature
            short_walk_applied_count += 1
            short_walk_added_area_total += difference_metrics[
                "candidate_minus_reference_area"
            ]

        if short_walk_applied_count > 0:
            operator_applied[short_walk_operator] = (
                operator_applied.get(short_walk_operator, 0)
                + short_walk_applied_count
            )

        fill_chain_applied_count = 0
        fill_chain_added_area_total = 0.0
        fill_chain_added_area_budget = max(
            0.005 * float(current_polygon.area),
            32.0 * grid * grid,
            1.0e-9,
        )
        fill_chain_operator = "final_shape_fill_chain_collapse"
        while True:
            operator_attempts[fill_chain_operator] = (
                operator_attempts.get(fill_chain_operator, 0) + 1
            )
            candidate = _try_polygon_fill_chain_collapse(
                current_polygon,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                break
            normalized_candidate = _normalize_final_shape_candidate(candidate.polygon)
            if normalized_candidate is None:
                break
            candidate_polygon, candidate_signature = normalized_candidate
            if not _signature_not_worse(
                current_signature,
                candidate_signature,
                grid=grid,
            ):
                break
            if (
                candidate_signature.clearance is None
                or current_signature.clearance is None
                or candidate_signature.clearance + grid < current_signature.clearance
            ):
                break
            if (
                candidate_signature.min_edge_length is None
                or current_signature.min_edge_length is None
                or candidate_signature.min_edge_length + grid
                < current_signature.min_edge_length
            ):
                break
            if candidate_signature.vertex_count >= current_signature.vertex_count:
                break

            difference_metrics = _difference_area_metrics(
                current_polygon,
                candidate_polygon,
            )
            missing_tolerance = max(16.0 * grid * grid, 1.0e-9)
            area_growth_budget = max(
                0.004 * float(current_polygon.area),
                0.75 * float(candidate.edit_zone.area),
                16.0 * grid * grid,
                1.0e-9,
            )
            if difference_metrics["reference_minus_candidate_area"] > missing_tolerance:
                break
            if difference_metrics["candidate_minus_reference_area"] > area_growth_budget:
                break
            if (
                fill_chain_added_area_total
                + difference_metrics["candidate_minus_reference_area"]
                > fill_chain_added_area_budget
            ):
                break

            current_polygon = candidate_polygon
            current_signature = candidate_signature
            fill_chain_applied_count += 1
            fill_chain_added_area_total += difference_metrics[
                "candidate_minus_reference_area"
            ]

        if fill_chain_applied_count > 0:
            operator_applied[fill_chain_operator] = (
                operator_applied.get(fill_chain_operator, 0)
                + fill_chain_applied_count
            )

        bevel_corner_applied_count = 0
        bevel_corner_operator = "final_shape_bevel_corner_collapse"
        while True:
            operator_attempts[bevel_corner_operator] = (
                operator_attempts.get(bevel_corner_operator, 0) + 1
            )
            candidate = _try_polygon_bevel_corner_collapse(
                current_polygon,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                break
            normalized_candidate = _normalize_final_shape_candidate(candidate.polygon)
            if normalized_candidate is None:
                break
            candidate_polygon, candidate_signature = normalized_candidate
            if not _signature_not_worse(
                current_signature,
                candidate_signature,
                grid=grid,
            ):
                break
            if (
                candidate_signature.clearance is None
                or current_signature.clearance is None
                or candidate_signature.clearance + grid < current_signature.clearance
            ):
                break
            if (
                candidate_signature.min_edge_length is None
                or current_signature.min_edge_length is None
                or candidate_signature.min_edge_length + grid
                < current_signature.min_edge_length
            ):
                break
            if candidate_signature.vertex_count >= current_signature.vertex_count:
                break

            difference_metrics = _difference_area_metrics(
                current_polygon,
                candidate_polygon,
            )
            missing_tolerance = max(16.0 * grid * grid, 1.0e-9)
            area_growth_budget = max(
                0.001 * float(current_polygon.area),
                0.25 * float(candidate.edit_zone.area),
                16.0 * grid * grid,
                1.0e-9,
            )
            if difference_metrics["reference_minus_candidate_area"] > missing_tolerance:
                break
            if difference_metrics["candidate_minus_reference_area"] > area_growth_budget:
                break

            current_polygon = candidate_polygon
            current_signature = candidate_signature
            bevel_corner_applied_count += 1

        if bevel_corner_applied_count > 0:
            operator_applied[bevel_corner_operator] = (
                operator_applied.get(bevel_corner_operator, 0)
                + bevel_corner_applied_count
            )

        shape_polygon = current_polygon
        shape_signature = current_signature

        best_polygon = shape_polygon
        best_signature = shape_signature
        best_score: tuple[float, int, float] | None = None
        best_operator: str | None = None

        def consider_final_shape_candidate(
            candidate_polygon: Polygon,
            candidate_signature: _PolygonDefectSignature,
            *,
            operator: str,
            max_relative_change: float,
            max_area_growth_ratio: float,
        ) -> None:
            nonlocal best_polygon, best_signature, best_score, best_operator
            if not _signature_not_worse(
                shape_signature,
                candidate_signature,
                grid=grid,
            ):
                return
            if (
                candidate_signature.min_edge_length is None
                or best_signature.min_edge_length is None
                or candidate_signature.min_edge_length
                <= best_signature.min_edge_length + grid
            ):
                return
            if (
                candidate_signature.clearance is None
                or best_signature.clearance is None
                or candidate_signature.clearance + grid < best_signature.clearance
            ):
                return
            if candidate_signature.vertex_count >= best_signature.vertex_count:
                return

            difference_metrics = _difference_area_metrics(
                shape_polygon,
                candidate_polygon,
            )
            relative_change_budget = max(
                max_relative_change * float(shape_polygon.area),
                32.0 * grid * grid,
                1.0e-9,
            )
            area_growth_budget = max(
                max_area_growth_ratio * float(shape_polygon.area),
                16.0 * grid * grid,
                1.0e-9,
            )
            if difference_metrics["symmetric_difference_area"] > relative_change_budget:
                return
            if difference_metrics["union_area_delta"] > area_growth_budget:
                return
            if difference_metrics["candidate_minus_reference_area"] > max(
                area_growth_budget,
                0.1 * difference_metrics["reference_minus_candidate_area"] + grid * grid,
            ):
                return

            score = (
                difference_metrics["symmetric_difference_area"],
                candidate_signature.vertex_count,
                -(candidate_signature.min_edge_length or 0.0),
            )
            if best_score is None or score < best_score:
                best_polygon = candidate_polygon
                best_signature = candidate_signature
                best_score = score
                best_operator = operator

        should_try_local_shape_simplify = (
            micro_detour_applied_count == 0
            and short_walk_applied_count == 0
            and fill_chain_applied_count == 0
            and bevel_corner_applied_count == 0
            and shape_signature.vertex_count >= 12
            and (
                shape_signature.clearance <= near_threshold_limit
                or shape_signature.min_edge_length <= local_shape_edge_limit
            )
        )
        if should_try_local_shape_simplify:
            for tolerance in (
                min_segment_length,
                min_segment_length * 1.25,
                min_segment_length * 1.5,
            ):
                operator = f"final_shape_local_simplify_{tolerance:.3f}"
                operator_attempts[operator] = operator_attempts.get(operator, 0) + 1
                candidate = _try_polygon_local_simplify(
                    shape_polygon,
                    tolerance=tolerance,
                    grid=grid,
                    diagnostics=diagnostics,
                )
                if candidate is None:
                    continue
                normalized_candidate = _normalize_final_shape_candidate(candidate.polygon)
                if normalized_candidate is None:
                    continue
                candidate_polygon, candidate_signature = normalized_candidate
                consider_final_shape_candidate(
                    candidate_polygon,
                    candidate_signature,
                    operator=operator,
                    max_relative_change=0.01,
                    max_area_growth_ratio=0.01,
                )

        refined_polygons.append(best_polygon)
        refined_sources.append(list(indices))
        if best_operator is not None:
            operator_applied[best_operator] = operator_applied.get(best_operator, 0) + 1

    diagnostics["final_shape_regularization_applied"] = bool(operator_applied)
    diagnostics["final_shape_regularization_operator_attempts"] = operator_attempts
    diagnostics["final_shape_regularization_operator_applied"] = operator_applied
    return _stable_sort(refined_polygons, refined_sources)


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
    def polygon_sort_key(polygon: Polygon) -> tuple[Any, ...]:
        return (
            polygon.bounds[0],
            polygon.bounds[1],
            -polygon.area,
            len(polygon.exterior.coords),
            len(polygon.interiors),
            bytes(polygon.wkb),
        )

    paired = list(zip(polygons, source_map))
    paired.sort(key=lambda item: polygon_sort_key(item[0]))
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


def _regularize_coverage_ring_contacts(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    min_segment_length: float,
    grid: float,
) -> tuple[tuple[list[Polygon], list[list[int]]], dict[str, float | int]] | None:
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

    for polygon, indices in zip(polygons, source_map):
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

    if changed_count == 0:
        return None

    stats: dict[str, float | int] = {
        "changed_count": changed_count,
        "component_count": component_count,
        "failed_count": failed_count,
        "area_delta": area_delta,
    }
    return _stable_sort(candidate_polygons, candidate_sources), stats


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
        ring_contact_count=_polygon_ring_boundary_contact_count(polygon),
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
    if candidate.ring_contact_count > reference.ring_contact_count:
        return False
    if candidate.ring_contact_count < reference.ring_contact_count:
        if candidate.clearance_deficit > max(reference.clearance_deficit, tolerance):
            return False
        if candidate.short_edge_count > reference.short_edge_count:
            return False
        reference_min_edge = reference.min_edge_length or 0.0
        candidate_min_edge = candidate.min_edge_length or 0.0
        if candidate_min_edge + tolerance < reference_min_edge:
            return False
        return True

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
    if candidate.close_pair_count < reference.close_pair_count:
        return True
    if candidate.close_pair_count > reference.close_pair_count:
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
    if candidate.ring_contact_count > reference.ring_contact_count:
        return False
    if candidate.ring_contact_count < reference.ring_contact_count:
        if candidate.clearance_deficit > max(reference.clearance_deficit, tolerance):
            return False
        if candidate.short_edge_count > reference.short_edge_count:
            return False
        reference_min_edge = reference.min_edge_length or 0.0
        candidate_min_edge = candidate.min_edge_length or 0.0
        if candidate_min_edge + tolerance < reference_min_edge:
            return False
        return True

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
    if signature.ring_contact_count > 0:
        return False
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
        candidate.pair_issue_count == reference.pair_issue_count
        and candidate.close_pair_count > reference.close_pair_count
    ):
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


def _record_polygon_repair_noop(
    diagnostics: dict[str, Any],
    *,
    stage_prefix: str,
    tolerance: float,
    before_stats: dict[str, Any],
) -> None:
    diagnostics[_candidate_metric_key(stage_prefix, "segment_length_before")] = (
        before_stats
    )
    diagnostics[_candidate_metric_key(stage_prefix, "short_edge_count_before")] = (
        before_stats["short_edge_count"]
    )
    diagnostics[_candidate_metric_key(stage_prefix, "target_min_edge_length")] = (
        tolerance
    )
    diagnostics[_candidate_metric_key(stage_prefix, "tolerance")] = tolerance
    diagnostics[_candidate_metric_key(stage_prefix, "candidate_count")] = 0
    diagnostics[_candidate_metric_key(stage_prefix, "applied_count")] = 0
    diagnostics[_candidate_metric_key(stage_prefix, "edit_zone_area")] = 0.0
    diagnostics[_candidate_metric_key(stage_prefix, "change_outside_edit_zone")] = 0.0
    diagnostics[_candidate_metric_key(stage_prefix, "rejected_nonlocal_count")] = 0
    diagnostics[_candidate_metric_key(stage_prefix, "rejected_area_imbalance_count")] = 0
    diagnostics[_candidate_metric_key(stage_prefix, "rejected_non_improving_count")] = 0
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
    diagnostics[_candidate_metric_key(stage_prefix, "reference_minus_candidate_area")] = (
        0.0
    )
    diagnostics[_candidate_metric_key(stage_prefix, "candidate_minus_reference_area")] = (
        0.0
    )
    diagnostics[_candidate_metric_key(stage_prefix, "signed_area_delta")] = 0.0
    diagnostics[_candidate_metric_key(stage_prefix, "operator_attempts")] = {}
    diagnostics[_candidate_metric_key(stage_prefix, "operator_applied")] = {}


def _coverage_signature_score(
    signature: _CoverageDefectSignature,
    *,
    target_scale: float,
) -> tuple[float, int, int, int, int, float, int]:
    return (
        max(target_scale - (signature.min_clearance or 0.0), 0.0),
        signature.ring_contact_count,
        signature.pair_issue_count,
        signature.close_pair_count,
        signature.short_edge_count,
        -(signature.min_edge_length or 0.0),
        signature.vertex_count,
    )


def _polygon_signature_score(
    signature: _PolygonDefectSignature,
    *,
    target_scale: float,
) -> tuple[int, float, int, float, int]:
    return (
        signature.ring_contact_count,
        max(target_scale - (signature.clearance or 0.0), 0.0),
        signature.short_edge_count,
        -(signature.min_edge_length or 0.0),
        signature.vertex_count,
    )


def _polygon_acute_tip_metrics(
    polygon: Polygon,
    *,
    max_angle_degrees: float = 10.0,
    min_tip_span: float = 0.0,
) -> tuple[int, float]:
    coords = list(polygon.exterior.coords[:-1])
    count = len(coords)
    if count < 3:
        return 0, 0.0

    acute_tip_count = 0
    acute_tip_span = 0.0
    for index in range(count):
        a = coords[(index - 1) % count]
        b = coords[index]
        c = coords[(index + 1) % count]
        prev_length = float(np.hypot(a[0] - b[0], a[1] - b[1]))
        next_length = float(np.hypot(c[0] - b[0], c[1] - b[1]))
        tip_span = min(prev_length, next_length)
        if tip_span + 1.0e-12 < min_tip_span:
            continue
        if prev_length <= 1.0e-12 or next_length <= 1.0e-12:
            continue
        dot = ((a[0] - b[0]) * (c[0] - b[0])) + ((a[1] - b[1]) * (c[1] - b[1]))
        dot /= prev_length * next_length
        dot = max(-1.0, min(1.0, dot))
        angle = float(np.degrees(np.arccos(dot)))
        if angle > max_angle_degrees:
            continue
        acute_tip_count += 1
        acute_tip_span += tip_span

    return acute_tip_count, acute_tip_span


def _contact_resolution_signature_score(
    signature: _CoverageDefectSignature,
    *,
    target_scale: float,
) -> tuple[int, int, int, float, int, float, int]:
    return (
        signature.ring_contact_count,
        signature.pair_issue_count,
        signature.short_edge_count,
        max(target_scale - (signature.min_clearance or 0.0), 0.0),
        signature.close_pair_count,
        -(signature.min_edge_length or 0.0),
        signature.vertex_count,
    )


def _contact_resolution_candidate_score(
    signature: _CoverageDefectSignature,
    difference_metrics: dict[str, float],
    *,
    target_scale: float,
) -> tuple[int, int, int, float, int, float, int, float, float, float, float]:
    clearance_deficit = max(target_scale - (signature.min_clearance or 0.0), 0.0)
    clearance_contract_unmet = int(clearance_deficit > max(target_scale * 1e-3, 1e-9))
    return (
        signature.ring_contact_count,
        signature.pair_issue_count,
        signature.short_edge_count,
        clearance_contract_unmet,
        signature.close_pair_count,
        difference_metrics["candidate_minus_reference_area"],
        difference_metrics["reference_minus_candidate_area"],
        difference_metrics["symmetric_difference_area"],
        abs(difference_metrics["union_area_delta"]),
        clearance_deficit,
        -(signature.min_edge_length or 0.0),
        signature.vertex_count,
    )


def _prefer_fill_only_point_contact_candidate(
    operator_name: str,
    difference_metrics: dict[str, float],
    *,
    target_scale: float,
    grid: float,
) -> int:
    if not operator_name.startswith(
        ("coverage_pair_issue_point_", "coverage_pair_issue_point_cluster_")
    ):
        return 1

    small_fill_budget = max(target_scale * target_scale, 16.0 * grid * grid, 1.0e-9)
    tiny_loss_budget = max(
        0.05 * target_scale * target_scale,
        16.0 * grid * grid,
        1.0e-9,
    )
    if difference_metrics["candidate_minus_reference_area"] > small_fill_budget:
        return 1
    if difference_metrics["reference_minus_candidate_area"] > tiny_loss_budget:
        return 1
    return 0


def _should_accept_post_contact_local_repair(
    reference_signature: _CoverageDefectSignature,
    reference_difference_metrics: dict[str, float],
    candidate_signature: _CoverageDefectSignature,
    candidate_difference_metrics: dict[str, float],
    *,
    target_scale: float,
    grid: float,
) -> bool:
    if not _coverage_signature_not_worse(
        reference_signature,
        candidate_signature,
        grid=grid,
        target_scale=target_scale,
    ):
        return False
    candidate_score = _contact_resolution_candidate_score(
        candidate_signature,
        candidate_difference_metrics,
        target_scale=target_scale,
    )
    reference_score = _contact_resolution_candidate_score(
        reference_signature,
        reference_difference_metrics,
        target_scale=target_scale,
    )
    if candidate_score < reference_score:
        return True

    reference_clearance = reference_signature.min_clearance or 0.0
    candidate_clearance = candidate_signature.min_clearance or 0.0
    crossed_subgrid_clearance_floor = (
        reference_clearance + 1.0e-9 < grid
        and candidate_clearance + 1.0e-9 >= grid
    )
    if not crossed_subgrid_clearance_floor:
        return False

    local_area_budget = max(5.0 * target_scale * target_scale, 256.0 * grid * grid, 1.0e-9)
    if (
        candidate_difference_metrics["symmetric_difference_area"]
        - reference_difference_metrics["symmetric_difference_area"]
        > local_area_budget
    ):
        return False
    if (
        candidate_difference_metrics["candidate_minus_reference_area"]
        - reference_difference_metrics["candidate_minus_reference_area"]
        > local_area_budget
    ):
        return False
    if (
        candidate_difference_metrics["reference_minus_candidate_area"]
        - reference_difference_metrics["reference_minus_candidate_area"]
        > local_area_budget
    ):
        return False

    return True


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


def _coverage_signature_contract_priority(
    signature: _CoverageDefectSignature,
    *,
    target_scale: float,
    grid: float,
) -> int:
    return 0 if _coverage_signature_satisfies_scale_contract(
        signature,
        target_scale=target_scale,
        grid=grid,
    ) else 1


def _coverage_signature_requires_polygon_regularization(
    signature: _CoverageDefectSignature,
    *,
    target_scale: float,
    grid: float,
) -> bool:
    return (
        signature.short_edge_count > 0
        or signature.ring_contact_count > 0
        or max(target_scale - (signature.min_clearance or 0.0), 0.0) > max(grid, 1e-9)
    )


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


def _nearest_ring_classification(
    polygon: Polygon,
    point: Point,
) -> tuple[str | None, int | None, float]:
    best_kind: str | None = None
    best_ring_index: int | None = None
    best_distance = float("inf")

    candidates: list[tuple[float, str, int | None]] = [
        (float(polygon.exterior.distance(point)), "exterior", None)
    ]
    candidates.extend(
        (
            float(LineString(ring.coords).distance(point)),
            "hole",
            index,
        )
        for index, ring in enumerate(polygon.interiors)
    )

    for distance, kind, ring_index in candidates:
        if distance + 1e-12 < best_distance:
            best_kind = kind
            best_ring_index = ring_index
            best_distance = distance
            continue
        if abs(distance - best_distance) > 1e-12:
            continue
        if best_kind == "exterior":
            continue
        if kind == "exterior":
            best_kind = kind
            best_ring_index = ring_index
            best_distance = distance
            continue
        if best_ring_index is None or (
            ring_index is not None and ring_index < best_ring_index
        ):
            best_kind = kind
            best_ring_index = ring_index
            best_distance = distance

    return best_kind, best_ring_index, best_distance


def _nearest_vertex_index(coords: np.ndarray, point: np.ndarray) -> int:
    distances = np.linalg.norm(coords - point, axis=1)
    return int(np.argmin(distances))


def _nearest_vertex_index_and_distance(
    coords: np.ndarray,
    point: np.ndarray,
) -> tuple[int, float]:
    distances = np.linalg.norm(coords - point, axis=1)
    index = int(np.argmin(distances))
    return index, float(distances[index])


def _ring_coords(
    polygon: Polygon,
    *,
    ring_kind: str,
    ring_index: int | None,
) -> np.ndarray | None:
    if ring_kind == "exterior":
        return np.asarray(polygon.exterior.coords[:-1], dtype=float)
    if ring_kind == "hole" and ring_index is not None:
        return np.asarray(polygon.interiors[ring_index].coords[:-1], dtype=float)
    return None


def _nearest_nonadjacent_ring_segment(
    coords: np.ndarray,
    point: np.ndarray,
    *,
    excluded_segments: set[int],
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray] | None:
    if len(coords) < 2:
        return None

    best: tuple[float, int, np.ndarray, np.ndarray, np.ndarray] | None = None
    wrapped = np.vstack([coords[1:], coords[:1]])
    for segment_index, (start, end) in enumerate(zip(coords, wrapped)):
        if segment_index in excluded_segments:
            continue
        segment = LineString([tuple(start), tuple(end)])
        try:
            projected = segment.interpolate(segment.project(Point(point)))
        except GEOSException:
            continue
        projected_xy = np.asarray(projected.coords[0], dtype=float)
        distance = float(np.hypot(*(projected_xy - point)))
        if best is None or distance < best[0]:
            best = (
                distance,
                segment_index,
                start,
                end,
                projected_xy,
            )

    if best is None:
        return None

    return best[1], best[2], best[3], best[4]


def _step_from_vertex_toward_neighbor(
    vertex: np.ndarray,
    neighbor: np.ndarray,
    *,
    distance: float,
) -> np.ndarray | None:
    edge = neighbor - vertex
    edge_length = float(np.hypot(edge[0], edge[1]))
    if edge_length <= 1e-9:
        return None
    step = min(float(distance), 0.45 * edge_length)
    if step <= 0.0:
        return None
    return vertex + (edge / edge_length) * step


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


def _iter_self_clearance_connector_half_widths(
    *,
    target_scale: float,
    grid: float,
) -> tuple[float, ...]:
    if target_scale <= 0:
        return ()

    widths: list[float] = []
    for factor in (0.25, 0.375, 0.5, 0.75, 1.0):
        width = max(float(target_scale) * factor, float(grid), 1e-9)
        if any(abs(width - existing) <= 1e-12 for existing in widths):
            continue
        widths.append(width)
    return tuple(widths)


def _iter_hole_pair_merge_distances(
    *,
    min_clearance: float,
    current_clearance: float,
    grid: float,
) -> tuple[float, ...]:
    base = max(float(min_clearance), 4.0 * float(grid), 1e-6)
    lower_bound = max(
        float(current_clearance) + max(0.05 * float(min_clearance), 4.0 * float(grid)),
        0.66 * float(min_clearance),
        4.0 * float(grid),
    )
    distances: list[float] = []
    for distance in (
        lower_bound,
        0.75 * base,
        0.9 * base,
        1.0 * base,
    ):
        distance = max(distance, 4.0 * float(grid), 1e-6)
        if any(abs(distance - existing) <= 1e-12 for existing in distances):
            continue
        distances.append(distance)
    return tuple(distances)


def _try_hole_pair_self_clearance_merge(
    polygon: Polygon,
    *,
    start_ring_index: int,
    end_ring_index: int,
    min_clearance: float,
    current_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
    require_signature_improvement: bool = True,
) -> _RepairCandidate | None:
    if (
        start_ring_index < 0
        or end_ring_index < 0
        or start_ring_index >= len(polygon.interiors)
        or end_ring_index >= len(polygon.interiors)
        or start_ring_index == end_ring_index
    ):
        return None

    try:
        hole_a = orient(
            Polygon(list(polygon.interiors[start_ring_index].coords)),
            sign=1.0,
        )
        hole_b = orient(
            Polygon(list(polygon.interiors[end_ring_index].coords)),
            sign=1.0,
        )
    except (GEOSException, ValueError):
        return None
    if hole_a.is_empty or hole_b.is_empty:
        return None

    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=min_clearance,
    )
    other_holes = [
        list(ring.coords)
        for ring_index, ring in enumerate(polygon.interiors)
        if ring_index not in {start_ring_index, end_ring_index}
    ]
    best: tuple[
        tuple[float, int, int, float, int, float, float, int],
        Polygon,
        BaseGeometry,
    ] | None = None

    def consider_merged_hole(merged_hole: BaseGeometry) -> None:
        nonlocal best
        if merged_hole.is_empty or merged_hole.geom_type != "Polygon":
            return

        candidate_polygon = _normalize_single_polygon_candidate(
            Polygon(
                list(polygon.exterior.coords),
                other_holes + [list(merged_hole.exterior.coords)],
            ),
            grid=grid,
            min_area=0.0,
            min_hole_area=0.0,
            diagnostics=diagnostics,
        )
        if candidate_polygon is None:
            return

        candidate_signature = _polygon_defect_signature(
            candidate_polygon,
            target_scale=min_clearance,
        )
        if require_signature_improvement and not _signature_improves(
            reference_signature,
            candidate_signature,
            grid=grid,
        ):
            return

        difference_metrics = _difference_area_metrics(
            polygon,
            candidate_polygon,
        )
        edit_zone = polygon.symmetric_difference(candidate_polygon)
        edge_deficit = max(
            min_clearance - (candidate_signature.min_edge_length or 0.0),
            0.0,
        )
        score = (
            candidate_signature.clearance_deficit,
            candidate_signature.short_edge_count,
            candidate_signature.ring_contact_count,
            edge_deficit,
            len(candidate_polygon.interiors),
            difference_metrics["symmetric_difference_area"],
            abs(difference_metrics["union_area_delta"]),
            candidate_signature.vertex_count,
        )
        if best is None or score < best[0]:
            best = (
                score,
                candidate_polygon,
                edit_zone,
            )

    for distance in _iter_hole_pair_merge_distances(
        min_clearance=min_clearance,
        current_clearance=current_clearance,
        grid=grid,
    ):
        try:
            merged_hole = unary_union(
                [
                    hole_a.buffer(
                        distance,
                        quad_segs=1,
                        join_style=BufferJoinStyle.mitre,
                        mitre_limit=1000.0,
                    ),
                    hole_b.buffer(
                        distance,
                        quad_segs=1,
                        join_style=BufferJoinStyle.mitre,
                        mitre_limit=1000.0,
                    ),
                ]
            ).buffer(
                -distance,
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
        except GEOSException as exc:
            _record_geos_exception(diagnostics, "hole_pair_clearance_merge", exc)
            continue
        consider_merged_hole(merged_hole)

    try:
        consider_merged_hole(unary_union([hole_a, hole_b]).convex_hull)
    except GEOSException as exc:
        _record_geos_exception(diagnostics, "hole_pair_clearance_merge", exc)

    if best is None:
        return None

    area_budget_override = max(float(best[2].area), 16.0 * grid * grid, 1e-9)
    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator="hole_pair_clearance_merge",
        area_balance_budget_override=area_budget_override,
    )


def _try_same_ring_self_clearance_connector_fill(
    polygon: Polygon,
    clearance_coords: np.ndarray,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
    ring_kind: str,
    ring_index: int | None,
    require_signature_improvement: bool = True,
) -> _RepairCandidate | None:
    coords = _ring_coords(
        polygon,
        ring_kind=ring_kind,
        ring_index=ring_index,
    )
    if coords is None or len(coords) < 4:
        return None

    start_vertex_index, start_vertex_distance = _nearest_vertex_index_and_distance(
        coords,
        clearance_coords[0],
    )
    end_vertex_index, end_vertex_distance = _nearest_vertex_index_and_distance(
        coords,
        clearance_coords[-1],
    )
    vertex_tolerance = max(2.0 * grid, 1e-6)
    if (
        start_vertex_distance > vertex_tolerance
        and end_vertex_distance > vertex_tolerance
    ):
        return None

    if start_vertex_distance <= end_vertex_distance:
        vertex_index = start_vertex_index
        vertex_xy = coords[vertex_index]
        other_xy = clearance_coords[-1]
    else:
        vertex_index = end_vertex_index
        vertex_xy = coords[vertex_index]
        other_xy = clearance_coords[0]

    previous_vertex = coords[(vertex_index - 1) % len(coords)]
    next_vertex = coords[(vertex_index + 1) % len(coords)]
    excluded_segments = {
        (vertex_index - 1) % len(coords),
        vertex_index % len(coords),
    }
    other_segment = _nearest_nonadjacent_ring_segment(
        coords,
        other_xy,
        excluded_segments=excluded_segments,
    )
    if other_segment is None:
        return None

    _, segment_start, segment_end, projected_xy = other_segment
    segment = segment_end - segment_start
    segment_length = float(np.hypot(segment[0], segment[1]))
    if segment_length <= 1e-9:
        return None
    segment_tangent = segment / segment_length

    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=min_clearance,
    )
    best: tuple[
        tuple[float, int, int, float, float, float, int],
        Polygon,
        BaseGeometry,
        str,
    ] | None = None

    for distance in _iter_ring_contact_connector_distances(
        min_clearance=min_clearance,
        grid=grid,
    ):
        previous_support = _step_from_vertex_toward_neighbor(
            vertex_xy,
            previous_vertex,
            distance=distance,
        )
        next_support = _step_from_vertex_toward_neighbor(
            vertex_xy,
            next_vertex,
            distance=distance,
        )
        if previous_support is None or next_support is None:
            continue

        local_envelope = polygon.buffer(
            distance,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        )

        for half_width in _iter_self_clearance_connector_half_widths(
            target_scale=min_clearance,
            grid=grid,
        ):
            segment_offset = min(float(half_width), 0.45 * segment_length)
            if segment_offset <= 0.0:
                continue

            left_support = projected_xy - segment_tangent * segment_offset
            right_support = projected_xy + segment_tangent * segment_offset
            patch = shapely.convex_hull(
                MultiPoint(
                    [
                        tuple(previous_support),
                        tuple(next_support),
                        tuple(left_support),
                        tuple(right_support),
                    ]
                )
            ).intersection(local_envelope)
            if patch.is_empty:
                continue

            candidate_polygon = _normalize_single_polygon_candidate(
                unary_union([polygon, patch]),
                grid=grid,
                min_area=0.0,
                min_hole_area=0.0,
                diagnostics=diagnostics,
            )
            if candidate_polygon is None:
                continue

            candidate_signature = _polygon_defect_signature(
                candidate_polygon,
                target_scale=min_clearance,
            )
            if require_signature_improvement and not _signature_improves(
                reference_signature,
                candidate_signature,
                grid=grid,
            ):
                continue

            difference_metrics = _difference_area_metrics(
                polygon,
                candidate_polygon,
            )
            edge_deficit = max(
                min_clearance - (candidate_signature.min_edge_length or 0.0),
                0.0,
            )
            score = (
                candidate_signature.clearance_deficit,
                candidate_signature.short_edge_count,
                candidate_signature.ring_contact_count,
                edge_deficit,
                difference_metrics["symmetric_difference_area"],
                abs(difference_metrics["union_area_delta"]),
                candidate_signature.vertex_count,
            )
            if best is None or score < best[0]:
                best = (
                    score,
                    candidate_polygon,
                    patch,
                    "self_clearance_connector",
                )

    if best is None:
        start_point = Point(
            float(clearance_coords[0][0]),
            float(clearance_coords[0][1]),
        )
        end_point = Point(
            float(clearance_coords[-1][0]),
            float(clearance_coords[-1][1]),
        )
        tolerance = max(grid, min_clearance * 1.0e-6, 1.0e-9)
        start_direction = _polygon_inward_direction_at_contact(
            polygon,
            start_point,
            tolerance=tolerance,
        )
        end_direction = _polygon_inward_direction_at_contact(
            polygon,
            end_point,
            tolerance=tolerance,
        )
        if start_direction is None or end_direction is None:
            return None

        for distance in _iter_ring_contact_connector_distances(
            min_clearance=min_clearance,
            grid=grid,
        ):
            start_support = _inset_contact_support_point(
                polygon,
                start_point,
                start_direction,
                distance=distance,
                tolerance=grid,
            )
            end_support = _inset_contact_support_point(
                polygon,
                end_point,
                end_direction,
                distance=distance,
                tolerance=grid,
            )
            if start_support is None or end_support is None:
                continue

            local_envelope = polygon.buffer(
                distance,
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
            if local_envelope.is_empty:
                continue

            for width_factor in _iter_ring_contact_connector_width_factors():
                corridor_half_width = _connector_half_width(
                    radius=distance,
                    grid=grid,
                    width_factor=width_factor,
                )
                try:
                    patch = LineString(
                        [tuple(start_support), tuple(end_support)]
                    ).buffer(
                        corridor_half_width,
                        quad_segs=1,
                        cap_style=BufferCapStyle.flat,
                        join_style=BufferJoinStyle.mitre,
                        mitre_limit=1000.0,
                    )
                except GEOSException as exc:
                    _record_geos_exception(
                        diagnostics,
                        "self_clearance_connector_segment",
                        exc,
                    )
                    continue
                if patch.is_empty:
                    continue

                patch = patch.intersection(local_envelope)
                if patch.is_empty:
                    continue

                candidate_polygon = _normalize_single_polygon_candidate(
                    unary_union([polygon, patch]),
                    grid=grid,
                    min_area=0.0,
                    min_hole_area=0.0,
                    diagnostics=diagnostics,
                )
                if candidate_polygon is None:
                    continue

                candidate_signature = _polygon_defect_signature(
                    candidate_polygon,
                    target_scale=min_clearance,
                )
                if require_signature_improvement and not _signature_improves(
                    reference_signature,
                    candidate_signature,
                    grid=grid,
                ):
                    continue

                difference_metrics = _difference_area_metrics(
                    polygon,
                    candidate_polygon,
                )
                edge_deficit = max(
                    min_clearance - (candidate_signature.min_edge_length or 0.0),
                    0.0,
                )
                score = (
                    candidate_signature.clearance_deficit,
                    candidate_signature.short_edge_count,
                    candidate_signature.ring_contact_count,
                    edge_deficit,
                    difference_metrics["symmetric_difference_area"],
                    abs(difference_metrics["union_area_delta"]),
                    candidate_signature.vertex_count,
                )
                if best is None or score < best[0]:
                    best = (
                        score,
                        candidate_polygon,
                        patch,
                        "self_clearance_connector_segment",
                    )

    if best is None:
        return None

    area_budget_override = max(float(best[2].area), 16.0 * grid * grid, 1e-9)
    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator=best[3],
        area_balance_budget_override=area_budget_override,
    )


def _try_cross_ring_self_clearance_connector_cut(
    polygon: Polygon,
    clearance_coords: np.ndarray,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
    start_kind: str,
    start_ring_index: int | None,
    end_kind: str,
    end_ring_index: int | None,
    require_signature_improvement: bool = True,
) -> _RepairCandidate | None:
    if (
        start_kind == end_kind
        and start_ring_index == end_ring_index
    ) or (start_kind == "exterior" and end_kind == "exterior"):
        return None

    segment_start = np.asarray(clearance_coords[0], dtype=float)
    segment_end = np.asarray(clearance_coords[-1], dtype=float)
    segment_vector = segment_end - segment_start
    segment_length = float(np.hypot(segment_vector[0], segment_vector[1]))
    if segment_length <= 1e-9:
        return None
    segment_direction = segment_vector / segment_length

    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=min_clearance,
    )
    best: tuple[
        tuple[float, int, int, float, int, float, float, int],
        Polygon,
        BaseGeometry,
    ] | None = None

    for half_width in _iter_self_clearance_connector_half_widths(
        target_scale=min_clearance,
        grid=grid,
    ):
        extension = max(float(half_width), float(grid), 1e-9)
        cutter_axis = LineString(
            [
                tuple(segment_start - segment_direction * extension),
                tuple(segment_end + segment_direction * extension),
            ]
        )
        cutter = cutter_axis.buffer(
            half_width,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        )
        if cutter.is_empty:
            continue

        candidate_polygon = _normalize_single_polygon_candidate(
            polygon.difference(cutter),
            grid=grid,
            min_area=0.0,
            min_hole_area=0.0,
            diagnostics=diagnostics,
        )
        if candidate_polygon is None:
            continue

        candidate_signature = _polygon_defect_signature(
            candidate_polygon,
            target_scale=min_clearance,
        )
        if not _signature_improves(
            reference_signature,
            candidate_signature,
            grid=grid,
        ):
            continue

        difference_metrics = _difference_area_metrics(
            polygon,
            candidate_polygon,
        )
        edge_deficit = max(
            min_clearance - (candidate_signature.min_edge_length or 0.0),
            0.0,
        )
        score = (
            candidate_signature.clearance_deficit,
            candidate_signature.short_edge_count,
            candidate_signature.ring_contact_count,
            edge_deficit,
            len(candidate_polygon.interiors),
            difference_metrics["symmetric_difference_area"],
            abs(difference_metrics["union_area_delta"]),
            candidate_signature.vertex_count,
        )
        if best is None or score < best[0]:
            best = (
                score,
                candidate_polygon,
                cutter,
            )

    if best is None:
        return None

    area_budget_override = max(float(best[2].area), 16.0 * grid * grid, 1e-9)
    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator="self_clearance_connector_cut",
        area_balance_budget_override=area_budget_override,
    )


def _try_exterior_hole_self_clearance_connector_fill(
    polygon: Polygon,
    *,
    exterior_point: Point,
    hole_point: Point,
    hole_ring_index: int,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
    require_signature_improvement: bool = True,
) -> _RepairCandidate | None:
    if hole_ring_index < 0 or hole_ring_index >= len(polygon.interiors):
        return None

    try:
        shell = Polygon(list(polygon.exterior.coords))
        hole_polygon = orient(
            Polygon(list(polygon.interiors[hole_ring_index].coords)),
            sign=1.0,
        )
    except (GEOSException, ValueError):
        return None
    if shell.is_empty or hole_polygon.is_empty:
        return None

    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=min_clearance,
    )
    hole_coords = np.asarray(hole_polygon.exterior.coords[:-1], dtype=float)
    best: tuple[
        tuple[float, int, int, float, float, float, int],
        Polygon,
        BaseGeometry,
    ] | None = None
    tolerance = max(min_clearance * 1.0e-3, grid, 1.0e-9)

    def consider_candidate(
        updated_hole_coords: np.ndarray,
        edit_zone: BaseGeometry,
    ) -> None:
        nonlocal best
        updated_ring = [
            (float(x), float(y))
            for x, y in np.vstack([updated_hole_coords, updated_hole_coords[0]])
        ]
        holes: list[Sequence[Sequence[float]]] = []
        for candidate_hole_index, ring in enumerate(polygon.interiors):
            if candidate_hole_index == hole_ring_index:
                holes.append(updated_ring)
            else:
                holes.append(list(ring.coords))
        candidate_polygon = _normalize_single_polygon_candidate(
            Polygon(list(polygon.exterior.coords), holes),
            grid=grid,
            min_area=0.0,
            min_hole_area=0.0,
            diagnostics=diagnostics,
        )
        if candidate_polygon is None:
            return

        candidate_signature = _polygon_defect_signature(
            candidate_polygon,
            target_scale=min_clearance,
        )
        if require_signature_improvement and not _signature_improves(
            reference_signature,
            candidate_signature,
            grid=grid,
        ):
            return

        difference_metrics = _difference_area_metrics(
            polygon,
            candidate_polygon,
        )
        candidate_edit_zone = polygon.symmetric_difference(candidate_polygon)
        if not candidate_edit_zone.is_empty:
            edit_zone = unary_union([edit_zone, candidate_edit_zone])
        edge_deficit = max(
            min_clearance - (candidate_signature.min_edge_length or 0.0),
            0.0,
        )
        score = (
            candidate_signature.clearance_deficit,
            candidate_signature.short_edge_count,
            candidate_signature.ring_contact_count,
            edge_deficit,
            difference_metrics["symmetric_difference_area"],
            abs(difference_metrics["union_area_delta"]),
            candidate_signature.vertex_count,
        )
        if best is None or score < best[0]:
            best = (
                score,
                candidate_polygon,
                edit_zone,
            )

    exterior_direction = _polygon_inward_direction_at_contact(
        polygon,
        exterior_point,
        tolerance=tolerance,
    )
    hole_direction = _polygon_inward_direction_at_contact(
        hole_polygon,
        hole_point,
        tolerance=tolerance,
    )
    if exterior_direction is None or hole_direction is None:
        return None

    hole_vertex_index, hole_vertex_distance = _nearest_vertex_index_and_distance(
        hole_coords,
        np.asarray(hole_point.coords[0], dtype=float),
    )
    hole_segment = _nearest_nonadjacent_ring_segment(
        hole_coords,
        np.asarray(hole_point.coords[0], dtype=float),
        excluded_segments=set(),
    )

    for distance in _iter_ring_contact_fill_distances(
        min_clearance=min_clearance,
        grid=grid,
    ):
        if hole_segment is not None:
            segment_index, _, _, projected_xy = hole_segment
            if float(np.hypot(*(projected_xy - np.asarray(hole_point.coords[0], dtype=float)))) <= tolerance:
                inward = _boundary_edge_inward_direction(
                    hole_coords,
                    projected_xy,
                    tolerance=max(tolerance, 1.0e-6),
                )
                if inward is not None:
                    updated_coords = hole_coords.copy()
                    updated_coords[segment_index] = updated_coords[segment_index] + inward * distance
                    updated_coords[(segment_index + 1) % len(updated_coords)] = (
                        updated_coords[(segment_index + 1) % len(updated_coords)] + inward * distance
                    )
                    edit_zone = LineString(
                        [
                            tuple(hole_coords[segment_index]),
                            tuple(projected_xy),
                            tuple(hole_coords[(segment_index + 1) % len(hole_coords)]),
                        ]
                    ).buffer(
                        max(distance, grid),
                        quad_segs=1,
                        join_style=BufferJoinStyle.mitre,
                        mitre_limit=1000.0,
                    )
                    consider_candidate(updated_coords, edit_zone)

        if hole_vertex_distance <= tolerance:
            updated_coords = hole_coords.copy()
            updated_coords[hole_vertex_index] = (
                updated_coords[hole_vertex_index] + hole_direction * distance
            )
            edit_zone = LineString(
                [
                    tuple(hole_coords[hole_vertex_index]),
                    tuple(updated_coords[hole_vertex_index]),
                ]
            ).buffer(
                max(distance, grid),
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
            consider_candidate(updated_coords, edit_zone)

            for segment_index in (
                (hole_vertex_index - 1) % len(hole_coords),
                hole_vertex_index % len(hole_coords),
            ):
                segment_start = hole_coords[segment_index]
                segment_end = hole_coords[(segment_index + 1) % len(hole_coords)]
                segment_midpoint = 0.5 * (segment_start + segment_end)
                segment_inward = _boundary_edge_inward_direction(
                    hole_coords,
                    segment_midpoint,
                    tolerance=max(tolerance, 1.0e-6),
                )
                if segment_inward is None:
                    continue
                updated_coords = hole_coords.copy()
                updated_coords[segment_index] = (
                    updated_coords[segment_index] + segment_inward * distance
                )
                updated_coords[(segment_index + 1) % len(updated_coords)] = (
                    updated_coords[(segment_index + 1) % len(updated_coords)]
                    + segment_inward * distance
                )
                edit_zone = LineString(
                    [
                        tuple(segment_start),
                        tuple(segment_end),
                    ]
                ).buffer(
                    max(distance, grid),
                    quad_segs=1,
                    join_style=BufferJoinStyle.mitre,
                    mitre_limit=1000.0,
                )
                consider_candidate(updated_coords, edit_zone)

    for distance in _iter_ring_contact_connector_distances(
        min_clearance=min_clearance,
        grid=grid,
    ):
        exterior_support = _inset_contact_support_point(
            polygon,
            exterior_point,
            exterior_direction,
            distance=distance,
            tolerance=grid,
        )
        hole_support = _inset_contact_support_point(
            hole_polygon,
            hole_point,
            hole_direction,
            distance=distance,
            tolerance=grid,
        )
        if exterior_support is None or hole_support is None:
            continue

        for width_factor in _iter_ring_contact_connector_width_factors():
            corridor_half_width = _connector_half_width(
                radius=distance,
                grid=grid,
                width_factor=width_factor,
            )
            try:
                patch = LineString(
                    [tuple(exterior_support), tuple(hole_support)]
                ).buffer(
                    corridor_half_width,
                    quad_segs=1,
                    cap_style=BufferCapStyle.flat,
                    join_style=BufferJoinStyle.mitre,
                    mitre_limit=1000.0,
                )
            except GEOSException as exc:
                _record_geos_exception(
                    diagnostics,
                    "self_clearance_connector_fill",
                    exc,
                )
                continue
            if patch.is_empty:
                continue

            patch = patch.intersection(shell)
            if patch.is_empty:
                continue

            candidate_polygon = _normalize_single_polygon_candidate(
                unary_union([polygon, patch]),
                grid=grid,
                min_area=0.0,
                min_hole_area=0.0,
                diagnostics=diagnostics,
            )
            if candidate_polygon is None:
                continue

            candidate_signature = _polygon_defect_signature(
                candidate_polygon,
                target_scale=min_clearance,
            )
            if require_signature_improvement and not _signature_improves(
                reference_signature,
                candidate_signature,
                grid=grid,
            ):
                continue

            difference_metrics = _difference_area_metrics(
                polygon,
                candidate_polygon,
            )
            edge_deficit = max(
                min_clearance - (candidate_signature.min_edge_length or 0.0),
                0.0,
            )
            score = (
                candidate_signature.clearance_deficit,
                candidate_signature.short_edge_count,
                candidate_signature.ring_contact_count,
                edge_deficit,
                difference_metrics["symmetric_difference_area"],
                abs(difference_metrics["union_area_delta"]),
                candidate_signature.vertex_count,
            )
            if best is None or score < best[0]:
                best = (
                    score,
                    candidate_polygon,
                    patch,
                )

    if best is None:
        return None

    area_budget_override = max(float(best[2].area), 16.0 * grid * grid, 1e-9)
    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator="self_clearance_connector_fill",
        area_balance_budget_override=area_budget_override,
    )


def _try_polygon_self_clearance_connector_fill(
    polygon: Polygon,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
    require_signature_improvement: bool = True,
) -> _RepairCandidate | None:
    if min_clearance <= 0:
        return None

    try:
        clearance_line = shapely.minimum_clearance_line(polygon)
    except GEOSException as exc:
        _record_geos_exception(diagnostics, "minimum_clearance_line", exc)
        return None
    if clearance_line is None or clearance_line.is_empty:
        return None

    clearance_coords = np.asarray(clearance_line.coords, dtype=float)
    if len(clearance_coords) < 2:
        return None

    start_point = Point(float(clearance_coords[0][0]), float(clearance_coords[0][1]))
    end_point = Point(float(clearance_coords[-1][0]), float(clearance_coords[-1][1]))
    start_kind, start_ring_index, _ = _nearest_ring_classification(
        polygon,
        start_point,
    )
    end_kind, end_ring_index, _ = _nearest_ring_classification(
        polygon,
        end_point,
    )
    if start_kind is None or end_kind is None:
        return None
    if start_kind == end_kind and start_ring_index == end_ring_index:
        return _try_same_ring_self_clearance_connector_fill(
            polygon,
            clearance_coords,
            min_clearance=min_clearance,
            grid=grid,
            diagnostics=diagnostics,
            ring_kind=start_kind,
            ring_index=start_ring_index,
            require_signature_improvement=require_signature_improvement,
        )
    if (
        start_kind == "hole"
        and end_kind == "hole"
        and start_ring_index is not None
        and end_ring_index is not None
    ):
        hole_pair_candidate = _try_hole_pair_self_clearance_merge(
            polygon,
            start_ring_index=start_ring_index,
            end_ring_index=end_ring_index,
            min_clearance=min_clearance,
            current_clearance=float(np.linalg.norm(clearance_coords[-1] - clearance_coords[0])),
            grid=grid,
            diagnostics=diagnostics,
            require_signature_improvement=require_signature_improvement,
        )
        if hole_pair_candidate is not None:
            return hole_pair_candidate
    if {start_kind, end_kind} == {"exterior", "hole"}:
        exterior_point = start_point if start_kind == "exterior" else end_point
        hole_point = start_point if start_kind == "hole" else end_point
        hole_ring_index = (
            start_ring_index if start_kind == "hole" else end_ring_index
        )
        if hole_ring_index is not None:
            return _try_exterior_hole_self_clearance_connector_fill(
                polygon,
                exterior_point=exterior_point,
                hole_point=hole_point,
                hole_ring_index=hole_ring_index,
                min_clearance=min_clearance,
                grid=grid,
                diagnostics=diagnostics,
                require_signature_improvement=require_signature_improvement,
            )
    return _try_cross_ring_self_clearance_connector_cut(
        polygon,
        clearance_coords,
        min_clearance=min_clearance,
        grid=grid,
        diagnostics=diagnostics,
        start_kind=start_kind,
        start_ring_index=start_ring_index,
        end_kind=end_kind,
        end_ring_index=end_ring_index,
        require_signature_improvement=require_signature_improvement,
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


def _hole_contact_vertex_indices(
    polygon: Polygon,
    *,
    tolerance: float,
) -> dict[int, set[int]]:
    if not polygon.interiors:
        return {}

    exterior = LineString(polygon.exterior.coords)
    hole_lines = [LineString(ring.coords) for ring in polygon.interiors]
    contact_indices: dict[int, set[int]] = {}

    for hole_index, ring in enumerate(polygon.interiors):
        coords = np.asarray(ring.coords[:-1], dtype=float)
        if len(coords) < 3:
            continue
        for vertex_index, coord in enumerate(coords):
            point = Point(float(coord[0]), float(coord[1]))
            if exterior.distance(point) <= tolerance:
                contact_indices.setdefault(hole_index, set()).add(vertex_index)
                continue
            for other_index, other in enumerate(hole_lines):
                if other_index == hole_index:
                    continue
                if other.distance(point) <= tolerance:
                    contact_indices.setdefault(hole_index, set()).add(vertex_index)
                    break

    return contact_indices


def _hole_contact_inward_direction(
    hole_coords: np.ndarray,
    index: int,
) -> np.ndarray | None:
    if len(hole_coords) < 3:
        return None

    contact = hole_coords[index]
    centroid = np.asarray(Polygon(hole_coords).centroid.coords[0], dtype=float)
    direction = centroid - contact
    norm = float(np.hypot(direction[0], direction[1]))
    if norm > 1e-9:
        return direction / norm

    midpoint = 0.5 * (
        hole_coords[index - 1] + hole_coords[(index + 1) % len(hole_coords)]
    )
    direction = midpoint - contact
    norm = float(np.hypot(direction[0], direction[1]))
    if norm > 1e-9:
        return direction / norm

    return None


def _ring_contact_groups(
    polygon: Polygon,
    *,
    tolerance: float,
) -> tuple[dict[str, Any], ...]:
    contact_vertices = _hole_contact_vertex_indices(
        polygon,
        tolerance=tolerance,
    )
    if not contact_vertices:
        return ()

    exterior = LineString(polygon.exterior.coords)
    groups: dict[tuple[int, int], dict[str, Any]] = {}
    scale = max(tolerance, 1e-9)

    for hole_index, vertex_indices in contact_vertices.items():
        coords = np.asarray(polygon.interiors[hole_index].coords[:-1], dtype=float)
        for vertex_index in sorted(vertex_indices):
            coord = coords[vertex_index]
            point = Point(float(coord[0]), float(coord[1]))
            key = (
                int(round(float(coord[0]) / scale)),
                int(round(float(coord[1]) / scale)),
            )
            group = groups.setdefault(
                key,
                {
                    "point": point,
                    "touches_exterior": False,
                    "contacts": [],
                },
            )
            if exterior.distance(point) <= tolerance:
                group["touches_exterior"] = True
            contact = (int(hole_index), int(vertex_index))
            if contact not in group["contacts"]:
                group["contacts"].append(contact)

    return tuple(groups[key] for key in sorted(groups))


def _iter_ring_contact_connector_distances(
    *,
    min_clearance: float,
    grid: float,
) -> tuple[float, ...]:
    base = max(float(min_clearance), 4.0 * float(grid), 1e-6)
    distances: list[float] = []
    for factor in (1.0, 1.25, 1.5, 1.75, 2.0):
        distance = max(base * factor, 4.0 * float(grid))
        if any(abs(distance - existing) <= 1e-12 for existing in distances):
            continue
        distances.append(distance)
    return tuple(distances)


def _iter_ring_contact_connector_width_factors() -> tuple[float, ...]:
    return (0.55, 0.7, 0.85, 1.0)


def _try_ring_contact_connector_fill(
    polygon: Polygon,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    tolerance = max(grid, min_clearance * 1e-6, 1e-9)
    contact_groups = _ring_contact_groups(
        polygon,
        tolerance=tolerance,
    )
    if not contact_groups:
        return None

    shell = Polygon(list(polygon.exterior.coords))
    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=min_clearance,
    )
    best: tuple[
        tuple[int, int, float, float, float, float, int],
        Polygon,
        BaseGeometry,
    ] | None = None

    for distance in _iter_ring_contact_connector_distances(
        min_clearance=min_clearance,
        grid=grid,
    ):
        for width_factor in _iter_ring_contact_connector_width_factors():
            corridor_half_width = _connector_half_width(
                radius=distance,
                grid=grid,
                width_factor=width_factor,
            )
            patches: list[BaseGeometry] = []
            for group in contact_groups:
                contact = group["point"]
                support_points: list[tuple[float, float]] = []
                seen_support_keys: set[tuple[int, int]] = set()
                for hole_index, vertex_index in group["contacts"]:
                    hole_coords = np.asarray(
                        polygon.interiors[hole_index].coords[:-1],
                        dtype=float,
                    )
                    direction = _hole_contact_inward_direction(
                        hole_coords,
                        vertex_index,
                    )
                    if direction is None:
                        continue
                    support = _inset_contact_support_point(
                        Polygon(hole_coords),
                        contact,
                        direction,
                        distance=distance,
                        tolerance=grid,
                    )
                    if support is None:
                        continue
                    support_key = (
                        int(round(float(support[0]) / max(grid, 1e-9))),
                        int(round(float(support[1]) / max(grid, 1e-9))),
                    )
                    if support_key in seen_support_keys:
                        continue
                    seen_support_keys.add(support_key)
                    support_points.append((float(support[0]), float(support[1])))

                if group["touches_exterior"]:
                    for support in support_points:
                        try:
                            patch = LineString(
                                [contact.coords[0], support]
                            ).buffer(
                                corridor_half_width,
                                quad_segs=1,
                                cap_style=BufferCapStyle.flat,
                                join_style=BufferJoinStyle.mitre,
                                mitre_limit=1000.0,
                            )
                        except GEOSException as exc:
                            _record_geos_exception(
                                diagnostics,
                                "ring_contact_connector",
                                exc,
                            )
                            continue
                        if not patch.is_empty:
                            patches.append(patch)
                    continue

                for left_index, left_support in enumerate(support_points):
                    for right_support in support_points[left_index + 1 :]:
                        try:
                            patch = LineString(
                                [left_support, right_support]
                            ).buffer(
                                corridor_half_width,
                                quad_segs=1,
                                cap_style=BufferCapStyle.flat,
                                join_style=BufferJoinStyle.mitre,
                                mitre_limit=1000.0,
                            )
                        except GEOSException as exc:
                            _record_geos_exception(
                                diagnostics,
                                "ring_contact_connector",
                                exc,
                            )
                            continue
                        if not patch.is_empty:
                            patches.append(patch)

            if not patches:
                continue

            patch_geometry = unary_union(patches).intersection(shell)
            if patch_geometry.is_empty:
                continue

            candidate_polygon = _normalize_single_polygon_candidate(
                unary_union([polygon, patch_geometry]),
                grid=grid,
                min_area=0.0,
                min_hole_area=0.0,
                diagnostics=diagnostics,
            )
            if candidate_polygon is None:
                continue

            candidate_signature = _polygon_defect_signature(
                candidate_polygon,
                target_scale=min_clearance,
            )
            if candidate_signature.ring_contact_count >= reference_signature.ring_contact_count:
                continue

            difference_metrics = _difference_area_metrics(
                polygon,
                candidate_polygon,
            )
            score = (
                candidate_signature.ring_contact_count,
                candidate_signature.short_edge_count,
                candidate_signature.clearance_deficit,
                -(candidate_signature.min_edge_length or 0.0),
                difference_metrics["symmetric_difference_area"],
                abs(difference_metrics["union_area_delta"]),
                candidate_signature.vertex_count,
            )
            if best is None or score < best[0]:
                best = (
                    score,
                    candidate_polygon,
                    patch_geometry,
                )

    if best is None:
        return None

    area_budget_override = max(float(best[2].area), 16.0 * grid * grid, 1e-9)
    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator="ring_contact_connector",
        area_balance_budget_override=area_budget_override,
    )


def _iter_ring_contact_fill_distances(
    *,
    min_clearance: float,
    grid: float,
) -> tuple[float, ...]:
    base = max(float(min_clearance), 2.0 * float(grid), 1e-6)
    distances = [
        base,
        max(base + 2.0 * float(grid), 1.5 * base),
        max(base + 4.0 * float(grid), 2.0 * base),
    ]
    ordered: list[float] = []
    seen: set[float] = set()
    for value in distances:
        candidate = float(value)
        if candidate in seen:
            continue
        seen.add(candidate)
        ordered.append(candidate)
    return tuple(ordered)


def _try_fill_ring_contact_vertices(
    polygon: Polygon,
    *,
    min_clearance: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    tolerance = max(grid, min_clearance * 1e-6, 1e-9)
    contact_vertices = _hole_contact_vertex_indices(
        polygon,
        tolerance=tolerance,
    )
    if not contact_vertices:
        return None

    best: tuple[tuple[float, ...], Polygon, BaseGeometry, float] | None = None
    for distance in _iter_ring_contact_fill_distances(
        min_clearance=min_clearance,
        grid=grid,
    ):
        moved_holes: list[list[tuple[float, float]]] = []
        edit_geometries: list[BaseGeometry] = []
        valid = True

        for hole_index, ring in enumerate(polygon.interiors):
            coords = np.asarray(ring.coords[:-1], dtype=float)
            moved = coords.copy()
            for vertex_index in sorted(contact_vertices.get(hole_index, ())):
                direction = _hole_contact_inward_direction(moved, vertex_index)
                if direction is None:
                    valid = False
                    break
                source = moved[vertex_index].copy()
                moved[vertex_index] = source + direction * distance
                edit_geometries.append(
                    LineString(
                        [
                            (float(source[0]), float(source[1])),
                            (
                                float(moved[vertex_index][0]),
                                float(moved[vertex_index][1]),
                            ),
                        ]
                    ).buffer(
                        max(distance, grid),
                        quad_segs=1,
                        join_style=BufferJoinStyle.mitre,
                        mitre_limit=1000.0,
                    )
                )
            if not valid:
                break
            moved_holes.append(
                [(float(x), float(y)) for x, y in np.vstack([moved, moved[0]])]
            )

        if not valid:
            continue

        try:
            candidate_geometry = make_valid(
                Polygon(list(polygon.exterior.coords), moved_holes)
            )
        except (GEOSException, ValueError) as exc:
            _record_geos_exception(diagnostics, "ring_contact_fill", exc)
            continue

        candidate_parts = _extract_polygon_parts(candidate_geometry)
        if len(candidate_parts) != 1:
            continue

        candidate_polygon = orient(candidate_parts[0], sign=1.0)
        candidate_signature = _polygon_defect_signature(
            candidate_polygon,
            target_scale=min_clearance,
        )
        if not _signature_satisfies_scale_contract(
            candidate_signature,
            target_scale=min_clearance,
            grid=grid,
        ):
            continue

        difference_metrics = _difference_area_metrics(polygon, candidate_polygon)
        score = (
            difference_metrics["symmetric_difference_area"],
            difference_metrics["candidate_minus_reference_area"],
            float(candidate_signature.vertex_count),
        )
        edit_zone = _as_geometry(edit_geometries) if edit_geometries else Point()
        area_budget_override = max(
            8.0 * len(edit_geometries) * min_clearance * min_clearance,
            16.0 * grid * grid,
        )
        if best is None or score < best[0]:
            best = (score, candidate_polygon, edit_zone, area_budget_override)

    if best is None:
        return None

    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator="ring_contact_fill",
        area_balance_budget_override=best[3],
    )


def _iter_short_edge_angle_open_amounts(
    *,
    target_scale: float,
    grid: float,
) -> tuple[float, ...]:
    if target_scale <= 0:
        return ()

    amounts: list[float] = []
    for factor in (0.25, 0.5, 0.75, 1.0, 1.25):
        amount = max(float(target_scale) * factor, 4.0 * float(grid), 1e-9)
        if any(abs(amount - existing) <= 1e-12 for existing in amounts):
            continue
        amounts.append(amount)
    return tuple(amounts)


def _open_short_edge_ring_coords(
    coords: Sequence[Sequence[float]],
    *,
    segment_index: int,
    amount: float,
) -> tuple[list[tuple[float, float]], BaseGeometry] | None:
    if amount <= 0:
        return None

    unique = [tuple(map(float, coord[:2])) for coord in coords]
    if len(unique) >= 2 and unique[0] == unique[-1]:
        unique.pop()
    if len(unique) < 4:
        return None

    index = segment_index % len(unique)
    updated = np.asarray(unique, dtype=float).copy()
    previous = updated[index - 1].copy()
    start = updated[index].copy()
    end = updated[(index + 1) % len(updated)].copy()
    next_coord = updated[(index + 2) % len(updated)].copy()

    start_direction = previous - start
    end_direction = next_coord - end
    start_norm = float(np.hypot(start_direction[0], start_direction[1]))
    end_norm = float(np.hypot(end_direction[0], end_direction[1]))
    if start_norm <= 1e-9 or end_norm <= 1e-9:
        return None

    start_move = min(float(amount), 0.45 * start_norm)
    end_move = min(float(amount), 0.45 * end_norm)
    if start_move <= 0.0 or end_move <= 0.0:
        return None

    updated[index] = start + (start_direction / start_norm) * start_move
    updated[(index + 1) % len(updated)] = end + (end_direction / end_norm) * end_move

    edit_zone = unary_union(
        [
            LineString([tuple(start), tuple(updated[index])]),
            LineString(
                [
                    tuple(end),
                    tuple(updated[(index + 1) % len(updated)]),
                ]
            ),
            LineString([tuple(start), tuple(end)]),
        ]
    ).buffer(
        max(float(amount), 1e-9),
        quad_segs=1,
        join_style=BufferJoinStyle.mitre,
        mitre_limit=1000.0,
    )
    ring = [(float(x), float(y)) for x, y in updated]
    ring.append(ring[0])
    return ring, edit_zone


def _try_polygon_short_edge_angle_open(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=target_scale,
    )
    if reference_signature.short_edge_count == 0:
        return None

    best: tuple[
        tuple[int, float, float, float, float, int],
        Polygon,
        BaseGeometry,
        float,
    ] | None = None
    rings: list[tuple[str, int | None, Sequence[Sequence[float]]]] = [
        ("exterior", None, polygon.exterior.coords),
    ]
    rings.extend(
        ("hole", hole_index, ring.coords)
        for hole_index, ring in enumerate(polygon.interiors)
    )

    for ring_kind, hole_index, ring_coords in rings:
        coords = list(ring_coords)
        for segment_index, (start, end) in enumerate(zip(coords, coords[1:])):
            length = float(np.hypot(float(end[0] - start[0]), float(end[1] - start[1])))
            if length + 1e-12 >= target_scale:
                continue

            for amount in _iter_short_edge_angle_open_amounts(
                target_scale=target_scale,
                grid=grid,
            ):
                updated = _open_short_edge_ring_coords(
                    ring_coords,
                    segment_index=segment_index,
                    amount=amount,
                )
                if updated is None:
                    continue
                updated_ring, edit_zone = updated
                if ring_kind == "exterior":
                    candidate_geometry = Polygon(
                        updated_ring,
                        [list(ring.coords) for ring in polygon.interiors],
                    )
                else:
                    holes: list[Sequence[Sequence[float]]] = []
                    for candidate_hole_index, ring in enumerate(polygon.interiors):
                        if candidate_hole_index == hole_index:
                            holes.append(updated_ring)
                        else:
                            holes.append(list(ring.coords))
                    candidate_geometry = Polygon(list(polygon.exterior.coords), holes)

                candidate_polygon = _normalize_single_polygon_candidate(
                    candidate_geometry,
                    grid=grid,
                    min_area=0.0,
                    min_hole_area=0.0,
                    diagnostics=diagnostics,
                )
                if candidate_polygon is None:
                    continue

                candidate_signature = _polygon_defect_signature(
                    candidate_polygon,
                    target_scale=target_scale,
                )
                if not _signature_improves(
                    reference_signature,
                    candidate_signature,
                    grid=grid,
                ):
                    continue

                difference_metrics = _difference_area_metrics(
                    polygon,
                    candidate_polygon,
                )
                score = (
                    candidate_signature.short_edge_count,
                    candidate_signature.clearance_deficit,
                    -(candidate_signature.min_edge_length or 0.0),
                    difference_metrics["symmetric_difference_area"],
                    abs(difference_metrics["union_area_delta"]),
                    candidate_signature.vertex_count,
                )
                if best is None or score < best[0]:
                    best = (
                        score,
                        candidate_polygon,
                        edit_zone,
                        amount,
                    )

    if best is None:
        return None

    return _RepairCandidate(
        polygon=best[1],
        edit_zone=best[2],
        operator=f"short_edge_angle_open_{best[3]:.3f}",
        area_balance_budget_override=max(float(best[2].area), 16.0 * grid * grid, 1e-9),
    )


def _iteratively_open_polygon_short_edges(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
    operator_prefix: str,
    initial_edit_zone: BaseGeometry | None = None,
    area_balance_budget_override: float | None = None,
    max_steps: int = 4,
) -> _RepairCandidate | None:
    current_polygon = polygon
    current_signature = _polygon_defect_signature(
        current_polygon,
        target_scale=target_scale,
    )
    if current_signature.short_edge_count == 0:
        return None

    edit_zones: list[BaseGeometry] = []
    if initial_edit_zone is not None and not initial_edit_zone.is_empty:
        edit_zones.append(initial_edit_zone)

    applied = 0
    while applied < max_steps and current_signature.short_edge_count > 0:
        candidate = _try_polygon_short_edge_angle_open(
            current_polygon,
            target_scale=target_scale,
            grid=grid,
            diagnostics=diagnostics,
        )
        if candidate is None:
            break
        candidate_signature = _polygon_defect_signature(
            candidate.polygon,
            target_scale=target_scale,
        )
        if not _signature_improves(
            current_signature,
            candidate_signature,
            grid=grid,
        ):
            break

        current_polygon = candidate.polygon
        current_signature = candidate_signature
        edit_zones.append(candidate.edit_zone)
        applied += 1

    if applied == 0:
        return None

    edit_zone = (
        unary_union(edit_zones)
        if edit_zones
        else GeometryCollection()
    )
    area_budget = max(float(edit_zone.area), 16.0 * grid * grid, 1e-9)
    if area_balance_budget_override is not None:
        area_budget = max(area_budget, float(area_balance_budget_override))
    return _RepairCandidate(
        polygon=current_polygon,
        edit_zone=edit_zone,
        operator=f"{operator_prefix}_chain",
        area_balance_budget_override=area_budget,
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


def _try_polygon_micro_detour_chain_simplify(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    coords = list(polygon.exterior.coords[:-1])
    count = len(coords)
    if count < 4:
        return None

    mouth_width_limit = max(3.0 * target_scale, 8.0 * grid, 1.0e-9)
    area_limit = max(128.0 * target_scale * target_scale, 128.0 * grid * grid, 1.0e-9)
    compactness_threshold = 0.15

    best_candidate: tuple[tuple[float, float, float], Polygon, BaseGeometry] | None = None

    for index in range(count):
        a = coords[(index - 1) % count]
        b = coords[index]
        c = coords[(index + 1) % count]
        d = coords[(index + 2) % count]

        mouth_width = float(np.hypot(d[0] - a[0], d[1] - a[1]))
        if mouth_width <= max(grid, 1.0e-9) or mouth_width > mouth_width_limit:
            continue

        cross_b = ((d[0] - a[0]) * (b[1] - a[1])) - ((d[1] - a[1]) * (b[0] - a[0]))
        cross_c = ((d[0] - a[0]) * (c[1] - a[1])) - ((d[1] - a[1]) * (c[0] - a[0]))
        if abs(cross_b) <= 1.0e-9 or abs(cross_c) <= 1.0e-9:
            continue
        if cross_b * cross_c <= 0.0:
            continue

        notch_polygon = Polygon([a, b, c, d])
        local_area = notch_polygon.area
        if local_area > area_limit:
            continue
        if _compactness_ratio(notch_polygon) >= compactness_threshold:
            continue

        candidate_shell = [point for j, point in enumerate(coords) if j not in {index, (index + 1) % count}]
        if len(candidate_shell) < 3:
            continue

        try:
            candidate_polygon = Polygon(
                candidate_shell,
                [list(ring.coords) for ring in polygon.interiors],
            )
        except (GEOSException, ValueError) as exc:
            _record_geos_exception(diagnostics, "micro_detour_chain", exc)
            continue

        if (
            candidate_polygon.is_empty
            or candidate_polygon.area <= 0.0
            or not candidate_polygon.is_valid
        ):
            continue

        edit_zone = LineString([a, b, c, d]).buffer(
            max(target_scale, grid),
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        )
        score = (
            _compactness_ratio(notch_polygon),
            mouth_width,
            abs(float(candidate_polygon.area - polygon.area)),
        )
        if best_candidate is None or score < best_candidate[0]:
            best_candidate = (score, orient(candidate_polygon, sign=1.0), edit_zone)

    if best_candidate is None:
        return None

    return _RepairCandidate(
        polygon=best_candidate[1],
        edit_zone=best_candidate[2],
        operator="micro_detour_chain",
    )


def _try_polygon_same_turn_short_walk_collapse(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    coords = list(polygon.exterior.coords[:-1])
    count = len(coords)
    if count < 4:
        return None

    short_walk_limit = max(12.0 * target_scale, 16.0 * grid, 1.0e-9)
    missing_tolerance = max(16.0 * grid * grid, 1.0e-9)
    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=target_scale,
    )

    best_candidate: tuple[tuple[float, ...], Polygon, BaseGeometry] | None = None

    for index in range(count):
        a = coords[index]
        b = coords[(index + 1) % count]
        c = coords[(index + 2) % count]
        d = coords[(index + 3) % count]

        cross_ab_bc = ((b[0] - a[0]) * (c[1] - b[1])) - ((b[1] - a[1]) * (c[0] - b[0]))
        cross_bc_cd = ((c[0] - b[0]) * (d[1] - c[1])) - ((c[1] - b[1]) * (d[0] - c[0]))
        if abs(cross_ab_bc) <= 1.0e-9 or abs(cross_bc_cd) <= 1.0e-9:
            continue
        if cross_ab_bc * cross_bc_cd <= 0.0:
            continue

        short_walk = float(np.hypot(c[0] - b[0], c[1] - b[1]))
        if short_walk <= max(grid, 1.0e-9) or short_walk > short_walk_limit:
            continue

        candidate_shell = [
            point
            for j, point in enumerate(coords)
            if j not in {(index + 1) % count, (index + 2) % count}
        ]
        if len(candidate_shell) < 3:
            continue

        try:
            candidate_polygon = Polygon(
                candidate_shell,
                [list(ring.coords) for ring in polygon.interiors],
            )
        except (GEOSException, ValueError) as exc:
            _record_geos_exception(diagnostics, "same_turn_short_walk", exc)
            continue

        if (
            candidate_polygon.is_empty
            or candidate_polygon.area <= 0.0
            or not candidate_polygon.is_valid
        ):
            continue

        candidate_signature = _polygon_defect_signature(
            candidate_polygon,
            target_scale=target_scale,
        )
        if not _signature_not_worse(
            reference_signature,
            candidate_signature,
            grid=grid,
        ):
            continue
        if (
            candidate_signature.clearance is None
            or reference_signature.clearance is None
            or candidate_signature.clearance + grid < reference_signature.clearance
        ):
            continue
        if (
            candidate_signature.min_edge_length is None
            or reference_signature.min_edge_length is None
            or candidate_signature.min_edge_length + grid
            < reference_signature.min_edge_length
        ):
            continue
        if candidate_signature.vertex_count >= reference_signature.vertex_count:
            continue

        difference_metrics = _difference_area_metrics(polygon, candidate_polygon)
        if difference_metrics["reference_minus_candidate_area"] > missing_tolerance:
            continue

        edit_zone_scale = 1.5 * max(short_walk, target_scale, grid, 1.0e-9)
        edit_zone = LineString([a, b, c, d]).buffer(
            edit_zone_scale,
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        )
        area_growth_budget = max(
            0.003 * float(polygon.area),
            0.5 * float(edit_zone.area),
            16.0 * grid * grid,
            1.0e-9,
        )
        if difference_metrics["candidate_minus_reference_area"] > area_growth_budget:
            continue
        score = (
            difference_metrics["candidate_minus_reference_area"],
            difference_metrics["symmetric_difference_area"],
            short_walk,
            -(candidate_signature.min_edge_length or 0.0),
            candidate_signature.vertex_count,
        )
        if best_candidate is None or score < best_candidate[0]:
            best_candidate = (score, orient(candidate_polygon, sign=1.0), edit_zone)

    if best_candidate is None:
        return None

    return _RepairCandidate(
        polygon=best_candidate[1],
        edit_zone=best_candidate[2],
        operator="same_turn_short_walk",
    )


def _try_polygon_fill_chain_collapse(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    coords = list(polygon.exterior.coords[:-1])
    count = len(coords)
    if count < 6:
        return None

    min_internal_vertices = 2
    max_internal_vertices = min(6, count - 3)
    missing_tolerance = max(16.0 * grid * grid, 1.0e-9)
    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=target_scale,
    )

    best_candidate: tuple[tuple[float, ...], Polygon, BaseGeometry] | None = None

    for start in range(count):
        for internal_vertices in range(min_internal_vertices, max_internal_vertices + 1):
            span = internal_vertices + 1
            if span >= count - 1:
                continue
            end = (start + span) % count
            internal_indices = {
                (start + step) % count for step in range(1, internal_vertices + 1)
            }
            if len(internal_indices) != internal_vertices:
                continue

            chain = [coords[(start + step) % count] for step in range(span + 1)]
            chain_length = sum(
                float(np.hypot(b[0] - a[0], b[1] - a[1]))
                for a, b in zip(chain, chain[1:])
            )
            chord_length = float(
                np.hypot(
                    chain[-1][0] - chain[0][0],
                    chain[-1][1] - chain[0][1],
                )
            )
            if chord_length <= max(grid, 1.0e-9):
                continue
            path_excess = chain_length - chord_length
            if path_excess <= max(target_scale, 8.0 * grid, 1.0e-9):
                continue

            candidate_shell = [
                point for j, point in enumerate(coords) if j not in internal_indices
            ]
            if len(candidate_shell) < 3:
                continue

            try:
                candidate_polygon = Polygon(
                    candidate_shell,
                    [list(ring.coords) for ring in polygon.interiors],
                )
            except (GEOSException, ValueError) as exc:
                _record_geos_exception(diagnostics, "fill_chain_collapse", exc)
                continue

            if (
                candidate_polygon.is_empty
                or candidate_polygon.area <= 0.0
                or not candidate_polygon.is_valid
            ):
                continue

            candidate_signature = _polygon_defect_signature(
                candidate_polygon,
                target_scale=target_scale,
            )
            if not _signature_not_worse(
                reference_signature,
                candidate_signature,
                grid=grid,
            ):
                continue
            if (
                candidate_signature.clearance is None
                or reference_signature.clearance is None
                or candidate_signature.clearance + grid < reference_signature.clearance
            ):
                continue
            if (
                candidate_signature.min_edge_length is None
                or reference_signature.min_edge_length is None
                or candidate_signature.min_edge_length + grid
                < reference_signature.min_edge_length
            ):
                continue
            if candidate_signature.vertex_count >= reference_signature.vertex_count:
                continue
            if (
                candidate_signature.clearance
                <= reference_signature.clearance + grid
                and candidate_signature.min_edge_length
                <= reference_signature.min_edge_length + grid
            ):
                continue

            difference_metrics = _difference_area_metrics(polygon, candidate_polygon)
            if difference_metrics["reference_minus_candidate_area"] > missing_tolerance:
                continue

            fill_area_budget = max(
                0.75 * chord_length * target_scale,
                32.0 * grid * grid,
                1.0e-9,
            )
            if difference_metrics["candidate_minus_reference_area"] > fill_area_budget:
                continue

            edit_zone = candidate_polygon.difference(polygon)
            if edit_zone.is_empty:
                edit_zone = Polygon(chain)

            clearance_gain = (candidate_signature.clearance or 0.0) - (
                reference_signature.clearance or 0.0
            )
            min_edge_gain = (candidate_signature.min_edge_length or 0.0) - (
                reference_signature.min_edge_length or 0.0
            )
            score = (
                -clearance_gain,
                -min_edge_gain,
                difference_metrics["candidate_minus_reference_area"],
                path_excess,
                candidate_signature.vertex_count,
            )
            if best_candidate is None or score < best_candidate[0]:
                best_candidate = (score, orient(candidate_polygon, sign=1.0), edit_zone)

    if best_candidate is None:
        return None

    return _RepairCandidate(
        polygon=best_candidate[1],
        edit_zone=best_candidate[2],
        operator="fill_chain_collapse",
    )


def _try_polygon_bevel_corner_collapse(
    polygon: Polygon,
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> _RepairCandidate | None:
    coords = list(polygon.exterior.coords[:-1])
    count = len(coords)
    if count < 4:
        return None

    bevel_length_limit = max(2.5 * target_scale, 16.0 * grid, 1.0e-9)
    missing_tolerance = max(16.0 * grid * grid, 1.0e-9)
    reference_signature = _polygon_defect_signature(
        polygon,
        target_scale=target_scale,
    )

    best_candidate: tuple[tuple[float, ...], Polygon, BaseGeometry] | None = None

    for index in range(count):
        a = coords[(index - 1) % count]
        b = coords[index]
        c = coords[(index + 1) % count]
        d = coords[(index + 2) % count]

        bevel_length = float(np.hypot(c[0] - b[0], c[1] - b[1]))
        if bevel_length <= max(grid, 1.0e-9) or bevel_length > bevel_length_limit:
            continue

        leg_ab = float(np.hypot(b[0] - a[0], b[1] - a[1]))
        leg_cd = float(np.hypot(d[0] - c[0], d[1] - c[1]))
        if leg_ab < max(2.0 * bevel_length, target_scale) or leg_cd < max(
            2.0 * bevel_length,
            target_scale,
        ):
            continue

        x1, y1 = a
        x2, y2 = b
        x3, y3 = c
        x4, y4 = d
        denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if abs(denominator) <= 1.0e-12:
            continue

        intersection_x = (
            ((x1 * y2) - (y1 * x2)) * (x3 - x4)
            - (x1 - x2) * ((x3 * y4) - (y3 * x4))
        ) / denominator
        intersection_y = (
            ((x1 * y2) - (y1 * x2)) * (y3 - y4)
            - (y1 - y2) * ((x3 * y4) - (y3 * x4))
        ) / denominator
        intersection = (intersection_x, intersection_y)

        if (
            float(np.hypot(intersection_x - b[0], intersection_y - b[1]))
            > max(3.0 * target_scale, 2.0 * bevel_length, 1.0e-9)
            or float(np.hypot(intersection_x - c[0], intersection_y - c[1]))
            > max(3.0 * target_scale, 2.0 * bevel_length, 1.0e-9)
        ):
            continue

        candidate_shell: list[tuple[float, float]] = []
        for j, point in enumerate(coords):
            if j == index:
                candidate_shell.append(intersection)
            elif j == (index + 1) % count:
                continue
            else:
                candidate_shell.append(point)
        if len(candidate_shell) < 3:
            continue

        try:
            candidate_polygon = Polygon(
                candidate_shell,
                [list(ring.coords) for ring in polygon.interiors],
            )
        except (GEOSException, ValueError) as exc:
            _record_geos_exception(diagnostics, "bevel_corner_collapse", exc)
            continue

        if (
            candidate_polygon.is_empty
            or candidate_polygon.area <= 0.0
            or not candidate_polygon.is_valid
        ):
            continue

        candidate_signature = _polygon_defect_signature(
            candidate_polygon,
            target_scale=target_scale,
        )
        if not _signature_not_worse(
            reference_signature,
            candidate_signature,
            grid=grid,
        ):
            continue
        if (
            candidate_signature.clearance is None
            or reference_signature.clearance is None
            or candidate_signature.clearance + grid < reference_signature.clearance
        ):
            continue
        if (
            candidate_signature.min_edge_length is None
            or reference_signature.min_edge_length is None
            or candidate_signature.min_edge_length + grid
            < reference_signature.min_edge_length
        ):
            continue
        if candidate_signature.vertex_count >= reference_signature.vertex_count:
            continue

        difference_metrics = _difference_area_metrics(polygon, candidate_polygon)
        if difference_metrics["reference_minus_candidate_area"] > missing_tolerance:
            continue
        area_growth_budget = max(
            0.5 * bevel_length * target_scale,
            16.0 * grid * grid,
            1.0e-9,
        )
        if difference_metrics["candidate_minus_reference_area"] > area_growth_budget:
            continue

        edit_zone = candidate_polygon.symmetric_difference(polygon)
        score = (
            difference_metrics["symmetric_difference_area"],
            difference_metrics["candidate_minus_reference_area"],
            candidate_signature.vertex_count,
        )
        if best_candidate is None or score < best_candidate[0]:
            best_candidate = (score, orient(candidate_polygon, sign=1.0), edit_zone)

    if best_candidate is None:
        return None

    return _RepairCandidate(
        polygon=best_candidate[1],
        edit_zone=best_candidate[2],
        operator="bevel_corner_collapse",
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
    if candidate.area_balance_budget_override is not None:
        area_balance_budget = max(
            area_balance_budget,
            float(candidate.area_balance_budget_override),
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
    recovery_signature_cache = _CoverageEvalCache()
    before_signature = _cached_coverage_defect_signature(
        recovery_signature_cache,
        polygons,
        target_scale=min_segment_length,
    )
    diagnostics["source_coordinate_recovery_short_edge_count_before"] = before_stats[
        "short_edge_count"
    ]
    diagnostics["source_coordinate_recovery_pair_issue_count_before"] = (
        before_signature.pair_issue_count
    )
    diagnostics["source_coordinate_recovery_ring_contact_count_before"] = (
        before_signature.ring_contact_count
    )
    diagnostics["source_coordinate_recovery_min_clearance_before"] = (
        before_signature.min_clearance
    )
    if min_segment_length <= 0 or not polygons:
        diagnostics["source_coordinate_recovery_applied"] = False
        diagnostics["source_coordinate_recovery_short_edge_count_after"] = before_stats[
            "short_edge_count"
        ]
        diagnostics["source_coordinate_recovery_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["source_coordinate_recovery_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        diagnostics["source_coordinate_recovery_min_clearance_after"] = (
            before_signature.min_clearance
        )
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
        diagnostics["source_coordinate_recovery_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["source_coordinate_recovery_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        diagnostics["source_coordinate_recovery_min_clearance_after"] = (
            before_signature.min_clearance
        )
        diagnostics["source_coordinate_recovery_applied"] = False
        diagnostics["source_coordinate_recovery_reference_minus_candidate_area"] = 0.0
        diagnostics["source_coordinate_recovery_candidate_minus_reference_area"] = 0.0
        diagnostics["source_coordinate_recovery_signed_area_delta"] = 0.0
        return polygons, source_map

    current_polygons, current_sources = _stable_sort(current_polygons, current_sources)
    after_signature = _cached_coverage_defect_signature(
        recovery_signature_cache,
        current_polygons,
        target_scale=min_segment_length,
    )
    if not _coverage_signature_not_worse(
        before_signature,
        after_signature,
        grid=grid,
        target_scale=min_segment_length,
    ):
        diagnostics["source_coordinate_recovery_contract_reject_count"] += 1
        diagnostics["source_coordinate_recovery_short_edge_count_after"] = before_stats[
            "short_edge_count"
        ]
        diagnostics["source_coordinate_recovery_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["source_coordinate_recovery_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        diagnostics["source_coordinate_recovery_min_clearance_after"] = (
            before_signature.min_clearance
        )
        diagnostics["source_coordinate_recovery_applied"] = False
        diagnostics["source_coordinate_recovery_applied_count"] = 0
        diagnostics["source_coordinate_recovery_exact_count"] = 0
        diagnostics["source_coordinate_recovery_vertex_count"] = 0
        diagnostics["source_coordinate_recovery_operator_applied"] = {}
        diagnostics["source_coordinate_recovery_reference_minus_candidate_area"] = 0.0
        diagnostics["source_coordinate_recovery_candidate_minus_reference_area"] = 0.0
        diagnostics["source_coordinate_recovery_signed_area_delta"] = 0.0
        diagnostics["source_coordinate_recovery_support_reference_minus_candidate_area"] = (
            0.0
        )
        diagnostics["source_coordinate_recovery_support_candidate_minus_reference_area"] = (
            0.0
        )
        diagnostics["source_coordinate_recovery_support_signed_area_delta"] = 0.0
        return polygons, source_map

    after_stats = _segment_length_stats(
        current_polygons,
        short_edge_threshold=min_segment_length,
    )
    diagnostics["source_coordinate_recovery_short_edge_count_after"] = after_stats[
        "short_edge_count"
    ]
    diagnostics["source_coordinate_recovery_pair_issue_count_after"] = (
        after_signature.pair_issue_count
    )
    diagnostics["source_coordinate_recovery_ring_contact_count_after"] = (
        after_signature.ring_contact_count
    )
    diagnostics["source_coordinate_recovery_min_clearance_after"] = (
        after_signature.min_clearance
    )
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
    def _meshing_regularization_score(
        signature: _CoverageDefectSignature,
        difference_metrics: dict[str, float],
    ) -> tuple[float, ...]:
        return (
            float(
                _coverage_signature_contract_priority(
                    signature,
                    target_scale=min_segment_length,
                    grid=grid,
                )
            ),
            *_coverage_signature_score(
                signature,
                target_scale=min_segment_length,
            ),
            difference_metrics["reference_minus_candidate_area"],
            difference_metrics["candidate_minus_reference_area"],
            difference_metrics["symmetric_difference_area"],
            abs(difference_metrics["union_area_delta"]),
        )

    def _regularize_ring_contacts(
        input_polygons: list[Polygon],
        input_sources: list[list[int]],
    ) -> tuple[list[Polygon], list[list[int]]] | None:
        candidate = _regularize_coverage_ring_contacts(
            input_polygons,
            input_sources,
            min_segment_length=min_segment_length,
            grid=grid,
        )
        if candidate is None:
            return None
        (candidate_polygons, candidate_sources), stats = candidate
        diagnostics["coverage_meshing_regularization_ring_contact_polygon_count"] = (
            stats["changed_count"]
        )
        diagnostics["coverage_meshing_regularization_ring_contact_component_count"] = (
            stats["component_count"]
        )
        diagnostics["coverage_meshing_regularization_ring_contact_failed_count"] = (
            stats["failed_count"]
        )
        diagnostics["coverage_meshing_regularization_ring_contact_area_delta"] = (
            stats["area_delta"]
        )
        return candidate_polygons, candidate_sources

    def _repair_residual_self_clearance_for_meshing(
        input_polygons: list[Polygon],
        input_sources: list[list[int]],
    ) -> tuple[
        list[Polygon],
        list[list[int]],
        _CoverageDefectSignature,
        dict[str, float],
        dict[str, int],
        dict[str, int],
    ] | None:
        tolerance = max(grid, 1e-9)

        def _residual_clearance_repair_score(
            signature: _CoverageDefectSignature,
            difference_metrics: dict[str, float],
        ) -> tuple[float, ...]:
            # Favor the smallest overall support drift once the contract metrics
            # tie. For residual self-clearance hotspots like case 54, this keeps
            # us from trading one local slit for large unsupported wedges.
            return (
                float(
                    _coverage_signature_contract_priority(
                        signature,
                        target_scale=min_segment_length,
                        grid=grid,
                    )
                ),
                *_coverage_signature_score(
                    signature,
                    target_scale=min_segment_length,
                ),
                difference_metrics["symmetric_difference_area"],
                max(
                    difference_metrics["reference_minus_candidate_area"],
                    difference_metrics["candidate_minus_reference_area"],
                ),
                difference_metrics["candidate_minus_reference_area"],
                difference_metrics["reference_minus_candidate_area"],
                abs(difference_metrics["union_area_delta"]),
            )

        reference_score = _residual_clearance_repair_score(
            best_signature,
            best_difference_metrics,
        )

        best_variant: tuple[
            list[Polygon],
            list[list[int]],
            _CoverageDefectSignature,
            dict[str, float],
            tuple[float, ...],
            dict[str, int],
            dict[str, int],
        ] | None = None

        def consider_variant(
            candidate_polygons: list[Polygon],
            candidate_sources: list[list[int]],
            *,
            operator_attempts: dict[str, int],
            operator_applied: dict[str, int],
        ) -> None:
            nonlocal best_variant
            candidate_signature = _cached_coverage_defect_signature(
                cache,
                candidate_polygons,
                target_scale=min_segment_length,
            )
            candidate_difference_metrics = _cached_difference_area_metrics(
                cache,
                polygons,
                candidate_polygons,
            )
            reference_clearance_deficit = max(
                min_segment_length - (best_signature.min_clearance or 0.0),
                0.0,
            )
            candidate_clearance_deficit = max(
                min_segment_length - (candidate_signature.min_clearance or 0.0),
                0.0,
            )
            candidate_near_threshold_short_edges_ok = (
                candidate_signature.short_edge_count > best_signature.short_edge_count
                and (candidate_signature.min_edge_length or 0.0) + tolerance
                >= min_segment_length
            )
            if candidate_signature.ring_contact_count > best_signature.ring_contact_count:
                return
            if candidate_signature.pair_issue_count > best_signature.pair_issue_count:
                return
            if (
                candidate_signature.pair_issue_count == best_signature.pair_issue_count
                and candidate_signature.close_pair_count > best_signature.close_pair_count
            ):
                return
            if candidate_clearance_deficit > reference_clearance_deficit + tolerance:
                return
            if (
                best_signature.min_clearance is not None
                and candidate_signature.min_clearance is not None
                and candidate_signature.min_clearance + tolerance
                < best_signature.min_clearance
            ):
                return
            if (
                candidate_signature.short_edge_count > best_signature.short_edge_count
                and not candidate_near_threshold_short_edges_ok
            ):
                return

            candidate_score = _residual_clearance_repair_score(
                candidate_signature,
                candidate_difference_metrics,
            )
            if candidate_score >= reference_score:
                return
            if best_variant is None or candidate_score < best_variant[4]:
                best_variant = (
                    candidate_polygons,
                    candidate_sources,
                    candidate_signature,
                    candidate_difference_metrics,
                    candidate_score,
                    operator_attempts,
                    operator_applied,
                )

        direct_diagnostics = _empty_diagnostics(len(input_polygons))
        direct_diagnostics["collect_stage_metrics"] = False
        direct_diagnostics["enable_logging"] = False
        direct_polygons, direct_sources = _regularize_low_clearance_polygons(
            input_polygons,
            input_sources,
            min_clearance=min_segment_length,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=direct_diagnostics,
        )
        if _polygon_sequence_key(direct_polygons) != _polygon_sequence_key(input_polygons):
            consider_variant(
                direct_polygons,
                direct_sources,
                operator_attempts={"meshing_contract_direct_clearance_regularization": 1},
                operator_applied={"meshing_contract_direct_clearance_regularization": 1},
            )

        local_diagnostics = _empty_diagnostics(len(input_polygons))
        local_diagnostics["collect_stage_metrics"] = False
        local_diagnostics["enable_logging"] = False
        local_polygons, local_sources = _apply_local_polygon_repairs(
            input_polygons,
            input_sources,
            min_segment_length=min_segment_length,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=local_diagnostics,
            stage_prefix="meshing_contract_repair",
            enable_defect_operators=True,
            enable_simplify_operators=True,
        )
        local_operator_attempts = dict(
            local_diagnostics.get("meshing_contract_repair_operator_attempts", {})
        )
        local_operator_applied = dict(
            local_diagnostics.get("meshing_contract_repair_operator_applied", {})
        )
        if _polygon_sequence_key(local_polygons) != _polygon_sequence_key(input_polygons):
            consider_variant(
                local_polygons,
                local_sources,
                operator_attempts=local_operator_attempts,
                operator_applied=local_operator_applied,
            )

        local_regularized_diagnostics = _empty_diagnostics(len(local_polygons))
        local_regularized_diagnostics["collect_stage_metrics"] = False
        local_regularized_diagnostics["enable_logging"] = False
        local_regularized_polygons, local_regularized_sources = (
            _regularize_low_clearance_polygons(
                local_polygons,
                local_sources,
                min_clearance=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=local_regularized_diagnostics,
            )
        )
        local_regularized_operator_attempts = dict(local_operator_attempts)
        local_regularized_operator_applied = dict(local_operator_applied)
        if _polygon_sequence_key(local_regularized_polygons) != _polygon_sequence_key(
            local_polygons
        ):
            local_regularized_operator_attempts[
                "meshing_contract_clearance_regularization"
            ] = (
                local_regularized_operator_attempts.get(
                    "meshing_contract_clearance_regularization",
                    0,
                )
                + 1
            )
            local_regularized_operator_applied[
                "meshing_contract_clearance_regularization"
            ] = (
                local_regularized_operator_applied.get(
                    "meshing_contract_clearance_regularization",
                    0,
                )
                + 1
            )
        if _polygon_sequence_key(local_regularized_polygons) != _polygon_sequence_key(
            input_polygons
        ):
            consider_variant(
                local_regularized_polygons,
                local_regularized_sources,
                operator_attempts=local_regularized_operator_attempts,
                operator_applied=local_regularized_operator_applied,
            )

        if best_variant is None:
            return None

        return (
            best_variant[0],
            best_variant[1],
            best_variant[2],
            best_variant[3],
            best_variant[5],
            best_variant[6],
        )

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
        *_meshing_regularization_score(
            best_signature,
            best_difference_metrics,
        ),
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
        score = _meshing_regularization_score(
            candidate_signature,
            difference_metrics,
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

    def _repair_residual_pair_issues_for_meshing(
        input_polygons: list[Polygon],
        input_sources: list[list[int]],
    ) -> tuple[list[Polygon], list[list[int]], dict[str, int]] | None:
        return _repair_residual_pair_issues(
            input_polygons,
            input_sources,
            target_scale=min_segment_length,
            grid=grid,
            diagnostics=diagnostics,
            cache=cache,
        )

    if best_signature.pair_issue_count > 0:
        rescue_candidate = _repair_residual_pair_issues_for_meshing(
            best_polygons,
            best_sources,
        )
        if rescue_candidate is not None:
            rescue_polygons, rescue_sources, rescue_operator_counts = rescue_candidate
            rescue_baseline_signature = _cached_coverage_defect_signature(
                cache,
                rescue_polygons,
                target_scale=min_segment_length,
            )
            rescue_baseline_difference_metrics = _cached_difference_area_metrics(
                cache,
                polygons,
                rescue_polygons,
            )
            rescue_variants: list[tuple[list[Polygon], list[list[int]]]] = [
                (rescue_polygons, rescue_sources)
            ]
            if _coverage_signature_requires_polygon_regularization(
                rescue_baseline_signature,
                target_scale=min_segment_length,
                grid=grid,
            ):
                rescue_diagnostics = _empty_diagnostics(len(rescue_polygons))
                rescue_diagnostics["collect_stage_metrics"] = False
                rescue_diagnostics["enable_logging"] = False
                postprocessed_rescue_polygons, postprocessed_rescue_sources = (
                    _simplify_polygons_for_meshing(
                        rescue_polygons,
                        rescue_sources,
                        min_segment_length=min_segment_length,
                        grid=grid,
                        min_area=min_area,
                        min_hole_area=min_hole_area,
                        diagnostics=rescue_diagnostics,
                    )
                )
                (
                    postprocessed_rescue_polygons,
                    postprocessed_rescue_sources,
                ) = _regularize_low_clearance_polygons(
                    postprocessed_rescue_polygons,
                    postprocessed_rescue_sources,
                    min_clearance=min_segment_length,
                    grid=grid,
                    min_area=min_area,
                    min_hole_area=min_hole_area,
                    diagnostics=rescue_diagnostics,
                )
                rescue_variants.append(
                    (
                        postprocessed_rescue_polygons,
                        postprocessed_rescue_sources,
                    )
                )

            best_rescue_variant: tuple[
                list[Polygon],
                list[list[int]],
                _CoverageDefectSignature,
                dict[str, float],
                tuple[float, int, int, int, float, int, float, float, float, float],
            ] | None = None
            for candidate_polygons, candidate_sources in rescue_variants:
                candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    candidate_polygons,
                    target_scale=min_segment_length,
                )
                candidate_difference_metrics = _cached_difference_area_metrics(
                    cache,
                    polygons,
                    candidate_polygons,
                )
                rescue_refinement_tolerance = max(
                    min_segment_length * min_segment_length,
                    16.0 * grid * grid,
                    1e-9,
                )
                if (
                    candidate_signature.pair_issue_count
                    == rescue_baseline_signature.pair_issue_count
                    and candidate_signature.short_edge_count
                    == rescue_baseline_signature.short_edge_count
                    and candidate_difference_metrics["reference_minus_candidate_area"]
                    > rescue_baseline_difference_metrics[
                        "reference_minus_candidate_area"
                    ]
                    + rescue_refinement_tolerance
                ):
                    continue
                if (
                    candidate_signature.pair_issue_count
                    == rescue_baseline_signature.pair_issue_count
                    and candidate_signature.short_edge_count
                    == rescue_baseline_signature.short_edge_count
                    and candidate_difference_metrics["candidate_minus_reference_area"]
                    > rescue_baseline_difference_metrics[
                        "candidate_minus_reference_area"
                    ]
                    + rescue_refinement_tolerance
                ):
                    continue
                if (
                    candidate_signature.pair_issue_count
                    == rescue_baseline_signature.pair_issue_count
                    and candidate_signature.short_edge_count
                    == rescue_baseline_signature.short_edge_count
                    and candidate_difference_metrics["symmetric_difference_area"]
                    > rescue_baseline_difference_metrics["symmetric_difference_area"]
                    + 2.0 * rescue_refinement_tolerance
                ):
                    continue
                if not _coverage_signature_improves(
                    best_signature,
                    candidate_signature,
                    grid=grid,
                    target_scale=min_segment_length,
                ):
                    continue
                candidate_score = (
                    *_contact_resolution_candidate_score(
                        candidate_signature,
                        candidate_difference_metrics,
                        target_scale=min_segment_length,
                    ),
                )
                if (
                    best_rescue_variant is None
                    or candidate_score < best_rescue_variant[4]
                ):
                    best_rescue_variant = (
                        candidate_polygons,
                        candidate_sources,
                        candidate_signature,
                        candidate_difference_metrics,
                        candidate_score,
                    )

            if best_rescue_variant is not None:
                (
                    rescue_polygons,
                    rescue_sources,
                    rescue_signature,
                    rescue_difference_metrics,
                    _,
                ) = best_rescue_variant
                best_polygons = rescue_polygons
                best_sources = rescue_sources
                best_signature = rescue_signature
                best_difference_metrics = rescue_difference_metrics
                best_label = "residual_pair_issue_rescue"
                best_score = (
                    *_meshing_regularization_score(
                        best_signature,
                        best_difference_metrics,
                    ),
                )
                best_operator_attempts = dict(best_operator_attempts)
                for operator_name, count in rescue_operator_counts.items():
                    best_operator_attempts[operator_name] = (
                        best_operator_attempts.get(operator_name, 0) + count
                    )
                best_operator_applied = dict(best_operator_applied)
                for operator_name, count in rescue_operator_counts.items():
                    best_operator_applied[operator_name] = (
                        best_operator_applied.get(operator_name, 0) + count
                    )

    if best_signature.pair_issue_count > 0 or before_signature.pair_issue_count > 0:
        direct_rescue_diagnostics = _empty_diagnostics(len(best_polygons))
        direct_rescue_diagnostics["collect_stage_metrics"] = False
        direct_rescue_diagnostics["enable_logging"] = False
        direct_rescue_candidate: tuple[
            list[Polygon],
            list[list[int]],
            _CoverageDefectSignature,
            dict[str, float],
            tuple[float, int, int, int, float, int, float, float, float, float],
            str,
        ] | None = None
        candidate_streams: list[tuple[list[Polygon], list[list[int]], bool]] = [
            (best_polygons, best_sources, False)
        ]
        if before_signature.pair_issue_count > 0:
            candidate_streams.append((polygons, source_map, True))

        for (
            candidate_input_polygons,
            candidate_input_sources,
            from_original_coverage,
        ) in candidate_streams:
            for (
                _affected_indices,
                operator_name,
                candidate_polygons,
                candidate_sources,
            ) in _direct_pair_issue_cluster_candidates(
                candidate_input_polygons,
                candidate_input_sources,
                tolerance=min_segment_length,
                grid=grid,
                diagnostics=direct_rescue_diagnostics,
                cache=cache,
            ):
                candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    candidate_polygons,
                    target_scale=min_segment_length,
                )
                if from_original_coverage and candidate_signature.pair_issue_count > 0:
                    continue
                if not _coverage_signature_improves(
                    best_signature,
                    candidate_signature,
                    grid=grid,
                    target_scale=min_segment_length,
                ):
                    continue
                candidate_difference_metrics = _cached_difference_area_metrics(
                    cache,
                    polygons,
                    candidate_polygons,
                )
                candidate_score = (
                    *_meshing_regularization_score(
                        candidate_signature,
                        candidate_difference_metrics,
                    ),
                )
                if (
                    direct_rescue_candidate is None
                    or candidate_score < direct_rescue_candidate[4]
                ):
                    direct_rescue_candidate = (
                        candidate_polygons,
                        candidate_sources,
                        candidate_signature,
                        candidate_difference_metrics,
                        candidate_score,
                        operator_name,
                    )

        if direct_rescue_candidate is not None:
            (
                best_polygons,
                best_sources,
                best_signature,
                best_difference_metrics,
                _,
                applied_operator_name,
            ) = direct_rescue_candidate
            best_label = "residual_pair_issue_rescue"
            best_score = (
                *_meshing_regularization_score(
                    best_signature,
                    best_difference_metrics,
                ),
            )
            best_operator_attempts = dict(best_operator_attempts)
            best_operator_attempts[applied_operator_name] = (
                best_operator_attempts.get(applied_operator_name, 0) + 1
            )
            best_operator_applied = dict(best_operator_applied)
            best_operator_applied[applied_operator_name] = (
                best_operator_applied.get(applied_operator_name, 0) + 1
            )

    if best_signature.ring_contact_count > 0:
        final_ring_contact_candidate = _regularize_ring_contacts(
            best_polygons,
            best_sources,
        )
        if final_ring_contact_candidate is not None:
            candidate_polygons, candidate_sources = final_ring_contact_candidate
            candidate_signature = _cached_coverage_defect_signature(
                cache,
                candidate_polygons,
                target_scale=min_segment_length,
            )
            if _coverage_signature_improves(
                best_signature,
                candidate_signature,
                grid=grid,
                target_scale=min_segment_length,
            ):
                best_polygons = candidate_polygons
                best_sources = candidate_sources
                best_signature = candidate_signature
                best_difference_metrics = _cached_difference_area_metrics(
                    cache,
                    polygons,
                    candidate_polygons,
                )
                best_score = (
                    *_meshing_regularization_score(
                        best_signature,
                        best_difference_metrics,
                    ),
                )
                ring_contact_applied = True

    residual_clearance_deficit = max(
        min_segment_length - (best_signature.min_clearance or 0.0),
        0.0,
    )
    clearance_repair_iterations = 0
    while (
        best_signature.pair_issue_count == 0
        and best_signature.ring_contact_count == 0
        and residual_clearance_deficit > max(grid, 1e-9)
        and clearance_repair_iterations < 3
    ):
        clearance_repair_candidate = _repair_residual_self_clearance_for_meshing(
            best_polygons,
            best_sources,
        )
        if clearance_repair_candidate is None:
            break

        (
            candidate_polygons,
            candidate_sources,
            candidate_signature,
            candidate_difference_metrics,
            repair_operator_attempts,
            repair_operator_applied,
        ) = clearance_repair_candidate
        if _polygon_sequence_key(candidate_polygons) == _polygon_sequence_key(best_polygons):
            break

        best_polygons = candidate_polygons
        best_sources = candidate_sources
        best_signature = candidate_signature
        best_difference_metrics = candidate_difference_metrics
        best_label = "residual_self_clearance_repair"
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
        best_operator_attempts = dict(best_operator_attempts)
        for operator_name, count in repair_operator_attempts.items():
            best_operator_attempts[operator_name] = (
                best_operator_attempts.get(operator_name, 0) + count
            )
        best_operator_applied = dict(best_operator_applied)
        for operator_name, count in repair_operator_applied.items():
            best_operator_applied[operator_name] = (
                best_operator_applied.get(operator_name, 0) + count
            )
        clearance_repair_iterations += 1
        residual_clearance_deficit = max(
            min_segment_length - (best_signature.min_clearance or 0.0),
            0.0,
        )

    diagnostics["coverage_meshing_regularization_clearance_repair_iterations"] = (
        clearance_repair_iterations
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


def _regularize_coverage_contacts(
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
    if cache is None:
        cache = _CoverageEvalCache()

    before_signature = _cached_coverage_defect_signature(
        cache,
        polygons,
        target_scale=min_segment_length,
    )
    diagnostics["coverage_contact_regularization_short_edge_count_before"] = (
        before_signature.short_edge_count
    )
    diagnostics["coverage_contact_regularization_pair_issue_count_before"] = (
        before_signature.pair_issue_count
    )
    diagnostics["coverage_contact_regularization_ring_contact_count_before"] = (
        before_signature.ring_contact_count
    )

    if min_segment_length <= 0 or not polygons or before_signature.pair_issue_count == 0:
        diagnostics["coverage_contact_regularization_applied"] = False
        diagnostics["coverage_contact_regularization_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        diagnostics["coverage_contact_regularization_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["coverage_contact_regularization_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        return polygons, source_map

    current_polygons = list(polygons)
    current_sources = [list(indices) for indices in source_map]
    operator_attempts: dict[str, int] = {}
    operator_applied: dict[str, int] = {}
    progress = False
    max_contact_shrink_area_loss = max(
        32.0 * min_segment_length * min_segment_length,
        128.0 * grid * grid,
        1e-9,
    )

    while True:
        current_signature = _cached_coverage_defect_signature(
            cache,
            current_polygons,
            target_scale=min_segment_length,
        )
        if current_signature.pair_issue_count == 0:
            break

        current_union = _cached_union(cache, current_polygons)
        local_radius = _derive_local_edit_radius(min_segment_length)

        best_step: tuple[
            list[Polygon],
            list[list[int]],
            _CoverageDefectSignature,
            dict[str, float],
            tuple[float, int, int, int, float, int, float, float, float, float],
            str,
        ] | None = None

        for (
            affected_indices,
            operator_name,
            candidate_polygons,
            candidate_sources,
        ) in _direct_pair_issue_cluster_candidates(
            current_polygons,
            current_sources,
            tolerance=min_segment_length,
            grid=grid,
            diagnostics=diagnostics,
            cache=cache,
        ):
            operator_attempts[operator_name] = operator_attempts.get(operator_name, 0) + 1

            candidate_signature = _cached_coverage_defect_signature(
                cache,
                candidate_polygons,
                target_scale=min_segment_length,
            )
            if candidate_signature.pair_issue_count >= current_signature.pair_issue_count:
                continue
            if candidate_signature.ring_contact_count > current_signature.ring_contact_count:
                continue

            candidate_difference_metrics = _cached_difference_area_metrics(
                cache,
                current_polygons,
                candidate_polygons,
            )
            candidate_edit_zone = _cached_coverage_edit_zone(
                cache,
                [current_polygons[index] for index in affected_indices],
                target_scale=min_segment_length,
                radius=local_radius,
            )
            if (
                "shrink" in operator_name
                and candidate_difference_metrics["reference_minus_candidate_area"]
                > max_contact_shrink_area_loss
            ):
                continue
            if "shrink" in operator_name and len(affected_indices) == 2:
                shrink_budget_context = _pair_issue_shrink_budget_context(
                    [current_polygons[index] for index in affected_indices],
                    tolerance=min_segment_length,
                    grid=grid,
                )
                if shrink_budget_context is not None:
                    candidate_edit_zone = shrink_budget_context[0]
            candidate_union = _cached_union(cache, candidate_polygons)
            change_outside_edit_zone = _change_outside_edit_zone(
                current_union,
                candidate_union,
                edit_zone=candidate_edit_zone,
            )
            if change_outside_edit_zone > max(
                float(candidate_edit_zone.area),
                grid * grid,
                1e-9,
            ):
                continue

            candidate_score = (
                _prefer_fill_only_point_contact_candidate(
                    operator_name,
                    candidate_difference_metrics,
                    target_scale=min_segment_length,
                    grid=grid,
                ),
                *_contact_resolution_candidate_score(
                    candidate_signature,
                    candidate_difference_metrics,
                    target_scale=min_segment_length,
                ),
            )
            if best_step is None or candidate_score < best_step[4]:
                best_step = (
                    candidate_polygons,
                    candidate_sources,
                    candidate_signature,
                    candidate_difference_metrics,
                    candidate_score,
                    operator_name,
                )

        if best_step is None:
            break

        (
            current_polygons,
            current_sources,
            _,
            _,
            _,
            applied_operator_name,
        ) = best_step
        operator_applied[applied_operator_name] = (
            operator_applied.get(applied_operator_name, 0) + 1
        )
        progress = True

    if not progress:
        diagnostics["coverage_contact_regularization_applied"] = False
        diagnostics["coverage_contact_regularization_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        diagnostics["coverage_contact_regularization_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["coverage_contact_regularization_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        return polygons, source_map

    current_signature = _cached_coverage_defect_signature(
        cache,
        current_polygons,
        target_scale=min_segment_length,
    )
    candidate_variants: list[tuple[list[Polygon], list[list[int]]]] = [
        (current_polygons, current_sources)
    ]
    if _coverage_signature_requires_polygon_regularization(
        current_signature,
        target_scale=min_segment_length,
        grid=grid,
    ):
        local_diagnostics = _empty_diagnostics(len(current_polygons))
        local_diagnostics["collect_stage_metrics"] = False
        local_diagnostics["enable_logging"] = False
        postprocessed_polygons, postprocessed_sources = _simplify_polygons_for_meshing(
            current_polygons,
            current_sources,
            min_segment_length=min_segment_length,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=local_diagnostics,
        )
        postprocessed_polygons, postprocessed_sources = _regularize_low_clearance_polygons(
            postprocessed_polygons,
            postprocessed_sources,
            min_clearance=min_segment_length,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=local_diagnostics,
        )
        candidate_variants.append((postprocessed_polygons, postprocessed_sources))

    best_variant: tuple[
        list[Polygon],
        list[list[int]],
        _CoverageDefectSignature,
        dict[str, float],
        tuple[float, int, int, int, float, int, float, float, float, float],
    ] | None = None
    for candidate_polygons, candidate_sources in candidate_variants:
        candidate_signature = _cached_coverage_defect_signature(
            cache,
            candidate_polygons,
            target_scale=min_segment_length,
        )
        if not _coverage_signature_improves(
            before_signature,
            candidate_signature,
            grid=grid,
            target_scale=min_segment_length,
        ):
            continue
        candidate_difference_metrics = _cached_difference_area_metrics(
            cache,
            polygons,
            candidate_polygons,
        )
        candidate_score = _contact_resolution_candidate_score(
            candidate_signature,
            candidate_difference_metrics,
            target_scale=min_segment_length,
        )
        if best_variant is None or candidate_score < best_variant[4]:
            best_variant = (
                candidate_polygons,
                candidate_sources,
                candidate_signature,
                candidate_difference_metrics,
                candidate_score,
            )

    if best_variant is None:
        diagnostics["coverage_contact_regularization_applied"] = False
        diagnostics["coverage_contact_regularization_short_edge_count_after"] = (
            before_signature.short_edge_count
        )
        diagnostics["coverage_contact_regularization_pair_issue_count_after"] = (
            before_signature.pair_issue_count
        )
        diagnostics["coverage_contact_regularization_ring_contact_count_after"] = (
            before_signature.ring_contact_count
        )
        diagnostics["coverage_contact_regularization_operator_attempts"] = operator_attempts
        diagnostics["coverage_contact_regularization_operator_applied"] = {}
        return polygons, source_map

    (
        best_polygons,
        best_sources,
        best_signature,
        best_difference_metrics,
        _,
    ) = best_variant
    if best_signature.pair_issue_count > 0:
        rescue_candidate = _repair_residual_pair_issues(
            best_polygons,
            best_sources,
            target_scale=min_segment_length,
            grid=grid,
            diagnostics=diagnostics,
            cache=cache,
        )
        if rescue_candidate is not None:
            rescue_polygons, rescue_sources, rescue_operator_counts = rescue_candidate
            rescue_baseline_signature = _cached_coverage_defect_signature(
                cache,
                rescue_polygons,
                target_scale=min_segment_length,
            )
            rescue_baseline_difference_metrics = _cached_difference_area_metrics(
                cache,
                polygons,
                rescue_polygons,
            )
            rescue_variants: list[tuple[list[Polygon], list[list[int]]]] = [
                (rescue_polygons, rescue_sources)
            ]
            rescue_diagnostics = _empty_diagnostics(len(rescue_polygons))
            rescue_diagnostics["collect_stage_metrics"] = False
            rescue_diagnostics["enable_logging"] = False
            postprocessed_rescue_polygons, postprocessed_rescue_sources = (
                _simplify_polygons_for_meshing(
                    rescue_polygons,
                    rescue_sources,
                    min_segment_length=min_segment_length,
                    grid=grid,
                    min_area=min_area,
                    min_hole_area=min_hole_area,
                    diagnostics=rescue_diagnostics,
                )
            )
            (
                postprocessed_rescue_polygons,
                postprocessed_rescue_sources,
            ) = _regularize_low_clearance_polygons(
                postprocessed_rescue_polygons,
                postprocessed_rescue_sources,
                min_clearance=min_segment_length,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=rescue_diagnostics,
            )
            rescue_variants.append(
                (
                    postprocessed_rescue_polygons,
                    postprocessed_rescue_sources,
                )
            )

            best_rescue_variant: tuple[
                list[Polygon],
                list[list[int]],
                _CoverageDefectSignature,
                dict[str, float],
                tuple[float, int, int, int, float, int, float, float, float, float],
            ] | None = None
            for candidate_polygons, candidate_sources in rescue_variants:
                candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    candidate_polygons,
                    target_scale=min_segment_length,
                )
                candidate_difference_metrics = _cached_difference_area_metrics(
                    cache,
                    polygons,
                    candidate_polygons,
                )
                rescue_refinement_tolerance = max(
                    min_segment_length * min_segment_length,
                    16.0 * grid * grid,
                    1e-9,
                )
                if (
                    candidate_signature.pair_issue_count
                    == rescue_baseline_signature.pair_issue_count
                    and candidate_signature.short_edge_count
                    == rescue_baseline_signature.short_edge_count
                    and candidate_difference_metrics["reference_minus_candidate_area"]
                    > rescue_baseline_difference_metrics[
                        "reference_minus_candidate_area"
                    ]
                    + rescue_refinement_tolerance
                ):
                    continue
                if (
                    candidate_signature.pair_issue_count
                    == rescue_baseline_signature.pair_issue_count
                    and candidate_signature.short_edge_count
                    == rescue_baseline_signature.short_edge_count
                    and candidate_difference_metrics["candidate_minus_reference_area"]
                    > rescue_baseline_difference_metrics[
                        "candidate_minus_reference_area"
                    ]
                    + rescue_refinement_tolerance
                ):
                    continue
                if (
                    candidate_signature.pair_issue_count
                    == rescue_baseline_signature.pair_issue_count
                    and candidate_signature.short_edge_count
                    == rescue_baseline_signature.short_edge_count
                    and candidate_difference_metrics["symmetric_difference_area"]
                    > rescue_baseline_difference_metrics["symmetric_difference_area"]
                    + 2.0 * rescue_refinement_tolerance
                ):
                    continue
                if not _coverage_signature_improves(
                    best_signature,
                    candidate_signature,
                    grid=grid,
                    target_scale=min_segment_length,
                ):
                    continue
                candidate_score = (
                    *_contact_resolution_candidate_score(
                        candidate_signature,
                        candidate_difference_metrics,
                        target_scale=min_segment_length,
                    ),
                )
                if (
                    best_rescue_variant is None
                    or candidate_score < best_rescue_variant[4]
                ):
                    best_rescue_variant = (
                        candidate_polygons,
                        candidate_sources,
                        candidate_signature,
                        candidate_difference_metrics,
                        candidate_score,
                    )

            if best_rescue_variant is not None:
                (
                    best_polygons,
                    best_sources,
                    best_signature,
                    best_difference_metrics,
                    _,
                ) = best_rescue_variant
                for operator_name, count in rescue_operator_counts.items():
                    operator_attempts[operator_name] = (
                        operator_attempts.get(operator_name, 0) + count
                    )
                    operator_applied[operator_name] = (
                        operator_applied.get(operator_name, 0) + count
                    )
    if (
        best_signature.pair_issue_count > 0
        or _coverage_signature_requires_polygon_regularization(
            best_signature,
            target_scale=min_segment_length,
            grid=grid,
        )
    ):
        post_contact_polygons, post_contact_sources = _apply_local_polygon_repairs(
            best_polygons,
            best_sources,
            min_segment_length=min_segment_length,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=diagnostics,
            stage_prefix="post_contact_local_defect_repair",
            enable_defect_operators=True,
            enable_simplify_operators=True,
        )
        post_contact_signature = _cached_coverage_defect_signature(
            cache,
            post_contact_polygons,
            target_scale=min_segment_length,
        )
        if _should_accept_post_contact_local_repair(
            best_signature,
            best_difference_metrics,
            post_contact_signature,
            post_contact_difference_metrics := _cached_difference_area_metrics(
                cache,
                polygons,
                post_contact_polygons,
            ),
            target_scale=min_segment_length,
            grid=grid,
        ):
            best_polygons = post_contact_polygons
            best_sources = post_contact_sources
            best_signature = post_contact_signature
            best_difference_metrics = post_contact_difference_metrics
        elif _coverage_signature_improves(
            best_signature,
            post_contact_signature,
            grid=grid,
            target_scale=min_segment_length,
        ):
            post_contact_difference_metrics = _cached_difference_area_metrics(
                cache,
                polygons,
                post_contact_polygons,
            )
            post_contact_score = (
                *_contact_resolution_candidate_score(
                    post_contact_signature,
                    post_contact_difference_metrics,
                    target_scale=min_segment_length,
                ),
            )
            current_best_score = (
                *_contact_resolution_candidate_score(
                    best_signature,
                    best_difference_metrics,
                    target_scale=min_segment_length,
                ),
            )
            if post_contact_score < current_best_score:
                best_polygons = post_contact_polygons
                best_sources = post_contact_sources
                best_signature = post_contact_signature
                best_difference_metrics = post_contact_difference_metrics
    else:
        _record_polygon_repair_noop(
            diagnostics,
            stage_prefix="post_contact_local_defect_repair",
            tolerance=min_segment_length,
            before_stats=_segment_length_stats(
                best_polygons,
                short_edge_threshold=min_segment_length,
            ),
        )

    diagnostics["coverage_contact_regularization_applied"] = True
    diagnostics["coverage_contact_regularization_selected_branch"] = "pair_contacts"
    diagnostics["coverage_contact_regularization_short_edge_count_after"] = (
        best_signature.short_edge_count
    )
    diagnostics["coverage_contact_regularization_pair_issue_count_after"] = (
        best_signature.pair_issue_count
    )
    diagnostics["coverage_contact_regularization_ring_contact_count_after"] = (
        best_signature.ring_contact_count
    )
    diagnostics["coverage_contact_regularization_reference_minus_candidate_area"] = (
        best_difference_metrics["reference_minus_candidate_area"]
    )
    diagnostics["coverage_contact_regularization_candidate_minus_reference_area"] = (
        best_difference_metrics["candidate_minus_reference_area"]
    )
    diagnostics["coverage_contact_regularization_signed_area_delta"] = (
        best_difference_metrics["union_area_delta"]
    )
    diagnostics["coverage_contact_regularization_operator_attempts"] = operator_attempts
    diagnostics["coverage_contact_regularization_operator_applied"] = operator_applied
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


_STAGE_LABEL_TITLES = {
    "atomic_input": "input coverage",
    "opened": "scale opening",
    "regularized_groups": "group regularization",
    "reconstructed": "coverage reconstruction",
    "presimplify": "pre-simplification",
    "local_defect_repaired": "local defect repair",
    "coverage_simplified": "coverage simplification",
    "source_reclaimed": "source reclaim",
    "small_component_absorbed": "small-component absorption",
    "boundary_regularized": "boundary regularization",
    "clearance_regularized": "final clearance repair",
    "source_coordinate_recovered": "source coordinate recovery",
    "post_recovery_regularized": "post-recovery regularization",
    "coverage_contact_regularized": "contact regularization",
    "coverage_meshing_regularized": "mesher-ready regularization",
    "final_output": "output coverage",
}


def _format_stage_label(stage: str) -> str:
    return _STAGE_LABEL_TITLES.get(
        stage,
        stage.replace("_", " ").strip().lower(),
    )


def _format_stage_metrics(
    stage: str,
    metrics: dict[str, Any],
    *,
    short_edge_threshold: float,
) -> str:
    line = (
        f"{_format_stage_label(stage)}: "
        f"polys={metrics['polygon_count']} verts={metrics['vertex_count']} "
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
        "Footprint cleaning | "
        + _format_stage_metrics(
            label,
            metrics,
            short_edge_threshold=short_edge_threshold,
        )
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
    summary_lines = ["Footprint cleaning complete"]
    if diagnostics.get("collect_stage_metrics", True):
        summary_lines.extend(
            [
                (
                    f"  Input  | {_format_stage_metrics('atomic_input', atomic, short_edge_threshold=short_edge_threshold)}"
                    if atomic
                    else "  Input  | n/a"
                ),
                (
                    f"  Output | {_format_stage_metrics('final_output', final, short_edge_threshold=short_edge_threshold)}"
                    if final
                    else "  Output | n/a"
                ),
            ]
        )
    else:
        summary_lines.append("  Stage metrics disabled")

    summary_lines.extend(
        [
        (
            f"  Coverage | overlap={_fmt_metric(diagnostics['overlap_area_before'])} -> "
            f"{_fmt_metric(diagnostics['overlap_area_after'])} m2 "
            f"clearance={_fmt_metric(diagnostics['min_clearance_before'])} -> "
            f"{_fmt_metric(diagnostics['min_clearance_after'])} m"
        ),
        (
            f"  Output   | polygons={diagnostics['output_count']} "
            f"repaired_invalid={diagnostics['repaired_invalid_count']} "
            f"collapsed={diagnostics['collapsed_count']} "
            f"dropped_small={diagnostics['dropped_small_count']} "
            f"geos_exceptions={diagnostics['geos_exception_count']}"
        ),
        (
            f"  Actions  | local_defect={diagnostics['local_defect_repair_applied_count']}/"
            f"{diagnostics['local_defect_repair_candidate_count']} "
            f"coverage_simplify={diagnostics['coverage_simplify_patch_applied_count']}/"
            f"{diagnostics['coverage_simplify_patch_count']} "
            f"({diagnostics['coverage_simplify_selected_branch']}) "
            f"source_recovery={diagnostics['source_coordinate_recovery_applied_count']}/"
            f"{diagnostics['source_coordinate_recovery_candidate_count']} "
            f"boundary_regularization={diagnostics['polygon_simplify_applied_count']}/"
            f"{diagnostics['polygon_simplify_candidate_count']}"
        ),
        ]
    )
    info("\n".join(summary_lines))

    detail_lines = [
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
            f"contract_rejected={diagnostics['source_coordinate_recovery_contract_reject_count']} "
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
            f"  coverage_contact_regularization branch="
            f"{diagnostics['coverage_contact_regularization_selected_branch']} "
            f"short_edges={diagnostics['coverage_contact_regularization_short_edge_count_before']} -> "
            f"{diagnostics['coverage_contact_regularization_short_edge_count_after']} "
            f"pair_issues={diagnostics['coverage_contact_regularization_pair_issue_count_before']} -> "
            f"{diagnostics['coverage_contact_regularization_pair_issue_count_after']} "
            f"missing={_fmt_metric(diagnostics['coverage_contact_regularization_reference_minus_candidate_area'])} "
            f"extra={_fmt_metric(diagnostics['coverage_contact_regularization_candidate_minus_reference_area'])} "
            f"signed_delta={_fmt_metric(diagnostics['coverage_contact_regularization_signed_area_delta'])}"
        ),
        (
            f"  coverage_contact_regularization operators "
            f"attempted={diagnostics.get('coverage_contact_regularization_operator_attempts', {})} "
            f"applied={diagnostics.get('coverage_contact_regularization_operator_applied', {})}"
        ),
        (
            f"  coverage_void_regularization applied="
            f"{diagnostics['coverage_void_regularization_applied']} "
            f"holes={diagnostics['coverage_void_regularization_hole_cleanup_count']} "
            f"gap_patches={diagnostics['coverage_void_regularization_gap_patch_applied_count']} / "
            f"{diagnostics['coverage_void_regularization_gap_patch_count']} "
            f"notch_simplify={diagnostics['coverage_void_regularization_notch_simplify_count']}"
        ),
        (
            f"  coverage_void_regularization operators "
            f"attempted={diagnostics.get('coverage_void_regularization_operator_attempts', {})} "
            f"applied={diagnostics.get('coverage_void_regularization_operator_applied', {})}"
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
        *(
            [
                (
                    f"  coverage_residual_polygon_simplify targets="
                    f"{diagnostics.get('coverage_residual_polygon_target_count', 0)} "
                    f"attempts={diagnostics.get('coverage_residual_polygon_attempt_count', 0)} "
                    f"viable={diagnostics.get('coverage_residual_polygon_viable_candidate_count', 0)} "
                    f"changed={diagnostics.get('coverage_residual_polygon_changed_count', 0)} "
                    f"no_candidate={diagnostics.get('coverage_residual_polygon_no_candidate_count', 0)} "
                    f"reject_contract={diagnostics.get('coverage_residual_polygon_contract_reject_count', 0)} "
                    f"reject_drift={diagnostics.get('coverage_residual_polygon_drift_reject_count', 0)} "
                    f"reject_normalize={diagnostics.get('coverage_residual_polygon_normalization_failed_count', 0)}"
                )
            ]
            if diagnostics.get("coverage_residual_polygon_attempt_count", 0) > 0
            else []
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
    ]
    if diagnostics["geos_exception_messages"]:
        detail_lines.append(
            "  messages=" + " | ".join(diagnostics["geos_exception_messages"][:3])
        )
    debug("Footprint Cleaning Details\n" + "\n".join(detail_lines))


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


def _courtyard_pair_indices(
    polygons: Sequence[Polygon],
) -> tuple[int, int] | None:
    if len(polygons) != 2:
        return None

    larger_index = max(
        range(len(polygons)),
        key=lambda index: float(polygons[index].area),
    )
    smaller_index = 1 - larger_index
    larger = polygons[larger_index]
    smaller = polygons[smaller_index]
    if not larger.interiors:
        return None

    representative = smaller.representative_point()
    if not any(Polygon(ring).covers(representative) for ring in larger.interiors):
        return None

    return larger_index, smaller_index


def _apply_local_courtyard_pair_merge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or len(subset_polygons) != 2:
        return None
    if _courtyard_pair_indices(subset_polygons) is None:
        return None

    try:
        shortest = shapely.shortest_line(subset_polygons[0], subset_polygons[1])
    except GEOSException:
        return None
    if shortest.is_empty or shortest.length <= 1e-9:
        return None

    try:
        passage = shortest.buffer(
            max(radius, grid),
            quad_segs=1,
            join_style=BufferJoinStyle.mitre,
            mitre_limit=1000.0,
        )
    except GEOSException:
        return None
    if passage.is_empty:
        return None

    patch_geometry = unary_union([*subset_polygons, passage])
    candidate = _apply_local_pair_merge_operator(
        subset_polygons,
        subset_sources,
        radius=radius,
        grid=grid,
        diagnostics=diagnostics,
        patch_geometry=patch_geometry,
    )
    if candidate is None:
        return None

    candidate_polygons, candidate_sources = candidate
    if len(candidate_polygons) != 1:
        return None
    return candidate_polygons, candidate_sources


def _local_offset_overlap_envelope(
    subset_polygons: Sequence[Polygon],
    *,
    radius: float,
    grid: float,
) -> BaseGeometry | None:
    if radius <= 0 or len(subset_polygons) < 2:
        return None

    buffered_polygons: list[BaseGeometry] = []
    for polygon in subset_polygons:
        try:
            buffered = polygon.buffer(
                max(radius, grid),
                quad_segs=1,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
        except GEOSException:
            return None
        if buffered.is_empty:
            return None
        buffered_polygons.append(buffered)

    overlap_parts: list[BaseGeometry] = []
    for index, buffered in enumerate(buffered_polygons):
        for other_index in range(index + 1, len(buffered_polygons)):
            overlap = buffered.intersection(buffered_polygons[other_index])
            if overlap.is_empty or overlap.area <= max(grid * grid, 1e-9):
                continue
            overlap_parts.append(overlap)

    if not overlap_parts:
        return None

    return unary_union(overlap_parts)


def _extract_point_contacts(geometry: BaseGeometry) -> list[Point]:
    if geometry.is_empty:
        return []
    if isinstance(geometry, Point):
        return [geometry]

    contacts: list[Point] = []
    for part in shapely.get_parts(geometry):
        if isinstance(part, Point):
            contacts.append(part)
        elif not part.is_empty:
            contacts.extend(_extract_point_contacts(part))
    return contacts


def _nearest_boundary_contact_point(
    polygon: Polygon,
    geometry: BaseGeometry,
    *,
    tolerance: float,
) -> Point | None:
    try:
        closest = shapely.shortest_line(geometry, polygon.boundary)
    except GEOSException:
        return None
    coords = list(closest.coords)
    if len(coords) < 2:
        return None

    first = Point(coords[0])
    last = Point(coords[-1])
    boundary = polygon.boundary
    if boundary.distance(last) <= tolerance:
        return last
    if boundary.distance(first) <= tolerance:
        return first
    return None


def _boundary_edge_inward_direction(
    coords: np.ndarray,
    point_xy: np.ndarray,
    *,
    tolerance: float,
) -> np.ndarray | None:
    if len(coords) < 2:
        return None

    for start, end in zip(coords, np.vstack([coords[1:], coords[:1]])):
        edge = end - start
        edge_length = float(np.hypot(edge[0], edge[1]))
        if edge_length <= 1e-9:
            continue
        projection = np.dot(point_xy - start, edge) / (edge_length * edge_length)
        if projection <= 1e-6 or projection >= 1.0 - 1e-6:
            continue
        nearest = start + projection * edge
        if float(np.hypot(*(point_xy - nearest))) > tolerance:
            continue

        inward = np.array([-edge[1], edge[0]], dtype=float) / edge_length
        inward_norm = float(np.hypot(inward[0], inward[1]))
        if inward_norm <= 1e-9:
            return None
        return inward / inward_norm

    return None


def _ring_vertex_bisector_direction(
    coords: np.ndarray,
    index: int,
) -> np.ndarray | None:
    if len(coords) < 3:
        return None

    current = coords[index]
    previous = coords[index - 1]
    next_coord = coords[(index + 1) % len(coords)]
    previous_vector = previous - current
    next_vector = next_coord - current

    previous_norm = float(np.hypot(previous_vector[0], previous_vector[1]))
    next_norm = float(np.hypot(next_vector[0], next_vector[1]))
    if previous_norm <= 1e-9 or next_norm <= 1e-9:
        return None

    direction = previous_vector / previous_norm + next_vector / next_norm
    direction_norm = float(np.hypot(direction[0], direction[1]))
    if direction_norm <= 1e-9:
        midpoint = 0.5 * (previous + next_coord)
        direction = midpoint - current
        direction_norm = float(np.hypot(direction[0], direction[1]))
        if direction_norm <= 1e-9:
            return None

    return direction / direction_norm


def _polygon_inward_direction_at_contact(
    polygon: Polygon,
    contact: Point,
    *,
    tolerance: float,
) -> np.ndarray | None:
    oriented = orient(polygon, sign=1.0)
    contact_xy = np.asarray(contact.coords[0], dtype=float)
    ring_coords = [np.asarray(oriented.exterior.coords[:-1], dtype=float)]
    ring_coords.extend(
        np.asarray(interior.coords[:-1], dtype=float) for interior in oriented.interiors
    )
    for coords in ring_coords:
        if len(coords) < 3:
            continue
        distances = np.hypot(
            coords[:, 0] - contact_xy[0],
            coords[:, 1] - contact_xy[1],
        )
        matching = np.flatnonzero(distances <= tolerance)
        for index in matching:
            direction = _ring_vertex_bisector_direction(coords, int(index))
            if direction is not None:
                return direction
        direction = _boundary_edge_inward_direction(
            coords,
            contact_xy,
            tolerance=tolerance,
        )
        if direction is not None:
            return direction

    centroid = np.asarray(oriented.centroid.coords[0], dtype=float)
    direction = centroid - contact_xy
    direction_norm = float(np.hypot(direction[0], direction[1]))
    if direction_norm <= 1e-9:
        return None
    return direction / direction_norm


def _connector_half_width(
    *,
    radius: float,
    grid: float,
    width_factor: float = 0.4,
) -> float:
    return max(float(width_factor) * float(radius), float(grid), 1e-9)


def _iter_point_touch_connector_variants(
    *,
    radius: float,
    grid: float,
) -> tuple[tuple[float, bool], ...]:
    if radius <= 0:
        return ()

    variants: list[tuple[float, bool]] = [
        (0.4, False),
        (0.4, True),
        (0.55, False),
    ]
    if radius >= max(4.0 * grid, 1e-9):
        variants.append((1.0, False))

    deduped: list[tuple[float, bool]] = []
    for variant in variants:
        if variant in deduped:
            continue
        deduped.append(variant)
    return tuple(deduped)


def _select_best_connector_bridge_candidate(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    anchor_pairs: Sequence[tuple[int, Point, int, Point]],
    radius: float,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
    variants: Sequence[tuple[float, bool]],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or target_scale <= 0 or not variants or not anchor_pairs:
        return None

    reference_signature = _coverage_defect_signature(
        list(subset_polygons),
        target_scale=target_scale,
    )
    reference_union = unary_union(subset_polygons)
    best_candidate: tuple[
        tuple[int, int, int, float, int, float, int, float, float, float, float],
        list[Polygon],
        list[list[int]],
    ] | None = None

    for width_factor, allow_patch_union_cleanup in variants:
        candidate = _apply_local_connector_bridge_operator(
            subset_polygons,
            subset_sources,
            anchor_pairs=anchor_pairs,
            radius=radius,
            target_scale=target_scale,
            grid=grid,
            diagnostics=diagnostics,
            corridor_width_factor=width_factor,
            allow_patch_union_cleanup=allow_patch_union_cleanup,
        )
        if candidate is None:
            continue

        candidate_polygons, candidate_sources = candidate
        candidate_signature = _coverage_defect_signature(
            candidate_polygons,
            target_scale=target_scale,
        )
        if not _coverage_signature_improves(
            reference_signature,
            candidate_signature,
            grid=grid,
            target_scale=target_scale,
        ):
            continue

        difference_metrics = _difference_area_metrics(
            reference_union,
            unary_union(candidate_polygons),
        )
        candidate_score = _contact_resolution_candidate_score(
            candidate_signature,
            difference_metrics,
            target_scale=target_scale,
        )
        if best_candidate is None or candidate_score < best_candidate[0]:
            best_candidate = (
                candidate_score,
                candidate_polygons,
                candidate_sources,
            )

    if best_candidate is None:
        return None

    return best_candidate[1], best_candidate[2]


def _iter_connector_cleanup_radii(
    *,
    radius: float,
    grid: float,
) -> tuple[float, ...]:
    if radius <= 0:
        return ()

    cleanup_radii: list[float] = []
    for factor in (0.25,):
        cleanup_radius = max(float(grid), float(radius) * factor, 1e-9)
        if any(abs(cleanup_radius - existing) <= 1e-12 for existing in cleanup_radii):
            continue
        cleanup_radii.append(cleanup_radius)
    return tuple(cleanup_radii)


def _inset_contact_support_point(
    polygon: Polygon,
    contact: Point,
    direction: np.ndarray | None,
    *,
    distance: float,
    tolerance: float,
) -> np.ndarray | None:
    if direction is None or distance <= 0.0:
        return None

    contact_xy = np.asarray(contact.coords[0], dtype=float)
    point_tolerance = max(tolerance, 1e-9)
    for scale in (1.0, 0.75, 0.5, 0.25):
        candidate_xy = contact_xy + direction * (distance * scale)
        candidate_point = Point(float(candidate_xy[0]), float(candidate_xy[1]))
        if polygon.buffer(point_tolerance).covers(candidate_point):
            return candidate_xy
    return None


def _apply_local_connector_bridge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    anchor_pairs: Sequence[tuple[int, Point, int, Point]],
    radius: float,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
    corridor_width_factor: float = 0.4,
    allow_patch_union_cleanup: bool = False,
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or target_scale <= 0 or len(subset_polygons) < 2 or not anchor_pairs:
        return None

    overlap_envelope = _local_offset_overlap_envelope(
        subset_polygons,
        radius=radius,
        grid=grid,
    )
    if overlap_envelope is None or overlap_envelope.is_empty:
        return None

    corridor_half_width = _connector_half_width(
        radius=radius,
        grid=grid,
        width_factor=corridor_width_factor,
    )
    corridors: list[BaseGeometry] = []
    for left_index, left_contact, right_index, right_contact in anchor_pairs:
        if left_index < 0 or right_index < 0:
            continue
        if left_index >= len(subset_polygons) or right_index >= len(subset_polygons):
            continue
        left_polygon = subset_polygons[left_index]
        right_polygon = subset_polygons[right_index]
        left_direction = _polygon_inward_direction_at_contact(
            left_polygon,
            left_contact,
            tolerance=max(radius * 1e-3, grid, 1e-9),
        )
        right_direction = _polygon_inward_direction_at_contact(
            right_polygon,
            right_contact,
            tolerance=max(radius * 1e-3, grid, 1e-9),
        )
        left_support = _inset_contact_support_point(
            left_polygon,
            left_contact,
            left_direction,
            distance=radius,
            tolerance=grid,
        )
        right_support = _inset_contact_support_point(
            right_polygon,
            right_contact,
            right_direction,
            distance=radius,
            tolerance=grid,
        )
        if left_support is None or right_support is None:
            continue
        try:
            corridor = LineString(
                [tuple(left_support), tuple(right_support)]
            ).buffer(
                corridor_half_width,
                quad_segs=1,
                cap_style=BufferCapStyle.flat,
                join_style=BufferJoinStyle.mitre,
                mitre_limit=1000.0,
            )
        except GEOSException:
            continue
        if corridor.is_empty:
            continue
        corridors.append(corridor)

    if not corridors:
        return None

    patch_geometry = overlap_envelope.intersection(unary_union(corridors))
    if patch_geometry.is_empty:
        return None

    reference_signature = _coverage_defect_signature(
        list(subset_polygons),
        target_scale=target_scale,
    )
    reference_union = unary_union(subset_polygons)
    patch_union = unary_union([*subset_polygons, patch_geometry])
    best_candidate: tuple[
        tuple[float, ...],
        list[Polygon],
        list[list[int]],
    ] | None = None

    def consider_candidate(
        candidate: tuple[list[Polygon], list[list[int]]] | None,
    ) -> None:
        nonlocal best_candidate
        if candidate is None:
            return

        candidate_polygons, candidate_sources = candidate
        candidate_signature = _coverage_defect_signature(
            candidate_polygons,
            target_scale=target_scale,
        )
        if not _coverage_signature_improves(
            reference_signature,
            candidate_signature,
            grid=grid,
            target_scale=target_scale,
        ):
            return

        difference_metrics = _difference_area_metrics(
            reference_union,
            unary_union(candidate_polygons),
        )
        score = _contact_resolution_candidate_score(
            candidate_signature,
            difference_metrics,
            target_scale=target_scale,
        )
        if best_candidate is None or score < best_candidate[0]:
            best_candidate = (score, candidate_polygons, candidate_sources)

    consider_candidate(
        _apply_local_pair_merge_operator(
            subset_polygons,
            subset_sources,
            radius=radius,
            grid=grid,
            diagnostics=diagnostics,
            patch_geometry=patch_union,
        )
    )

    if allow_patch_union_cleanup:
        for cleanup_radius in _iter_connector_cleanup_radii(radius=radius, grid=grid):
            for operator in ("open_close", "close_open"):
                consider_candidate(
                    _apply_local_patch_union_operator(
                        subset_polygons,
                        subset_sources,
                        radius=cleanup_radius,
                        operator=operator,
                        grid=grid,
                        diagnostics=diagnostics,
                        merge_components=True,
                        patch_geometry=patch_union,
                    )
                )

    if best_candidate is None:
        return None

    return best_candidate[1], best_candidate[2]


def _apply_local_point_touch_bridge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or len(subset_polygons) < 2:
        return None

    point_tolerance = max(grid, radius) * 1e-3
    tree = STRtree(subset_polygons)
    anchor_pairs: list[tuple[int, Point, int, Point]] = []
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

            for contact in _extract_point_contacts(boundary_intersection):
                anchor_pairs.append((index, contact, other_index, contact))

    if not anchor_pairs:
        return None

    return _select_best_connector_bridge_candidate(
        subset_polygons,
        subset_sources,
        radius=radius,
        target_scale=target_scale,
        grid=grid,
        diagnostics=diagnostics,
        anchor_pairs=anchor_pairs,
        variants=_iter_point_touch_connector_variants(
            radius=radius,
            grid=grid,
        ),
    )


def _apply_local_close_pair_bridge_operator(
    subset_polygons: Sequence[Polygon],
    subset_sources: Sequence[Sequence[int]],
    *,
    radius: float,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
) -> tuple[list[Polygon], list[list[int]]] | None:
    if radius <= 0 or len(subset_polygons) < 2:
        return None

    anchor_pairs: list[tuple[int, Point, int, Point]] = []
    point_tolerance = max(grid, radius) * 1e-3
    if len(subset_polygons) == 2:
        overlap_envelope = _local_offset_overlap_envelope(
            subset_polygons,
            radius=radius,
            grid=grid,
        )
        if overlap_envelope is not None and not overlap_envelope.is_empty:
            seen_anchor_keys: set[
                tuple[tuple[int, int], tuple[int, int]]
            ] = set()
            for component in shapely.get_parts(overlap_envelope):
                if component.is_empty:
                    continue
                left_contact = _nearest_boundary_contact_point(
                    subset_polygons[0],
                    component,
                    tolerance=max(grid, 1e-9),
                )
                right_contact = _nearest_boundary_contact_point(
                    subset_polygons[1],
                    component,
                    tolerance=max(grid, 1e-9),
                )
                if left_contact is None or right_contact is None:
                    continue
                left_key = tuple(
                    int(round(value / max(grid, 1e-9))) for value in left_contact.coords[0]
                )
                right_key = tuple(
                    int(round(value / max(grid, 1e-9))) for value in right_contact.coords[0]
                )
                anchor_key = (left_key, right_key)
                if anchor_key in seen_anchor_keys:
                    continue
                seen_anchor_keys.add(anchor_key)
                anchor_pairs.append((0, left_contact, 1, right_contact))

    if not anchor_pairs:
        tree = STRtree(subset_polygons)
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
                    closest = shapely.shortest_line(polygon, other)
                except GEOSException:
                    continue
                coords = list(closest.coords)
                if len(coords) < 2:
                    continue
                anchor_pairs.append(
                    (index, Point(coords[0]), other_index, Point(coords[-1]))
                )

    if not anchor_pairs:
        return None

    return _select_best_connector_bridge_candidate(
        subset_polygons,
        subset_sources,
        radius=radius,
        target_scale=target_scale,
        grid=grid,
        diagnostics=diagnostics,
        anchor_pairs=anchor_pairs,
        variants=((0.55, True),),
    )


def _iter_pair_issue_shrink_radii(
    *,
    tolerance: float,
    grid: float,
) -> tuple[float, ...]:
    if tolerance <= 0:
        return ()

    radii: list[float] = []
    for factor in (0.5, 0.625, 0.75, 1.0):
        radius = max(grid, tolerance * factor)
        if any(abs(radius - existing) <= 1e-12 for existing in radii):
            continue
        radii.append(radius)
    return tuple(radii)


def _iter_pair_issue_bridge_radii(
    *,
    issue_kind: Literal["point", "close"],
    distance: float,
    tolerance: float,
    grid: float,
) -> tuple[float, ...]:
    if tolerance <= 0:
        return ()

    radii: list[float] = []
    if issue_kind == "point":
        factors = (0.0625, 0.125, 0.25, 0.5)
        for factor in factors:
            radius = max(grid, tolerance * factor)
            if any(abs(radius - existing) <= 1e-12 for existing in radii):
                continue
            radii.append(radius)
        return tuple(radii)

    for factor in (0.5, 0.625, 0.75, 1.0):
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

    pair_indices = _courtyard_pair_indices(polygons)
    if pair_indices is None:
        return None

    larger_index, smaller_index = pair_indices
    larger = polygons[larger_index]
    smaller = polygons[smaller_index]

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


def _pair_issue_shrink_budget_context(
    polygons: Sequence[Polygon],
    *,
    tolerance: float,
    grid: float,
) -> tuple[BaseGeometry, float] | None:
    if len(polygons) != 2 or tolerance <= 0:
        return None

    smaller = min(polygons, key=lambda polygon: float(polygon.area))
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
) -> list[tuple[tuple[int, ...], str, list[Polygon], list[list[int]]]]:
    candidates: list[tuple[tuple[int, ...], str, list[Polygon], list[list[int]]]] = []
    if tolerance <= 0 or len(subset_polygons) < 2:
        return candidates

    seen_signatures: set[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]] = set()

    def append_candidate(
        affected_indices: Sequence[int],
        operator: str,
        replacement_polygons: Sequence[Polygon],
        replacement_sources: Sequence[Sequence[int]],
    ) -> None:
        candidate_polygons, candidate_sources = _replace_polygon_subset(
            subset_polygons,
            subset_sources,
            affected_indices=affected_indices,
            replacement_polygons=replacement_polygons,
            replacement_sources=replacement_sources,
        )
        signature = (
            _polygon_sequence_key(candidate_polygons),
            tuple(tuple(indices) for indices in candidate_sources),
        )
        if signature in seen_signatures:
            return
        seen_signatures.add(signature)
        candidates.append(
            (
                tuple(sorted(int(index) for index in affected_indices)),
                operator,
                candidate_polygons,
                candidate_sources,
            )
        )

    for cluster_indices in _point_touch_clusters(
        subset_polygons,
        target_scale=tolerance,
        cache=cache,
    ):
        if len(cluster_indices) < 2:
            continue
        cluster_polygons = [subset_polygons[index] for index in cluster_indices]
        cluster_sources = [subset_sources[index] for index in cluster_indices]
        for bridge_radius in _iter_pair_issue_bridge_radii(
            issue_kind="point",
            distance=0.0,
            tolerance=tolerance,
            grid=grid,
        ):
            candidate = _apply_local_point_touch_bridge_operator(
                cluster_polygons,
                cluster_sources,
                radius=bridge_radius,
                target_scale=tolerance,
                grid=grid,
                diagnostics=diagnostics,
            )
            if candidate is None:
                continue
            repaired_cluster_polygons, repaired_cluster_sources = candidate
            append_candidate(
                cluster_indices,
                f"coverage_pair_issue_point_cluster_{bridge_radius:.3f}",
                repaired_cluster_polygons,
                repaired_cluster_sources,
            )

    for left_index, right_index, distance, issue_kind in _cached_pair_issue_candidates(
        cache,
        subset_polygons,
        target_scale=tolerance,
    ):
        pair_polygons = [subset_polygons[left_index], subset_polygons[right_index]]
        pair_sources = [subset_sources[left_index], subset_sources[right_index]]
        operator_candidates: list[
            tuple[str, tuple[list[Polygon], list[list[int]]] | None]
        ] = []
        if issue_kind == "point":
            for bridge_radius in _iter_pair_issue_bridge_radii(
                issue_kind=issue_kind,
                distance=distance,
                tolerance=tolerance,
                grid=grid,
            ):
                operator_candidates.append(
                    (
                        f"coverage_pair_issue_point_{bridge_radius:.3f}",
                        _apply_local_point_touch_bridge_operator(
                            pair_polygons,
                            pair_sources,
                            radius=bridge_radius,
                            target_scale=tolerance,
                            grid=grid,
                            diagnostics=diagnostics,
                        ),
                    )
                )
        else:
            for bridge_radius in _iter_pair_issue_bridge_radii(
                issue_kind=issue_kind,
                distance=distance,
                tolerance=tolerance,
                grid=grid,
            ):
                operator_candidates.append(
                    (
                        f"coverage_pair_issue_bridge_{bridge_radius:.3f}",
                        _apply_local_close_pair_bridge_operator(
                            pair_polygons,
                            pair_sources,
                            radius=bridge_radius,
                            target_scale=tolerance,
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
            append_candidate(
                (left_index, right_index),
                operator,
                repaired_pair_polygons,
                repaired_pair_sources,
            )

    return candidates


def _repair_residual_pair_issues(
    input_polygons: list[Polygon],
    input_sources: list[list[int]],
    *,
    target_scale: float,
    grid: float,
    diagnostics: dict[str, Any],
    cache: _CoverageEvalCache | None = None,
) -> tuple[list[Polygon], list[list[int]], dict[str, int]] | None:
    if target_scale <= 0 or len(input_polygons) < 2:
        return None
    if cache is None:
        cache = _CoverageEvalCache()

    current_polygons = list(input_polygons)
    current_sources = [list(indices) for indices in input_sources]
    progress = False
    applied_counts: dict[str, int] = {}
    max_contact_shrink_area_loss = max(
        32.0 * target_scale * target_scale,
        128.0 * grid * grid,
        1e-9,
    )

    while True:
        pair_candidates = _cached_pair_issue_candidates(
            cache,
            current_polygons,
            target_scale=target_scale,
        )
        if not pair_candidates:
            break
        current_signature = _cached_coverage_defect_signature(
            cache,
            current_polygons,
            target_scale=target_scale,
        )

        repaired_this_round = False
        for left_index, right_index, distance, issue_kind in pair_candidates:
            pair_polygons = [current_polygons[left_index], current_polygons[right_index]]
            pair_sources = [current_sources[left_index], current_sources[right_index]]
            budget_context = _pair_issue_shrink_budget_context(
                pair_polygons,
                tolerance=target_scale,
                grid=grid,
            )
            if budget_context is None:
                continue

            reference_signature = _cached_coverage_defect_signature(
                cache,
                pair_polygons,
                target_scale=target_scale,
            )
            reference_union = _cached_union(cache, pair_polygons)
            edit_zone, shrink_area_budget = budget_context
            operator_candidates: list[
                tuple[str, tuple[list[Polygon], list[list[int]]] | None]
            ] = []
            if issue_kind == "point":
                for bridge_radius in _iter_pair_issue_bridge_radii(
                    issue_kind=issue_kind,
                    distance=distance,
                    tolerance=target_scale,
                    grid=grid,
                ):
                    operator_candidates.append(
                        (
                            f"coverage_pair_issue_point_residual_{bridge_radius:.3f}",
                            _apply_local_point_touch_bridge_operator(
                                pair_polygons,
                                pair_sources,
                                radius=bridge_radius,
                                target_scale=target_scale,
                                grid=grid,
                                diagnostics=diagnostics,
                            ),
                        )
                    )
            else:
                for bridge_radius in _iter_pair_issue_bridge_radii(
                    issue_kind=issue_kind,
                    distance=distance,
                    tolerance=target_scale,
                    grid=grid,
                ):
                    operator_candidates.append(
                        (
                            f"coverage_pair_issue_bridge_residual_{bridge_radius:.3f}",
                            _apply_local_close_pair_bridge_operator(
                                pair_polygons,
                                pair_sources,
                                radius=bridge_radius,
                                target_scale=target_scale,
                                grid=grid,
                                diagnostics=diagnostics,
                            ),
                        )
                    )
            for shrink_radius in _iter_pair_issue_shrink_radii(
                tolerance=target_scale,
                grid=grid,
            ):
                operator_candidates.append(
                    (
                        f"coverage_pair_issue_shrink_residual_{shrink_radius:.3f}",
                        _apply_local_smaller_polygon_shrink_operator(
                            pair_polygons,
                            pair_sources,
                            radius=shrink_radius,
                            grid=grid,
                            diagnostics=diagnostics,
                        ),
                    )
                )
            best_pair_candidate: tuple[
                list[Polygon],
                list[list[int]],
                _CoverageDefectSignature,
                dict[str, float],
                str,
            ] | None = None
            best_pair_score: tuple[float, ...] | None = None
            for operator_name, candidate in operator_candidates:
                if candidate is None:
                    continue
                candidate_polygons, candidate_sources = candidate
                candidate_pair_signature = _cached_coverage_defect_signature(
                    cache,
                    candidate_polygons,
                    target_scale=target_scale,
                )
                if (
                    candidate_pair_signature.pair_issue_count
                    >= reference_signature.pair_issue_count
                ):
                    continue
                candidate_difference = _cached_difference_area_metrics(
                    cache,
                    pair_polygons,
                    candidate_polygons,
                )
                if (
                    "shrink" in operator_name
                    and candidate_difference["reference_minus_candidate_area"]
                    > max_contact_shrink_area_loss
                ):
                    continue
                candidate_union = _cached_union(cache, candidate_polygons)
                change_outside_edit_zone = _change_outside_edit_zone(
                    reference_union,
                    candidate_union,
                    edit_zone=edit_zone,
                )
                if change_outside_edit_zone > max(
                    float(edit_zone.area),
                    grid * grid,
                    1e-9,
                ):
                    continue
                extra_area_budget = max(
                    shrink_area_budget,
                    target_scale * target_scale,
                    0.1 * float(edit_zone.area),
                    16.0 * grid * grid,
                )
                if candidate_difference["reference_minus_candidate_area"] > shrink_area_budget:
                    continue
                if candidate_difference["candidate_minus_reference_area"] > extra_area_budget:
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
                full_candidate_polygons, full_candidate_sources = _stable_sort(
                    unaffected_polygons + candidate_polygons,
                    unaffected_sources + candidate_sources,
                )
                full_candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    full_candidate_polygons,
                    target_scale=target_scale,
                )
                if not _coverage_signature_improves(
                    current_signature,
                    full_candidate_signature,
                    grid=grid,
                    target_scale=target_scale,
                ):
                    continue
                full_difference_metrics = _cached_difference_area_metrics(
                    cache,
                    input_polygons,
                    full_candidate_polygons,
                )
                candidate_score = (
                    _prefer_fill_only_point_contact_candidate(
                        operator_name,
                        full_difference_metrics,
                        target_scale=target_scale,
                        grid=grid,
                    ),
                    *_contact_resolution_candidate_score(
                        full_candidate_signature,
                        full_difference_metrics,
                        target_scale=target_scale,
                    ),
                )
                if (
                    best_pair_candidate is None
                    or best_pair_score is None
                    or candidate_score < best_pair_score
                ):
                    best_pair_candidate = (
                        full_candidate_polygons,
                        full_candidate_sources,
                        full_candidate_signature,
                        full_difference_metrics,
                        operator_name,
                    )
                    best_pair_score = candidate_score

            if best_pair_candidate is None:
                continue

            (
                replacement_polygons,
                replacement_sources,
                _,
                _,
                applied_operator,
            ) = best_pair_candidate
            current_polygons = replacement_polygons
            current_sources = replacement_sources
            repaired_this_round = True
            progress = True
            applied_counts[applied_operator] = applied_counts.get(applied_operator, 0) + 1
            break

        if not repaired_this_round:
            break

    if not progress:
        return None
    return current_polygons, current_sources, applied_counts


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
    diagnostics.setdefault("coverage_residual_polygon_target_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_attempt_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_normalization_failed_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_contract_reject_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_drift_reject_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_viable_candidate_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_no_candidate_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_changed_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_overlap_reject_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_improvement_reject_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_final_contract_reject_count", 0)
    diagnostics.setdefault("coverage_residual_polygon_final_drift_reject_count", 0)

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
        diagnostics["coverage_residual_polygon_target_count"] += 1

        best_candidate: Polygon | None = None
        best_metrics: dict[str, float] | None = None
        best_score: tuple[float, float, float, float] | None = None
        for factor in (0.75, 1.0, 1.25):
            diagnostics["coverage_residual_polygon_attempt_count"] += 1
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
                diagnostics["coverage_residual_polygon_normalization_failed_count"] += 1
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
                diagnostics["coverage_residual_polygon_contract_reject_count"] += 1
                continue

            difference_metrics = _difference_area_metrics(polygon, normalized)
            if not _difference_metrics_have_small_fidelity_drift(
                difference_metrics,
                edit_zone_area=float(polygon.area),
                target_scale=tolerance,
                grid=grid,
            ):
                diagnostics["coverage_residual_polygon_drift_reject_count"] += 1
                continue
            diagnostics["coverage_residual_polygon_viable_candidate_count"] += 1

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
            diagnostics["coverage_residual_polygon_no_candidate_count"] += 1
            continue

        candidate_polygons[index] = best_candidate
        changed_count += 1
        changed_area += float(polygon.area)

    if changed_count == 0:
        return None
    diagnostics["coverage_residual_polygon_changed_count"] += changed_count

    if _cached_overlap_area(cache, candidate_polygons) > max(grid * grid, 1e-9):
        diagnostics["coverage_residual_polygon_overlap_reject_count"] += 1
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
        diagnostics["coverage_residual_polygon_improvement_reject_count"] += 1
        return None
    if not _coverage_signature_satisfies_scale_contract(
        candidate_signature,
        target_scale=tolerance,
        grid=grid,
    ):
        diagnostics["coverage_residual_polygon_final_contract_reject_count"] += 1
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
        diagnostics["coverage_residual_polygon_final_drift_reject_count"] += 1
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


def _point_touch_clusters(
    polygons: Sequence[Polygon],
    *,
    target_scale: float,
    cache: _CoverageEvalCache | None = None,
) -> list[list[int]]:
    point_touch_components: list[set[int]] = []
    for left_index, right_index, _, issue_kind in _cached_pair_issue_candidates(
        cache,
        polygons,
        target_scale=target_scale,
    ):
        if issue_kind != "point":
            continue
        point_touch_components.append({left_index, right_index})
    return _merge_index_clusters(point_touch_components)


def _replace_polygon_subset(
    polygons: Sequence[Polygon],
    source_map: Sequence[Sequence[int]],
    *,
    affected_indices: Sequence[int],
    replacement_polygons: Sequence[Polygon],
    replacement_sources: Sequence[Sequence[int]],
) -> tuple[list[Polygon], list[list[int]]]:
    affected_index_set = set(int(index) for index in affected_indices)
    unaffected_polygons = [
        polygon
        for index, polygon in enumerate(polygons)
        if index not in affected_index_set
    ]
    unaffected_sources = [
        indices
        for index, indices in enumerate(source_map)
        if index not in affected_index_set
    ]
    return _stable_sort(
        unaffected_polygons + list(replacement_polygons),
        unaffected_sources + [list(indices) for indices in replacement_sources],
    )


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
        current_signature = _cached_coverage_defect_signature(
            cache,
            current_polygons,
            target_scale=tolerance,
        )
        accepted = False
        best_candidate: tuple[
            list[Polygon],
            list[list[int]],
            str,
            tuple[float, ...],
        ] | None = None
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
            shrink_budget_context = _courtyard_pair_shrink_budget_context(
                subset_polygons,
                tolerance=tolerance,
                grid=grid,
            )
            operator_candidates: list[
                tuple[str, tuple[list[Polygon], list[list[int]]] | None]
            ] = []
            if issue_kind == "point":
                for bridge_radius in _iter_pair_issue_bridge_radii(
                    issue_kind=issue_kind,
                    distance=distance,
                    tolerance=tolerance,
                    grid=grid,
                ):
                    operator_candidates.append(
                        (
                            f"coverage_pair_issue_point_{bridge_radius:.3f}",
                            _apply_local_point_touch_bridge_operator(
                                subset_polygons,
                                subset_sources,
                                radius=bridge_radius,
                                target_scale=tolerance,
                                grid=grid,
                                diagnostics=diagnostics,
                            ),
                        )
                    )
            else:
                for bridge_radius in _iter_pair_issue_bridge_radii(
                    issue_kind=issue_kind,
                    distance=distance,
                    tolerance=tolerance,
                    grid=grid,
                ):
                    operator_candidates.append(
                        (
                            f"coverage_pair_issue_bridge_{bridge_radius:.3f}",
                            _apply_local_close_pair_bridge_operator(
                                subset_polygons,
                                subset_sources,
                                radius=bridge_radius,
                                target_scale=tolerance,
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
                full_candidate_polygons, full_candidate_sources = _stable_sort(
                    unaffected_polygons + replacement_polygons,
                    unaffected_sources + replacement_sources,
                )
                full_candidate_signature = _cached_coverage_defect_signature(
                    cache,
                    full_candidate_polygons,
                    target_scale=tolerance,
                )
                if not _coverage_signature_improves(
                    current_signature,
                    full_candidate_signature,
                    grid=grid,
                    target_scale=tolerance,
                ):
                    continue
                full_difference_metrics = _cached_difference_area_metrics(
                    cache,
                    polygons,
                    full_candidate_polygons,
                )
                candidate_score = _contact_resolution_candidate_score(
                    full_candidate_signature,
                    full_difference_metrics,
                    target_scale=tolerance,
                )
                if best_candidate is None or candidate_score < best_candidate[3]:
                    best_candidate = (
                        full_candidate_polygons,
                        full_candidate_sources,
                        operator,
                        candidate_score,
                    )

        if best_candidate is not None:
            current_polygons, current_sources, applied_operator, _ = best_candidate
            operator_applied[applied_operator] = (
                operator_applied.get(applied_operator, 0) + 1
            )
            patch_applied_count += 1
            accepted = True

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
                    _affected_indices,
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
                        _affected_indices,
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
    def _accept_optional_post_stage_candidate(
        reference_polygons: list[Polygon],
        candidate_polygons: list[Polygon],
    ) -> bool:
        overlap_tolerance = max(grid * grid, 1e-9)
        if _cached_overlap_area(cache, candidate_polygons) > overlap_tolerance:
            return False
        reference_signature = _cached_coverage_defect_signature(
            cache,
            reference_polygons,
            target_scale=min_feature_size,
        )
        candidate_signature = _cached_coverage_defect_signature(
            cache,
            candidate_polygons,
            target_scale=min_feature_size,
        )
        return _coverage_signature_not_worse(
            reference_signature,
            candidate_signature,
            grid=grid,
            target_scale=min_feature_size,
        )

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
    boundary_ring_contact_candidate = None
    if post_boundary_signature.ring_contact_count > 0:
        boundary_ring_contact_candidate = _regularize_coverage_ring_contacts(
            boundary_regularized_polygons,
            boundary_regularized_sources,
            min_segment_length=min_feature_size,
            grid=grid,
        )
    if boundary_ring_contact_candidate is not None:
        (
            (boundary_regularized_polygons, boundary_regularized_sources),
            boundary_ring_contact_stats,
        ) = boundary_ring_contact_candidate
        diagnostics["boundary_regularization_ring_contact_polygon_count"] = (
            boundary_ring_contact_stats["changed_count"]
        )
        diagnostics["boundary_regularization_ring_contact_component_count"] = (
            boundary_ring_contact_stats["component_count"]
        )
        diagnostics["boundary_regularization_ring_contact_failed_count"] = (
            boundary_ring_contact_stats["failed_count"]
        )
        diagnostics["boundary_regularization_ring_contact_area_delta"] = (
            boundary_ring_contact_stats["area_delta"]
        )
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
    source_coordinate_repaired_polygons, source_coordinate_repaired_sources = (
        _apply_local_polygon_repairs(
            source_coordinate_recovered_polygons,
            source_coordinate_recovered_sources,
            min_segment_length=min_feature_size,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=diagnostics,
            stage_prefix="post_recovery_local_defect_repair",
            enable_defect_operators=True,
            enable_simplify_operators=False,
        )
    )
    post_recovery_signature = _cached_coverage_defect_signature(
        cache,
        source_coordinate_repaired_polygons,
        target_scale=min_feature_size,
    )
    post_recovery_clearance_deficit = max(
        min_feature_size - (post_recovery_signature.min_clearance or 0.0),
        0.0,
    )
    should_post_recovery_regularize_clearance = (
        post_recovery_signature.short_edge_count > 0
        or post_recovery_signature.ring_contact_count > 0
        or post_recovery_signature.pair_issue_count > 0
        or post_recovery_clearance_deficit > max(grid, 1e-9)
    )
    diagnostics["post_recovery_clearance_regularization_skipped"] = (
        not should_post_recovery_regularize_clearance
    )
    if should_post_recovery_regularize_clearance:
        post_recovery_diagnostics = _empty_diagnostics(
            len(source_coordinate_repaired_polygons)
        )
        post_recovery_diagnostics["collect_stage_metrics"] = False
        post_recovery_diagnostics["enable_logging"] = False
        (
            post_recovery_regularized_polygons,
            post_recovery_regularized_sources,
        ) = _regularize_low_clearance_polygons(
            source_coordinate_repaired_polygons,
            source_coordinate_repaired_sources,
            min_clearance=min_feature_size,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=post_recovery_diagnostics,
        )
    else:
        post_recovery_regularized_polygons = source_coordinate_repaired_polygons
        post_recovery_regularized_sources = source_coordinate_repaired_sources
        post_recovery_diagnostics = _empty_diagnostics(
            len(source_coordinate_repaired_polygons)
        )
        post_recovery_diagnostics["collect_stage_metrics"] = False
        post_recovery_diagnostics["enable_logging"] = False
    diagnostics["geos_exception_count"] += post_recovery_diagnostics[
        "geos_exception_count"
    ]
    diagnostics["geos_exception_messages"].extend(
        post_recovery_diagnostics["geos_exception_messages"]
    )
    diagnostics["post_recovery_clearance_regularization_applied"] = (
        post_recovery_diagnostics["clearance_regularization_applied"]
    )
    diagnostics["post_recovery_clearance_regularization_candidate_count"] = (
        post_recovery_diagnostics["clearance_regularization_candidate_count"]
    )
    diagnostics["post_recovery_clearance_regularization_improved_count"] = (
        post_recovery_diagnostics["clearance_regularization_improved_count"]
    )
    diagnostics["post_recovery_clearance_regularization_failed_count"] = (
        post_recovery_diagnostics["clearance_regularization_failed_count"]
    )
    diagnostics["post_recovery_dropped_meshing_hole_count"] = (
        post_recovery_diagnostics.get("dropped_meshing_hole_count", 0)
    )
    (
        coverage_contact_regularized_polygons,
        coverage_contact_regularized_sources,
    ) = _regularize_coverage_contacts(
        post_recovery_regularized_polygons,
        post_recovery_regularized_sources,
        min_segment_length=min_feature_size,
        grid=grid,
        min_area=min_area,
        min_hole_area=min_hole_area,
        diagnostics=diagnostics,
        cache=cache,
    )
    (
        coverage_void_regularized_polygons,
        coverage_void_regularized_sources,
    ) = _regularize_meshing_hostile_voids(
        coverage_contact_regularized_polygons,
        coverage_contact_regularized_sources,
        min_segment_length=min_feature_size,
        grid=grid,
        min_hole_area=min_hole_area,
        diagnostics=diagnostics,
        cache=cache,
    )
    diagnostics["coverage_void_regularization_reverted"] = False
    if not _accept_optional_post_stage_candidate(
        coverage_contact_regularized_polygons,
        coverage_void_regularized_polygons,
    ):
        coverage_void_regularized_polygons = coverage_contact_regularized_polygons
        coverage_void_regularized_sources = coverage_contact_regularized_sources
        diagnostics["coverage_void_regularization_reverted"] = (
            diagnostics["coverage_void_regularization_applied"]
        )
        diagnostics["coverage_void_regularization_applied"] = False
        diagnostics["coverage_void_regularization_hole_cleanup_count"] = 0
        diagnostics["coverage_void_regularization_gap_patch_applied_count"] = 0
        diagnostics["coverage_void_regularization_notch_simplify_count"] = 0
        diagnostics["coverage_void_regularization_operator_applied"] = {}
    if apply_coverage_meshing_regularization:
        coverage_meshing_regularized_polygons, coverage_meshing_regularized_sources = (
            _regularize_coverage_for_meshing(
                coverage_void_regularized_polygons,
                coverage_void_regularized_sources,
                min_segment_length=min_feature_size,
                grid=grid,
                min_area=min_area,
                min_hole_area=min_hole_area,
                diagnostics=diagnostics,
                cache=cache,
            )
        )
    else:
        coverage_meshing_regularized_polygons = coverage_void_regularized_polygons
        coverage_meshing_regularized_sources = coverage_void_regularized_sources

    (
        final_shape_regularized_polygons,
        final_shape_regularized_sources,
    ) = _regularize_final_polygon_shapes(
        coverage_meshing_regularized_polygons,
        coverage_meshing_regularized_sources,
        min_segment_length=min_feature_size,
        grid=grid,
        min_hole_area=min_hole_area,
        diagnostics=diagnostics,
    )
    diagnostics["final_shape_regularization_reverted"] = False
    if not _accept_optional_post_stage_candidate(
        coverage_meshing_regularized_polygons,
        final_shape_regularized_polygons,
    ):
        final_shape_regularized_polygons = coverage_meshing_regularized_polygons
        final_shape_regularized_sources = coverage_meshing_regularized_sources
        diagnostics["final_shape_regularization_reverted"] = (
            diagnostics.get("final_shape_regularization_applied", False)
        )
        diagnostics["final_shape_regularization_applied"] = False
        diagnostics["final_shape_regularization_operator_applied"] = {}

    final_output_polygons, final_output_sources = _filter_small_output_polygons(
        final_shape_regularized_polygons,
        final_shape_regularized_sources,
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
        post_recovery_regularized_polygons=post_recovery_regularized_polygons,
        post_recovery_regularized_source_map=post_recovery_regularized_sources,
        coverage_contact_regularized_polygons=coverage_contact_regularized_polygons,
        coverage_contact_regularized_source_map=coverage_contact_regularized_sources,
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
        elif key.startswith("post_recovery_local_defect_repair_"):
            diagnostics[key] = value
        elif key.startswith("post_contact_local_defect_repair_"):
            diagnostics[key] = value
        elif key.startswith("coverage_contact_regularization_"):
            diagnostics[key] = value
        elif key.startswith("coverage_void_regularization_"):
            diagnostics[key] = value
        elif key.startswith("coverage_meshing_regularization_"):
            diagnostics[key] = value
        elif key.startswith("coverage_residual_polygon_"):
            diagnostics[key] = value
        elif key.startswith("final_shape_regularization_"):
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
    if identity_branch is None and global_branch is None and local_branch is None:
        raise ValueError("at least one post-coverage branch must be available")
    identity_score = (
        _post_coverage_branch_score(
            identity_branch,
            target_scale=target_scale,
        )
        if identity_branch is not None
        else None
    )
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

    scored_candidates: list[tuple[tuple[float, ...], _PostCoverageBranch]] = []
    if identity_branch is not None and identity_score is not None:
        scored_candidates.append(
            (
                (
                    float(
                        _coverage_signature_contract_priority(
                            identity_branch.final_signature,
                            target_scale=target_scale,
                            grid=grid,
                        )
                    ),
                    *identity_score,
                ),
                identity_branch,
            )
        )
    if global_branch is not None and global_score is not None:
        scored_candidates.append(
            (
                (
                    float(
                        _coverage_signature_contract_priority(
                            global_branch.final_signature,
                            target_scale=target_scale,
                            grid=grid,
                        )
                    ),
                    *global_score,
                ),
                global_branch,
            )
        )
    if local_branch is not None and local_score is not None:
        scored_candidates.append(
            (
                (
                    float(
                        _coverage_signature_contract_priority(
                            local_branch.final_signature,
                            target_scale=target_scale,
                            grid=grid,
                        )
                    ),
                    *local_score,
                ),
                local_branch,
            )
        )

    chosen_branch = min(scored_candidates, key=lambda item: item[0])[1]
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
        needs_ring_contact_repair = enable_defect_operators and (
            polygon_signature.ring_contact_count > 0
        )
        needs_clearance_repair = enable_defect_operators and (
            polygon_signature.clearance_deficit > max(grid, 1e-9)
        )
        needs_short_edge_repair = enable_simplify_operators and (
            polygon_signature.short_edge_count > 0
        )
        if not (
            needs_ring_contact_repair
            or needs_clearance_repair
            or needs_short_edge_repair
        ):
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        candidate_count += 1
        accepted_candidate: Polygon | None = None
        local_candidates: list[_RepairCandidate] = []

        def append_candidate(candidate: _RepairCandidate | None) -> None:
            if candidate is None:
                return
            candidate_signature = _polygon_defect_signature(
                candidate.polygon,
                target_scale=min_segment_length,
            )
            if candidate_signature.short_edge_count == 0:
                local_candidates.append(candidate)
                return
            if enable_simplify_operators:
                refined = _iteratively_open_polygon_short_edges(
                    candidate.polygon,
                    target_scale=min_segment_length,
                    grid=grid,
                    diagnostics=diagnostics,
                    operator_prefix=f"{candidate.operator}_short_edge_angle_open",
                    initial_edit_zone=candidate.edit_zone,
                    area_balance_budget_override=candidate.area_balance_budget_override,
                )
                if refined is not None:
                    local_candidates.append(refined)
            local_candidates.append(candidate)

        if needs_ring_contact_repair:
            ring_contact_connector_candidate = _try_ring_contact_connector_fill(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            append_candidate(ring_contact_connector_candidate)

            ring_contact_candidate = _try_fill_ring_contact_vertices(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            append_candidate(ring_contact_candidate)

        if needs_clearance_repair:
            self_clearance_candidate = _try_polygon_self_clearance_connector_fill(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            append_candidate(self_clearance_candidate)

            courtyard_candidate = _try_close_courtyard_passage(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            append_candidate(courtyard_candidate)

            clearance_candidate = _try_polygon_clearance_opening(
                polygon,
                min_clearance=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
            )
            append_candidate(clearance_candidate)

            if enable_simplify_operators:
                clearance_simplify_candidate = _try_polygon_local_simplify(
                    polygon,
                    tolerance=min_segment_length * 3.0,
                    grid=grid,
                    diagnostics=diagnostics,
                )
                append_candidate(clearance_simplify_candidate)

        if needs_short_edge_repair:
            angle_open_candidate = _iteratively_open_polygon_short_edges(
                polygon,
                target_scale=min_segment_length,
                grid=grid,
                diagnostics=diagnostics,
                operator_prefix="short_edge_angle_open",
            )
            append_candidate(angle_open_candidate)

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
                append_candidate(simplify_candidate)

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
    applied_count = 0
    failed_count = 0
    candidate_polygons: list[Polygon] = []
    candidate_sources: list[list[int]] = []
    simplify_tolerances = (
        min_clearance * 0.5,
        min_clearance * 0.75,
        min_clearance,
    )
    acute_tip_min_span = max(2.0 * min_clearance, 8.0 * grid, 1.0e-9)

    def _normalize_clearance_candidate(
        candidate_polygon: Polygon,
    ) -> Polygon | None:
        cleaned = _remove_meshing_hostile_holes(
            candidate_polygon,
            min_hole_area=min_hole_area,
            min_clearance=min_clearance,
            diagnostics=diagnostics,
        )
        cleaned_parts = _canonicalize(cleaned, grid, diagnostics)
        if len(cleaned_parts) != 1 or cleaned_parts[0].area < min_area:
            return None
        return orient(cleaned_parts[0], sign=1.0)

    for polygon, indices in zip(polygons, source_map):
        reference_signature = _polygon_defect_signature(
            polygon,
            target_scale=min_clearance,
        )
        try:
            clearance_before = float(shapely.minimum_clearance(polygon))
        except GEOSException as exc:
            _record_geos_exception(diagnostics, "minimum_clearance", exc)
            clearance_before = np.inf

        needs_ring_contact_repair = reference_signature.ring_contact_count > 0
        if (
            (not np.isfinite(clearance_before) or clearance_before + 1e-12 >= min_clearance)
            and not needs_ring_contact_repair
        ):
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        candidate_count += 1
        best_candidate: Polygon | None = None
        best_score: (
            tuple[int, int, float, int, float, int, float, int, float, float] | None
        ) = None
        working_polygon = polygon

        normalized_identity = _normalize_clearance_candidate(polygon)
        if (
            normalized_identity is not None
            and not normalized_identity.equals_exact(polygon, tolerance=0.0)
        ):
            working_polygon = normalized_identity

        reference_acute_tip_count, reference_acute_tip_span = _polygon_acute_tip_metrics(
            working_polygon,
            min_tip_span=acute_tip_min_span,
        )

        def consider_candidate(candidate: Polygon | None) -> None:
            nonlocal best_candidate, best_score
            if candidate is None:
                return
            candidate_signature = _polygon_defect_signature(
                candidate,
                target_scale=min_clearance,
            )
            if not _signature_not_worse(
                reference_signature,
                candidate_signature,
                grid=grid,
            ):
                return
            clearance_contract_unmet = int(
                not _signature_satisfies_scale_contract(
                    candidate_signature,
                    target_scale=min_clearance,
                    grid=grid,
                )
            )
            acute_tip_growth = 0.0
            acute_tip_count_growth = 0
            if clearance_contract_unmet:
                (
                    candidate_acute_tip_count,
                    candidate_acute_tip_span,
                ) = _polygon_acute_tip_metrics(
                    candidate,
                    min_tip_span=acute_tip_min_span,
                )
                acute_tip_growth = max(
                    candidate_acute_tip_span - reference_acute_tip_span,
                    0.0,
                )
                acute_tip_count_growth = max(
                    candidate_acute_tip_count - reference_acute_tip_count,
                    0,
                )
            score = (
                candidate_signature.ring_contact_count,
                clearance_contract_unmet,
                acute_tip_growth,
                acute_tip_count_growth,
                candidate_signature.clearance_deficit,
                candidate_signature.short_edge_count,
                -(candidate_signature.min_edge_length or 0.0),
                candidate_signature.vertex_count,
                abs(float(candidate.area - polygon.area)),
                candidate.area,
            )
            if best_score is None or score < best_score:
                best_score = score
                best_candidate = candidate

        def consider_repair_candidate(repair_candidate: _RepairCandidate | None) -> None:
            if repair_candidate is None:
                return
            normalized_candidate = _normalize_clearance_candidate(
                repair_candidate.polygon
            )
            consider_candidate(normalized_candidate)
            for tolerance in simplify_tolerances:
                simplify_candidate = _try_polygon_local_simplify(
                    repair_candidate.polygon,
                    tolerance=tolerance,
                    grid=grid,
                    diagnostics=diagnostics,
                )
                if simplify_candidate is None:
                    continue
                consider_candidate(
                    _normalize_clearance_candidate(simplify_candidate.polygon)
                )

        consider_candidate(working_polygon)

        opened_parts = _apply_opening(working_polygon, radius, grid, diagnostics)
        if len(opened_parts) == 1:
            closed_parts = _apply_closing(opened_parts[0], radius, grid, diagnostics)
            if len(closed_parts) == 1:
                consider_candidate(
                    _normalize_clearance_candidate(closed_parts[0])
                )

        relaxed_self_clearance_candidate = _try_polygon_self_clearance_connector_fill(
            working_polygon,
            min_clearance=min_clearance,
            grid=grid,
            diagnostics=diagnostics,
            require_signature_improvement=False,
        )
        consider_repair_candidate(relaxed_self_clearance_candidate)

        local_repair_diagnostics = _empty_diagnostics(1)
        local_repair_diagnostics["collect_stage_metrics"] = False
        local_repair_diagnostics["enable_logging"] = False
        repaired_polygons, _ = _apply_local_polygon_repairs(
            [working_polygon],
            [indices],
            min_segment_length=min_clearance,
            grid=grid,
            min_area=min_area,
            min_hole_area=min_hole_area,
            diagnostics=local_repair_diagnostics,
            stage_prefix="clearance_regularization_local_defect_repair",
            enable_defect_operators=True,
            enable_simplify_operators=True,
        )
        if len(repaired_polygons) == 1 and not repaired_polygons[0].equals_exact(
            working_polygon,
            tolerance=0.0,
        ):
            consider_candidate(
                _normalize_clearance_candidate(repaired_polygons[0])
            )

        if best_candidate is None:
            failed_count += 1
            candidate_polygons.append(polygon)
            candidate_sources.append(indices)
            continue

        try:
            clearance_after = float(shapely.minimum_clearance(best_candidate))
        except GEOSException:
            clearance_after = clearance_before
        if not best_candidate.equals_exact(polygon, tolerance=0.0):
            applied_count += 1
        if clearance_after > clearance_before + 1e-12:
            improved_count += 1
        candidate_polygons.append(best_candidate)
        candidate_sources.append(indices)

    overlap_area = _coverage_overlap_area(candidate_polygons)
    diagnostics["clearance_regularization_overlap_area"] = overlap_area
    diagnostics["clearance_regularization_candidate_count"] = candidate_count
    diagnostics["clearance_regularization_improved_count"] = improved_count
    diagnostics["clearance_regularization_applied_count"] = applied_count
    diagnostics["clearance_regularization_failed_count"] = failed_count
    if applied_count == 0:
        diagnostics["clearance_regularization_overlap_area"] = 0.0
        diagnostics["clearance_regularization_applied"] = False
        return polygons, source_map

    if overlap_area > overlap_tolerance:
        diagnostics["clearance_regularization_applied"] = False
        return polygons, source_map

    diagnostics["clearance_regularization_applied"] = True
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
    11. regularize residual polygon-pair contacts coverage-wide at the
        declared scale using local bridge operators
    12. reserve the final mesher-specific stage for any residual short-edge
        cleanup that still remains

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
    diagnostics["stage_seconds"] = {}

    _log_conditioning_start(
        len(polygons),
        options,
        grid=grid,
        output_grid=output_grid,
        opening_radius=r_open,
        closing_radius=r_close,
    )

    stage_started_at = time.perf_counter()
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

    atomic_polygons, atomic_sources = _stable_sort(atomic_polygons, atomic_sources)

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
    _record_stage_seconds(diagnostics, "atomic_input", stage_started_at)

    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "opened", stage_started_at)

    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "regularized_groups", stage_started_at)

    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "reconstructed", stage_started_at)

    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "presimplify", stage_started_at)

    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "local_defect_repaired", stage_started_at)

    stage_started_at = time.perf_counter()
    coverage_eval_cache = _CoverageEvalCache()
    coverage_reference_union = _cached_union(coverage_eval_cache, final_polygons)
    coverage_simplify_tolerance = _derive_coverage_simplify_tolerance(options)
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
    diagnostics["coverage_simplify_contract_fallback_triggered"] = False
    diagnostics["coverage_simplify_contract_fallback_evaluated"] = []
    finalized_branches = {
        candidate.label: _evaluate_post_coverage_branch(
            candidate,
            reference_union=coverage_reference_union,
            source_lookup=source_lookup,
            raw_support_union=raw_support_union,
            min_feature_size=options.min_feature_size,
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
    if not _coverage_signature_satisfies_scale_contract(
        chosen_branch.final_signature,
        target_scale=meshing_scale,
        grid=output_grid,
    ):
        if local_candidate is None:
            diagnostics["coverage_simplify_local_candidate_attempted"] = True
            local_candidate = _simplify_coverage_locally(
                final_polygons,
                final_sources,
                tolerance=coverage_simplify_tolerance,
                grid=output_grid,
                diagnostics=diagnostics,
                # The fast local pass intentionally skips the more expensive
                # cluster rescues. When the selected branch still misses the
                # declared meshing contract, evaluate the full deterministic
                # local candidate set before failing strict mode.
                enable_patch_union_fallback=True,
                enable_pair_cluster_rescue=True,
                cache=coverage_eval_cache,
            )
            if local_candidate is not None:
                coverage_candidates.append(local_candidate)
        fallback_candidates = [
            candidate
            for candidate in coverage_candidates
            if candidate.label not in finalized_branches
        ]
        diagnostics["coverage_simplify_contract_fallback_triggered"] = True
        diagnostics["coverage_simplify_contract_fallback_evaluated"] = [
            candidate.label for candidate in fallback_candidates
        ]
        for candidate in fallback_candidates:
            finalized_branches[candidate.label] = _evaluate_post_coverage_branch(
                candidate,
                reference_union=coverage_reference_union,
                source_lookup=source_lookup,
                raw_support_union=raw_support_union,
                min_feature_size=options.min_feature_size,
                source_recovery_scale=options.min_feature_size,
                grid=output_grid,
                min_area=0.0,
                output_min_area=options.min_area,
                min_hole_area=options.min_hole_area,
                cache=coverage_eval_cache,
            )
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
    _record_stage_seconds(diagnostics, "coverage_simplified", stage_started_at)

    final_polygons = chosen_branch.source_reclaimed_polygons
    final_sources = chosen_branch.source_reclaimed_source_map
    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "source_reclaimed", stage_started_at)

    final_polygons = chosen_branch.small_component_absorbed_polygons
    final_sources = chosen_branch.small_component_absorbed_source_map
    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "small_component_absorbed", stage_started_at)

    final_polygons = chosen_branch.boundary_regularized_polygons
    final_sources = chosen_branch.boundary_regularized_source_map
    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "boundary_regularized", stage_started_at)

    final_polygons = chosen_branch.clearance_regularized_polygons
    final_sources = chosen_branch.clearance_regularized_source_map
    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(diagnostics, "clearance_regularized", stage_started_at)

    final_polygons = chosen_branch.source_coordinate_recovered_polygons
    final_sources = chosen_branch.source_coordinate_recovered_source_map
    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(
        diagnostics,
        "source_coordinate_recovered",
        stage_started_at,
    )

    final_polygons = chosen_branch.post_recovery_regularized_polygons
    final_sources = chosen_branch.post_recovery_regularized_source_map
    stage_started_at = time.perf_counter()
    post_recovery_regularized_metrics = _record_stage_metrics(
        diagnostics,
        "post_recovery_regularized",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if post_recovery_regularized_metrics is not None:
        _log_conditioning_stage(
            "post_recovery_regularized",
            post_recovery_regularized_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )
    _record_stage_seconds(diagnostics, "post_recovery_regularized", stage_started_at)

    final_polygons = chosen_branch.coverage_contact_regularized_polygons
    final_sources = chosen_branch.coverage_contact_regularized_source_map
    stage_started_at = time.perf_counter()
    coverage_contact_regularized_metrics = _record_stage_metrics(
        diagnostics,
        "coverage_contact_regularized",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    if coverage_contact_regularized_metrics is not None:
        _log_conditioning_stage(
            "coverage_contact_regularized",
            coverage_contact_regularized_metrics,
            short_edge_threshold=meshing_scale,
            enabled=options.enable_logging,
        )
    _record_stage_seconds(
        diagnostics,
        "coverage_contact_regularized",
        stage_started_at,
    )

    final_polygons = chosen_branch.coverage_meshing_regularized_polygons
    final_sources = chosen_branch.coverage_meshing_regularized_source_map
    stage_started_at = time.perf_counter()
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
    _record_stage_seconds(
        diagnostics,
        "coverage_meshing_regularized",
        stage_started_at,
    )

    final_polygons = chosen_branch.final_output_polygons
    final_sources = chosen_branch.final_output_source_map
    final_polygons, final_sources = _stable_sort(final_polygons, final_sources)
    diagnostics["output_count"] = len(final_polygons)
    diagnostics["overlap_area_after"] = _coverage_overlap_area(final_polygons)
    diagnostics["min_clearance_after"] = _minimum_clearance(final_polygons)
    stage_started_at = time.perf_counter()
    _record_stage_metrics(
        diagnostics,
        "final_output",
        final_polygons,
        short_edge_threshold=meshing_scale,
    )
    _record_stage_seconds(diagnostics, "final_output", stage_started_at)
    _log_conditioning_summary(
        diagnostics,
        short_edge_threshold=meshing_scale,
    )

    return ConditioningResult(
        polygons=final_polygons,
        source_map=final_sources,
        diagnostics=diagnostics,
    )
