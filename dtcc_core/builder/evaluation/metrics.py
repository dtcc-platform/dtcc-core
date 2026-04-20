"""Per-building evaluation metrics. Each function is a pure, standalone measurement."""
from __future__ import annotations

from typing import Optional

from dtcc_core.model.object.building import Building
from dtcc_core.builder.evaluation.dataset import GroundTruth


STAGE_OUTCOMES = (
    "success",
    "insufficient_points",
    "classification_failed",
    "complex_footprint",
    "geometry_construction_failed",
    "validation_failed",
    "unknown",
)


def stage_outcome(building: Building) -> str:
    """Return the terminal pipeline outcome for one building.

    Reads building.attributes['fallback_reason'] (written by build_lod2_buildings).
    Returns 'success' if an LoD2 geometry is present and no fallback_reason set,
    'unknown' if neither is present.
    """
    reason = building.attributes.get("fallback_reason")
    if reason:
        return str(reason)
    if building.lod2 is not None:
        return "success"
    return "unknown"


def roof_type_prediction(building: Building) -> Optional[str]:
    """Return the predicted roof type name written to building.attributes, or None."""
    pred = building.attributes.get("roof_type")
    return str(pred) if pred is not None else None


def roof_type_correct(building: Building, gt: GroundTruth) -> bool:
    """True if predicted roof type matches ground-truth label (case-sensitive)."""
    pred = roof_type_prediction(building)
    if pred is None:
        return False
    return pred == gt.roof_type


_DEFAULT_PLANE_COUNT = {
    "FLAT": 1,
    "GABLED": 2,
    "HIPPED": 4,
}


def detected_plane_count(diagnostics: dict) -> int:
    """Number of planes detected for this building, from the diagnostics record."""
    planes = diagnostics.get("planes") if diagnostics else None
    return 0 if planes is None else len(planes)


def expected_plane_count(gt: GroundTruth) -> Optional[int]:
    """Expected plane count given the ground-truth roof type, or None if unknown."""
    if gt.expected_plane_count is not None:
        return int(gt.expected_plane_count)
    return _DEFAULT_PLANE_COUNT.get(gt.roof_type)


def plane_count_delta(diagnostics: dict, gt: GroundTruth) -> Optional[int]:
    """Absolute difference between detected and expected plane count, or None."""
    expected = expected_plane_count(gt)
    if expected is None:
        return None
    return abs(detected_plane_count(diagnostics) - expected)


def point_coverage_ratio(diagnostics: dict, filtered_point_count: int) -> Optional[float]:
    """Fraction of filtered roof points assigned to a detected plane.

    Mirrors the `coverage_factor` term used in the rule-based classifier's
    confidences. Returns None if filtered_point_count is zero.
    """
    if filtered_point_count <= 0:
        return None
    planes = diagnostics.get("planes") or []
    inlier_total = sum(len(p.inliers) for p in planes)
    return min(1.0, inlier_total / filtered_point_count)


def plane_areas(planes) -> list[float]:
    """Return the list of plane areas in the order given."""
    return [float(p.area) for p in planes]


def top2_area_share(diagnostics: dict) -> Optional[float]:
    """Area of the two largest detected planes / total detected plane area.

    Parallels the hipped-classification area-ratio rule. Returns None when no
    planes are present or when all areas are zero.
    """
    planes = diagnostics.get("planes") or []
    areas = plane_areas(planes)
    total = sum(areas)
    if not areas or total <= 0.0:
        return None
    top2 = sum(sorted(areas, reverse=True)[:2])
    return top2 / total


def slope_symmetry_error(diagnostics: dict) -> Optional[float]:
    """Absolute tilt difference (degrees) between the two largest detected planes.

    Parallels the symmetric-pair test in the rule-based gabled classifier.
    Returns None if fewer than two planes. Uses RoofPlane.slope_deg.
    """
    planes = diagnostics.get("planes") or []
    if len(planes) < 2:
        return None
    sorted_by_area = sorted(planes, key=lambda p: p.area, reverse=True)
    top2 = sorted_by_area[:2]
    return abs(top2[0].slope_deg - top2[1].slope_deg)


def is_watertight(
    building: Building,
    tolerance: float = 0.01,
    diagnostics: Optional[dict] = None,
) -> bool:
    """True iff the building's LoD2 is watertight.

    Prefers `diagnostics["validated"]` when available (set by
    build_lod2_buildings), avoiding a redundant shell validation.
    Falls back to calling validate_shell on building.lod2 otherwise.
    """
    if diagnostics is not None and diagnostics.get("validated") is not None:
        return bool(diagnostics["validated"])

    from dtcc_core.builder.geometry.shell_validation import validate_shell
    if building.lod2 is None:
        return False
    valid, _ = validate_shell(building.lod2, tolerance)
    return bool(valid)


def _predicted_height(diagnostics: dict, attr: str) -> Optional[float]:
    clsr = diagnostics.get("classification") if diagnostics else None
    if clsr is None:
        return None
    v = getattr(clsr, attr, None)
    return None if v is None else float(v)


def ridge_height_error(diagnostics: dict, gt: GroundTruth) -> Optional[float]:
    """Absolute error (m) between predicted and ground-truth ridge height, or None."""
    if gt.ridge_height is None:
        return None
    pred = _predicted_height(diagnostics, "ridge_height")
    if pred is None:
        return None
    return abs(pred - gt.ridge_height)


def eave_height_error(diagnostics: dict, gt: GroundTruth) -> Optional[float]:
    """Absolute error (m) between predicted and ground-truth eave height, or None."""
    if gt.eave_height is None:
        return None
    pred = _predicted_height(diagnostics, "eave_height")
    if pred is None:
        return None
    return abs(pred - gt.eave_height)


def _project_xy(surface):
    from shapely.geometry import Polygon as ShapelyPolygon

    verts = surface.vertices
    if verts is None or len(verts) < 3:
        return None
    ring = [(float(x), float(y)) for x, y, *_ in verts]
    poly = ShapelyPolygon(ring)
    if not poly.is_valid or poly.area <= 0.0:
        poly = poly.buffer(0)
    return poly if poly.is_valid and poly.area > 0.0 else None


def _union_by_semantic(ms):
    from collections import defaultdict
    from shapely.ops import unary_union

    if ms is None or not ms.surfaces:
        return {}
    semantics = ms.semantics if ms.semantics is not None else [None] * len(ms.surfaces)
    bag = defaultdict(list)
    for surf, sem in zip(ms.surfaces, semantics):
        poly = _project_xy(surf)
        if poly is None:
            continue
        bag[sem].append(poly)
    return {k: unary_union(v) for k, v in bag.items() if v}


def semantic_overlap(predicted, truth) -> dict:
    """Per-semantic intersection-over-union of XY-projected surfaces.

    Returns a dict keyed by SurfaceSemantic (drawn from predicted ∪ truth) with
    values in [0, 1]. Missing semantic in either side → 0.0.
    """
    pred_unions = _union_by_semantic(predicted)
    truth_unions = _union_by_semantic(truth)

    keys = set(pred_unions) | set(truth_unions)
    out = {}
    for k in keys:
        if k is None:
            continue
        p = pred_unions.get(k)
        t = truth_unions.get(k)
        if p is None or t is None:
            out[k] = 0.0
            continue
        inter = p.intersection(t).area
        union = p.union(t).area
        out[k] = float(inter / union) if union > 0 else 0.0
    return out
