from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy.spatial import ConvexHull

from dtcc_core.model.enums import RoofType
from dtcc_core.builder.geometry_builders.roof_config import (
    RoofDetectionConfig,
    RoofPlane,
)


@dataclass
class ClassificationResult:
    roof_type: RoofType
    confidence: float
    roof_planes: List[RoofPlane]
    ridge_line: Optional[np.ndarray] = None  # (2, 3) two endpoints, or None
    eave_height: float = 0.0
    ridge_height: float = 0.0


def merge_planes(
    planes: List[RoofPlane],
    config: RoofDetectionConfig,
) -> List[RoofPlane]:
    """Merge over-segmented planes. Repeat-until-stable."""
    result = list(planes)
    threshold_cos = np.cos(np.radians(config.merge_normal_threshold))
    changed = True

    while changed:
        changed = False
        for i in range(len(result)):
            for j in range(i + 1, len(result)):
                cos_angle = np.dot(result[i].normal, result[j].normal)
                if cos_angle < threshold_cos:
                    continue
                dist = np.linalg.norm(result[i].centroid - result[j].centroid)
                if dist > config.merge_distance:
                    continue

                # Merge j into i
                merged_inliers = np.union1d(result[i].inliers, result[j].inliers)
                merged_boundary = np.vstack(
                    [result[i].boundary_3d, result[j].boundary_3d]
                )

                # Refit plane via SVD
                centroid = merged_boundary.mean(axis=0)
                _, _, vh = np.linalg.svd(merged_boundary - centroid)
                new_normal = vh[-1]
                if new_normal[2] < 0:
                    new_normal *= -1
                new_normal = new_normal / np.linalg.norm(new_normal)
                new_offset = np.dot(new_normal, centroid)

                try:
                    hull = ConvexHull(merged_boundary[:, :2])
                    new_boundary = merged_boundary[hull.vertices]
                except Exception:
                    new_boundary = merged_boundary

                result[i] = RoofPlane(
                    normal=new_normal,
                    offset=new_offset,
                    inliers=merged_inliers,
                    boundary_3d=new_boundary,
                )
                result.pop(j)
                changed = True
                break
            if changed:
                break

    return result


def classify_roof(
    planes: List[RoofPlane],
    footprint_area: float,
    total_points: int,
    config: RoofDetectionConfig,
    footprint_vertices: Optional[np.ndarray] = None,
    all_points: Optional[np.ndarray] = None,
) -> ClassificationResult:
    """Classify roof type from dominant planes. Returns best match.

    Parameters
    ----------
    footprint_vertices : np.ndarray, optional
        (N, 2) or (N, 3) footprint vertices. Used for hipped residual analysis.
    all_points : np.ndarray, optional
        (M, 3) all roof points. Used for hipped residual analysis.
    """
    # Filter by min area ratio
    dominant = [
        p
        for p in planes
        if p.area >= config.min_plane_area_ratio * footprint_area
    ]
    if not dominant:
        dominant = planes[:1] if planes else []

    candidates: List[ClassificationResult] = []

    flat_result = _try_flat(dominant, total_points, config)
    if flat_result:
        candidates.append(flat_result)

    gabled_result = _try_gabled(dominant, total_points, config)
    if gabled_result:
        candidates.append(gabled_result)
        # Try upgrading gabled to hipped via residual analysis
        if footprint_vertices is not None and all_points is not None:
            hipped_upgrade = _try_upgrade_gabled_to_hipped(
                gabled_result, all_points, footprint_vertices, total_points, config,
            )
            if hipped_upgrade:
                candidates.append(hipped_upgrade)

    hipped_result = _try_hipped(dominant, total_points, config)
    if hipped_result:
        candidates.append(hipped_result)

    if not candidates:
        return ClassificationResult(
            roof_type=RoofType.UNKNOWN,
            confidence=0.0,
            roof_planes=planes,
        )

    best = max(candidates, key=lambda r: r.confidence)
    if best.confidence < config.min_roof_confidence:
        return ClassificationResult(
            roof_type=RoofType.UNKNOWN,
            confidence=best.confidence,
            roof_planes=planes,
        )
    return best


def _try_flat(
    planes: List[RoofPlane],
    total_points: int,
    config: RoofDetectionConfig,
) -> Optional[ClassificationResult]:
    if not planes:
        return None

    max_slope = max(p.slope_deg for p in planes)
    if max_slope >= config.flat_angle_threshold:
        return None

    total_inliers = sum(len(p.inliers) for p in planes)
    coverage = total_inliers / max(total_points, 1)
    flatness = 1.0 - (max_slope / config.flat_angle_threshold)
    confidence = min(coverage, max(0.0, flatness))

    for p in planes:
        p.role = "flat_roof"

    eave_height = min(p.boundary_3d[:, 2].min() for p in planes)

    return ClassificationResult(
        roof_type=RoofType.FLAT,
        confidence=confidence,
        roof_planes=planes,
        eave_height=eave_height,
        ridge_height=eave_height,
    )


def _try_gabled(
    planes: List[RoofPlane],
    total_points: int,
    config: RoofDetectionConfig,
) -> Optional[ClassificationResult]:
    if len(planes) != 2:
        return None

    p1, p2 = planes
    slope_diff = abs(p1.slope_deg - p2.slope_deg)
    if slope_diff > config.slope_symmetry_tolerance:
        return None

    az_diff = abs(p1.azimuth_deg - p2.azimuth_deg)
    az_diff = min(az_diff, 360 - az_diff)
    if abs(az_diff - 180) > config.azimuth_opposition_tolerance:
        return None

    ridge_line = _intersect_planes(p1, p2)

    symmetry = max(0.0, 1.0 - (slope_diff / config.slope_symmetry_tolerance))
    total_inliers = len(p1.inliers) + len(p2.inliers)
    coverage = total_inliers / max(total_points, 1)
    pattern = 1.0
    confidence = min(symmetry, coverage, pattern)

    p1.role = "ridge_left"
    p2.role = "ridge_right"

    eave_height = min(p1.boundary_3d[:, 2].min(), p2.boundary_3d[:, 2].min())
    ridge_height = max(p1.boundary_3d[:, 2].max(), p2.boundary_3d[:, 2].max())

    return ClassificationResult(
        roof_type=RoofType.GABLED,
        confidence=confidence,
        roof_planes=[p1, p2],
        ridge_line=ridge_line,
        eave_height=eave_height,
        ridge_height=ridge_height,
    )


def _try_hipped(
    planes: List[RoofPlane],
    total_points: int,
    config: RoofDetectionConfig,
) -> Optional[ClassificationResult]:
    if len(planes) != 4:
        return None

    sorted_planes = sorted(planes, key=lambda p: p.area, reverse=True)
    total_roof_area = sum(p.area for p in sorted_planes)
    if total_roof_area == 0:
        return None

    main1, main2, end1, end2 = sorted_planes

    if main1.area / total_roof_area < config.hip_main_area_ratio:
        return None
    if main2.area / total_roof_area < config.hip_main_area_ratio:
        return None

    main_az_diff = abs(main1.azimuth_deg - main2.azimuth_deg)
    main_az_diff = min(main_az_diff, 360 - main_az_diff)
    if abs(main_az_diff - 180) > config.azimuth_opposition_tolerance:
        return None

    # Check hip ends are perpendicular to the main slope axis.
    # Use main1's azimuth as reference (not the average, which is meaningless
    # for opposing normals like 90 and 270 averaging to 180).
    main_ref_az = main1.azimuth_deg
    for end in [end1, end2]:
        end_diff = abs(end.azimuth_deg - main_ref_az)
        end_diff = min(end_diff, 360 - end_diff)
        # Hip ends should be ~90 degrees from main slope direction
        if abs(end_diff - 90) > config.hip_perpendicular_tolerance:
            return None

    slope_diff = abs(main1.slope_deg - main2.slope_deg)
    symmetry = max(0.0, 1.0 - (slope_diff / config.slope_symmetry_tolerance))
    total_inliers = sum(len(p.inliers) for p in sorted_planes)
    coverage = total_inliers / max(total_points, 1)
    pattern = 1.0
    confidence = min(symmetry, coverage, pattern)

    main1.role = "main_left"
    main2.role = "main_right"
    end1.role = "hip_end_1"
    end2.role = "hip_end_2"

    eave_height = min(p.boundary_3d[:, 2].min() for p in sorted_planes)
    ridge_height = max(p.boundary_3d[:, 2].max() for p in sorted_planes)

    return ClassificationResult(
        roof_type=RoofType.HIPPED,
        confidence=confidence,
        roof_planes=sorted_planes,
        ridge_line=_intersect_planes(main1, main2),
        eave_height=eave_height,
        ridge_height=ridge_height,
    )


def _try_upgrade_gabled_to_hipped(
    gabled_result: ClassificationResult,
    all_points: np.ndarray,
    footprint_vertices: np.ndarray,
    total_points: int,
    config: RoofDetectionConfig,
) -> Optional[ClassificationResult]:
    """Check if a gabled classification should be upgraded to hipped.

    Strategy: project the ridge line and footprint onto the ridge direction axis.
    If the ridge is significantly shorter than the footprint along that axis,
    the building likely has hip ends. Then verify that residual points at the
    building ends (beyond the ridge) slope inward rather than being vertical.
    """
    ridge_line = gabled_result.ridge_line
    if ridge_line is None:
        return None

    fp_2d = footprint_vertices[:, :2] if footprint_vertices.shape[1] >= 2 else footprint_vertices

    # Ridge direction in 2D
    ridge_dir = ridge_line[1, :2] - ridge_line[0, :2]
    ridge_len = np.linalg.norm(ridge_dir)
    if ridge_len < 1e-6:
        return None
    ridge_dir_unit = ridge_dir / ridge_len

    # Project footprint onto ridge axis to get building extent
    ridge_mid = (ridge_line[0, :2] + ridge_line[1, :2]) / 2
    fp_projections = np.dot(fp_2d - ridge_mid, ridge_dir_unit)
    fp_extent = fp_projections.max() - fp_projections.min()

    # Determine actual ridge extent from the inlier point distribution
    # The ridge endpoints are artificially extended -- instead measure how far
    # the two main planes' high points (near ridge height) actually span
    p1, p2 = gabled_result.roof_planes[:2]
    ridge_height = gabled_result.ridge_height
    eave_height = gabled_result.eave_height
    height_threshold = eave_height + (ridge_height - eave_height) * 0.7

    # Points near the ridge (in the top 30% of roof height)
    high_mask = all_points[:, 2] > height_threshold
    if high_mask.sum() < 5:
        return None
    high_pts = all_points[high_mask]
    high_projections = np.dot(high_pts[:, :2] - ridge_mid, ridge_dir_unit)
    ridge_extent = high_projections.max() - high_projections.min()

    # Hipped criterion: the high ridge zone should be noticeably shorter
    # than the building along the ridge axis
    ridge_ratio = ridge_extent / max(fp_extent, 1e-6)
    if ridge_ratio > 0.85:
        # Ridge spans most of the building -- pure gabled, not hipped
        return None

    # Check residual points at building ends
    p1, p2 = gabled_result.roof_planes[:2]

    # Get points beyond the high-ridge zone (the hip end zones)
    point_projections = np.dot(all_points[:, :2] - ridge_mid, ridge_dir_unit)
    ridge_min_proj = high_projections.min()
    ridge_max_proj = high_projections.max()

    end_mask_1 = point_projections < ridge_min_proj
    end_mask_2 = point_projections > ridge_max_proj
    end_point_count = end_mask_1.sum() + end_mask_2.sum()

    # If there are substantial points beyond the ridge, it's hipped
    if end_point_count < max(10, total_points * 0.05):
        return None

    # Upgrade: reuse gabled planes as main slopes, infer hip ends
    p1_copy = RoofPlane(
        normal=p1.normal.copy(), offset=p1.offset,
        inliers=p1.inliers.copy(), boundary_3d=p1.boundary_3d.copy(),
        role="main_left",
    )
    p2_copy = RoofPlane(
        normal=p2.normal.copy(), offset=p2.offset,
        inliers=p2.inliers.copy(), boundary_3d=p2.boundary_3d.copy(),
        role="main_right",
    )

    # Create synthetic hip end planes from the end-zone points
    hip_planes = []
    for end_mask, label in [(end_mask_1, "hip_end_1"), (end_mask_2, "hip_end_2")]:
        end_pts = all_points[end_mask]
        if len(end_pts) < 5:
            continue
        # Fit a plane to the end points
        centroid = end_pts.mean(axis=0)
        _, _, vh = np.linalg.svd(end_pts - centroid)
        normal = vh[-1]
        if normal[2] < 0:
            normal *= -1
        normal = normal / np.linalg.norm(normal)
        offset = np.dot(normal, centroid)

        hip_plane = RoofPlane(
            normal=normal, offset=offset,
            inliers=np.where(end_mask)[0],
            boundary_3d=end_pts[ConvexHull(end_pts[:, :2]).vertices]
                if len(end_pts) >= 3 else end_pts,
            role=label,
        )
        hip_planes.append(hip_plane)

    if len(hip_planes) < 2:
        return None

    all_planes = [p1_copy, p2_copy] + hip_planes

    # Confidence: slightly lower than gabled since this is inferred
    confidence = gabled_result.confidence * 0.95 * (1.0 - ridge_ratio)
    # But boost above gabled if ridge_ratio is clearly hipped
    if ridge_ratio < 0.7:
        confidence = max(confidence, gabled_result.confidence * 1.02)

    # Clip ridge to footprint
    ridge_clipped = _clip_ridge_to_footprint(ridge_line, fp_2d)
    if ridge_clipped is not None:
        ridge_line = ridge_clipped

    return ClassificationResult(
        roof_type=RoofType.HIPPED,
        confidence=min(confidence, 1.0),
        roof_planes=all_planes,
        ridge_line=ridge_line,
        eave_height=gabled_result.eave_height,
        ridge_height=gabled_result.ridge_height,
    )


def _clip_ridge_to_footprint(
    ridge_line: np.ndarray,
    fp_2d: np.ndarray,
) -> Optional[np.ndarray]:
    """Clip a ridge line to the footprint boundary, returning shortened endpoints."""
    from shapely.geometry import LineString as SLineString, Polygon as SPolygon
    fp_poly = SPolygon(fp_2d)
    ridge_dir = ridge_line[1, :2] - ridge_line[0, :2]
    ridge_dir = ridge_dir / max(np.linalg.norm(ridge_dir), 1e-10)
    mid = (ridge_line[0, :2] + ridge_line[1, :2]) / 2
    extent = max(fp_poly.bounds[2] - fp_poly.bounds[0],
                 fp_poly.bounds[3] - fp_poly.bounds[1]) * 2
    long_line = SLineString([mid - ridge_dir * extent, mid + ridge_dir * extent])
    clipped = long_line.intersection(fp_poly)
    if clipped.is_empty or clipped.geom_type != "LineString":
        return None
    coords = list(clipped.coords)
    if len(coords) < 2:
        return None
    ridge_z = (ridge_line[0, 2] + ridge_line[1, 2]) / 2
    return np.array([
        [coords[0][0], coords[0][1], ridge_z],
        [coords[-1][0], coords[-1][1], ridge_z],
    ])


def _intersect_planes(p1: RoofPlane, p2: RoofPlane) -> Optional[np.ndarray]:
    """Compute the intersection line of two planes as two 3D points.

    Uses the centroid of both plane boundaries as an anchor point so the
    resulting line segment is positioned near the actual roof geometry.
    """
    direction = np.cross(p1.normal, p2.normal)
    if np.linalg.norm(direction) < 1e-10:
        return None
    direction = direction / np.linalg.norm(direction)

    # Find a point on the intersection line near the actual roof
    # Use the midpoint of both plane centroids as anchor
    anchor = (p1.centroid + p2.centroid) / 2

    # Project anchor onto the intersection line by solving
    # the two plane equations with minimum distance from anchor
    A = np.array([p1.normal, p2.normal, direction])
    b = np.array([p1.offset, p2.offset, np.dot(direction, anchor)])
    try:
        point = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return None

    return np.array([point - direction * 20, point + direction * 20])
