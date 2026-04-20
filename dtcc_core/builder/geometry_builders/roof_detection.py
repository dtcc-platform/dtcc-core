from __future__ import annotations

import time
from typing import List, Optional, Tuple

import numpy as np
import open3d as o3d
from scipy.spatial import ConvexHull, cKDTree

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.builder.geometry_builders.roof_config import (
    RoofDetectionConfig,
    RoofPlane,
)


def filter_roof_points(
    pc: PointCloud,
    config: RoofDetectionConfig,
    stage_timings: Optional[dict] = None,
) -> Tuple[Optional[PointCloud], Optional[np.ndarray]]:
    """Filter a per-building point cloud to keep likely roof points.

    Estimates normals first (reusable by detect_roof_planes), then removes
    facade hits and optionally filters by classification.

    Returns (filtered_pointcloud, normals) or (None, None) if too few remain.

    When stage_timings is provided, records 'normals' and 'filter' keys in ms.
    """
    if pc is None or len(pc.points) < config.min_roof_points:
        return None, None

    # Step 1: Estimate normals using Open3D
    t0 = time.perf_counter()
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pc.points)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamKNN(knn=config.normal_estimation_k)
    )
    normals = np.asarray(pcd.normals)
    # Orient normals upward
    flip_mask = normals[:, 2] < 0
    normals[flip_mask] *= -1
    if stage_timings is not None:
        stage_timings["normals"] = (time.perf_counter() - t0) * 1000.0

    # Remaining steps: facade filter + classification filter + mask apply
    t_filter = time.perf_counter()
    keep_mask = np.ones(len(pc.points), dtype=bool)

    # Step 2: Remove facade points (near-vertical normals)
    angle_to_z = np.degrees(np.arccos(np.clip(normals[:, 2], -1.0, 1.0)))
    keep_mask &= angle_to_z < config.facade_angle_threshold

    # Step 3: Use classification if available
    if pc.classification is not None and len(pc.classification) == len(pc.points):
        has_building_class = np.any(pc.classification == 6)
        if has_building_class:
            keep_mask &= pc.classification == 6

    # Apply mask
    filtered_points = pc.points[keep_mask]
    filtered_normals = normals[keep_mask]
    filtered_cls = None
    if pc.classification is not None and len(pc.classification) == len(pc.points):
        filtered_cls = pc.classification[keep_mask]

    # Step 4: Check minimum count
    if len(filtered_points) < config.min_roof_points:
        if stage_timings is not None:
            stage_timings["filter"] = (time.perf_counter() - t_filter) * 1000.0
        return None, None

    filtered_pc = PointCloud(
        points=filtered_points,
        classification=filtered_cls if filtered_cls is not None else np.empty(0),
    )
    if stage_timings is not None:
        stage_timings["filter"] = (time.perf_counter() - t_filter) * 1000.0
    return filtered_pc, filtered_normals


def detect_roof_planes(
    pc: PointCloud,
    normals: np.ndarray,
    config: RoofDetectionConfig,
    stage_timings: Optional[dict] = None,
) -> List[RoofPlane]:
    """Detect roof planes via iterative RANSAC + region-growing refinement.

    When stage_timings is provided, records 'ransac' and 'region_growing'
    keys in ms, summed across all plane iterations.
    """
    if len(pc.points) < config.min_plane_points:
        if stage_timings is not None:
            stage_timings.setdefault("ransac", 0.0)
            stage_timings.setdefault("region_growing", 0.0)
        return []

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pc.points)
    pcd.normals = o3d.utility.Vector3dVector(normals)

    total_points = len(pc.points)
    remaining_indices = set(range(total_points))
    all_points = np.asarray(pcd.points)
    all_normals = normals.copy()
    planes: List[RoofPlane] = []

    ransac_ms = 0.0
    region_ms = 0.0

    for _ in range(config.max_roof_planes):
        if len(remaining_indices) < config.min_plane_points:
            break
        if len(remaining_indices) / total_points < config.min_remaining_ratio:
            break

        # Build sub-cloud from remaining points
        remaining_list = sorted(remaining_indices)
        sub_pcd = o3d.geometry.PointCloud()
        sub_pcd.points = o3d.utility.Vector3dVector(all_points[remaining_list])

        # RANSAC plane fit
        t0 = time.perf_counter()
        plane_model, inlier_sub_indices = sub_pcd.segment_plane(
            distance_threshold=config.ransac_distance_threshold,
            ransac_n=3,
            num_iterations=config.ransac_iterations,
        )
        ransac_ms += (time.perf_counter() - t0) * 1000.0

        if len(inlier_sub_indices) < config.min_plane_points:
            break

        # Map sub-indices back to global indices
        inlier_global = {remaining_list[i] for i in inlier_sub_indices}
        plane_normal = np.array(plane_model[:3])
        plane_offset = -plane_model[3]

        # Orient normal upward
        if plane_normal[2] < 0:
            plane_normal *= -1
            plane_offset *= -1

        # Region-growing refinement
        t1 = time.perf_counter()
        inlier_global = _region_grow(
            all_points, all_normals, plane_normal, plane_offset,
            inlier_global, remaining_indices, config,
        )
        region_ms += (time.perf_counter() - t1) * 1000.0

        if len(inlier_global) < config.min_plane_points:
            remaining_indices -= inlier_global
            continue

        # Build RoofPlane
        inlier_arr = np.array(sorted(inlier_global))
        inlier_pts = all_points[inlier_arr]

        # Convex hull boundary in 3D world coordinates
        try:
            hull = ConvexHull(inlier_pts[:, :2])
            boundary_3d = inlier_pts[hull.vertices]
        except Exception:
            boundary_3d = inlier_pts

        plane = RoofPlane(
            normal=plane_normal / np.linalg.norm(plane_normal),
            offset=plane_offset,
            inliers=inlier_arr,
            boundary_3d=boundary_3d,
        )
        planes.append(plane)
        remaining_indices -= inlier_global

    planes.sort(key=lambda p: len(p.inliers), reverse=True)
    if stage_timings is not None:
        stage_timings["ransac"] = ransac_ms
        stage_timings["region_growing"] = region_ms
    return planes


def _region_grow(
    points: np.ndarray,
    normals: np.ndarray,
    plane_normal: np.ndarray,
    plane_offset: float,
    seed_indices: set,
    candidate_indices: set,
    config: RoofDetectionConfig,
) -> set:
    """Grow a plane region from RANSAC seeds using normal consistency."""
    candidate_list = sorted(candidate_indices)
    if not candidate_list:
        return seed_indices

    tree = cKDTree(points[candidate_list])
    local_to_global = {i: candidate_list[i] for i in range(len(candidate_list))}
    global_to_local = {v: k for k, v in local_to_global.items()}

    region = set(seed_indices)
    frontier = list(seed_indices)
    threshold_cos = np.cos(np.radians(config.region_normal_threshold))
    dist_threshold = config.ransac_distance_threshold * 2

    while frontier:
        current = frontier.pop()
        if current not in global_to_local:
            continue
        local_idx = global_to_local[current]

        neighbor_locals = tree.query_ball_point(
            points[current], r=config.merge_distance, p=2
        )

        for nl in neighbor_locals:
            g = local_to_global[nl]
            if g in region:
                continue
            if g not in candidate_indices:
                continue

            # Check normal consistency
            cos_angle = np.dot(normals[g], plane_normal)
            if cos_angle < threshold_cos:
                continue

            # Check distance to plane
            dist = abs(np.dot(points[g], plane_normal) - plane_offset)
            if dist > dist_threshold:
                continue

            region.add(g)
            frontier.append(g)

    return region
