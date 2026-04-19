from __future__ import annotations

from typing import Optional

import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import split

from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.roof_config import RoofDetectionConfig
from dtcc_core.builder.geometry_builders.roof_classification import ClassificationResult


def build_flat_geometry(lod1: MultiSurface) -> MultiSurface:
    """Build LoD2 flat-roof geometry from LoD1. Reorders to ground-walls-roof."""
    return _relabel_lod1_as_lod2(lod1)


def build_fallback_geometry(lod1: MultiSurface) -> MultiSurface:
    """Build fallback LoD2 geometry (same as flat)."""
    return _relabel_lod1_as_lod2(lod1)


def build_gabled_geometry(
    lod1: MultiSurface,
    footprint: Surface,
    ridge_line: np.ndarray,
    eave_height: float,
    ridge_height: float,
) -> MultiSurface:
    """Build LoD2 gabled roof geometry from LoD1 + classification results.

    Constructs geometry directly from footprint corners and ridge points
    using shared vertices to ensure watertightness. No Shapely split needed.
    """
    surfaces = []
    semantics = []

    # Ground surface (lowest avg Z from LoD1)
    avg_z = [s.vertices[:, 2].mean() for s in lod1.surfaces]
    bottom_idx = int(np.argmin(avg_z))
    surfaces.append(lod1.surfaces[bottom_idx])
    semantics.append(SurfaceSemantic.GROUND)

    fp_verts_2d = footprint.vertices[:, :2]
    ground_z = lod1.surfaces[bottom_idx].vertices[:, 2].min()
    n_fp = len(fp_verts_2d)

    # Ridge direction and perpendicular for classifying corners
    ridge_dir = ridge_line[1, :2] - ridge_line[0, :2]
    ridge_dir_norm = ridge_dir / max(np.linalg.norm(ridge_dir), 1e-10)
    ridge_perp = np.array([-ridge_dir_norm[1], ridge_dir_norm[0]])
    ridge_mid_2d = (ridge_line[0, :2] + ridge_line[1, :2]) / 2

    # Classify each footprint corner as left or right of the ridge
    left_corners_2d = []   # corners on one side of ridge
    right_corners_2d = []  # corners on the other side
    for v in fp_verts_2d:
        side = np.dot(v - ridge_mid_2d, ridge_perp)
        if side >= 0:
            left_corners_2d.append(v)
        else:
            right_corners_2d.append(v)

    if len(left_corners_2d) < 2 or len(right_corners_2d) < 2:
        return build_flat_geometry(lod1)

    # Find where the ridge line intersects the footprint boundary
    # to get the two ridge endpoints on the footprint edge
    fp_poly = Polygon(fp_verts_2d)
    extent = max(fp_poly.bounds[2] - fp_poly.bounds[0],
                 fp_poly.bounds[3] - fp_poly.bounds[1]) * 2
    extended_line = LineString([
        ridge_mid_2d - ridge_dir_norm * extent,
        ridge_mid_2d + ridge_dir_norm * extent,
    ])
    intersection = extended_line.intersection(fp_poly.boundary)

    # Extract the two ridge-footprint intersection points
    ridge_pts_2d = []
    if intersection.geom_type == "MultiPoint":
        ridge_pts_2d = [np.array(p.coords[0]) for p in intersection.geoms]
    elif intersection.geom_type == "Point":
        ridge_pts_2d = [np.array(intersection.coords[0])]

    if len(ridge_pts_2d) < 2:
        return build_flat_geometry(lod1)

    # Use the two intersection points as ridge endpoints at ridge_height
    r1_2d = ridge_pts_2d[0]
    r2_2d = ridge_pts_2d[1]
    r1_3d = np.array([r1_2d[0], r1_2d[1], ridge_height])
    r2_3d = np.array([r2_2d[0], r2_2d[1], ridge_height])

    # Build shared 3D corner vertices at eave height
    left_eave = [np.array([v[0], v[1], eave_height]) for v in left_corners_2d]
    right_eave = [np.array([v[0], v[1], eave_height]) for v in right_corners_2d]

    # Walls: ground to eave for each footprint edge
    for i in range(n_fp):
        v1 = fp_verts_2d[i]
        v2 = fp_verts_2d[(i + 1) % n_fp]
        wall = Surface(vertices=np.array([
            [v1[0], v1[1], ground_z],
            [v2[0], v2[1], ground_z],
            [v2[0], v2[1], eave_height],
            [v1[0], v1[1], eave_height],
        ], dtype=float))
        surfaces.append(wall)
        semantics.append(SurfaceSemantic.WALL)

    # Gable triangles: on edges where ridge meets the footprint boundary
    # These are the two edges that contain the ridge intersection points
    for r_pt_2d, r_pt_3d in [(r1_2d, r1_3d), (r2_2d, r2_3d)]:
        # Find which footprint edge contains this ridge point
        for i in range(n_fp):
            v1 = fp_verts_2d[i]
            v2 = fp_verts_2d[(i + 1) % n_fp]
            edge_line = LineString([v1, v2])
            if edge_line.distance(Point(r_pt_2d)) < 0.1:
                gable = Surface(vertices=np.array([
                    [v1[0], v1[1], eave_height],
                    [v2[0], v2[1], eave_height],
                    r_pt_3d,
                ], dtype=float))
                surfaces.append(gable)
                semantics.append(SurfaceSemantic.WALL)
                break

    # Roof surfaces: two slopes, each using footprint corners + ridge endpoints
    # Left roof: left_eave corners + r2 + r1 (sharing ridge and eave vertices with walls)
    left_roof_verts = np.array(left_eave + [r2_3d, r1_3d], dtype=float)
    surfaces.append(Surface(vertices=left_roof_verts))
    semantics.append(SurfaceSemantic.ROOF)

    # Right roof: right_eave corners + r1 + r2
    right_roof_verts = np.array(right_eave + [r1_3d, r2_3d], dtype=float)
    surfaces.append(Surface(vertices=right_roof_verts))
    semantics.append(SurfaceSemantic.ROOF)

    assert len(surfaces) == len(semantics)
    ms = MultiSurface(surfaces=surfaces, semantics=semantics)
    return _snap_multisurface_vertices(ms)


def build_hipped_geometry(
    lod1: MultiSurface,
    footprint: Surface,
    ridge_line: np.ndarray,
    eave_height: float,
    ridge_height: float,
) -> MultiSurface:
    """Build LoD2 hipped roof geometry from LoD1 + classification results.

    Constructs geometry from shared footprint corner vertices and ridge points.
    A hipped roof has: 4 walls (ground to eave), 2 trapezoidal main slopes,
    2 triangular hip ends. All surfaces share vertices at their junctions.
    """
    surfaces = []
    semantics = []

    # Ground surface
    avg_z = [s.vertices[:, 2].mean() for s in lod1.surfaces]
    bottom_idx = int(np.argmin(avg_z))
    surfaces.append(lod1.surfaces[bottom_idx])
    semantics.append(SurfaceSemantic.GROUND)

    fp_verts_2d = footprint.vertices[:, :2]
    ground_z = lod1.surfaces[bottom_idx].vertices[:, 2].min()
    n_fp = len(fp_verts_2d)

    # Ridge direction and perpendicular
    ridge_dir = ridge_line[1, :2] - ridge_line[0, :2]
    ridge_dir_norm = ridge_dir / max(np.linalg.norm(ridge_dir), 1e-10)
    ridge_perp = np.array([-ridge_dir_norm[1], ridge_dir_norm[0]])
    ridge_mid_2d = (ridge_line[0, :2] + ridge_line[1, :2]) / 2

    # Clip the ridge line to the footprint interior
    # For hipped roofs, the ridge is SHORTER than the building -- find where
    # it should end by projecting footprint corners onto the ridge axis
    fp_proj = np.dot(fp_verts_2d - ridge_mid_2d, ridge_dir_norm)
    fp_min_proj = fp_proj.min()
    fp_max_proj = fp_proj.max()

    # Ridge endpoints are inset from the footprint ends
    # Use the classified ridge line direction but clip to ~70% of footprint extent
    # (the actual ridge extent comes from the detection, but we ensure it's inside)
    r1_2d = ridge_mid_2d + ridge_dir_norm * max(fp_min_proj * 0.6, fp_min_proj)
    r2_2d = ridge_mid_2d + ridge_dir_norm * min(fp_max_proj * 0.6, fp_max_proj)

    # If the detection gave us ridge points closer to center, use those
    det_r1_proj = np.dot(ridge_line[0, :2] - ridge_mid_2d, ridge_dir_norm)
    det_r2_proj = np.dot(ridge_line[1, :2] - ridge_mid_2d, ridge_dir_norm)
    if abs(det_r1_proj) < abs(fp_min_proj) and abs(det_r2_proj) < abs(fp_max_proj):
        r1_2d = ridge_line[0, :2]
        r2_2d = ridge_line[1, :2]

    # Shared 3D ridge endpoints at ridge height
    r1_3d = np.array([r1_2d[0], r1_2d[1], ridge_height])
    r2_3d = np.array([r2_2d[0], r2_2d[1], ridge_height])

    # Classify footprint corners for the 4 hipped roof faces.
    # For a rectangular footprint with ridge R1-R2:
    #   - Each corner belongs to the nearest ridge endpoint (R1 or R2)
    #   - Within each endpoint group, corners are split left/right of the ridge
    # This gives: left_r1, right_r1 (hip end 1), left_r2, right_r2 (hip end 2)
    # Main slopes: left_r1 + left_r2 + ridge (left trapezoid)
    #              right_r1 + right_r2 + ridge (right trapezoid)
    # Hip ends:    left_r1 + right_r1 + r1 (triangle at end 1)
    #              left_r2 + right_r2 + r2 (triangle at end 2)

    r1_proj = np.dot(r1_2d - ridge_mid_2d, ridge_dir_norm)
    r2_proj = np.dot(r2_2d - ridge_mid_2d, ridge_dir_norm)

    # Assign each corner to nearest ridge endpoint AND side
    left_r1 = []   # left of ridge, closer to r1
    right_r1 = []  # right of ridge, closer to r1
    left_r2 = []   # left of ridge, closer to r2
    right_r2 = []  # right of ridge, closer to r2

    for v in fp_verts_2d:
        side = np.dot(v - ridge_mid_2d, ridge_perp)
        proj = np.dot(v - ridge_mid_2d, ridge_dir_norm)
        corner_3d = np.array([v[0], v[1], eave_height])

        # Closer to r1 or r2?
        near_r1 = abs(proj - r1_proj) < abs(proj - r2_proj)

        if near_r1:
            if side >= 0:
                left_r1.append(corner_3d)
            else:
                right_r1.append(corner_3d)
        else:
            if side >= 0:
                left_r2.append(corner_3d)
            else:
                right_r2.append(corner_3d)

    # Walls: ground to eave for each footprint edge (shared eave vertices)
    for i in range(n_fp):
        v1 = fp_verts_2d[i]
        v2 = fp_verts_2d[(i + 1) % n_fp]
        wall = Surface(vertices=np.array([
            [v1[0], v1[1], ground_z],
            [v2[0], v2[1], ground_z],
            [v2[0], v2[1], eave_height],
            [v1[0], v1[1], eave_height],
        ], dtype=float))
        surfaces.append(wall)
        semantics.append(SurfaceSemantic.WALL)

    # Main slope left: left_r1 + left_r2 corners + both ridge points (trapezoid)
    if left_r1 and left_r2:
        main_left = Surface(vertices=np.array(
            left_r1 + left_r2 + [r2_3d, r1_3d], dtype=float
        ))
        surfaces.append(main_left)
        semantics.append(SurfaceSemantic.ROOF)

    # Main slope right: right_r1 + right_r2 corners + both ridge points (trapezoid)
    if right_r1 and right_r2:
        main_right = Surface(vertices=np.array(
            right_r2 + right_r1 + [r1_3d, r2_3d], dtype=float
        ))
        surfaces.append(main_right)
        semantics.append(SurfaceSemantic.ROOF)

    # Hip end 1: left_r1 + right_r1 corners + r1 ridge point (triangle)
    if left_r1 and right_r1:
        hip1 = Surface(vertices=np.array(
            [right_r1[0], left_r1[0], r1_3d], dtype=float
        ))
        surfaces.append(hip1)
        semantics.append(SurfaceSemantic.ROOF)

    # Hip end 2: left_r2 + right_r2 corners + r2 ridge point (triangle)
    if left_r2 and right_r2:
        hip2 = Surface(vertices=np.array(
            [left_r2[0], right_r2[0], r2_3d], dtype=float
        ))
        surfaces.append(hip2)
        semantics.append(SurfaceSemantic.ROOF)

    assert len(surfaces) == len(semantics)
    ms = MultiSurface(surfaces=surfaces, semantics=semantics)
    return _snap_multisurface_vertices(ms)


def _snap_multisurface_vertices(ms: MultiSurface, tolerance: float = 0.05) -> MultiSurface:
    """Snap nearby vertices across all surfaces to shared coordinates.

    Collects all vertices, clusters those within tolerance, replaces each
    cluster with its centroid. This closes small gaps at surface junctions
    (e.g., where walls meet roof edges) caused by independent construction.
    """
    # Collect all vertices with (surface_idx, vertex_idx) tracking
    all_verts = []
    index_map = []  # (surface_idx, vertex_idx)
    for si, surface in enumerate(ms.surfaces):
        for vi, v in enumerate(surface.vertices):
            all_verts.append(v)
            index_map.append((si, vi))

    if not all_verts:
        return ms

    all_verts = np.array(all_verts)
    n = len(all_verts)

    # Build KDTree and find clusters within tolerance
    from scipy.spatial import cKDTree
    tree = cKDTree(all_verts)
    pairs = tree.query_pairs(r=tolerance)

    # Union-find to group vertices into clusters
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in pairs:
        union(a, b)

    # Compute centroid for each cluster
    clusters = {}
    for i in range(n):
        root = find(i)
        if root not in clusters:
            clusters[root] = []
        clusters[root].append(i)

    snapped = all_verts.copy()
    for root, members in clusters.items():
        if len(members) > 1:
            centroid = all_verts[members].mean(axis=0)
            for m in members:
                snapped[m] = centroid

    # Write back to surfaces
    new_surfaces = [Surface(vertices=s.vertices.copy()) for s in ms.surfaces]
    for idx, (si, vi) in enumerate(index_map):
        new_surfaces[si].vertices[vi] = snapped[idx]

    return MultiSurface(surfaces=new_surfaces, semantics=list(ms.semantics) if ms.semantics else None)


def _relabel_lod1_as_lod2(lod1: MultiSurface) -> MultiSurface:
    """Reorder LoD1 surfaces to ground-walls-roof convention and assign semantics."""
    if len(lod1.surfaces) < 2:
        return lod1

    surfaces = list(lod1.surfaces)

    # Identify bottom (lowest avg Z) and top (highest avg Z)
    avg_z = [s.vertices[:, 2].mean() for s in surfaces]
    bottom_idx = int(np.argmin(avg_z))
    top_idx = int(np.argmax(avg_z))

    ordered_surfaces = []
    ordered_semantics = []

    # Ground first
    ordered_surfaces.append(surfaces[bottom_idx])
    ordered_semantics.append(SurfaceSemantic.GROUND)

    # Walls (everything that's not top or bottom)
    for i, s in enumerate(surfaces):
        if i != bottom_idx and i != top_idx:
            ordered_surfaces.append(s)
            ordered_semantics.append(SurfaceSemantic.WALL)

    # Roof last
    ordered_surfaces.append(surfaces[top_idx])
    ordered_semantics.append(SurfaceSemantic.ROOF)

    assert len(ordered_surfaces) == len(ordered_semantics)
    return MultiSurface(surfaces=ordered_surfaces, semantics=ordered_semantics)
