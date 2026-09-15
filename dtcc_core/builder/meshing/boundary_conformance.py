# Copyright(C) 2026 DTCC
# Licensed under the MIT License

"""Boundary conformance for tiled surface meshes.

When 3D-printing tiled urban surface meshes and expanding the exhibit with
new neighboring tiles, terrain heights at shared boundaries may not match
already-printed tiles. This module implements constrained Laplacian smoothing
with pinned boundary vertices.

For axis-aligned rectangular tiles, the algorithm:
1. Auto-detects the shared boundary edge (x_min, x_max, y_min, or y_max)
2. Projects new boundary vertices onto the exact edge line
3. Samples z-values from old mesh at projected positions via barycentric interpolation
4. Pins those vertices and applies Laplacian smoothing to propagate adjustments

Building vertices are never modified during smoothing.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np

from dtcc_core.model.geometry.mesh import Mesh
from dtcc_core.builder.logging import info, warning

if TYPE_CHECKING:
    from scipy import sparse


def conform_boundary(
    new_mesh: Mesh,
    old_meshes: List[Mesh],
    smoothing_iterations: int = 10,
    tolerance: Optional[float] = None,
) -> Mesh:
    """Conform the boundary of *new_mesh* to match adjacent *old_meshes*.

    Terrain vertices along shared boundaries are pinned to the z-values
    sampled from the old meshes, then Laplacian smoothing propagates the
    adjustment into the interior while building vertices remain fixed.

    Parameters
    ----------
    new_mesh : Mesh
        The mesh whose boundary should be adjusted.
    old_meshes : list of Mesh
        Previously-printed meshes that the new mesh should conform to.
    smoothing_iterations : int, optional
        Number of Jacobi smoothing iterations (default 10).
    tolerance : float or None, optional
        Maximum XY distance for two boundary vertices to be considered in
        contact.  If *None*, computed as 0.5 * mean edge length of *new_mesh*.

    Returns
    -------
    Mesh
        A copy of *new_mesh* with adjusted z-values.

    Raises
    ------
    ValueError
        If *new_mesh* has no vertices or faces.
    """
    if len(new_mesh.vertices) == 0 or len(new_mesh.faces) == 0:
        raise ValueError("new_mesh must have vertices and faces")

    # Return copy if no old meshes
    if not old_meshes:
        info("No old meshes provided — returning unchanged copy")
        return new_mesh.copy(geometry_only=True)

    # Never mutate input
    mesh = new_mesh.copy(geometry_only=True)

    # Auto-compute tolerance from mean edge length
    if tolerance is None:
        tolerance = 0.5 * _compute_mean_edge_length(mesh.faces, mesh.vertices)

    # Classify vertices
    is_building = _classify_vertices(mesh.faces, mesh.markers, mesh.num_vertices)

    # Find boundary vertices
    boundary_mask = _find_boundary_vertices(mesh.faces, mesh.num_vertices)

    # Detect contact with old meshes
    pinned_mask, target_z = _detect_contact_vertices(
        mesh.vertices, boundary_mask, is_building, old_meshes, tolerance
    )

    num_pinned = int(pinned_mask.sum())
    num_candidates = int((boundary_mask & ~is_building).sum())

    if num_pinned == 0:
        warning(
            f"No contact vertices found (searched {num_candidates} boundary candidates "
            f"with tolerance {tolerance:.3f}m) — returning unchanged copy"
        )
        return mesh

    # Apply pinned z-values
    mesh.vertices[pinned_mask, 2] = target_z[pinned_mask]

    # Build adjacency and smooth
    adj = _build_adjacency(mesh.faces, mesh.num_vertices)
    freeze = is_building | pinned_mask
    mesh.vertices[:, 2] = _laplacian_smooth_z(
        mesh.vertices, adj, freeze, smoothing_iterations
    )

    info(
        f"Boundary conformance: {num_pinned}/{num_candidates} vertices pinned "
        f"({100.0*num_pinned/num_candidates:.1f}% of boundary), "
        f"tolerance={tolerance:.3f}m, {smoothing_iterations} iterations"
    )
    return mesh


def _classify_vertices(
    faces: np.ndarray,
    markers: np.ndarray,
    num_vertices: int,
) -> np.ndarray:
    """Classify mesh vertices as building or terrain.

    Parameters
    ----------
    faces : np.ndarray
        (M, 3) face connectivity.
    markers : np.ndarray
        (M,) per-face markers.  marker >= 0 means building.
    num_vertices : int
        Total number of vertices.

    Returns
    -------
    np.ndarray
        Boolean array of length *num_vertices*; True = building vertex.
    """
    is_building = np.zeros(num_vertices, dtype=bool)
    if len(markers) == 0 or len(markers) != len(faces):
        return is_building
    building_faces = markers >= 0
    if building_faces.any():
        building_vertex_ids = np.unique(faces[building_faces].ravel())
        is_building[building_vertex_ids] = True
    return is_building


def _find_boundary_vertices(
    faces: np.ndarray,
    num_vertices: int,
) -> np.ndarray:
    """Identify boundary vertices of an open mesh.

    A boundary vertex belongs to at least one edge that is shared by exactly
    one face.

    Parameters
    ----------
    faces : np.ndarray
        (M, 3) face connectivity.
    num_vertices : int
        Total number of vertices.

    Returns
    -------
    np.ndarray
        Boolean array of length *num_vertices*; True = boundary vertex.
    """
    is_boundary = np.zeros(num_vertices, dtype=bool)
    if len(faces) == 0:
        return is_boundary
    # Build all edges as sorted (min, max) pairs
    edges = np.vstack(
        [
            np.sort(faces[:, [0, 1]], axis=1),
            np.sort(faces[:, [1, 2]], axis=1),
            np.sort(faces[:, [2, 0]], axis=1),
        ]
    )
    edges = np.ascontiguousarray(edges)
    # Use structured array for unique counting
    edge_dtype = np.dtype([("a", edges.dtype), ("b", edges.dtype)])
    structured = edges.view(edge_dtype).ravel()
    unique_edges, inverse, counts = np.unique(
        structured, return_inverse=True, return_counts=True
    )
    boundary_unique = counts == 1
    boundary_edges_mask = boundary_unique[inverse]
    boundary_edge_verts = edges[boundary_edges_mask]
    if len(boundary_edge_verts) > 0:
        is_boundary[np.unique(boundary_edge_verts.ravel())] = True
    return is_boundary


def _find_shared_boundary_edge(
    new_vertices: np.ndarray,
    old_mesh: Mesh,
    tolerance: float,
) -> Tuple[Optional[str], Optional[float]]:
    """Find which bounding box edge of new mesh is shared with old mesh.

    Parameters
    ----------
    new_vertices : np.ndarray
        (N, 3) vertices of the new mesh.
    old_mesh : Mesh
        The old mesh to test against.
    tolerance : float
        Maximum distance for an edge to be considered shared.

    Returns
    -------
    edge_type : str or None
        One of 'x_min', 'x_max', 'y_min', 'y_max' if a shared edge is found.
    edge_value : float or None
        The coordinate value of the shared edge (e.g., x=4.0 or y=6398887).
    """
    if len(old_mesh.vertices) == 0:
        return None, None

    # Compute new mesh bounding box edges
    new_bbox = {
        'x_min': new_vertices[:, 0].min(),
        'x_max': new_vertices[:, 0].max(),
        'y_min': new_vertices[:, 1].min(),
        'y_max': new_vertices[:, 1].max(),
    }

    # Compute old mesh bounding box
    old_bbox = {
        'x_min': old_mesh.vertices[:, 0].min(),
        'x_max': old_mesh.vertices[:, 0].max(),
        'y_min': old_mesh.vertices[:, 1].min(),
        'y_max': old_mesh.vertices[:, 1].max(),
    }

    # Check each edge of new mesh against old mesh extent
    candidates = []

    # Check if new x_min edge is near old x_max edge (old mesh to the left)
    if abs(new_bbox['x_min'] - old_bbox['x_max']) < tolerance:
        # Verify y-ranges overlap
        y_overlap = (
            max(new_bbox['y_min'], old_bbox['y_min']) <=
            min(new_bbox['y_max'], old_bbox['y_max'])
        )
        if y_overlap:
            candidates.append(('x_min', new_bbox['x_min'],
                             abs(new_bbox['x_min'] - old_bbox['x_max'])))

    # Check if new x_max edge is near old x_min edge (old mesh to the right)
    if abs(new_bbox['x_max'] - old_bbox['x_min']) < tolerance:
        y_overlap = (
            max(new_bbox['y_min'], old_bbox['y_min']) <=
            min(new_bbox['y_max'], old_bbox['y_max'])
        )
        if y_overlap:
            candidates.append(('x_max', new_bbox['x_max'],
                             abs(new_bbox['x_max'] - old_bbox['x_min'])))

    # Check if new y_min edge is near old y_max edge (old mesh below)
    if abs(new_bbox['y_min'] - old_bbox['y_max']) < tolerance:
        x_overlap = (
            max(new_bbox['x_min'], old_bbox['x_min']) <=
            min(new_bbox['x_max'], old_bbox['x_max'])
        )
        if x_overlap:
            candidates.append(('y_min', new_bbox['y_min'],
                             abs(new_bbox['y_min'] - old_bbox['y_max'])))

    # Check if new y_max edge is near old y_min edge (old mesh above)
    if abs(new_bbox['y_max'] - old_bbox['y_min']) < tolerance:
        x_overlap = (
            max(new_bbox['x_min'], old_bbox['x_min']) <=
            min(new_bbox['x_max'], old_bbox['x_max'])
        )
        if x_overlap:
            candidates.append(('y_max', new_bbox['y_max'],
                             abs(new_bbox['y_max'] - old_bbox['y_min'])))

    if not candidates:
        return None, None

    # Return the closest edge match
    candidates.sort(key=lambda x: x[2])  # sort by distance
    return candidates[0][0], candidates[0][1]


def _detect_contact_vertices(
    new_vertices: np.ndarray,
    new_boundary_mask: np.ndarray,
    is_building: np.ndarray,
    old_meshes: List[Mesh],
    tolerance: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Detect which boundary vertices are in contact with old meshes.

    For axis-aligned rectangular tiles, identifies the shared boundary edge
    (e.g., x_min, x_max, y_min, or y_max) and projects new boundary vertices
    onto that edge line, then samples z-values from the old mesh.

    Parameters
    ----------
    new_vertices : np.ndarray
        (N, 3) vertex array of the new mesh.
    new_boundary_mask : np.ndarray
        Boolean array marking boundary vertices of the new mesh.
    is_building : np.ndarray
        Boolean array marking building vertices.
    old_meshes : list of Mesh
        Old meshes to test contact against.
    tolerance : float
        Maximum distance for edge detection (in XY).

    Returns
    -------
    pinned_mask : np.ndarray
        Boolean array of length N; True = vertex is pinned.
    target_z : np.ndarray
        Float array of length N; target z-value for pinned vertices.
    """
    n = len(new_vertices)
    pinned_mask = np.zeros(n, dtype=bool)
    target_z = np.zeros(n, dtype=np.float64)
    hit_count = np.zeros(n, dtype=np.int32)

    # Only terrain boundary vertices are candidates
    candidate_mask = new_boundary_mask & ~is_building
    candidate_indices = np.where(candidate_mask)[0]
    if len(candidate_indices) == 0:
        return pinned_mask, target_z

    for old_mesh in old_meshes:
        if len(old_mesh.vertices) == 0 or len(old_mesh.faces) == 0:
            continue

        # Find shared boundary edge
        edge_type, edge_value = _find_shared_boundary_edge(
            new_vertices, old_mesh, tolerance
        )
        if edge_type is None:
            continue

        # Identify vertices on the shared edge
        edge_tol = tolerance  # tolerance for "on the edge" check
        if edge_type == 'x_min':
            on_edge = np.abs(new_vertices[candidate_indices, 0] - edge_value) < edge_tol
            axis_idx = 0
        elif edge_type == 'x_max':
            on_edge = np.abs(new_vertices[candidate_indices, 0] - edge_value) < edge_tol
            axis_idx = 0
        elif edge_type == 'y_min':
            on_edge = np.abs(new_vertices[candidate_indices, 1] - edge_value) < edge_tol
            axis_idx = 1
        elif edge_type == 'y_max':
            on_edge = np.abs(new_vertices[candidate_indices, 1] - edge_value) < edge_tol
            axis_idx = 1
        else:
            continue

        edge_candidate_idx = candidate_indices[on_edge]
        if len(edge_candidate_idx) == 0:
            continue

        # Project vertices onto the exact edge line (snap to edge_value)
        query_xy = new_vertices[edge_candidate_idx, :2].copy()
        query_xy[:, axis_idx] = edge_value  # snap to exact edge coordinate

        # Sample z from old mesh at projected positions
        sampled_z, _ = _sample_old_mesh_heights(query_xy, old_mesh)

        # Accumulate (for averaging when near multiple old meshes)
        target_z[edge_candidate_idx] += sampled_z
        hit_count[edge_candidate_idx] += 1

    # Average and set pinned
    has_hits = hit_count > 0
    pinned_mask[has_hits] = True
    target_z[has_hits] /= hit_count[has_hits]

    return pinned_mask, target_z


def _sample_old_mesh_heights(
    query_points_xy: np.ndarray,
    old_mesh: Mesh,
) -> Tuple[np.ndarray, np.ndarray]:
    """Sample z-values from *old_mesh* at the given XY positions.

    Uses barycentric interpolation within the containing triangle.  Falls
    back to clamped barycentric coords on the nearest triangle, then to the
    nearest vertex z.

    Parameters
    ----------
    query_points_xy : np.ndarray
        (K, 2) array of XY query positions.
    old_mesh : Mesh
        The mesh to sample from.

    Returns
    -------
    z_values : np.ndarray
        (K,) sampled z-values.
    valid_mask : np.ndarray
        (K,) boolean; True if the point was inside a triangle (no fallback).
    """
    eps = 1e-8
    k = len(query_points_xy)
    z_values = np.zeros(k, dtype=np.float64)
    valid_mask = np.zeros(k, dtype=bool)

    if len(old_mesh.faces) == 0 or len(old_mesh.vertices) == 0:
        return z_values, valid_mask

    # Imported here rather than at module scope to keep scipy off the
    # `import dtcc_core` path. See issue #87.
    from scipy.spatial import cKDTree

    verts = old_mesh.vertices
    faces = old_mesh.faces

    # Build cKDTree on triangle centroids (XY)
    centroids_xy = (
        verts[faces[:, 0], :2] + verts[faces[:, 1], :2] + verts[faces[:, 2], :2]
    ) / 3.0
    tree = cKDTree(centroids_xy)

    # Also build vertex tree for last-resort fallback
    vertex_tree = cKDTree(verts[:, :2])

    # Search radius: use max edge length as heuristic
    v0 = verts[faces[:, 0], :2]
    v1 = verts[faces[:, 1], :2]
    v2 = verts[faces[:, 2], :2]
    max_edge = max(
        np.max(np.linalg.norm(v1 - v0, axis=1)),
        np.max(np.linalg.norm(v2 - v1, axis=1)),
        np.max(np.linalg.norm(v0 - v2, axis=1)),
    )
    search_radius = max_edge * 1.5

    for i in range(k):
        pt = query_points_xy[i]
        # Find candidate triangles within radius
        candidate_indices = tree.query_ball_point(pt, search_radius)

        found = False
        best_dist = np.inf
        best_z = 0.0

        for fi in candidate_indices:
            tri_v0 = verts[faces[fi, 0]]
            tri_v1 = verts[faces[fi, 1]]
            tri_v2 = verts[faces[fi, 2]]

            # Barycentric coordinates
            v01 = tri_v1[:2] - tri_v0[:2]
            v02 = tri_v2[:2] - tri_v0[:2]
            denom = v01[0] * v02[1] - v01[1] * v02[0]
            if abs(denom) < 1e-12:
                continue

            v0p = pt - tri_v0[:2]
            u = (v0p[0] * v02[1] - v0p[1] * v02[0]) / denom
            v = (v01[0] * v0p[1] - v01[1] * v0p[0]) / denom
            w = 1.0 - u - v

            if u >= -eps and v >= -eps and w >= -eps:
                # Point is inside this triangle
                z_values[i] = w * tri_v0[2] + u * tri_v1[2] + v * tri_v2[2]
                valid_mask[i] = True
                found = True
                break

            # Track nearest triangle for fallback (clamp barycentric)
            centroid_dist = np.linalg.norm(pt - centroids_xy[fi])
            if centroid_dist < best_dist:
                best_dist = centroid_dist
                uc = max(0.0, min(1.0, u))
                vc = max(0.0, min(1.0, v))
                wc = 1.0 - uc - vc
                if wc < 0:
                    total = uc + vc
                    uc /= total
                    vc /= total
                    wc = 0.0
                best_z = wc * tri_v0[2] + uc * tri_v1[2] + vc * tri_v2[2]

        if not found:
            if best_dist < np.inf:
                z_values[i] = best_z
            else:
                # Last resort: nearest vertex z
                _, nearest_vi = vertex_tree.query(pt)
                z_values[i] = verts[nearest_vi, 2]

    return z_values, valid_mask


def _build_adjacency(
    faces: np.ndarray,
    num_vertices: int,
) -> sparse.csr_matrix:
    """Build a sparse vertex adjacency matrix from faces.

    Parameters
    ----------
    faces : np.ndarray
        (M, 3) face connectivity.
    num_vertices : int
        Total number of vertices.

    Returns
    -------
    scipy.sparse.csr_matrix
        Binary symmetric adjacency matrix of shape (N, N).
    """
    # Imported here rather than at module scope to keep scipy off the
    # `import dtcc_core` path. See issue #87.
    from scipy import sparse

    edges = np.vstack(
        [
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        ]
    )
    data = np.ones(len(edges), dtype=np.float64)
    adj = sparse.coo_matrix(
        (data, (edges[:, 0], edges[:, 1])),
        shape=(num_vertices, num_vertices),
    )
    adj = adj + adj.T
    adj = adj.tocsr()
    # Binarize
    adj.data[:] = 1.0
    return adj


def _laplacian_smooth_z(
    vertices: np.ndarray,
    adjacency: sparse.csr_matrix,
    freeze_mask: np.ndarray,
    num_iterations: int,
) -> np.ndarray:
    """Jacobi Laplacian smoothing of z-values only.

    Parameters
    ----------
    vertices : np.ndarray
        (N, 3) vertex positions — only the z column is modified.
    adjacency : scipy.sparse.csr_matrix
        Binary symmetric adjacency matrix.
    freeze_mask : np.ndarray
        Boolean array; True = vertex z is held fixed.
    num_iterations : int
        Number of Jacobi iterations.

    Returns
    -------
    np.ndarray
        (N,) smoothed z-values.
    """
    z = vertices[:, 2].copy()
    if num_iterations == 0:
        return z

    degree = np.array(adjacency.sum(axis=1), dtype=np.float64).ravel()
    degree[degree == 0] = 1.0  # guard isolated vertices
    free = ~freeze_mask

    for _ in range(num_iterations):
        neighbor_sum = adjacency.dot(z)
        mean_z = neighbor_sum / degree
        z[free] = mean_z[free]

    return z


def _compute_mean_edge_length(faces: np.ndarray, vertices: np.ndarray) -> float:
    """Compute the mean edge length of a triangular mesh.

    Parameters
    ----------
    faces : np.ndarray
        (M, 3) face connectivity.
    vertices : np.ndarray
        (N, 3) vertex positions.

    Returns
    -------
    float
        Mean edge length.
    """
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]
    lengths = np.concatenate(
        [
            np.linalg.norm(v1 - v0, axis=1),
            np.linalg.norm(v2 - v1, axis=1),
            np.linalg.norm(v0 - v2, axis=1),
        ]
    )
    return float(np.mean(lengths))
