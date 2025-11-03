from ...model import Mesh, Surface, MultiSurface
from ..register import register_model_method

from .. import _dtcc_builder

from ..model_conversion import (
    create_builder_multisurface,
    create_builder_surface,
    builder_mesh_to_mesh,
    mesh_to_builder_mesh,
)

from dtcc_core.builder.polygons.surface import clean_surface, clean_multisurface

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from typing import List, Tuple
from dtcc_core.builder.logging import warning, info


def mesh_multisurface(
    ms: MultiSurface, triangle_size=None, weld=False, snap=0, clean=False
) -> Mesh:
    """
    Mesh a `MultiSurface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).
        snap (float): The snap distance for the mesh vertices (default 0).
        clean (bool): Whether to clean and attempt to fix errors in the multisurface before meshing (default True).
                      Warning! meshing a multisurface with a max triangle size that is invalid may crash the program
                      or produce unexpected results.

    Returns:
        Mesh: A `Mesh` object representing the meshed `MultiSurface`.
    """
    if clean:
        ms = clean_multisurface(ms)
    if ms is None:
        warning("Failed to clean multisurface.")
        return Mesh()
    builder_ms = create_builder_multisurface(ms)
    min_mesh_angle = 20.7
    if triangle_size is None or triangle_size < 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_multisurface(
        builder_ms, triangle_size, min_mesh_angle, weld, snap
    )
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh


def mesh_surface(s: Surface, triangle_size=None, clean=False) -> Mesh:
    """
    Mesh a `Surface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        clean (bool): Whether to clean the surface before meshing (default True). Warning! meshing a surface with a max
        triangle size that is not clean may crash the program or produce unexpected results.

    Returns:
        Mesh: A `Mesh` object representing the meshed `Surface`.
    """
    if clean:
        s = clean_surface(s)
        if s is None:
            warning("Failed to clean surface.")
            return Mesh()
    builder_surface = create_builder_surface(s)
    if triangle_size is None or triangle_size < 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_surface(builder_surface, triangle_size, 20.7)
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh


def mesh_multisurfaces(
    multisurfaces: [MultiSurface],
    max_mesh_edge_size=-1,
    min_mesh_angle=20.7,
    weld=False,
    clean=False,
) -> [Mesh]:
    """
    Mesh multiple MultiSurface objects into a list of Mesh objects.
    
    This function processes multiple MultiSurface objects simultaneously,
    creating triangular meshes for each one with consistent parameters.
    
    Parameters
    ----------
    multisurfaces : List[MultiSurface]
        List of MultiSurface objects to mesh.
    max_mesh_edge_size : float, default -1
        Maximum edge size for mesh triangles. If -1, no size limit.
    min_mesh_angle : float, default 20.7
        Minimum angle in degrees for mesh triangles.
    weld : bool, default False
        Whether to weld vertices during meshing.
    clean : bool, default False
        Whether to clean MultiSurfaces before meshing.
        
    Returns
    -------
    List[Mesh]
        List of meshes corresponding to input MultiSurfaces.
    """

    if clean:
        multisurfaces = [clean_multisurface(ms) for ms in multisurfaces]
        multisurfaces = [ms for ms in multisurfaces if ms is not None]

    if len(multisurfaces) == 0:
        return []
    builder_multisurfaces = [create_builder_multisurface(ms) for ms in multisurfaces]
    # print(f"create builder multisurfaces took {time() - start_time} seconds")
    # start_time = time()
    meshes = _dtcc_builder.mesh_multisurfaces(
        builder_multisurfaces, max_mesh_edge_size, min_mesh_angle, weld
    )
    # print(f"mesh multisurfaces took {time() - start_time} seconds")
    # start_time = time()
    meshes = [builder_mesh_to_mesh(mesh) for mesh in meshes]
    # print(f"convert builder mesh to mesh took {time() - start_time} seconds")
    return meshes


def merge_meshes(meshes: [Mesh], weld=False, snap=0) -> Mesh:
    """
    Merge multiple meshes into a single mesh.
    
    This function combines multiple mesh objects into one unified mesh,
    with optional vertex welding and snapping operations.
    
    Parameters
    ----------
    meshes : List[Mesh]
        List of meshes to merge.
    weld : bool, default False
        Whether to weld vertices during merge.
    snap : float, default 0
        Snap distance for vertex snapping.
        
    Returns
    -------
    Mesh
        Single merged mesh containing all input meshes.
    """
    builder_meshes = [mesh_to_builder_mesh(mesh) for mesh in meshes]
    merged_mesh = _dtcc_builder.merge_meshes(builder_meshes, weld, snap)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


def merge(mesh: Mesh, other: Mesh, weld=False, snap=0) -> Mesh:
    """
    Merge two meshes into a single mesh.
    
    This function combines two mesh objects into one unified mesh,
    with optional vertex welding and snapping operations.
    
    Parameters
    ----------
    mesh : Mesh
        First mesh to merge.
    other : Mesh
        Second mesh to merge.
    weld : bool, default False
        Whether to weld vertices during merge.
    snap : float, default 0
        Snap distance for vertex snapping.
        
    Returns
    -------
    Mesh
        Single merged mesh containing both input meshes.
    """
    builder_mesh = mesh_to_builder_mesh(mesh)
    builder_other = mesh_to_builder_mesh(other)
    merged_mesh = _dtcc_builder.merge_meshes([builder_mesh, builder_other], weld, snap)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


def snap_vertices(mesh: Mesh, snap_distance: float) -> Mesh:
    """
    Snap mesh vertices within a specified distance.
    
    This function merges vertices that are within the snap distance of each other,
    helping to clean up mesh topology and reduce duplicate vertices.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh to snap vertices for.
    snap_distance : float
        Maximum distance between vertices to be snapped together.
        
    Returns
    -------
    Mesh
        Mesh with snapped vertices.
    """
    builder_mesh = mesh_to_builder_mesh(mesh)
    snapped_mesh = _dtcc_builder.snap_mesh_vertices(builder_mesh, snap_distance)
    snapped_mesh = builder_mesh_to_mesh(snapped_mesh)
    return snapped_mesh


def disjoint_meshes(mesh: Mesh) -> List[Mesh]:
    """
    Separate a mesh into disconnected components.
    
    This function analyzes mesh connectivity and splits it into separate
    mesh objects for each disconnected component using graph analysis.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh to separate into components.
        
    Returns
    -------
    List[Mesh]
        List of meshes, each containing one connected component.
    """
    num_vertices = len(mesh.vertices)
    edges = np.vstack(
        [
            mesh.faces[:, [0, 1]],  # First edge of each face
            mesh.faces[:, [1, 2]],  # Second edge of each face
            mesh.faces[:, [2, 0]],  # Third edge of each face
        ]
    )
    # Create sparse adjacency matrix
    adj_matrix = sparse.coo_matrix(
        (np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
        shape=(num_vertices, num_vertices),
    )

    # Make matrix symmetric (undirected graph)
    adj_matrix = adj_matrix + adj_matrix.T
    n_components, labels = connected_components(
        csgraph=adj_matrix, directed=False, return_labels=True
    )
    disjointed_meshes = []
    for component_id in range(n_components):
        # Get vertices in this component
        component_vertex_mask = labels == component_id
        component_vertex_indices = np.where(component_vertex_mask)[0]

        # Create vertex index mapping
        vertex_map = {
            old_idx: new_idx for new_idx, old_idx in enumerate(component_vertex_indices)
        }

        # Get faces that use these vertices
        face_vertex_mask = np.isin(mesh.faces, component_vertex_indices)
        valid_faces_mask = np.all(face_vertex_mask, axis=1)
        component_faces = mesh.faces[valid_faces_mask]

        # Vectorized vertex index remapping
        new_faces = np.vectorize(vertex_map.get)(component_faces)

        # Create new mesh component
        new_vertices = mesh.vertices[component_vertex_indices]

        disjointed_meshes.append(Mesh(vertices=new_vertices, faces=new_faces))

    return disjointed_meshes


def remove_degenerate_faces(mesh: Mesh, min_area: float = 1e-8) -> Mesh:
    """
    Remove triangles with near-zero area (degenerate faces).

    These are created during snapping operations when vertices merge together,
    collapsing triangles to a line or point.

    Parameters
    ----------
    mesh : Mesh
        The mesh to clean.
    min_area : float, default 1e-8
        Minimum area threshold. Triangles with area < min_area are removed.

    Returns
    -------
    Mesh
        Mesh with degenerate faces removed.
    """
    if len(mesh.faces) == 0:
        return mesh

    # Compute face areas using cross product: area = |AB × AC| / 2
    v0 = mesh.vertices[mesh.faces[:, 0]]  # All first vertices
    v1 = mesh.vertices[mesh.faces[:, 1]]  # All second vertices
    v2 = mesh.vertices[mesh.faces[:, 2]]  # All third vertices

    # Edge vectors
    edge1 = v1 - v0  # AB
    edge2 = v2 - v0  # AC

    # Cross product (area normal)
    cross = np.cross(edge1, edge2)

    # Area = |cross| / 2
    areas = np.linalg.norm(cross, axis=1) / 2.0

    # Keep only non-degenerate faces
    valid_mask = areas >= min_area
    valid_faces = mesh.faces[valid_mask]

    if len(valid_faces) == 0:
        info(f"All faces removed as degenerate (area < {min_area})")
        return Mesh(vertices=mesh.vertices, faces=np.array([], dtype=np.int32))

    removed_count = len(mesh.faces) - len(valid_faces)
    if removed_count > 0:
        info(f"Removed {removed_count} degenerate faces (area < {min_area})")

    return Mesh(vertices=mesh.vertices, faces=valid_faces)


def remove_duplicate_faces(mesh: Mesh) -> Mesh:
    """
    Remove exact duplicate faces (same 3 vertices).

    Parameters
    ----------
    mesh : Mesh
        The mesh to clean.

    Returns
    -------
    Mesh
        Mesh with duplicate faces removed.
    """
    if len(mesh.faces) == 0:
        return mesh

    # Create sorted face representation for comparison
    # Sort vertex indices in each face to normalize (0,1,2) == (1,0,2) etc.
    sorted_faces = np.sort(mesh.faces, axis=1)

    # Convert to tuple for hashing
    face_tuples = [tuple(f) for f in sorted_faces]

    # Track unique faces
    unique_faces = []
    seen = set()

    for i, face_tuple in enumerate(face_tuples):
        if face_tuple not in seen:
            seen.add(face_tuple)
            unique_faces.append(mesh.faces[i])  # Keep original (unsorted) face

    if len(unique_faces) == len(mesh.faces):
        return mesh

    unique_faces = np.array(unique_faces, dtype=np.int32)
    removed_count = len(mesh.faces) - len(unique_faces)
    info(f"Removed {removed_count} exact duplicate faces")

    return Mesh(vertices=mesh.vertices, faces=unique_faces)


def remove_internal_faces(mesh: Mesh, angle_threshold: float = 170.0) -> Mesh:
    """
    Remove internal faces (back-to-back triangles) from merged meshes.

    When meshes are merged, duplicate geometry from building interiors creates
    back-to-back faces (opposite normal directions). This function removes them.

    Algorithm:
    1. Compute face normal for each triangle
    2. Build adjacency: which faces share edges
    3. For each face, check adjacent faces for opposite normals
    4. Remove faces that are back-to-back with adjacent face

    Parameters
    ----------
    mesh : Mesh
        The mesh to clean.
    angle_threshold : float, default 170.0
        Angle threshold (degrees) for detecting opposite normals.
        Default 170° catches nearly-opposite normals (internal faces).
        Use 179° for very strict, 160° for loose.

    Returns
    -------
    Mesh
        Mesh with internal faces removed.
    """
    if len(mesh.faces) == 0:
        return mesh

    # Step 1: Compute face normals
    v0 = mesh.vertices[mesh.faces[:, 0]]
    v1 = mesh.vertices[mesh.faces[:, 1]]
    v2 = mesh.vertices[mesh.faces[:, 2]]

    edge1 = v1 - v0
    edge2 = v2 - v0

    normals = np.cross(edge1, edge2)

    # Normalize (handle near-zero length normals)
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    # Avoid division by zero
    norms[norms < 1e-10] = 1.0
    normals = normals / norms

    # Step 2: Build edge-to-faces mapping
    edge_to_faces = {}

    for face_id, face in enumerate(mesh.faces):
        # Three edges of the triangle (normalized: smaller vertex index first)
        edges = [
            tuple(sorted([face[0], face[1]])),
            tuple(sorted([face[1], face[2]])),
            tuple(sorted([face[2], face[0]])),
        ]

        for edge in edges:
            if edge not in edge_to_faces:
                edge_to_faces[edge] = []
            edge_to_faces[edge].append(face_id)

    # Step 3: Find internal faces (back-to-back pairs)
    internal_faces = set()

    for edge, face_ids in edge_to_faces.items():
        if len(face_ids) == 2:
            # Two faces share this edge
            face1_id, face2_id = face_ids[0], face_ids[1]

            # Compute angle between normals
            dot_product = np.dot(normals[face1_id], normals[face2_id])
            # Clamp to [-1, 1] to handle numerical errors
            dot_product = np.clip(dot_product, -1.0, 1.0)
            angle_rad = np.arccos(dot_product)
            angle_deg = np.degrees(angle_rad)

            # If normals point opposite directions, these are back-to-back (internal)
            if angle_deg > angle_threshold:
                # Mark both as internal (remove pairs)
                internal_faces.add(face1_id)
                internal_faces.add(face2_id)
        elif len(face_ids) > 2:
            # More than 2 faces share an edge (non-manifold or merged duplicates)
            # Mark all as suspicious
            for face_id in face_ids:
                internal_faces.add(face_id)

    # Step 4: Keep only non-internal faces
    external_mask = np.array([i not in internal_faces for i in range(len(mesh.faces))])
    external_faces = mesh.faces[external_mask]

    if len(external_faces) == len(mesh.faces):
        return mesh

    removed_count = len(mesh.faces) - len(external_faces)
    info(f"Removed {removed_count} internal faces (angle_threshold={angle_threshold}°)")

    if len(external_faces) == 0:
        warning("All faces removed as internal - mesh is now empty!")
        return Mesh(vertices=mesh.vertices, faces=np.array([], dtype=np.int32))

    return Mesh(vertices=mesh.vertices, faces=external_faces)
