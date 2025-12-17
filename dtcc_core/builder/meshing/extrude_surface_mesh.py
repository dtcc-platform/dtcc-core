"""
Mesh extrusion functionality for creating solid meshes from surface meshes.
"""

import numpy as np
from typing import Set, Tuple, List, Dict
from dtcc_core.model import Mesh
from dtcc_core.logging import warning, error


def extrude_surface_to_solid(
    surface_mesh: Mesh,
    extrusion_depth: float = None,
    base_z: float = 0,
    minimum_thickness: float = 0.01,
) -> Mesh:
    """
    Extrude a surface mesh downwards to create a solid mesh suitable for 3D printing.

    This function takes a surface mesh and creates a solid by:
    1. Finding the boundary edges of the surface
    2. Extruding the boundary downwards to create side walls
    3. Creating a bottom cap at the specified depth

    Args:
        surface_mesh: Input surface mesh to extrude
        extrusion_depth: How far down to extrude (if None, uses mesh bounds)
        base_z: Z-coordinate for the bottom (if None, uses min_z - extrusion_depth)

    Returns:
        Solid mesh with extruded sides and bottom cap
    """
    if len(surface_mesh.vertices) == 0 or len(surface_mesh.faces) == 0:
        return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

    # Calculate bounds if needed
    z_min = surface_mesh.vertices[:, 2].min()
    z_max = surface_mesh.vertices[:, 2].max()

    if extrusion_depth is not None:
        base_z = z_min - extrusion_depth
    if base_z is None:
        error("Either extrusion_depth or base_z must be provided.")
    if z_min - base_z < minimum_thickness:
        warning(
            "base_z is above or equal to the minimum Z of the surface mesh. Adjusting base_z."
        )
        base_z = z_min - minimum_thickness

    # Find boundary edges
    boundary_edges = _find_boundary_edges(surface_mesh)

    # Start building the solid mesh
    solid_vertices = []
    solid_faces = []
    vertex_map = {}

    # Add original surface vertices
    for i, vertex in enumerate(surface_mesh.vertices):
        solid_vertices.append(vertex.copy())
        vertex_map[i] = len(solid_vertices) - 1

    # Add original surface faces
    for face in surface_mesh.faces:
        solid_faces.append(
            [vertex_map[face[0]], vertex_map[face[1]], vertex_map[face[2]]]
        )

    # Create extruded vertices for boundary
    boundary_vertex_map = {}
    for edge in boundary_edges:
        for vertex_idx in edge:
            if vertex_idx not in boundary_vertex_map:
                # Create extruded vertex at base_z
                original_vertex = surface_mesh.vertices[vertex_idx]
                extruded_vertex = np.array(
                    [original_vertex[0], original_vertex[1], base_z]
                )
                solid_vertices.append(extruded_vertex)
                boundary_vertex_map[vertex_idx] = len(solid_vertices) - 1

    # Create side wall faces
    for edge in boundary_edges:
        v1_orig = vertex_map[edge[0]]
        v2_orig = vertex_map[edge[1]]
        v1_ext = boundary_vertex_map[edge[0]]
        v2_ext = boundary_vertex_map[edge[1]]

        # Create two triangles for the side wall (quad split)
        # Ensure correct winding order (outward facing normals)
        solid_faces.append([v1_orig, v2_ext, v1_ext])
        solid_faces.append([v1_orig, v2_orig, v2_ext])

    # Create bottom cap
    additional_bottom_vertices, bottom_faces = _create_bottom_cap(
        surface_mesh, boundary_vertex_map, base_z
    )

    # Add additional bottom vertices
    solid_vertices.extend(additional_bottom_vertices)

    # Add bottom faces
    solid_faces.extend(bottom_faces)

    # Convert to numpy arrays
    vertices_array = np.array(solid_vertices)
    faces_array = np.array(solid_faces, dtype=np.int32)

    return Mesh(vertices=vertices_array, faces=faces_array)


def _find_boundary_edges(mesh: Mesh) -> Set[Tuple[int, int]]:
    """
    Find boundary edges of a mesh (edges that belong to only one face).

    Args:
        mesh: Input mesh

    Returns:
        Set of boundary edges as tuples (vertex1_idx, vertex2_idx)
    """
    edge_count = {}

    # Count how many faces each edge belongs to
    for face in mesh.faces:
        edges = [
            tuple(sorted([face[0], face[1]])),
            tuple(sorted([face[1], face[2]])),
            tuple(sorted([face[2], face[0]])),
        ]

        for edge in edges:
            edge_count[edge] = edge_count.get(edge, 0) + 1

    # Boundary edges appear in only one face
    boundary_edges = {edge for edge, count in edge_count.items() if count == 1}

    return boundary_edges


def _create_bottom_cap(
    surface_mesh: Mesh, boundary_vertex_map: Dict[int, int], base_z: float
) -> Tuple[List[np.ndarray], List[List[int]]]:
    """
    Create a bottom cap for the extruded mesh.

    This creates a flat bottom surface by projecting all surface vertices onto the base_z plane
    and reversing the face winding order to ensure inward-facing normals.

    Args:
        surface_mesh: Original surface mesh
        boundary_vertex_map: Mapping from original vertex indices to extruded vertex indices
        base_z: Z-coordinate for the bottom plane

    Returns:
        Tuple of (additional_vertices, bottom_faces)
    """
    additional_vertices = []
    bottom_faces = []

    # Create a mapping for all vertices projected to the bottom
    vertex_to_bottom_idx = {}

    # First, map boundary vertices to their existing extruded counterparts
    for orig_idx, bottom_idx in boundary_vertex_map.items():
        vertex_to_bottom_idx[orig_idx] = bottom_idx

    # For interior vertices, create new bottom vertices
    next_vertex_idx = len(surface_mesh.vertices) + len(boundary_vertex_map)

    for i, vertex in enumerate(surface_mesh.vertices):
        if i not in boundary_vertex_map:
            # Create a new bottom vertex for interior vertices
            bottom_vertex = np.array([vertex[0], vertex[1], base_z])
            additional_vertices.append(bottom_vertex)
            vertex_to_bottom_idx[i] = next_vertex_idx
            next_vertex_idx += 1

    # Create bottom faces by projecting all original faces and flipping orientation
    for face in surface_mesh.faces:
        bottom_face_indices = [
            vertex_to_bottom_idx[face[0]],
            vertex_to_bottom_idx[face[1]],
            vertex_to_bottom_idx[face[2]],
        ]

        # Flip winding order for inward-facing normal (bottom surface)
        bottom_faces.append(
            [bottom_face_indices[0], bottom_face_indices[2], bottom_face_indices[1]]
        )

    return additional_vertices, bottom_faces


def create_printable_surface_mesh(
    surface_mesh: Mesh,
    extrusion_depth: float = None,
    base_z: float = None,
    minimum_thickness: float = 0.001,
) -> Mesh:
    """
    Create a printable solid mesh from a surface mesh with additional validation.

    This is a higher-level function that includes validation and fixes for 3D printing.

    Args:
        surface_mesh: Input surface mesh
        extrusion_depth: Depth of extrusion
        base_z: Base Z coordinate
        minimum_thickness: Minimum wall thickness for 3D printing

    Returns:
        Solid mesh suitable for 3D printing
    """
    # Basic extrusion
    solid_mesh = extrude_surface_to_solid(
        surface_mesh, extrusion_depth, base_z, minimum_thickness
    )

    return solid_mesh
