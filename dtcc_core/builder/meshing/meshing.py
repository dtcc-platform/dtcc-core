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
    ms: MultiSurface, triangle_size=None, weld=False, clean=False
) -> Mesh:
    """
    Mesh a `MultiSurface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).
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
    min_mesh_angle = 25
    if triangle_size is None or triangle_size < 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_multisurface(
        builder_ms, triangle_size, min_mesh_angle, weld
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
    builder_mesh = _dtcc_builder.mesh_surface(builder_surface, triangle_size, 25)
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh


def mesh_multisurfaces(
    multisurfaces: [MultiSurface],
    max_mesh_edge_size=-1,
    min_mesh_angle=25,
    weld=False,
    clean=False,
) -> [Mesh]:

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


def merge_meshes(meshes: [Mesh], weld=False) -> Mesh:
    builder_meshes = [mesh_to_builder_mesh(mesh) for mesh in meshes]
    merged_mesh = _dtcc_builder.merge_meshes(builder_meshes, weld)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


@register_model_method
def merge(mesh: Mesh, other: Mesh, weld=False) -> Mesh:
    builder_mesh = mesh_to_builder_mesh(mesh)
    builder_other = mesh_to_builder_mesh(other)
    merged_mesh = _dtcc_builder.merge_meshes([builder_mesh, builder_other], weld)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


def disjoint_meshes(mesh: Mesh) -> List[Mesh]:
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
