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
    Mesh a MultiSurface into a triangular Mesh.

    Parameters
    ----------
    ms : MultiSurface
        MultiSurface to mesh.
    triangle_size : float, optional
        Maximum triangle size; ``None`` leaves it unconstrained.
    weld : bool, optional
        Whether to weld mesh vertices.
    snap : float, optional
        Snap distance for mesh vertices.
    clean : bool, optional
        Whether to clean the MultiSurface before meshing. Warning: meshing an
        invalid MultiSurface with a max triangle size may crash or produce
        unexpected results.

    Returns
    -------
    Mesh
        Triangular mesh representation of the MultiSurface.
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
    Mesh a Surface into a triangular Mesh.

    Parameters
    ----------
    s : Surface
        Surface to mesh.
    triangle_size : float, optional
        Maximum triangle size; ``None`` leaves it unconstrained.
    clean : bool, optional
        Whether to clean the surface before meshing. Warning: meshing an
        unclean surface with a max triangle size may crash or produce
        unexpected results.

    Returns
    -------
    Mesh
        Triangular mesh representation of the Surface.
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

    Parameters
    ----------
    multisurfaces : list[MultiSurface]
        MultiSurface objects to mesh.
    max_mesh_edge_size : float, optional
        Maximum edge size for mesh triangles; ``-1`` means no limit.
    min_mesh_angle : float, optional
        Minimum triangle angle in degrees.
    weld : bool, optional
        Whether to weld vertices during meshing.
    clean : bool, optional
        Whether to clean MultiSurfaces before meshing.

    Returns
    -------
    list[Mesh]
        Meshes corresponding to each input MultiSurface.
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

    Parameters
    ----------
    meshes : list[Mesh]
        Meshes to merge.
    weld : bool, optional
        Whether to weld vertices during merge.
    snap : float, optional
        Snap distance for vertex snapping.

    Returns
    -------
    Mesh
        Merged mesh containing all input meshes.
    """
    builder_meshes = [mesh_to_builder_mesh(mesh) for mesh in meshes]
    merged_mesh = _dtcc_builder.merge_meshes(builder_meshes, weld, snap)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


def merge(mesh: Mesh, other: Mesh, weld=False, snap=0) -> Mesh:
    """
    Merge two meshes into a single mesh.

    Parameters
    ----------
    mesh : Mesh
        First mesh to merge.
    other : Mesh
        Second mesh to merge.
    weld : bool, optional
        Whether to weld vertices during merge.
    snap : float, optional
        Snap distance for vertex snapping.

    Returns
    -------
    Mesh
        Merged mesh containing both input meshes.
    """
    builder_mesh = mesh_to_builder_mesh(mesh)
    builder_other = mesh_to_builder_mesh(other)
    merged_mesh = _dtcc_builder.merge_meshes([builder_mesh, builder_other], weld, snap)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


def snap_vertices(mesh: Mesh, snap_distance: float) -> Mesh:
    """
    Snap mesh vertices within a specified distance.

    Parameters
    ----------
    mesh : Mesh
        Mesh to snap vertices for.
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

    Parameters
    ----------
    mesh : Mesh
        Mesh to separate into components.

    Returns
    -------
    list[Mesh]
        Meshes, each containing one connected component.
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
