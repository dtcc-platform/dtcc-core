from ...model import Mesh, Surface, MultiSurface
from ..register import register_model_method

from .. import _dtcc_builder

from ..model_conversion import (
    create_builder_multisurface,
    create_builder_surface,
    builder_mesh_to_mesh,
    _check_mesh_metadata,
)

from dtcc_core.builder.polygons.surface import clean_surface, clean_multisurface

import numpy as np
from typing import List, Tuple
from copy import deepcopy
from dtcc_core.builder.logging import warning, info

from .backends import resolve_2d_mesher
from .dtcc_mesher_backend import mesh_surface_with_dtcc_mesher
from .orientation import orient_faces_consistently


def _mesh_surface_with_builder(
    surface: Surface,
    triangle_size: float | None = None,
    min_mesh_angle: float = 20.7,
    mesher: str = "auto",
) -> Mesh:
    builder_surface = create_builder_surface(surface)
    if triangle_size is None or triangle_size < 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_surface(
        builder_surface,
        triangle_size,
        min_mesh_angle,
        mesher,
    )
    return builder_mesh_to_mesh(builder_mesh)


def mesh_multisurface(
    ms: MultiSurface,
    triangle_size=None,
    weld=False,
    snap=0,
    clean=False,
    mesher: str | None = None,
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
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        Select the 2D meshing backend. ``"auto"`` prefers ``dtcc_mesher``
        when it is installed, then ``triangle``.

    Returns
    -------
    Mesh
        Triangular mesh representation of the MultiSurface.

    Notes
    -----
    Semantic regions are copied and reindexed onto output triangles. This path
    preserves the enclosing transform and leaves the input unchanged. It requires
    nonzero-area rings and surfaces in the enclosing coordinate frame; cleaning,
    welding/snapping and attached field interpolation are not supported.
    """
    if ms is not None and any(surface.regions for surface in ms.surfaces):
        raise NotImplementedError(
            "Semantic regions must belong to the MultiSurface, not its individual surfaces"
        )
    if ms is not None and ms.regions:
        return _mesh_surface_collection(ms, triangle_size, weld, snap, clean, mesher)
    if clean:
        ms = clean_multisurface(ms)
    if ms is None:
        warning("Failed to clean multisurface.")
        return Mesh()
    active_mesher = resolve_2d_mesher(mesher)

    if active_mesher == "dtcc_mesher":
        meshes = [
            mesh_surface(surface, triangle_size=triangle_size, mesher=active_mesher)
            for surface in ms.surfaces
        ]
        meshes = [mesh for mesh in meshes if len(mesh.faces) > 0]
        if not meshes:
            return Mesh()
        # Surfaces are meshed one by one, so the assembled solid needs a
        # single consistent winding before it reaches an exporter.
        return orient_faces_consistently(merge_meshes(meshes, weld=weld, snap=snap))

    builder_ms = create_builder_multisurface(ms)
    min_mesh_angle = 20.7
    if triangle_size is None or triangle_size < 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_multisurface(
        builder_ms, triangle_size, min_mesh_angle, weld, snap, active_mesher
    )
    mesh = builder_mesh_to_mesh(builder_mesh)
    return orient_faces_consistently(mesh)


def _mesh_surface_collection(ms, triangle_size, weld, snap, clean, mesher):
    """Triangulate a Solid or semantic MultiSurface, transferring membership."""
    from ...model import exchange

    # A public numerical operation accepts mutable native input. Reuse canonical
    # admission for this supported subset rather than duplicating region rules.
    exchange.validate(ms)
    if any(surface.regions for surface in ms.surfaces):
        raise NotImplementedError(
            "Semantic regions must belong to the enclosing geometry, not its individual surfaces"
        )
    if clean or weld or snap != 0:
        raise NotImplementedError(
            "Meshing semantic regions with cleaning, welding or snapping needs a face mapping"
        )
    if ms.fields or ms.dataset_context is not None:
        raise NotImplementedError(
            "Meshing semantic regions cannot transfer fields or Dataset Context"
        )
    for surface in ms.surfaces:
        if surface.fields:
            raise NotImplementedError(
                "Meshing semantic regions cannot interpolate surface fields"
            )
        if not np.array_equal(
            surface.transform.affine, np.eye(4)
        ) or surface.transform.srs not in ("", ms.transform.srs):
            raise NotImplementedError(
                "Meshing semantic regions requires surfaces in their enclosing geometry frame"
            )
        for ring in (
            ([surface.vertices] + surface.holes) if surface.vertices.size else []
        ):
            local = ring - ring[0]
            area_vector = np.cross(local, np.roll(local, -1, axis=0)).sum(axis=0)
            if not np.any(area_vector):
                raise ValueError(
                    "Meshing semantic regions requires nonzero-area polygon rings"
                )
    active_mesher = resolve_2d_mesher(mesher)
    meshes, offsets = [], [0]
    for surface_index, surface in enumerate(ms.surfaces):
        # Backends may calculate/store normals. Keep the source model untouched.
        try:
            mesh = mesh_surface(
                surface.copy(), triangle_size=triangle_size, mesher=active_mesher
            )
        except RuntimeError as exc:
            raise RuntimeError(
                f"Meshing surface {surface_index} failed: {exc}"
            ) from exc
        if surface.vertices.size and not len(mesh.faces):
            raise ValueError(
                "Meshing produced no triangles for a nonempty semantic surface"
            )
        meshes.append(mesh)
        offsets.append(offsets[-1] + len(mesh.faces))
    # Unwelded merge concatenates triangles in input order. Regions are attached
    # afterwards, so no semantic facts cross the geometry-only C++ adapter.
    result = (
        merge_meshes([mesh for mesh in meshes if len(mesh.faces)])
        if offsets[-1]
        else Mesh()
    )
    result.transform = deepcopy(ms.transform)
    result.regions = deepcopy(ms.regions)
    for region in result.regions:
        spans = [
            np.arange(offsets[i], offsets[i + 1], dtype=np.int64)
            for i in region.indices
        ]
        region.indices = np.concatenate(spans) if spans else np.empty(0, dtype=np.int64)
    return result


def mesh_surface(
    s: Surface,
    triangle_size=None,
    clean=False,
    mesher: str | None = None,
) -> Mesh:
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
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        Select the 2D meshing backend. ``"auto"`` prefers ``dtcc_mesher``
        when it is installed, then ``triangle``.

    Returns
    -------
    Mesh
        Triangular mesh representation of the Surface.
    """
    if s.regions:
        raise NotImplementedError(
            "Semantic regions must belong to a MultiSurface for triangulation"
        )
    if clean:
        s = clean_surface(s)
        if s is None:
            warning("Failed to clean surface.")
            return Mesh()
    active_mesher = resolve_2d_mesher(mesher)

    if active_mesher == "dtcc_mesher":
        return mesh_surface_with_dtcc_mesher(
            s,
            triangle_size=triangle_size,
            min_mesh_angle=20.7,
        )

    return _mesh_surface_with_builder(s, triangle_size, 20.7, active_mesher)


def mesh_multisurfaces(
    multisurfaces: [MultiSurface],
    max_mesh_edge_size=-1,
    min_mesh_angle=20.7,
    weld=False,
    clean=False,
    mesher: str | None = None,
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
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        Select the 2D meshing backend. ``"auto"`` prefers ``dtcc_mesher``
        when it is installed, then ``triangle``.

    Returns
    -------
    list[Mesh]
        Meshes corresponding to each input MultiSurface.

    Notes
    -----
    Inputs with semantic regions use the same preservation rules as
    ``mesh_multisurface`` and currently require ``min_mesh_angle=20.7``.
    """

    if any(ms.regions or any(s.regions for s in ms.surfaces) for ms in multisurfaces):
        if min_mesh_angle != 20.7:
            raise NotImplementedError(
                "Batch meshing semantic regions currently requires the default minimum mesh angle"
            )
        return [
            mesh_multisurface(
                ms,
                triangle_size=max_mesh_edge_size,
                weld=weld,
                clean=clean,
                mesher=mesher,
            )
            for ms in multisurfaces
        ]
    if clean:
        multisurfaces = [clean_multisurface(ms) for ms in multisurfaces]
        multisurfaces = [ms for ms in multisurfaces if ms is not None]

    if len(multisurfaces) == 0:
        return []

    active_mesher = resolve_2d_mesher(mesher)
    if active_mesher == "dtcc_mesher":
        return [
            mesh_multisurface(
                ms,
                triangle_size=max_mesh_edge_size,
                weld=weld,
                clean=False,
                mesher=active_mesher,
            )
            for ms in multisurfaces
        ]

    builder_multisurfaces = [create_builder_multisurface(ms) for ms in multisurfaces]
    # print(f"create builder multisurfaces took {time() - start_time} seconds")
    # start_time = time()
    meshes = _dtcc_builder.mesh_multisurfaces(
        builder_multisurfaces,
        max_mesh_edge_size,
        min_mesh_angle,
        weld,
        active_mesher,
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

    Notes
    -----
    Inputs must share the exact same transform and SRS; neither reprojection nor
    coordinate-frame conversion is implicit. The output retains a copy of that
    frame. Face normals are retained, or recomputed when snapping changes faces.
    Fields, semantic regions, Dataset Context and schema declarations require
    explicit transfer and are rejected. Snap distance uses local units.
    """
    meshes = list(meshes)
    for mesh in meshes:
        _check_mesh_metadata(mesh, allow_transform=True)
    frame = meshes[0].transform if meshes else None
    if any(
        mesh.transform.srs != frame.srs
        or not np.array_equal(mesh.transform.affine, frame.affine)
        for mesh in meshes[1:]
    ):
        raise ValueError(
            "Merging meshes requires the same transform and coordinate system"
        )
    builder_meshes = [
        _dtcc_builder.create_mesh(mesh.vertices, mesh.faces, mesh.markers, mesh.normals)
        for mesh in meshes
    ]
    merged_mesh = _dtcc_builder.merge_meshes(builder_meshes, weld, snap)
    result = builder_mesh_to_mesh(merged_mesh)
    if frame is not None:
        result.transform = deepcopy(frame)
    if any(mesh.normals.size for mesh in meshes):
        if snap > 0:
            result.normals = _face_normals(result)
        else:
            # Welding only identifies equal coordinates and retains face order.
            result.normals = np.concatenate(
                [
                    mesh.normals if mesh.normals.size else _face_normals(mesh)
                    for mesh in meshes
                ]
            )
    return result


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

    Notes
    -----
    Inputs must share the exact same transform and SRS; neither reprojection nor
    coordinate-frame conversion is implicit. The output retains a copy of that
    frame. Face normals are retained, or recomputed when snapping changes faces.
    Fields, semantic regions, Dataset Context and schema declarations require
    explicit transfer and are rejected. Snap distance uses local units.
    """
    return merge_meshes([mesh, other], weld=weld, snap=snap)


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

    Notes
    -----
    Preserves a copy of the transform/SRS; distance is measured in local units.
    Existing face normals are recomputed after snapping. Degenerate faces cannot
    receive a normal and raise ValueError. Fields, semantic regions, Dataset
    Context and schema declarations require explicit transfer and are rejected.
    """
    _check_mesh_metadata(mesh, allow_transform=True)
    builder_mesh = _dtcc_builder.create_mesh(
        mesh.vertices, mesh.faces, mesh.markers, mesh.normals
    )
    result = builder_mesh_to_mesh(
        _dtcc_builder.snap_mesh_vertices(builder_mesh, snap_distance)
    )
    result.transform = deepcopy(mesh.transform)
    if mesh.normals.size:
        result.normals = _face_normals(result)
    return result


def _face_normals(mesh: Mesh) -> np.ndarray:
    """Compute local face normals after geometry changes or for missing inputs."""
    if not len(mesh.faces):
        return np.empty((0, 3))
    triangles = np.asarray(mesh.vertices, dtype=np.float64)[mesh.faces]
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    lengths = np.linalg.norm(normals, axis=1)
    if not np.isfinite(lengths).all() or np.any(lengths == 0):
        raise ValueError(
            "Cannot compute normals for degenerate or nonfinite mesh faces"
        )
    return normals / lengths[:, None]


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
    # Imported here rather than at module scope to keep scipy off the
    # `import dtcc_core` path. See issue #87.
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components

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
