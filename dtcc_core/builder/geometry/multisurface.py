from ...model import MultiSurface, Mesh
from ..register import register_model_method
from ..meshing.meshing import mesh_multisurface
from ..polygons.polygons import merge_list_of_polygons
from ..polygons.surface import union_surfaces
from ..geometry.surface import are_coplanar

from collections import defaultdict

from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid

from ..model_conversion import create_builder_multisurface
import numpy as np

from ..logging import error, warning, info

from .. import _dtcc_builder


@register_model_method
def mesh(ms: MultiSurface, triangle_size=None, weld=False, snap=0, clean=False) -> Mesh:
    """
    Mesh a MultiSurface into a triangular Mesh.

    Parameters
    ----------
    ms : MultiSurface
        Surface collection to mesh.
    triangle_size : float, optional
        Maximum triangle size; ``None`` leaves it unconstrained.
    weld : bool, optional
        Whether to weld mesh vertices.
    snap : float, optional
        Snap distance for mesh vertices.
    clean : bool, optional
        Whether to clean the mesh after generation.

    Returns
    -------
    Mesh
        Triangular mesh representation of the MultiSurface.
    """

    return mesh_multisurface(ms, triangle_size, weld, snap, clean)


@register_model_method
def to_polygon(ms: MultiSurface, simplify=1e-2) -> Polygon:
    """
    Flatten a MultiSurface to a single 2D Polygon.

    Parameters
    ----------
    ms : MultiSurface
        The MultiSurface to flatten.
    simplify : float, optional
        Simplification tolerance passed to ``shapely.simplify``; default is 1e-2.

    Returns
    -------
    Polygon
        Flattened polygon.
    """
    polygons = [s.to_polygon() for s in ms.surfaces]
    polygons = [p for p in polygons if not p.is_empty and p.area > 1e-2]
    merged = merge_list_of_polygons(polygons)
    if simplify:
        merged = make_valid(merged)
        merged = merged.simplify(simplify, preserve_topology=True)
    return merged


def ray_intersection(
    ms: MultiSurface, origin: np.ndarray, direction: np.ndarray
) -> np.ndarray:
    """
    Intersect a ray with a MultiSurface.

    Parameters
    ----------
    ms : MultiSurface
        MultiSurface to intersect.
    origin : numpy.ndarray
        Origin of the ray as XYZ coordinates.
    direction : numpy.ndarray
        Direction vector of the ray.

    Returns
    -------
    numpy.ndarray
        Intersection points returned by the builder backend.
    """
    builder_multisurface = create_builder_multisurface(ms)
    origin = np.array(origin, dtype=np.float64)
    direction = np.array(direction, dtype=np.float64)
    return _dtcc_builder.ray_multisurface_intersection(
        builder_multisurface, origin, direction
    )


def find_edge_connections(ms, tol=1e-6):
    """
    Find shared edges and adjacency between surfaces in a MultiSurface.

    Builds a map from edges to the surfaces that contain them and derives
    adjacency relationships from shared edges.

    Parameters
    ----------
    ms : MultiSurface
        The MultiSurface object to analyze.
    tol : float, default 1e-6
        Tolerance for considering vertices equal when rounding coordinates.

    Returns
    -------
    tuple[dict, dict]
        A tuple containing (edge_map, adjacent) where:
        - edge_map: dict mapping edge tuples to lists of surface indices
        - adjacent: dict mapping surface indices to sets of adjacent surface indices
    """
    tol_decimals = round(np.log10(1 / tol))
    edge_map = defaultdict(list)
    for s_idx, s in enumerate(ms.surfaces):
        for i in range(len(s.vertices)):
            v0 = tuple(np.round(s.vertices[i], tol_decimals))
            v1 = tuple(np.round(s.vertices[(i + 1) % len(s.vertices)], tol_decimals))
            if v0 > v1:
                v0, v1 = v1, v0
            edge = (v0, v1)
            edge_map[edge].append(s_idx)
    adjacent = defaultdict(set)
    for edge, faces in edge_map.items():
        for i in range(len(faces)):
            for j in range(i + 1, len(faces)):
                adjacent[faces[i]].add(faces[j])
                adjacent[faces[j]].add(faces[i])
    return edge_map, adjacent


def group_coplanar_surfaces(ms, tol=1e-8):
    """
    Group coplanar surfaces that are connected by edges.

    Uses depth-first search to find connected components of coplanar surfaces
    within the MultiSurface, identifying surfaces that lie on the same plane
    and are connected through shared edges.

    Parameters
    ----------
    ms : MultiSurface
        The MultiSurface object containing surfaces to group.
    tol : float, default 1e-8
        Tolerance for determining if surfaces are coplanar.

    Returns
    -------
    list[set]
        List of sets, where each set contains indices of surfaces that are
        coplanar and connected.
    """
    visited = set()
    coplanar_groups = []

    edge_map, adjacent = find_edge_connections(ms, tol=tol)

    surfaces = ms.surfaces

    def dfs(surf_idx, component):
        """
        Depth-first search that collects connected coplanar surfaces.

        Parameters
        ----------
        surf_idx : int
            Index of the starting surface.
        component : set[int]
            Accumulator for surface indices that are coplanar and connected
            to ``surf_idx``.

        Returns
        -------
        None
            The function modifies ``component`` and ``visited`` in place.
        """
        visited.add(surf_idx)
        component.add(surf_idx)

        for neighbor in adjacent[surf_idx]:
            if neighbor not in visited and are_coplanar(
                surfaces[surf_idx], surfaces[neighbor], tol
            ):
                dfs(neighbor, component)

    for i in range(len(surfaces)):
        if i in visited:
            continue
        component = set()
        dfs(i, component)
        coplanar_groups.append(component)
    return coplanar_groups


def merge_coplanar(ms: MultiSurface, tol=1e-6) -> MultiSurface:
    """
    Merge coplanar surfaces in a MultiSurface.

    Parameters
    ----------
    ms : MultiSurface
        Input MultiSurface object.
    tol : float, optional
        Tolerance used to consider two vertices equal; default is 1e-6.

    Returns
    -------
    MultiSurface
        New MultiSurface with coplanar surfaces merged.
    """
    coplanar_groups = group_coplanar_surfaces(ms, tol=tol)
    union_ms = MultiSurface()
    for cpg in coplanar_groups:
        union_surf = union_surfaces([ms.surfaces[i] for i in cpg])
        union_ms.surfaces.append(union_surf)
    return union_ms
