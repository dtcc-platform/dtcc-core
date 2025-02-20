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
def mesh(ms: MultiSurface, triangle_size=None, weld=False, clean=False) -> Mesh:
    """
    Mesh a `MultiSurface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).

    Returns:
        Mesh: A `Mesh` object representing the meshed `MultiSurface`.
    """

    return mesh_multisurface(ms, triangle_size, weld, clean)


@register_model_method
def to_polygon(ms: MultiSurface, simplify=1e-2) -> Polygon:
    """Flatten a MultiSurface to a single 2D Polygon.

    Args:
        ms (MultiSurface): The MultiSurface to flatten.

    Returns:
        Polygon: The flattened Polygon.
    """
    polygons = [s.to_polygon() for s in ms.surfaces]
    polygons = [p for p in polygons if not p.is_empty and p.area > 1e-2]
    merged = merge_list_of_polygons(polygons)
    if simplify:
        merged = make_valid(merged)
        merged = merged.simplify(simplify, preserve_topology=True)
    return merged


@register_model_method
def ray_intersection(
    ms: MultiSurface, origin: np.ndarray, direction: np.ndarray
) -> np.ndarray:
    """
    Compute the intersection points of a ray with a MultiSurface.

    Args:
        ms (MultiSurface): The MultiSurface.
        origin (np.ndarray): The origin of the ray.
        direction (np.ndarray): The direction of the ray.

    Returns:
        np.ndarray: The intersection points.
    """
    builder_multisurface = create_builder_multisurface(ms)
    origin = np.array(origin, dtype=np.float64)
    direction = np.array(direction, dtype=np.float64)
    return _dtcc_builder.ray_multisurface_intersection(
        builder_multisurface, origin, direction
    )


def find_edge_connections(ms, tol=1e-6):
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
    visited = set()
    coplanar_groups = []

    edge_map, adjacent = find_edge_connections(ms, tol=tol)

    surfaces = ms.surfaces

    def dfs(surf_idx, component):
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
    Merge all coplanar surfaces in a MultiSurface.

    Args:
        ms (MultiSurface): The input MultiSurface object.
        tol (float): The tolerance to consider two vertices equal (default 1e-6)

    Returns:
        MultiSurface: A new MultiSurface object with coplanar surfaces merged.
    """
    coplanar_groups = group_coplanar_surfaces(ms, tol=tol)
    union_ms = MultiSurface()
    for cpg in coplanar_groups:
        union_surf = union_surfaces([ms.surfaces[i] for i in cpg])
        union_ms.surfaces.append(union_surf)
    return union_ms
