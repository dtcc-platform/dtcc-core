from ...model import MultiSurface, Mesh
from ..register import register_model_method
from ..meshing.meshing import mesh_multisurface
from ..polygons.polygons import merge_list_of_polygons

from shapely.geometry import Polygon
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid

from ..model_conversion import create_builder_multisurface
import numpy as np

from .. import _dtcc_builder

@register_model_method
def mesh(ms: MultiSurface, triangle_size=None, weld=False) -> Mesh:
    """
    Mesh a `MultiSurface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).

    Returns:
        Mesh: A `Mesh` object representing the meshed `MultiSurface`.
    """

    return mesh_multisurface(ms, triangle_size, weld)


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
