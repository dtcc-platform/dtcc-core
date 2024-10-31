from ...model import Surface, Mesh

from ..register import register_model_method

from ..meshing.meshing import mesh_surface
import numpy as np
from ..model_conversion import create_builder_surface

from .. import _dtcc_builder

@register_model_method
def mesh(s: Surface, triangle_size=None) -> Mesh:
    """
    Mesh a `Surface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).

    Returns:
        Mesh: A `Mesh` object representing the meshed `Surface`.
    """

    return mesh_surface(s, triangle_size)


@register_model_method
def ray_intersection(
    s: Surface, origin: np.ndarray, direction: np.ndarray
) -> np.ndarray:
    """
    Compute the intersection points of a ray with a surface.

    Args:
        s (Surface): The surface.
        origin (np.ndarray): The origin of the ray.
        direction (np.ndarray): The direction of the ray.

    Returns:
        np.ndarray: The intersection points.
    """
    builder_surface = create_builder_surface(s)
    origin = np.array(origin, dtype=np.float64)
    direction = np.array(direction, dtype=np.float64)
    return _dtcc_builder.ray_surface_intersection(builder_surface, origin, direction)
