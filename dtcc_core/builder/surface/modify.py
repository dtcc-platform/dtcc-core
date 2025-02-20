from ...model import Mesh, Surface, MultiSurface
import numpy as np
from ..logging import info, warning, error


def clean_surface(s: Surface, tol: float = 1e-6) -> Surface:
    pass


def merge_coplanar(ms: MultiSurface, angle: float = 0.1, tol=1e-6) -> MultiSurface:
    """
    Merge coplanar surfaces in a MultiSurface.

    Args:
        ms (MultiSurface): The input MultiSurface object.
        angle (float): The angle in radians to consider two surfaces coplanar (default 0.1).
        tol (float): The tolerance to consider two vertices equal (default 1e-6

    Returns:
        MultiSurface: A new MultiSurface object with coplanar surfaces merged.
    """
    tol_decimals = round(np.log10(1 / tol))
    if not ms.is_planar(tol=tol):
        error("MultiSurface has non-planar surfaces. Cannot merge coplanar surfaces.")
    vertex_map = {}
    normals = []
    for i, s in enumerate(ms.surfaces):
        normals.append(np.abs(s.calculate_normal()))
        for v in s.vertices:
            tv = tuple(np.round(v, tol_decimals))
            if tv not in vertex_map:
                vertex_map[tv] = len(vertex_map)
