from copy import deepcopy

import numpy as np

from ...model.geometry import MultiSurface, Surface


def _ring_normal(ring: np.ndarray) -> np.ndarray:
    """Newell normal of a closed ring, whose sign follows the winding."""
    local = np.asarray(ring, dtype=np.float64) - ring[0]
    return np.cross(local, np.roll(local, -1, axis=0)).sum(axis=0)


def _wind_ring(ring: np.ndarray, upwards: bool) -> np.ndarray:
    """Return the ring wound so its normal points up (or down)."""
    ring = np.asarray(ring, dtype=np.float64)
    normal_z = _ring_normal(ring)[2]
    if (normal_z < 0.0) if upwards else (normal_z > 0.0):
        return ring[::-1].copy()
    return ring.copy()


def _surface_like(source: Surface, vertices: np.ndarray, holes=None) -> Surface:
    """Build a surface from explicit rings, keeping the source frame."""
    surface = Surface(
        vertices=np.asarray(vertices, dtype=np.float64), holes=list(holes or [])
    )
    surface.transform = deepcopy(source.transform)
    return surface


def extrude_surface(surface: Surface, height: float) -> MultiSurface:
    """
    Extrude a surface to a given height. The height is the absolute height, not the height relative to the surface.

    The extrusion is built as a closed solid whose faces are wound so that
    every normal points outwards: the roof ring runs counter-clockwise seen
    from above, the ground cap is the reverse of it, and each wall follows
    from the roof winding. Callers therefore get consistent normals no matter
    how the source footprint happened to be wound.

    Parameters
    ----------
    `surface` : model.Surface
        The surface to extrude.
    `height` : float
        The height to extrude the surface to.

    Returns
    -------
    `model.MultiSurface`
        The extruded surface.
    """
    roof_ring = _wind_ring(surface.vertices, upwards=True)
    # Holes are wound against the exterior so that a cap keeps a single sense.
    roof_holes = [_wind_ring(hole, upwards=False) for hole in surface.holes]

    ground_ring = roof_ring.copy()
    ground_ring[:, 2] = height
    ground_holes = []
    for hole in roof_holes:
        ground_hole = hole.copy()
        ground_hole[:, 2] = height
        ground_holes.append(ground_hole)

    extrusion = MultiSurface()
    # Roof keeps the upward winding; the ground cap is reversed to face down.
    extrusion.surfaces.append(_surface_like(surface, roof_ring, roof_holes))
    extrusion.surfaces.append(
        _surface_like(
            surface,
            ground_ring[::-1].copy(),
            [hole[::-1].copy() for hole in ground_holes],
        )
    )

    for top, bottom in [(roof_ring, ground_ring)] + list(zip(roof_holes, ground_holes)):
        for i in range(top.shape[0]):
            j = (i + 1) % top.shape[0]
            wall = np.array([top[j], top[i], bottom[i], bottom[j]])
            extrusion.surfaces.append(_surface_like(surface, wall))

    return extrusion
