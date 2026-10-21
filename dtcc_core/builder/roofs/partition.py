"""Clip fitted plane domains to target cells, preserving holes and all pieces."""

import numpy as np
from shapely import affinity
from shapely.geometry import Polygon

from .types import RoofPatch


def polygon_parts(geometry):
    """Extract nonempty polygon components recursively from a geometry."""
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        return [geometry]
    return [p for g in getattr(geometry, "geoms", []) for p in polygon_parts(g)]


def clip_halfplane(coords, a, b, c):
    """Clip a convex polygon against ax + by + c <= 0."""
    result = []
    for p, q in zip(coords, np.roll(coords, -1, axis=0)):
        fp, fq = a * p[0] + b * p[1] + c, a * q[0] + b * q[1] + c
        if fp <= 1e-10:
            result.append(p)
        if (fp < 0) != (fq < 0):
            result.append(p + fp / (fp - fq) * (q - p))
    return np.asarray(result)


def planar_patches(model, target, *, block_id=0):
    """Clip a roof model's lower-envelope plane domains to the target polygon."""
    local = affinity.translate(target, -model.origin[0], -model.origin[1])
    xmin, ymin, xmax, ymax = local.bounds
    box = np.array(
        [
            [xmin - 1, ymin - 1],
            [xmax + 1, ymin - 1],
            [xmax + 1, ymax + 1],
            [xmin - 1, ymax + 1],
        ]
    )
    output = []
    for i, plane in enumerate(model.planes):
        cell = box.copy()
        for j, other in enumerate(model.planes):
            if i == j:
                continue
            cell = clip_halfplane(cell, *(plane - other))
            if len(cell) < 3:
                break
        if len(cell) < 3:
            continue
        for part in polygon_parts(local.intersection(Polygon(cell))):
            world = affinity.translate(part, *model.origin)
            output.append(
                RoofPatch(
                    world,
                    model.origin.copy(),
                    plane.copy(),
                    model.source_key,
                    i,
                    block_id,
                )
            )
    return output
