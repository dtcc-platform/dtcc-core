"""A shared numerical frame and precision model for the entire roof coverage."""
from dataclasses import replace

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import box, Point
from scipy.spatial import cKDTree

from .roof_interfaces import XY_PRECISION, RoofAssemblyError, node_interfaces
from ..roofs.partition import polygon_parts


def local_roof_coverage(patches, bounds):
    """Normalize roofs and terrain together, retaining every patch constraint.

    Independent floating-point differences of almost coincident roof edges
    can leave spurious nanometre-wide ground polygons. All operands use the
    same local fixed-precision model before union/difference. A collapsed roof
    patch is rejected, never silently omitted.
    """
    offset = np.asarray(bounds[:2], dtype=float)
    local_bounds = (0., 0., bounds[2]-bounds[0], bounds[3]-bounds[1])
    local, world = [], []
    for patch in patches:
        original = affinity.translate(patch.polygon, *-offset)
        polygon = original
        if polygon.is_empty or polygon.geom_type != 'Polygon':
            raise RoofAssemblyError(f'Patch {patch.patch_id} collapses or splits at numerical coverage precision',
                                    block_ids=[patch.block_id])
        if polygon.boundary.hausdorff_distance(original.boundary) > 2*XY_PRECISION:
            raise RoofAssemblyError('Coverage normalization changed a mandatory roof boundary', block_ids=[patch.block_id])
        local.append(replace(patch, polygon=polygon, origin=patch.origin-offset))
        world.append(replace(patch, polygon=affinity.translate(polygon, *offset)))
    # Node all translated boundaries jointly before constructing the terrain.
    # This also transfers roundoff faces using original patch ownership.
    local = node_interfaces(local)
    world = [replace(p, polygon=affinity.translate(p.polygon, *offset), origin=p.origin+offset) for p in local]
    # Both native backends need distinguishable constrained vertices. Pairs
    # closer than the native 1 micrometre coincidence scale are diagnosed
    # with their incident owners, rather than collapsed by a backend.
    vertices = []
    for patch in local:
        for ring in [patch.polygon.exterior, *patch.polygon.interiors]:
            vertices.extend(list(ring.coords)[:-1])
    if vertices:
        points = np.unique(np.asarray(vertices), axis=0)
        pairs = cKDTree(points).query_pairs(1e-6)
        if pairs:
            indices = sorted(set(v for pair in pairs for v in pair))
            probes = [Point(points[i]) for i in indices]
            affected = {p.block_id for p in local if any(p.polygon.boundary.distance(q) <= 1e-6 for q in probes)}
            distance = min(float(np.linalg.norm(points[a]-points[b])) for a,b in pairs)
            raise RoofAssemblyError(
                f'Mandatory coverage vertices are only {distance:.6g} m apart; '
                'the native 1e-6 m coincidence scale cannot retain both constraints', block_ids=affected)
    domain = shapely.set_precision(box(*local_bounds), XY_PRECISION)
    union = shapely.union_all([p.polygon for p in local], grid_size=XY_PRECISION)
    ground = shapely.difference(domain, union, grid_size=XY_PRECISION)
    return local, world, polygon_parts(ground), local_bounds, offset
