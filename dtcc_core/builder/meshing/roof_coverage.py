"""A shared numerical frame and precision model for the entire roof coverage."""

from dataclasses import replace

import numpy as np
import shapely
from scipy.spatial import cKDTree
from shapely import affinity
from shapely.geometry import Point, box

from ..geometry_policy import COINCIDENCE
from ..roofs.partition import polygon_parts
from .roof_interfaces import XY_PRECISION, RoofAssemblyError, node_interfaces


def local_roof_coverage(
    patches, bounds, *, height_crossings=True, domains=None, audit=None
):
    """Normalize roofs and terrain together, retaining every patch constraint.

    Independent floating-point differences of almost coincident roof edges
    can leave spurious nanometre-wide ground polygons. All operands use the
    same local fixed-precision model before union/difference. A collapsed roof
    patch is rejected, never silently omitted.
    """
    offset = np.floor(np.asarray(bounds[:2], dtype=float))
    local_bounds = (
        bounds[0] - offset[0],
        bounds[1] - offset[1],
        bounds[2] - offset[0],
        bounds[3] - offset[1],
    )
    local, world = [], []
    for patch in patches:
        original = affinity.translate(patch.polygon, *-offset)
        polygon = original
        if polygon.is_empty or polygon.geom_type != "Polygon":
            raise RoofAssemblyError(
                f"Patch {patch.patch_id} collapses or splits at numerical coverage precision",
                block_ids=[patch.block_id],
            )
        if polygon.boundary.hausdorff_distance(original.boundary) > 2 * XY_PRECISION:
            raise RoofAssemblyError(
                "Coverage normalization changed a mandatory roof boundary",
                block_ids=[patch.block_id],
            )
        local.append(replace(patch, polygon=polygon, origin=patch.origin - offset))
        world.append(replace(patch, polygon=affinity.translate(polygon, *offset)))
    # Node all translated boundaries jointly before constructing the terrain.
    # This also transfers roundoff faces using original patch ownership.
    local_domains = (
        {key: affinity.translate(domain, *-offset) for key, domain in domains.items()}
        if domains is not None
        else None
    )
    local = node_interfaces(
        local, height_crossings=height_crossings, domains=local_domains, audit=audit
    )
    world = [
        replace(
            p, polygon=affinity.translate(p.polygon, *offset), origin=p.origin + offset
        )
        for p in local
    ]
    # Both native backends need distinguishable constrained vertices. Pairs
    # closer than the shared 10 nm coincidence scale are diagnosed
    # with their incident owners, rather than collapsed by a backend.
    vertices = []
    for patch in local:
        for ring in [patch.polygon.exterior, *patch.polygon.interiors]:
            vertices.extend(list(ring.coords)[:-1])
    if vertices:
        points = np.unique(np.asarray(vertices), axis=0)
        pairs = cKDTree(points).query_pairs(COINCIDENCE)
        if pairs:
            indices = sorted(set(v for pair in pairs for v in pair))
            probes = [Point(points[i]) for i in indices]
            affected = {
                p.block_id
                for p in local
                if any(p.polygon.boundary.distance(q) <= COINCIDENCE for q in probes)
            }
            distance = min(
                float(np.linalg.norm(points[a] - points[b])) for a, b in pairs
            )
            raise RoofAssemblyError(
                f"Mandatory coverage vertices are only {distance:.6g} m apart; "
                "the numerical coincidence policy cannot retain both constraints",
                block_ids=affected,
            )
    domain = shapely.set_precision(box(*local_bounds), XY_PRECISION)
    union = shapely.union_all([p.polygon for p in local], grid_size=XY_PRECISION)
    ground = shapely.difference(domain, union, grid_size=XY_PRECISION)
    return local, world, polygon_parts(ground), local_bounds, offset
