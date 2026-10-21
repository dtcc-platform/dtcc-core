"""Coordinated XY subdivision of affine patches and their height crossings.

This module knows only polygons, planes and source ownership. All neighbours
enter the same segment graph before either triangulation backend is invoked.
"""

from dataclasses import replace

import numpy as np
import shapely
from shapely import affinity
from shapely.ops import polygonize

from ..cleaning.subdivision import canonical_boundary_graph
from ..geometry_policy import XY_PRECISION


class RoofAssemblyError(RuntimeError):
    """A bounded failure with exact block attribution when available."""

    def __init__(self, message, *, block_ids=()):
        super().__init__(message)
        self.block_ids = tuple(sorted(set(block_ids)))


def _lines(geometry):
    if geometry.geom_type in ("LineString", "LinearRing"):
        yield geometry
    else:
        for part in getattr(geometry, "geoms", []):
            yield from _lines(part)


def node_interfaces(
    patches, *, height_crossings=True, domains=None, audit=None, canonical=False
):
    """Retain creases, node T junctions, and split every affine height crossing.

    Coincident profiles have no isolated root. Canonical source cells stay
    fixed; candidate overlay uses the numerical grid independently of detail.
    """
    if not patches:
        return []
    # Normalize roundoff from independently clipped plane domains on a 100 nm
    # XY grid. This is numerical noding, unrelated to min_building_detail.
    # Never silently discard a collapsed patch.
    offset = np.floor(np.min([p.polygon.bounds[:2] for p in patches], axis=0))
    originals = [affinity.translate(p.polygon, *-offset) for p in patches]
    polygons = (
        originals
        if canonical
        else [shapely.set_precision(p, XY_PRECISION) for p in originals]
    )
    for patch, before, after in zip(patches, originals, polygons):
        displacement = after.boundary.hausdorff_distance(before.boundary)
        if after.is_empty or displacement > 2 * XY_PRECISION:
            raise RoofAssemblyError(
                f"Roof subdivision source {patch.source_key}, block {patch.block_id}, "
                f"patch {patch.patch_id} collapses at numerical noding precision "
                f"(boundary displacement {displacement:g} m; area {before.area:g} -> {after.area:g} m2)",
                block_ids=[patch.block_id],
            )
    original_tree = shapely.STRtree(originals)
    if any(p.is_empty or not p.is_valid for p in polygons):
        raise RoofAssemblyError("A roof patch collapses at numerical noding precision")
    tree = shapely.STRtree(polygons)
    graph = (
        None
        if canonical
        else shapely.union_all([p.boundary for p in polygons], grid_size=XY_PRECISION)
    )
    local_domains = (
        {key: affinity.translate(domain, *-offset) for key, domain in domains.items()}
        if domains is not None
        else None
    )
    occupied = (
        shapely.union_all(list(local_domains.values()))
        if local_domains is not None
        else None
    )
    if local_domains is not None and not canonical:
        graph = canonical_boundary_graph(graph, local_domains, audit)
    # Establish labelled overlay cells before querying incidence. Distance from
    # a tiny edge midpoint to independently rounded polygons is not adjacency.
    cells = list(zip(polygons, range(len(polygons)))) if canonical else []
    for polygon in [] if canonical else polygonize(list(_lines(graph))):
        if (
            occupied is not None
            and polygon.intersection(occupied).area <= 0.5 * polygon.area
        ):
            continue
        indices = tree.query(polygon.representative_point(), predicate="intersects")
        if len(indices) != 1:
            candidates = original_tree.query(polygon, predicate="intersects")
            areas = sorted(
                (
                    (originals[int(i)].intersection(polygon).area, int(i))
                    for i in candidates
                ),
                reverse=True,
            )
            if not areas or areas[0][0] <= 0:
                if occupied is None:
                    continue
                nearby = original_tree.query(polygon.buffer(2 * XY_PRECISION))
                candidates = sorted(
                    (originals[int(i)].distance(polygon), int(i)) for i in nearby
                )
                if not candidates or candidates[0][0] > 2 * XY_PRECISION:
                    raise RoofAssemblyError(
                        "Canonical occupied face has no source roof owner"
                    )
                areas = [(polygon.area, candidates[0][1])]
            uncertainty = 2 * XY_PRECISION * polygon.length
            if polygon.area > uncertainty and areas[0][0] < 0.5 * polygon.area:
                continue
            if polygon.area > uncertainty and len(indices) > 1:
                raise RoofAssemblyError(
                    "Overlapping roof patch ownership",
                    block_ids=[patches[int(i)].block_id for i in indices],
                )
            indices = [areas[0][1]]
        cells.append((polygon, int(indices[0])))
    incident = {}
    for polygon, owner in cells:
        for ring in [polygon.exterior, *polygon.interiors]:
            xy = np.asarray(ring.coords)
            for a, b in zip(xy[:-1], xy[1:]):
                key = tuple(sorted((tuple(a), tuple(b))))
                if key[0] == key[1]:
                    raise RoofAssemblyError("Zero-length mandatory roof edge")
                incident.setdefault(key, set()).add(owner)
    crossings = {}
    for key, adjacent in incident.items():
        if len(adjacent) > 2:
            raise RoofAssemblyError(
                "Nonmanifold labelled roof subdivision",
                block_ids=[patches[i].block_id for i in adjacent],
            )
        if len(adjacent) != 2 or not height_crossings:
            continue
        a, b = np.asarray(key)
        left, right = (patches[i] for i in sorted(adjacent))
        d0, d1 = left.height(np.array([a, b]) + offset) - right.height(
            np.array([a, b]) + offset
        )
        if abs(d0) > 1e-8 and abs(d1) > 1e-8 and d0 * d1 < 0:
            root = float(d0 / (d0 - d1))
            if np.linalg.norm(b - a) * min(root, 1 - root) > 2 * XY_PRECISION:
                crossings[key] = a + root * (b - a)
    owned = {}
    from shapely.geometry import Polygon

    for polygon, owner in cells:
        rings = []
        for ring in [polygon.exterior, *polygon.interiors]:
            xy = np.asarray(ring.coords)
            noded = []
            for a, b in zip(xy[:-1], xy[1:]):
                noded.append(a)
                crossing = crossings.get(tuple(sorted((tuple(a), tuple(b)))))
                if crossing is not None:
                    noded.append(crossing)
            rings.append(noded)
        owned.setdefault(owner, []).append(Polygon(rings[0], rings[1:]))
    ownership_audit = []
    for index, polygon in enumerate(originals):
        final = shapely.union_all(owned.get(index, []))
        error = polygon.symmetric_difference(final).area
        budget = 3 * XY_PRECISION * polygon.length + np.pi * (3 * XY_PRECISION) ** 2
        if error > budget:
            raise RoofAssemblyError(
                f"Canonical roof ownership changed source {patches[index].source_key}, "
                f"patch {patches[index].patch_id} by {error:g} m2 (budget {budget:g})",
                block_ids=[patches[index].block_id],
            )
        ownership_audit.append(
            dict(
                source_key=patches[index].source_key,
                patch_id=patches[index].patch_id,
                before_area_m2=float(polygon.area),
                after_area_m2=float(final.area),
                ownership_error_m2=float(error),
                ownership_budget_m2=float(budget),
            )
        )
    if {patches[i].source_key for i in owned} != {p.source_key for p in patches}:
        raise RoofAssemblyError("Canonical roof subdivision would erase a source")
    if audit is not None:
        audit["ownership"] = ownership_audit
        audit["height_crossings"] = len(crossings)
    result = []
    for index, pieces in sorted(owned.items()):
        # Retain overlay cells and their exact noded rings. Coplanar face merging
        # belongs to the shared 3D stage; dissolving here can remove a neighbour's
        # mandatory breakpoint and force a second geometric adjacency query.
        for polygon in pieces:
            result.append(
                replace(
                    patches[index],
                    polygon=affinity.translate(
                        shapely.geometry.polygon.orient(polygon, sign=1.0), *offset
                    ),
                )
            )
    return result
