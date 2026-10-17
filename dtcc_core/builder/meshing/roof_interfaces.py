"""Coordinated XY subdivision of affine patches and their height crossings.

This module knows only polygons, planes and source ownership. All neighbours
enter the same segment graph before either triangulation backend is invoked.
"""
from dataclasses import replace

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import LineString, Point
from shapely.ops import polygonize, unary_union


XY_PRECISION = 1e-7


class RoofAssemblyError(RuntimeError):
    """A bounded failure with exact block attribution when available."""

    def __init__(self, message, *, block_ids=()):
        super().__init__(message)
        self.block_ids = tuple(sorted(set(block_ids)))


def _lines(geometry):
    if geometry.geom_type in ('LineString', 'LinearRing'):
        yield geometry
    else:
        for part in getattr(geometry, 'geoms', []):
            yield from _lines(part)


def node_interfaces(patches):
    """Retain creases, node T junctions, and split every affine height crossing.

    Coincident profiles have no isolated root. Numerical coincidence uses
    1e-7 m independently of the physical footprint/detail scale.
    """
    if not patches:
        return []
    # Normalize roundoff from independently clipped plane domains on a 100 nm
    # XY grid. This is numerical noding, unrelated to min_building_detail.
    # Never silently discard a collapsed patch.
    offset = np.floor(np.min([p.polygon.bounds[:2] for p in patches], axis=0))
    originals = [affinity.translate(p.polygon, *-offset) for p in patches]
    polygons = [shapely.set_precision(p, XY_PRECISION) for p in originals]
    for patch, before, after in zip(patches, originals, polygons):
        if after.is_empty or after.boundary.hausdorff_distance(before.boundary) > 2*XY_PRECISION:
            raise RoofAssemblyError('A mandatory boundary collapses at numerical noding precision', block_ids=[patch.block_id])
    original_tree = shapely.STRtree(originals)
    if any(p.is_empty or not p.is_valid for p in polygons):
        raise RoofAssemblyError('A roof patch collapses at numerical noding precision')
    tree = shapely.STRtree(polygons)
    boundaries = shapely.STRtree([p.boundary for p in polygons])
    graph = shapely.union_all([p.boundary for p in polygons], grid_size=XY_PRECISION)
    segments = []
    for line in _lines(graph):
        coords = np.asarray(line.coords)
        for a, b in zip(coords[:-1], coords[1:]):
            length = np.linalg.norm(b - a)
            if length <= XY_PRECISION/4:
                raise RoofAssemblyError('Zero-length mandatory roof edge')
            middle = Point((a + b) / 2)
            candidates = boundaries.query(middle.buffer(XY_PRECISION/10))
            adjacent = [i for i in candidates if polygons[i].boundary.distance(middle) <= XY_PRECISION/10]
            roots = [0., 1.]
            if len(adjacent) > 2:
                raise RoofAssemblyError('More than two patches incident on a boundary interval')
            if len(adjacent) == 2:
                left, right = (patches[i] for i in adjacent)
                differences = left.height(np.array([a, b])+offset) - right.height(np.array([a, b])+offset)
                d0, d1 = differences
                if abs(d0) > 1e-8 and abs(d1) > 1e-8 and d0 * d1 < 0:
                    root = float(d0 / (d0 - d1))
                    # A crossing within numerical coincidence of an existing
                    # endpoint uses that endpoint; repeated noding stays stable.
                    if length*min(root, 1-root) > 2*XY_PRECISION:
                        roots.append(root)
            roots.sort()
            vertices = [a, *[a + root * (b-a) for root in roots[1:-1]], b]
            for start, end in zip(vertices[:-1], vertices[1:]):
                segments.append(LineString([start, end]))
    # Polygonization transfers all noded vertices to both sides of each edge.
    owned = {}
    for polygon in polygonize(segments):
        point = polygon.representative_point()
        indices = tree.query(point, predicate='intersects')
        if len(indices) != 1:
            # Ownership comes from the pre-rounding coverage. Snap rounding
            # can put a tiny face between two independently clipped edges;
            # it is still roof where the original labelled patch covers it.
            candidates = original_tree.query(polygon, predicate='intersects')
            areas = [(originals[int(i)].intersection(polygon).area, int(i)) for i in candidates]
            areas.sort(reverse=True)
            if not areas or areas[0][0] <= 0:
                continue  # courtyard / outside the original coverage
            if polygon.area > 2*XY_PRECISION*polygon.length and areas[0][0] < .5*polygon.area:
                continue  # true courtyard; boundary intersections have zero area
            if polygon.area > 2*XY_PRECISION*polygon.length and len(indices) > 1:
                raise RoofAssemblyError('Overlapping roof patch ownership',
                                        block_ids=[patches[int(i)].block_id for i in indices])
            indices = [areas[0][1]]
        owned.setdefault(int(indices[0]), []).append(polygon)
    result = []
    from ..roofs.partition import polygon_parts
    for index, pieces in owned.items():
        # Remove numerical face boundaries within a single original patch.
        merged = shapely.union_all(pieces, grid_size=XY_PRECISION)
        for polygon in polygon_parts(merged):
            result.append(replace(patches[index], polygon=affinity.translate(
                shapely.geometry.polygon.orient(polygon, sign=1.0), *offset)))
    return result
