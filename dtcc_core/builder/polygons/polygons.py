"""Polygon merging helpers used by surface and multisurface processing."""

from typing import List

import shapely
import shapely.affinity
import shapely.ops
from shapely.geometry import (
    JOIN_STYLE,
    LineString,
    MultiPolygon,
    Polygon,
)
from shapely.validation import make_valid
import numpy as np

from ..logging import info, warning


def merge_polygons_convexhull(p1, p2):
    """
    Merge two polygons by taking their convex hull.

    This function creates a convex hull around two polygons as a fallback
    method when other merge approaches fail.

    Parameters
    ----------
    p1 : Polygon
        First polygon to merge.
    p2 : Polygon
        Second polygon to merge.

    Returns
    -------
    Polygon
        Convex hull of the two polygons.
    """
    mp = MultiPolygon([p1, p2])
    return mp.convex_hull


def merge_polygon_hulls(p1, p2):
    """
    Merge two polygons by unioning their convex hulls.

    This function attempts to merge polygons by first creating convex hulls
    and then performing a union operation.

    Parameters
    ----------
    p1 : Polygon
        First polygon to merge.
    p2 : Polygon
        Second polygon to merge.

    Returns
    -------
    Polygon or None
        Merged polygon if successful, None if merge fails.
    """
    p1_hull = p1.convex_hull
    p2_hull = p2.convex_hull
    mp = p1.union(p2)
    mp = make_valid(mp)
    if mp.geom_type == "Polygon":
        return mp
    else:
        return None


def merge_polygons_buffering(p1, p2, tol):
    """
    Merge two polygons using buffering technique.

    This function buffers both polygons by a tolerance distance, unions them,
    then shrinks the result back to attempt a clean merge.

    Parameters
    ----------
    p1 : Polygon
        First polygon to merge.
    p2 : Polygon
        Second polygon to merge.
    tol : float
        Tolerance distance for buffering operations.

    Returns
    -------
    Polygon or None
        Merged polygon if successful, None if merge fails.
    """
    b1 = p1.buffer(tol, 1, join_style=JOIN_STYLE.mitre)
    b2 = p2.buffer(tol, 1, join_style=JOIN_STYLE.mitre)
    m = b1.union(b2)
    m = m.buffer(-tol, 1, join_style=JOIN_STYLE.mitre)
    m = shapely.make_valid(m)
    if m.geom_type == "Polygon":
        return m
    else:
        return None


def merge_polygons(p1, p2, tol):
    """
    Merge two polygons using multiple fallback strategies.

    This function attempts to merge polygons using various approaches:
    union, buffering, hull merging, and convex hull as fallbacks.

    Parameters
    ----------
    p1 : Polygon
        First polygon to merge.
    p2 : Polygon
        Second polygon to merge.
    tol : float
        Tolerance for merge operations.

    Returns
    -------
    Polygon
        Merged polygon using the most appropriate method.
    """
    mp = p1.union(p2)
    # if mp.geom_type != "Polygon":
    #     info("Failed to merge polygons. Trying snapping")
    #     mp = merge_polygons_snapping(p1, p2, tol)
    if mp.geom_type != "Polygon":
        # info("Failed to merge polygons. Trying buffering")
        mp = merge_polygons_buffering(p1, p2, tol)
        if mp is None:
            mp = merge_polygon_hulls(p1, p2)
            if mp is None:
                info("Failed to merge polygons. Falling back to merging convex hull")
                mp = merge_polygons_convexhull(p1, p2)
    return mp


def remove_slivers(p: Polygon, tol):
    """
    Remove thin sliver areas from a polygon.

    This function uses buffering operations to remove narrow features
    that are considered slivers based on the tolerance.

    Parameters
    ----------
    p : Polygon
        Polygon to process.
    tol : float
        Tolerance for sliver removal.

    Returns
    -------
    Polygon
        Polygon with slivers removed.
    """
    # info("Removing slivers")
    b_p = p.buffer(tol, 1, join_style=JOIN_STYLE.mitre).buffer(
        -tol, 1, join_style=JOIN_STYLE.mitre
    )
    b_p = shapely.make_valid(b_p)
    if b_p.geom_type == "Polygon":
        return b_p
    else:
        b_p = merge_multipolygon(b_p, tol)
    return b_p


def merge_multipolygon(multipolygon, tol=0.1):
    """
    Merge a multipolygon into a single polygon.

    This function attempts to merge multiple polygon geometries into
    a single polygon by finding intersecting parts and merging them.

    Parameters
    ----------
    multipolygon : MultiPolygon
        MultiPolygon to merge.
    tol : float, default 0.1
        Tolerance for merge operations.

    Returns
    -------
    Polygon
        Single merged polygon.
    """
    polygons = list(multipolygon.geoms)
    merged = []
    # orig_polygons = polygons.copy()
    if len(polygons) == 0:
        return Polygon()
    if len(polygons) == 1:
        m = polygons[0]
    elif len(polygons) == 2:
        m = merge_polygons(polygons[0], polygons[1], tol)
    else:
        while len(polygons) > 0:
            p = polygons.pop()
            for i, p2 in enumerate(polygons):
                try:
                    p.intersects(p2.buffer(tol))
                except shapely.errors.GEOSException:
                    print(f"Failed to merge polygons {p} and {p2}")
                    print(f"p1: {p.is_valid}, p2: {p2.is_valid}")
                if p.intersects(p2.buffer(tol)):
                    p = merge_polygons(p, p2, tol)
                    polygons.pop(i)
                    polygons.append(p)
        merged.append(p)
        m = shapely.ops.unary_union(merged)
    m = make_valid(m)
    return m


def merge_list_of_polygons(mcp: List[Polygon], tolerance=1e-2) -> Polygon:
    """
    Merge a list of polygons into a single polygon.

    This function takes a list of polygons and attempts to merge them
    into a single polygon using union operations and fallback strategies.

    Parameters
    ----------
    mcp : List[Polygon]
        List of polygons to merge.
    tolerance : float, default 1e-2
        Tolerance for merge operations.

    Returns
    -------
    Polygon
        Single merged polygon.
    """
    if len(mcp) == 0:
        return []
    if len(mcp) == 1:
        return mcp[0]
    else:
        m = shapely.ops.unary_union(mcp)
        m = shapely.make_valid(m)
        if m.geom_type == "Polygon":
            if tolerance > 0:
                simp_m = m.simplify(tolerance / 4, preserve_topology=True)
                if simp_m.is_valid and simp_m.geom_type == "Polygon":
                    m = simp_m
            return m
        else:
            m = MultiPolygon([p for p in m.geoms if p.is_valid and (not p.is_empty)])

            m = merge_multipolygon(m, tolerance)
            if m.geom_type != "Polygon":
                warning("Failed to merge polygon list. Falling back to convex hull")
                m = m.convex_hull
            else:
                simp_m = m.simplify(tolerance / 4, preserve_topology=True)
                if simp_m.is_valid and simp_m.geom_type == "Polygon":
                    m = simp_m

            return m


def _split_line(line, max_length):
    # Calculate the length of the line
    length = line.length
    if length <= max_length:
        return [line]

    # Number of parts to split into
    num_parts = int(np.ceil(length / max_length))
    new_length = length / num_parts
    points = []

    # Generate points along the line at intervals of new_length
    for i in range(num_parts + 1):
        point = line.interpolate(new_length * i)
        points.append(point)

    # Create new lines from these points
    new_lines = [LineString([points[i], points[i + 1]]) for i in range(len(points) - 1)]
    return new_lines


def _split_boundary(boundary, max_length):
    new_points = []
    for i, (p1, p2) in enumerate(zip(boundary.coords[:-1], boundary.coords[1:])):
        line = LineString([p1, p2])
        split_lines = _split_line(line, max_length)
        # Add the start point of the first segment
        if i == 0:
            new_points.append(split_lines[0].coords[0])
        # Add the end points of each segment
        for split_line in split_lines:
            new_points.append(split_line.coords[1])
    return new_points


def split_polygon_sides(polygon, max_length):
    """
    Split long polygon sides into shorter segments.

    This function divides polygon edges that are longer than the maximum
    length into multiple shorter segments while preserving the overall shape.

    Parameters
    ----------
    polygon : Polygon
        Polygon to split sides for.
    max_length : float
        Maximum allowed edge length.

    Returns
    -------
    Polygon
        Polygon with split sides.
    """
    # Handle exterior
    new_exterior_points = _split_boundary(polygon.exterior, max_length)
    new_exterior = Polygon(new_exterior_points)

    # Handle interiors (holes)
    new_interiors = []
    for interior in polygon.interiors:
        new_interior_points = _split_boundary(interior, max_length)
        new_interiors.append(new_interior_points)

    # Create new polygon with holes
    new_polygon = Polygon(new_exterior_points, [interior for interior in new_interiors])
    return new_polygon
