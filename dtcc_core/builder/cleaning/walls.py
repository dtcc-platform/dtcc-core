"""Walls between the buildings of a merged block (stage 2 of cleaning).

Stage 1, the occupancy construction in
:mod:`dtcc_core.builder.cleaning.construction`, merges footprints closer than
delta into regions. This module partitions a region supported by several
sources into cells, so that every source keeps a footprint of its own: the
walls between cells follow the boundaries between the sources, and the region
boundary from stage 1 never moves.

Labels and the input partition
    Every source is a label. Where sources overlap, the source that comes first
    in a precedence order wins (the caller orders by height, then area, then
    input order). The input partition V_c is the part of source c that no
    earlier source covers.

Contract
    The cells Q^(c) of each label satisfy the identity fidelity condition

        V_c^{-epsilon} ⊆ Q^(c) ⊆ V_c^{+epsilon},                       (2L)

    with the fidelity budget of the region's merge group. Labels whose V_c has
    an empty core V_c^{-epsilon} may vanish. The cells must also pass the
    incident-sector profile, and the walls must not add a feature pair closer
    than delta to the essential boundary graph: every vertex-edge pair closer
    than delta must already be a pair of the region boundary alone. For a
    region without residual separation this is mfs >= delta for the labelled
    subdivision. :func:`dtcc_core.builder.cleaning.contract.check_labelled_contract`
    checks the result independently.

Construction
    1. Exact labelling: inside the sources, each point gets the label of its
       input partition; material added by stage 1 gets the nearest label
       (Voronoi regions of sampled label boundaries).
    2. The exact labelling is kept if no material was added and it already
       satisfies the contract. Otherwise the walls are simplified (coverage
       simplification with the region boundary fixed, largest tolerance
       first) to the first start that satisfies (2L).
    3. Bounded repair, site by site and with whole-region wall operations:
       collapse short wall edges, snap or slide wall junctions on the region
       boundary, project wall vertices onto the boundary, drop wall vertices,
       and give a cell without a protected core to its neighbour. An edit is
       admitted only if (2L) holds for every label and the defect measure
       decreases: the number of sectors below the profile's minimum angle at
       vertices the walls add or change, then the number of wall pairs closer
       than delta, then their total deficit.
    4. If defects remain, the repair is retried once at the retry budget.

Outcomes are ``strict`` (no defects), ``warning`` (only wall pairs closer than
delta remain) and ``fallback`` (no start satisfies (2L), or sectors below the
minimum angle remain; the region keeps one footprint). A region with a single
label is ``single``.

All work is done in a local frame on a fixed precision grid, which keeps
near-coincident linework from separate overlays noding consistently. Results
are translated back: points at the region's vertices get their exact
coordinates, and wall ends on its edges are projected onto them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import time
from typing import Any, Mapping, Sequence

import numpy as np
import shapely
from shapely import STRtree
from shapely.affinity import translate
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.ops import polygonize

from .construction import polygon_parts
from .contract import (
    DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES,
    SEPARATION_RELATIVE_TOLERANCE,
    FidelityBudget,
    essential_boundary_graph,
    interpret_input,
)

# Precision grid of the local frame (metres).
GRID = 1.0e-6
# Distance at which a point counts as lying on the region boundary.
ON_BOUNDARY = 5.0e-6
# Coverage-simplification tolerances for the start, as multiples of delta,
# largest first; the exact labelling (tolerance 0) comes last.
WALL_LADDER = (1.0, 0.5, 0.25, 0.1, 0.0)
MAX_ROUNDS = 6
# Wall pairs closer than this multiple of delta are repaired together.
SITE_RADIUS = 2.0
# Wall vertices within this multiple of delta of a site are edited.
SITE_REACH = 1.5
# Sectors within this many degrees above the minimum count as defects, as in
# the contract's profile check.
_SECTOR_TOLERANCE_DEGREES = 1e-6

Cells = list[tuple[Polygon, int]]


@dataclass
class LabelledRegion:
    """Cells of one region, in world coordinates.

    Attributes
    ----------
    region : Polygon
        The region, unchanged.
    cells : list[tuple[Polygon, int]]
        Cells and their labels, ordered by label precedence. Empty for a
        fallback; one cell (the region) for a single label. The cells tile
        the region; walls that end on its edges add vertices there.
    outcome : str
        ``single``, ``strict``, ``warning`` or ``fallback``.
    report : dict
        Fidelity budget used, wall pairs, edits and work.
    """

    region: Polygon
    cells: Cells
    outcome: str
    report: dict[str, Any] = field(default_factory=dict)


def label_region(
    region: Polygon,
    sources: Mapping[int, Any],
    order: Sequence[int],
    *,
    delta: float,
    epsilon: float,
    retry_epsilon: float | None = None,
) -> LabelledRegion:
    """Partition a cleaned region into cells, one set per source.

    Parameters
    ----------
    region : Polygon
        A region from stage 1.
    sources : mapping of int to Polygon or MultiPolygon
        The raw footprints of the sources that support the region, by label.
    order : sequence of int
        All labels in precedence order; earlier labels win overlaps.
    delta : float
        Minimum feature size.
    epsilon : float
        Fidelity budget of the region's merge group.
    retry_epsilon : float or None, optional
        Budget for one more repair when defects remain at ``epsilon``.

    Returns
    -------
    LabelledRegion
    """
    started = time.perf_counter()
    rank = {label: position for position, label in enumerate(order)}
    if len(sources) < 2:
        label = next(iter(sources))
        return LabelledRegion(region, [(region, label)], "single", {"seconds": 0.0})

    origin = region.bounds[:2]
    local_region = _to_local(region, origin)
    local_sources = {}
    for label in sorted(sources, key=rank.__getitem__):
        occupied, _ = interpret_input([sources[label]])
        local_sources[label] = _to_local(occupied, origin)
    partition = _input_partition(local_sources, rank)
    exact, gap_filled = _exact_cells(local_region, partition, epsilon, rank)
    if len({label for _, label in exact}) <= 1:
        label = exact[0][1]
        return LabelledRegion(
            region,
            [(region, label)],
            "single",
            {"seconds": time.perf_counter() - started},
        )

    budgets = [epsilon]
    if retry_epsilon is not None and retry_epsilon > epsilon:
        budgets.append(retry_epsilon)
    best = None
    attempts = []
    for budget in budgets:
        judge = _Judge(local_region, partition, delta, budget)
        cells, report = _repair(local_region, exact, judge, delta, try_exact=not gap_filled)
        report["fidelity_budget"] = budget
        attempts.append(
            {
                key: report.get(key)
                for key in ("fidelity_budget", "outcome", "small_sectors", "wall_pairs", "evaluations")
            }
        )
        if cells is None or report["outcome"] == "fallback":
            continue
        key = (report["wall_pairs"], report["wall_deficit"])
        if best is None or key < best[0]:
            best = (key, cells, report)
        if report["outcome"] == "strict":
            break

    if best is None:
        return LabelledRegion(
            region,
            [],
            "fallback",
            {
                "reason": (
                    "small_sectors"
                    if any(attempt.get("small_sectors") for attempt in attempts)
                    else "no_faithful_start"
                ),
                "attempts": attempts,
                "seconds": time.perf_counter() - started,
            },
        )
    _, cells, report = best
    world = _to_world(cells, origin, region)
    if world is None:
        return LabelledRegion(
            region,
            [],
            "fallback",
            {
                "reason": "invalid_after_translation",
                "attempts": attempts,
                "seconds": time.perf_counter() - started,
            },
        )
    world.sort(key=lambda item: (rank[item[1]], item[0].wkb))
    report["attempts"] = attempts
    report["seconds"] = time.perf_counter() - started
    report["cells"] = len(world)
    return LabelledRegion(region, world, report["outcome"], report)


# ---------------------------------------------------------------------------
# Local frame
# ---------------------------------------------------------------------------


def _to_local(geometry, origin):
    return shapely.set_precision(translate(geometry, -origin[0], -origin[1]), GRID)


def _to_world(cells: Cells, origin, region: Polygon) -> Cells | None:
    """Translate cells back onto the exact region boundary.

    Points within a few grid steps of a region vertex become that vertex, and
    other points that close to the region boundary (wall junctions) are
    projected onto it.
    """
    vertices = np.vstack(
        [np.asarray(ring.coords)[:-1, :2] for ring in [region.exterior, *region.interiors]]
    )
    tree = STRtree(shapely.points(vertices))
    boundary = region.boundary
    tolerance = 4 * GRID

    def restore(coordinates):
        coordinates = coordinates + np.asarray(origin)
        points = shapely.points(coordinates)
        found, nearest = tree.query_nearest(points, max_distance=tolerance, all_matches=False)
        coordinates[found] = vertices[nearest]
        rest = np.setdiff1d(np.arange(len(coordinates)), found)
        if len(rest):
            distances = shapely.distance(points[rest], boundary)
            close = rest[(distances > 0) & (distances <= tolerance)]
            if len(close):
                projected = shapely.line_interpolate_point(
                    boundary, shapely.line_locate_point(boundary, points[close])
                )
                coordinates[close] = shapely.get_coordinates(projected)
        return coordinates

    out = []
    for polygon, label in cells:
        world = shapely.transform(polygon, restore).normalize()
        if world.is_empty or not world.is_valid:
            return None
        out.append((world, label))
    return out


# ---------------------------------------------------------------------------
# Labelled partitions of one region: region boundary plus walls
# ---------------------------------------------------------------------------


def _lines(geometry):
    if geometry is None or geometry.is_empty:
        return []
    if geometry.geom_type == "LineString":
        return [geometry]
    if geometry.geom_type == "LinearRing":
        return [LineString(geometry.coords)]
    if geometry.geom_type in ("Polygon", "MultiPolygon"):
        return _lines(geometry.boundary)
    if hasattr(geometry, "geoms"):
        return [line for part in geometry.geoms for line in _lines(part)]
    return []


def _union(geometries):
    return shapely.unary_union(list(geometries), grid_size=GRID)


def _dissolve(polygons):
    """Union of pieces of a polygonal coverage (pieces with matching edges)."""
    polygons = list(polygons)
    if len(polygons) == 1:
        return polygons[0]
    union = shapely.coverage_union_all(polygons)
    if not union.is_valid:
        union = _union(polygons)
    return union


def _polygonal(geometry):
    return shapely.unary_union(polygon_parts(geometry))


def _input_partition(sources, rank):
    """V_c: the part of each source that no earlier source covers."""
    partition = {}
    taken = None
    for label in sorted(sources, key=rank.__getitem__):
        geometry = sources[label]
        part = geometry if taken is None else geometry.difference(taken, grid_size=GRID)
        part = _polygonal(part)
        if not part.is_empty and part.area > 0:
            partition[label] = part
        taken = geometry if taken is None else shapely.union(taken, geometry, grid_size=GRID)
    return partition


def _walls(cells: Cells):
    """Shared boundaries between cells of different labels, as merged lines.

    Cells come from one noded arrangement, so a wall is a ring segment that
    two cells of different labels share.
    """
    owners = {}
    for polygon, label in cells:
        for ring in [polygon.exterior, *polygon.interiors]:
            coordinates = [tuple(c[:2]) for c in ring.coords]
            for a, b in zip(coordinates, coordinates[1:]):
                if a != b:
                    owners.setdefault((a, b) if a < b else (b, a), set()).add(label)
    segments = [key for key, labels in owners.items() if len(labels) > 1]
    if not segments:
        return []
    return _lines(shapely.line_merge(shapely.multilinestrings(segments)))


def _boundary_vertices(region):
    rings = [region.exterior, *region.interiors]
    return np.vstack([np.asarray(ring.coords)[:-1, :2] for ring in rings])


def _insert_junctions(region, walls):
    """Region rings with wall ends inserted as vertices, and the walls ending there."""
    rings = [np.asarray(region.exterior.coords)[:-1, :2]] + [
        np.asarray(ring.coords)[:-1, :2] for ring in region.interiors
    ]
    ring_lines = [LineString(np.vstack([ring, ring[:1]])) for ring in rings]
    inserts = [[] for _ in rings]
    fixed = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        for end in (0, len(coordinates) - 1):
            point = Point(coordinates[end])
            distances = [line.distance(point) for line in ring_lines]
            k = int(np.argmin(distances))
            if distances[k] > ON_BOUNDARY:
                continue  # an interior junction
            ring = rings[k]
            to_vertex = np.linalg.norm(ring - np.asarray(coordinates[end]), axis=1)
            v = int(np.argmin(to_vertex))
            if to_vertex[v] <= ON_BOUNDARY:
                coordinates[end] = tuple(ring[v])
                continue
            s = ring_lines[k].project(point)
            q = ring_lines[k].interpolate(s)
            coordinates[end] = (q.x, q.y)
            inserts[k].append((s, (q.x, q.y)))
        if len(set(coordinates)) >= 2:
            fixed.append(LineString(coordinates))
    new_rings = []
    for ring, extra in zip(rings, inserts):
        closed = np.vstack([ring, ring[:1]])
        if not extra:
            new_rings.append(LineString(closed))
            continue
        cumulative = np.concatenate(
            [[0.0], np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
        )
        extra = sorted(set(extra))
        out = []
        e = 0
        for i in range(len(ring)):
            out.append(tuple(ring[i]))
            while e < len(extra) and extra[e][0] < cumulative[i + 1]:
                if extra[e][0] > cumulative[i]:
                    out.append(extra[e][1])
                e += 1
        out.append(out[0])
        new_rings.append(LineString(out))
    return new_rings, fixed


def _rebuild(region, walls, reference: Cells) -> Cells:
    """Polygonize the region boundary and walls; label faces by largest overlap."""
    rings, walls = _insert_junctions(region, walls)
    noded = _union([*rings, *walls])
    faces = [face for face in polygonize(noded) if region.contains(face.representative_point())]
    tree = STRtree([polygon for polygon, _ in reference])
    labelled = {}
    for face in faces:
        overlaps = {}
        for k in tree.query(face, predicate="intersects"):
            polygon, label = reference[int(k)]
            overlaps[label] = overlaps.get(label, 0.0) + face.intersection(polygon).area
        if not overlaps:
            continue
        best = max(overlaps.values())
        label = min(label for label, area in overlaps.items() if area == best)
        labelled.setdefault(label, []).append(face)
    return [
        (part, label)
        for label, polygons in labelled.items()
        for part in polygon_parts(_dissolve(polygons))
    ]


def _canonical(region, cells: Cells) -> Cells:
    """Dissolve each label, keep only walls between labels, and rebuild exactly."""
    by_label = {}
    for polygon, label in cells:
        by_label.setdefault(label, []).append(polygon)
    if len(by_label) <= 1:
        return [(region, next(iter(by_label)))]
    dissolved = [
        (part, label)
        for label, polygons in by_label.items()
        for part in polygon_parts(_dissolve(polygons))
    ]
    return _rebuild(region, _walls(dissolved), dissolved)


def _nearest_regions(extent, near, generators, epsilon):
    """Nearest-label regions around ``near``, from sampled label boundaries."""
    coordinates, labels = [], []
    for label, geometry in generators.items():
        lines = geometry.boundary.intersection(near)
        if lines.is_empty:
            continue
        sampled = shapely.get_coordinates(shapely.segmentize(lines, epsilon / 8.0))
        coordinates.append(sampled)
        labels.append(np.full(len(sampled), label))
    if not coordinates:
        return {}
    coordinates = np.round(np.vstack(coordinates), 9)
    labels = np.concatenate(labels)
    coordinates, first = np.unique(coordinates, axis=0, return_index=True)
    labels = labels[first]
    cells = shapely.voronoi_polygons(
        shapely.multipoints(coordinates),
        extend_to=extent.buffer(4.0 * epsilon).envelope,
        ordered=True,
    )
    polygons = np.asarray(cells.geoms)
    return {
        int(label): shapely.coverage_union_all(polygons[labels == label])
        for label in np.unique(labels)
    }


def _exact_cells(region, partition, epsilon, rank) -> tuple[Cells, bool]:
    """Label the region by the input partition; added material by nearest label.

    Also returns whether material outside the sources was labelled. Its walls
    follow sampled Voronoi boundaries, which are never admissible as they are.
    """
    relevant = {
        label: geometry
        for label, geometry in partition.items()
        if shapely.dwithin(geometry, region, 2 * epsilon)
    }
    if len(relevant) <= 1:
        label = next(iter(relevant), min(partition, key=rank.__getitem__))
        return [(region, label)], False
    occupied = _union(relevant.values())
    gap = region.difference(occupied, grid_size=GRID)
    linework = [region.boundary]
    for geometry in relevant.values():
        linework.append(geometry.boundary.intersection(region, grid_size=GRID))
    nearest = {}
    if gap.area > 1e-12:
        nearest = _nearest_regions(region, gap.buffer(2 * epsilon), relevant, epsilon)
        for geometry in nearest.values():
            linework.append(geometry.boundary.intersection(gap.buffer(4 * GRID)))
    noded = _union([geometry for geometry in linework if not geometry.is_empty])
    faces = shapely.get_parts(shapely.polygonize(shapely.get_parts(noded)))
    points = shapely.point_on_surface(faces)
    inside_region = shapely.contains(region, points)
    faces, points = faces[inside_region], points[inside_region]

    def first_covering(generators):
        """For each point, the first label in precedence whose geometry covers it."""
        labels = sorted(generators, key=rank.__getitem__)
        found = np.full(len(points), -1)
        if labels:
            hits, targets = STRtree([generators[label] for label in labels]).query(
                points, predicate="covered_by"
            )
            for hit, target in zip(hits, targets):
                if found[hit] < 0 or target < found[hit]:
                    found[hit] = target
        return [None if index < 0 else labels[index] for index in found]

    inside = first_covering(relevant)
    near = first_covering(nearest) if nearest else [None] * len(points)
    labelled = {}
    for face, point, a, b in zip(faces, points, inside, near):
        label = a if a is not None else b
        if label is None:
            label = min(relevant, key=lambda c: (relevant[c].distance(point), rank[c]))
        labelled.setdefault(label, []).append(face)
    cells = [
        (part, label)
        for label, polygons in labelled.items()
        for part in polygon_parts(_dissolve(polygons))
    ]
    return _canonical(region, cells), bool(nearest)


# ---------------------------------------------------------------------------
# Wall operations; each returns new cells or None when nothing changes
# ---------------------------------------------------------------------------


def _simplify_walls(region, cells: Cells, tolerance) -> Cells | None:
    if tolerance <= 0:
        return cells
    simplified = shapely.coverage_simplify(
        [polygon for polygon, _ in cells], tolerance, simplify_boundary=False
    )
    out = [(polygon, label) for polygon, (_, label) in zip(simplified, cells)]
    if any(not polygon.is_valid or polygon.is_empty for polygon, _ in out):
        return None
    try:
        return _canonical(region, out)
    except shapely.errors.GEOSException:
        return None


def _move_junctions(region, cells: Cells, delta, mode) -> Cells | None:
    """Move wall ends within delta of a region vertex: onto it, or delta away."""
    walls = _walls(cells)
    vertices = _boundary_vertices(region)
    boundary = region.boundary
    changed = False
    new_walls = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        for end in (0, len(coordinates) - 1):
            p = np.asarray(coordinates[end])
            if boundary.distance(Point(p)) > ON_BOUNDARY:
                continue
            d = np.linalg.norm(vertices - p, axis=1)
            k = int(np.argmin(d))
            if not (ON_BOUNDARY < d[k] < delta):
                continue
            if mode == "snap":
                coordinates[end] = tuple(vertices[k])
            else:
                s_vertex = boundary.project(Point(vertices[k]))
                s_junction = boundary.project(Point(p))
                direction = 1.0 if s_junction >= s_vertex else -1.0
                target = boundary.interpolate(s_vertex + direction * delta * 1.0001)
                coordinates[end] = (target.x, target.y)
            changed = True
        if coordinates[0] != coordinates[-1] or len(coordinates) >= 4:
            new_walls.append(LineString(coordinates))
    if not changed:
        return None
    return _canonical(region, _rebuild(region, new_walls, cells))


def _line_intersection(p1, p2, p3, p4):
    """Intersection of the lines p1p2 and p3p4, or None if (nearly) parallel."""
    d1, d2 = p2 - p1, p4 - p3
    denominator = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denominator) < 1e-12 * max(1.0, np.linalg.norm(d1) * np.linalg.norm(d2)):
        return None
    t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / denominator
    return p1 + t * d1


def _collapse_short_edges(region, cells: Cells, delta) -> Cells | None:
    """Collapse wall edges shorter than delta into a corner point or onto an end."""
    walls = _walls(cells)
    changed = False
    new_walls = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        closed = coordinates[0] == coordinates[-1]
        points = [np.asarray(c, dtype=float) for c in (coordinates[:-1] if closed else coordinates)]
        progress = True
        while progress and len(points) > (3 if closed else 2):
            progress = False
            n = len(points)
            for e in range(n if closed else n - 1):
                a_i, b_i = e, (e + 1) % n
                a, b = points[a_i], points[b_i]
                if np.linalg.norm(b - a) >= delta:
                    continue
                a_end = not closed and a_i == 0
                b_end = not closed and b_i == n - 1
                if a_end and b_end:
                    continue
                if a_end:
                    del points[b_i]
                elif b_end:
                    del points[a_i]
                else:
                    q = _line_intersection(points[(a_i - 1) % n], a, b, points[(b_i + 1) % n])
                    middle = (a + b) / 2
                    if q is None or np.linalg.norm(q - middle) > delta:
                        q = middle
                    points[a_i] = q
                    del points[b_i]
                progress = changed = True
                break
        out = [tuple(p) for p in points]
        if closed:
            out.append(out[0])
        if len(set(out)) >= 2:
            new_walls.append(LineString(out))
    if not changed:
        return None
    return _canonical(region, _rebuild(region, new_walls, cells))


def _project_near_boundary(region, cells: Cells, delta) -> Cells | None:
    """Move wall vertices within delta of the region boundary onto it.

    The thin strip between such a wall and the boundary goes to the cell on
    the other side of the wall; wall pieces that end up on the boundary are
    dropped.
    """
    walls = _walls(cells)
    boundary = region.boundary
    vertices = _boundary_vertices(region)
    changed = False
    new_walls = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        closed = coordinates[0] == coordinates[-1]
        out, on = [], []
        for i, c in enumerate(coordinates):
            p = Point(c)
            d = boundary.distance(p)
            is_end = (not closed) and i in (0, len(coordinates) - 1)
            if is_end or d <= ON_BOUNDARY or d >= delta:
                out.append(c)
                on.append(d <= ON_BOUNDARY)
                continue
            to_vertex = np.linalg.norm(vertices - np.asarray(c), axis=1)
            k = int(np.argmin(to_vertex))
            if to_vertex[k] < delta:
                q = tuple(vertices[k])
            else:
                nearest = boundary.interpolate(boundary.project(p))
                q = (nearest.x, nearest.y)
            out.append(q)
            on.append(True)
            changed = True
        if closed:
            on[-1] = on[0]
            out[-1] = out[0]
        pieces, current = [], [out[0]]
        for i in range(1, len(out)):
            a, b = out[i - 1], out[i]
            middle = Point((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            if on[i - 1] and on[i] and boundary.distance(middle) <= 1e-6:
                if len(current) >= 2:
                    pieces.append(current)
                current = [b]
            else:
                current.append(b)
        if len(current) >= 2:
            pieces.append(current)
        if closed and len(pieces) >= 2 and pieces[0][0] == pieces[-1][-1] and not on[0]:
            pieces[0] = pieces[-1][:-1] + pieces[0]
            pieces.pop()
        for piece in pieces:
            if len(set(piece)) >= 2:
                new_walls.append(LineString(piece))
    if not changed:
        return None
    return _canonical(region, _rebuild(region, new_walls, cells))


def _absorb(region, cells: Cells, indices) -> Cells | None:
    """Give cells to the neighbour with the longest shared wall."""
    out = list(cells)
    changed = False
    for index in indices:
        polygon, label = out[index]
        shared = {}
        for other, (neighbour, other_label) in enumerate(out):
            if other == index or other_label == label:
                continue
            length = sum(
                line.length for line in _lines(polygon.boundary.intersection(neighbour.boundary))
            )
            if length > 0:
                shared[other_label] = shared.get(other_label, 0.0) + length
        if not shared:
            continue
        best = max(shared.values())
        out[index] = (polygon, min(label for label, length in shared.items() if length == best))
        changed = True
    if not changed:
        return None
    return _canonical(region, out)


# ---------------------------------------------------------------------------
# Judging candidates
# ---------------------------------------------------------------------------


def _close_pairs(vertices, edges, delta):
    """Nonincident vertex-edge pairs of a graph closer than delta."""
    if not vertices or not edges:
        return []
    tolerance = delta * SEPARATION_RELATIVE_TOLERANCE
    vertex_geometries = shapely.points(vertices)
    segments = shapely.linestrings(edges)
    index = {vertex: i for i, vertex in enumerate(vertices)}
    ends = np.array([[index[a], index[b]] for a, b in edges])
    source, target = STRtree(segments).query(
        vertex_geometries, predicate="dwithin", distance=delta
    )
    keep = np.all(ends[target] != source[:, None], axis=1)
    source, target = source[keep], target[keep]
    if not len(source):
        return []
    distances = shapely.distance(vertex_geometries[source], segments[target])
    return [
        (vertices[int(i)], edges[int(j)], float(d))
        for i, j, d in zip(source, target, distances)
        if d < delta - tolerance
    ]


def _neighbours(edges):
    neighbours = {}
    for a, b in edges:
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)
    return neighbours


def _small_sectors(neighbours, vertices):
    """Incident sectors below the profile's minimum angle at the given vertices."""
    limit = DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES + _SECTOR_TOLERANCE_DEGREES
    found = []
    for vertex in vertices:
        adjacent = neighbours.get(vertex, ())
        if len(adjacent) < 2:
            continue
        directions = np.asarray(adjacent, dtype=float) - np.asarray(vertex, dtype=float)
        angles = np.sort(np.arctan2(directions[:, 1], directions[:, 0]))
        sectors = np.degrees(np.diff(np.append(angles, angles[0] + 2 * math.pi)))
        smallest = float(sectors.min())
        if smallest < limit:
            found.append((vertex, smallest))
    return found


def _on_segment(point, a, b, tolerance):
    p, a, b = (np.asarray(value, dtype=float) for value in (point, a, b))
    d = b - a
    length2 = float(d @ d)
    t = 0.0 if length2 == 0 else float(np.clip((p - a) @ d / length2, 0.0, 1.0))
    return float(np.linalg.norm(p - (a + t * d))) <= tolerance


class _Judge:
    """Defects and fidelity of candidates against fixed per-label budgets.

    The defects of a candidate are the incident sectors below the profile's
    minimum angle at vertices that walls add or change, and the feature pairs
    closer than delta that the region alone does not have. Candidates are
    ranked by the measure (small sectors, wall pairs, total pair deficit).
    """

    def __init__(self, region, partition, delta, epsilon):
        self.region = region
        self.delta = delta
        self.area = region.area
        clip = region.buffer(2 * epsilon + 4 * GRID)
        self.budgets = {}
        for label, geometry in partition.items():
            reference = _polygonal(geometry.intersection(clip))
            if not reference.is_empty and reference.area > 0:
                self.budgets[label] = FidelityBudget([reference], epsilon)
        # Labels without a protected core may vanish.
        self.dispensable = {
            label for label, budget in self.budgets.items() if budget.protected_occupied_area == 0
        }
        vertices, edges = essential_boundary_graph([region.boundary])
        self.region_degree = {
            vertex: len(adjacent) for vertex, adjacent in _neighbours(edges).items()
        }
        self.region_pairs = {}
        for vertex, edge, _ in _close_pairs(vertices, edges, delta):
            self.region_pairs.setdefault(vertex, []).append(edge)
        self.evaluations = 0
        self._fidelity = {}

    def _is_region_pair(self, vertex, edge):
        for region_edge in self.region_pairs.get(vertex, ()):
            if _on_segment(edge[0], *region_edge, 4 * GRID) and _on_segment(
                edge[1], *region_edge, 4 * GRID
            ):
                return True
        return False

    def defects(self, cells: Cells):
        """Small sectors and wall pairs, or None for cells that are invalid or
        do not cover the region."""
        self.evaluations += 1
        polygons = [polygon for polygon, _ in cells]
        if any(polygon.is_empty or not polygon.is_valid for polygon in polygons):
            return None
        if abs(sum(polygon.area for polygon in polygons) - self.area) > 1e-6 * max(1.0, self.area):
            return None
        vertices, edges = essential_boundary_graph([polygon.boundary for polygon in polygons])
        neighbours = _neighbours(edges)
        changed = [
            vertex
            for vertex, adjacent in neighbours.items()
            if self.region_degree.get(vertex) != len(adjacent)
        ]
        sectors = _small_sectors(neighbours, changed)
        pairs = [
            pair
            for pair in _close_pairs(vertices, edges, self.delta)
            if not self._is_region_pair(pair[0], pair[1])
        ]
        return sectors, pairs

    def measure(self, sectors, pairs):
        return (
            len(sectors),
            len(pairs),
            sum(self.delta - distance for _, _, distance in pairs),
        )

    def faithful(self, cells: Cells) -> bool:
        """Whether (2L) holds for every label."""
        by_label = {}
        for polygon, label in cells:
            by_label.setdefault(label, []).append(polygon)
        if any(label not in self.budgets for label in by_label):
            return False
        for label, budget in self.budgets.items():
            own = by_label.get(label, [])
            key = (label, tuple(sorted(polygon.normalize().wkb for polygon in own)))
            passed = self._fidelity.get(key)
            if passed is None:
                passed = budget.measure(own)["status"] == "pass"
                self._fidelity[key] = passed
            if not passed:
                return False
        return True

    def evaluate(self, cells: Cells, *, below=None):
        """Defects and their measure for faithful cells measuring below ``below``."""
        found = self.defects(cells)
        if found is None:
            return None
        sectors, pairs = found
        measure = self.measure(sectors, pairs)
        if below is not None and measure >= below:
            return None
        if not self.faithful(cells):
            return None
        return {"sectors": sectors, "pairs": pairs, "measure": measure}


# ---------------------------------------------------------------------------
# Bounded repair
# ---------------------------------------------------------------------------


def _replace(walls, w_index, coordinates):
    new = list(walls)
    if len(set(map(tuple, coordinates))) >= 2:
        new[w_index] = LineString(coordinates)
    else:
        new.pop(w_index)
    return new


def _site_wall_offers(region, cells: Cells, site, delta):
    """Wall edits near one site, each as a list of walls."""
    walls = _walls(cells)
    boundary = region.boundary
    vertices = _boundary_vertices(region)
    site = np.asarray(site, dtype=float)
    for w_index, wall in enumerate(walls):
        coordinates = [tuple(c[:2]) for c in wall.coords]
        closed = coordinates[0] == coordinates[-1]
        ring = coordinates[:-1] if closed else coordinates
        for v_index, c in enumerate(ring):
            if np.min(np.hypot(site[:, 0] - c[0], site[:, 1] - c[1])) > SITE_REACH * delta:
                continue
            is_end = (not closed) and v_index in (0, len(ring) - 1)
            on_boundary = boundary.distance(Point(c)) <= ON_BOUNDARY
            point = np.asarray(c)
            if is_end and on_boundary:
                d = np.linalg.norm(vertices - point, axis=1)
                k = int(np.argmin(d))
                if ON_BOUNDARY < d[k] < 2 * delta:
                    moved = list(ring)
                    moved[v_index] = tuple(vertices[k])
                    yield "snap_junction", _replace(walls, w_index, moved)
                    s_vertex = boundary.project(Point(vertices[k]))
                    s_junction = boundary.project(Point(point))
                    direction = 1.0 if s_junction >= s_vertex else -1.0
                    target = boundary.interpolate(s_vertex + direction * delta * 1.0001)
                    moved = list(ring)
                    moved[v_index] = (target.x, target.y)
                    yield "slide_junction", _replace(walls, w_index, moved)
                neighbour = 1 if v_index == 0 else len(ring) - 2
                if 0 < neighbour < len(ring) - 1:
                    yield "drop_next", _replace(
                        walls, w_index, [p for i, p in enumerate(ring) if i != neighbour]
                    )
                continue
            if is_end:
                continue
            reduced = [p for i, p in enumerate(ring) if i != v_index]
            if closed:
                reduced = reduced + [reduced[0]]
            yield "drop_vertex", _replace(walls, w_index, reduced)
            if boundary.distance(Point(point)) < delta:
                d = np.linalg.norm(vertices - point, axis=1)
                k = int(np.argmin(d))
                if d[k] < delta:
                    q = tuple(vertices[k])
                else:
                    q = tuple(
                        np.asarray(boundary.interpolate(boundary.project(Point(point))).coords[0])
                    )
                new = list(walls)
                new.pop(w_index)
                if closed:
                    new.append(LineString([q] + ring[v_index + 1 :] + ring[:v_index] + [q]))
                else:
                    for piece in (ring[:v_index] + [q], [q] + ring[v_index + 1 :]):
                        if len(set(piece)) >= 2:
                            new.append(LineString(piece))
                yield "project_vertex", new
            n = len(ring)
            for other in (v_index - 1, v_index + 1):
                if not closed and not (0 < other < n - 1):
                    continue
                o = other % n
                a, b = np.asarray(ring[v_index]), np.asarray(ring[o])
                if np.linalg.norm(a - b) >= delta:
                    continue
                lo, hi = (v_index, o) if o == (v_index + 1) % n else (o, v_index)
                q = _line_intersection(
                    np.asarray(ring[(lo - 1) % n]),
                    np.asarray(ring[lo]),
                    np.asarray(ring[hi]),
                    np.asarray(ring[(hi + 1) % n]),
                )
                middle = (np.asarray(ring[lo]) + np.asarray(ring[hi])) / 2
                if q is None or np.linalg.norm(q - middle) > delta:
                    q = middle
                merged = list(ring)
                merged[lo] = tuple(q)
                del merged[hi]
                if closed:
                    merged = merged + [merged[0]]
                yield "collapse_edge", _replace(walls, w_index, merged)


def _site_offers(region, cells: Cells, site, delta, dispensable):
    """Candidate cells for one site."""
    for name, walls in _site_wall_offers(region, cells, site, delta):
        try:
            yield name, _canonical(region, _rebuild(region, walls, cells))
        except shapely.errors.GEOSException:
            continue
    site_geometry = MultiPoint([tuple(p) for p in site]).buffer(SITE_REACH * delta)
    near = [
        index
        for index, (polygon, label) in enumerate(cells)
        if label in dispensable and polygon.intersects(site_geometry)
    ]
    for index in near:
        moved = _absorb(region, cells, [index])
        if moved is not None:
            yield "absorb_cell", moved


def _global_offers(region, cells: Cells, delta, dispensable):
    """Whole-region wall operations, each followed by simplification at delta."""
    operations = (
        ("collapse_all", lambda c: _collapse_short_edges(region, c, delta)),
        ("project_all", lambda c: _project_near_boundary(region, c, delta)),
        ("snap_all", lambda c: _move_junctions(region, c, delta, "snap")),
        ("slide_all", lambda c: _move_junctions(region, c, delta, "slide")),
        (
            "absorb_all",
            lambda c: _absorb(
                region, c, [i for i, (_, label) in enumerate(c) if label in dispensable]
            ),
        ),
    )
    for name, operation in operations:
        try:
            out = operation(cells)
        except shapely.errors.GEOSException:
            out = None
        if out is None:
            continue
        again = _simplify_walls(region, out, delta)
        yield name, (out if again is None else again)


def _sites(evaluation, radius):
    """Group defect locations: small-sector vertices, and the vertex of each
    wall pair with its closest point on the edge."""
    groups = [[tuple(vertex)] for vertex, _ in evaluation["sectors"]]
    for vertex, edge, _ in evaluation["pairs"]:
        line = LineString(edge)
        foot = line.interpolate(line.project(Point(vertex)))
        groups.append([tuple(vertex), (foot.x, foot.y)])
    sites = []
    for group in groups:
        for site in sites:
            if min(math.hypot(p[0] - q[0], p[1] - q[1]) for p in group for q in site) <= radius:
                site.extend(group)
                break
        else:
            sites.append(list(group))
    return sites


def _repair(region, exact: Cells, judge: _Judge, delta, *, try_exact=True):
    """Repair the exact labelling, or its simplest faithful simplification.

    Returns the cells (None when no faithful start exists) and a report. The
    report's outcome is ``strict`` without defects, ``warning`` when only wall
    pairs remain, and ``fallback`` when small sectors remain.
    """
    report = {"edits": 0, "operations": []}
    evaluation, state = None, exact
    if try_exact:
        evaluation = judge.evaluate(exact)
        report["start_tolerance"] = 0.0
    if evaluation is None or evaluation["measure"] != (0, 0, 0):
        evaluation = None
        for factor in WALL_LADDER:
            if factor == 0 and not try_exact:
                continue
            candidate = (
                _simplify_walls(region, exact, factor * delta) if factor > 0 else exact
            )
            if candidate is None:
                continue
            found = judge.evaluate(candidate)
            if found is not None:
                evaluation, state = found, candidate
                report["start_tolerance"] = factor * delta
                break
    if evaluation is None:
        report.update(outcome="fallback", reason="no_faithful_start", evaluations=judge.evaluations)
        return None, report

    seen = {_key(state)}

    def best_of(offers):
        best = None
        for name, cells in offers:
            key = _key(cells)
            if key in seen:
                continue
            seen.add(key)
            below = evaluation["measure"] if best is None else best[0]["measure"]
            candidate = judge.evaluate(cells, below=below)
            if candidate is not None:
                best = (candidate, cells, name)
        return best

    for _ in range(MAX_ROUNDS):
        if evaluation["measure"] == (0, 0, 0):
            break
        improved = False
        best = best_of(_global_offers(region, state, delta, judge.dispensable))
        if best is not None:
            evaluation, state = best[0], best[1]
            report["edits"] += 1
            report["operations"].append(best[2])
            improved = True
        for site in _sites(evaluation, SITE_RADIUS * delta):
            if evaluation["measure"] == (0, 0, 0):
                break
            best = best_of(_site_offers(region, state, site, delta, judge.dispensable))
            if best is not None:
                evaluation, state = best[0], best[1]
                report["edits"] += 1
                report["operations"].append(best[2])
                improved = True
        if not improved:
            break
    sectors, pairs = evaluation["sectors"], evaluation["pairs"]
    report.update(
        outcome="fallback" if sectors else "warning" if pairs else "strict",
        small_sectors=len(sectors),
        wall_pairs=len(pairs),
        wall_deficit=evaluation["measure"][2],
        min_wall_distance=min((distance for _, _, distance in pairs), default=None),
        evaluations=judge.evaluations,
        walls=len(_walls(state)),
    )
    if sectors:
        report["reason"] = "small_sectors"
    return state, report


def _key(cells: Cells):
    return tuple(sorted((label, polygon.normalize().wkb) for polygon, label in cells))
