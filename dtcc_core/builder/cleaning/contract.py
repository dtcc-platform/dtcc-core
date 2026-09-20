"""Independent geometric checks for docs/design/footprint-cleaning-contract.md.

These observe conformance independently of construction and warning-continuation
policy. Distances are in the input's planar coordinate units. No repair operators
or mesher heuristics participate in these checks.
"""

from __future__ import annotations

import math

import numpy as np
from shapely import STRtree, distance, get_coordinates, linestrings, make_valid, points
from shapely.affinity import translate
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union

_SEPARATION_RELATIVE_TOLERANCE = 1e-6
_INCIDENT_SECTOR_DEGREE_TOLERANCE = 1e-6
DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES = 1.0
DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES = 3.0


def _validate_scale(value, name, *, positive):
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(
            f"{name} must be finite and {'positive' if positive else 'nonnegative'}"
        )


def _interpret_input(raw):
    polygons = []
    repaired_count = 0
    remnant_count = 0

    def collect(geometry):
        nonlocal remnant_count
        if geometry.is_empty:
            return
        if isinstance(geometry, Polygon):
            polygons.append(geometry)
        elif hasattr(geometry, "geoms"):
            for part in geometry.geoms:
                collect(part)
        else:
            remnant_count += 1

    for geometry in raw:
        if not isinstance(geometry, (Polygon, MultiPolygon)):
            raise ValueError(
                "Raw footprints must be Polygon or MultiPolygon geometries"
            )
        if not np.isfinite(get_coordinates(geometry)).all():
            raise ValueError("Raw footprints must have finite coordinates")
        repaired_count += not geometry.is_valid
        collect(make_valid(geometry))
    return unary_union(polygons), {
        "method": "GEOS make_valid; polygonal parts only",
        "repaired_input_count": repaired_count,
        "nonarea_remnant_count": remnant_count,
    }


def _lines(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type in {"LineString", "LinearRing"}:
        return [geometry]
    return [line for part in geometry.geoms for line in _lines(part)]


def _canonical_graph(boundaries):
    """Node exact shared boundaries; suppress only collinear degree-two nodes.

    No snapping or cleanup tolerance: a tiny genuine bend remains a feature.
    The exact-collinearity decision here uses floating-point arithmetic.
    """
    neighbors = {}
    for line in _lines(unary_union(boundaries)):
        coords = [tuple(c[:2]) for c in line.coords]
        for a, b in zip(coords, coords[1:]):
            if a != b:
                neighbors.setdefault(a, set()).add(b)
                neighbors.setdefault(b, set()).add(a)
    pending = list(neighbors)
    while pending:
        p = pending.pop()
        if p not in neighbors or len(neighbors[p]) != 2:
            continue
        a, b = sorted(neighbors[p])
        u, v = np.subtract(a, p), np.subtract(b, p)
        if u[0] * v[1] - u[1] * v[0] == 0 and np.dot(u, v) < 0:
            neighbors[a].remove(p)
            neighbors[b].remove(p)
            neighbors[a].add(b)
            neighbors[b].add(a)
            del neighbors[p]
            pending.extend([a, b])
    vertices = sorted(neighbors)
    edges = sorted((a, b) for a in vertices for b in neighbors[a] if a < b)
    return vertices, edges


def admissibility(polygons, delta):
    """Measure topology and feature separation without using cleaner operators.

    In a planar graph, distances between disjoint segments are attained at an
    endpoint, so vertex/vertex and nonincident vertex/edge pairs suffice.
    Spatial queries only seek pairs below delta: the reported minimum is capped.
    """
    _validate_scale(delta, "delta", positive=True)
    polygons = list(polygons)
    if any(not isinstance(p, Polygon) for p in polygons):
        return {"topology_ok": False, "resolved": False, "reason": "expected polygons"}
    if any(
        p.is_empty or not p.is_valid or not np.isfinite(np.asarray(r.coords)).all()
        for p in polygons
        for r in [p.exterior, *p.interiors]
    ):
        return {"topology_ok": False, "resolved": False, "reason": "invalid polygon"}
    if len(polygons) == 1:
        vertices, edges = _canonical_graph([polygons[0].boundary])
        occupancy_edges = edges
        overlap_pairs = 0
    else:
        occupied = unary_union(polygons)
        tree = STRtree(polygons)
        overlap_pairs = sum(
            1
            for i, p in enumerate(polygons)
            for j in tree.query(p, predicate="intersects")
            if j > i and p.relate_pattern(polygons[j], "T********")
        )
        _, occupancy_edges = _canonical_graph(
            [occupied.boundary] if polygons else []
        )
        vertices, edges = _canonical_graph([p.boundary for p in polygons])
    degree = {}
    for a, b in occupancy_edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    # Occupied/exterior boundaries must be disjoint simple cycles. Internal
    # parcel junctions may have higher valence in the full subdivision.
    nonmanifold = sum(d != 2 for d in degree.values())
    minimum = delta
    short_pairs = 0
    witness = None
    tolerance = delta * _SEPARATION_RELATIVE_TOLERANCE
    if vertices:
        point_geometries = points(vertices)
        segments = linestrings(edges)
        point_tree = STRtree(point_geometries)
        edge_tree = STRtree(segments)
        vertex_ids = {xy: i for i, xy in enumerate(vertices)}
        endpoints = np.array([[vertex_ids[a], vertex_ids[b]] for a, b in edges])
        closest_key = (delta, len(vertices), 2)
        # Batch GEOS queries/distances instead of crossing Python for each pair.
        # Fixed-size batches bound query-result memory on dense input graphs.
        for start in range(0, len(vertices), 256):
            batch = point_geometries[start : start + 256]
            for kind, tree, targets in (
                (0, point_tree, point_geometries),
                (1, edge_tree, segments),
            ):
                source, target = tree.query(batch, predicate="dwithin", distance=delta)
                source = source + start
                keep = (
                    target > source
                    if kind == 0
                    else np.all(endpoints[target] != source[:, None], axis=1)
                )
                source, target = source[keep], target[keep]
                if not len(source):
                    continue
                distances = distance(point_geometries[source], targets[target])
                short_pairs += int(np.count_nonzero(distances < delta - tolerance))
                k = int(np.argmin(distances))
                i, j = int(source[k]), int(target[k])
                key = (float(distances[k]), i, kind)
                # Preserve the old witness tie order: vertex order, then VV
                # before VE, then the tree's order within a vertex query.
                if key[0] < delta and key < closest_key:
                    closest_key = key
                    minimum = key[0]
                    witness = (
                        [vertices[i], vertices[j]]
                        if kind == 0
                        else [
                            vertices[i],
                            tuple(
                                segments[j]
                                .interpolate(segments[j].project(point_geometries[i]))
                                .coords[0]
                            ),
                        ]
                    )
    topology_ok = not overlap_pairs and nonmanifold == 0
    return {
        "topology_ok": topology_ok,
        "resolved": topology_ok and short_pairs == 0,
        "interior_overlap_pairs": overlap_pairs,
        "nonmanifold_boundary_vertices": nonmanifold,
        "canonical_vertex_count": len(vertices),
        "subscale_pairs": short_pairs,
        "separation_capped_at_delta": minimum,
        "closest_pair": witness,
    }


def incident_sectors(
    polygons,
    *,
    minimum_degrees=DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES,
):
    """Measure face sectors at vertices of the canonical labelled subdivision.

    Each polygon is a label-bearing region.  Exact shared walls are represented
    once in the canonical graph, while the point sampled inside each sector is
    tested against every region so occupied/occupied source junctions are not
    mistaken for open ground.  The angle is a property of consecutive incident
    graph rays, not an unsigned turn from one ring orientation.
    """
    _validate_scale(minimum_degrees, "minimum_degrees", positive=False)
    if minimum_degrees > 180:
        raise ValueError("minimum_degrees must not exceed 180")
    polygons = list(polygons)
    invalid_reason = None
    if any(not isinstance(polygon, Polygon) for polygon in polygons):
        invalid_reason = "expected polygons"
    elif any(
        polygon.is_empty
        or not polygon.is_valid
        or not np.isfinite(np.asarray(ring.coords)).all()
        for polygon in polygons
        for ring in [polygon.exterior, *polygon.interiors]
    ):
        invalid_reason = "invalid polygon"
    if invalid_reason is not None:
        return {
            "status": "fail",
            "reason": invalid_reason,
            "minimum_degrees": minimum_degrees,
            "angle_tolerance_degrees": _INCIDENT_SECTOR_DEGREE_TOLERANCE,
        }

    vertices, edges = _canonical_graph([polygon.boundary for polygon in polygons])
    polygon_tree = STRtree(polygons)
    neighbors = {vertex: [] for vertex in vertices}
    for a, b in edges:
        neighbors[a].append(b)
        neighbors[b].append(a)

    minimum = 360.0
    minimum_witness = None
    sector_count = 0
    occupied_count = 0
    open_count = 0
    below_count = 0
    below = []
    borderline = []
    for vertex in vertices:
        adjacent = neighbors[vertex]
        if len(adjacent) < 2:
            continue
        rays = sorted(
            (
                math.atan2(other[1] - vertex[1], other[0] - vertex[0]),
                other,
                math.hypot(other[0] - vertex[0], other[1] - vertex[1]),
            )
            for other in adjacent
        )
        for index, (angle, _, length) in enumerate(rays):
            next_angle, _, next_length = rays[(index + 1) % len(rays)]
            sector = (next_angle - angle) % (2 * math.pi)
            if sector <= 0:
                sector = 2 * math.pi
            degrees = math.degrees(sector)
            midpoint = angle + sector / 2
            sample_radius = min(0.1, max(1e-7, min(length, next_length) * 0.1))
            sample = Point(
                vertex[0] + sample_radius * math.cos(midpoint),
                vertex[1] + sample_radius * math.sin(midpoint),
            )
            labels = [
                int(i)
                for i in polygon_tree.query(sample)
                if polygons[int(i)].covers(sample)
            ]
            labels.sort()
            kind = "occupied" if labels else "open"
            sector_count += 1
            occupied_count += kind == "occupied"
            open_count += kind == "open"
            witness = {
                "vertex": list(vertex),
                "angle_degrees": degrees,
                "kind": kind,
                "region_indices": labels,
            }
            if degrees < minimum:
                minimum = degrees
                minimum_witness = witness
            if degrees < minimum_degrees - _INCIDENT_SECTOR_DEGREE_TOLERANCE:
                below_count += 1
                if len(below) < 8:
                    below.append(witness)
            elif abs(degrees - minimum_degrees) <= _INCIDENT_SECTOR_DEGREE_TOLERANCE:
                if len(borderline) < 8:
                    borderline.append(witness)

    status = "fail" if below else "borderline" if borderline else "pass"
    return {
        "status": status,
        "minimum_degrees": minimum_degrees,
        "angle_tolerance_degrees": _INCIDENT_SECTOR_DEGREE_TOLERANCE,
        "sector_count": sector_count,
        "occupied_sector_count": occupied_count,
        "open_sector_count": open_count,
        "minimum_angle_degrees": minimum if sector_count else None,
        "minimum_witness": minimum_witness,
        "below_minimum_count": below_count,
        "below_minimum_witnesses": below,
        "borderline_witnesses": borderline,
    }


def check_mesher_handoff_profile(
    polygons,
    *,
    minimum_sector_degrees=DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES,
):
    """Validate the independently declared city-meshing angle profile."""
    sectors = incident_sectors(
        polygons,
        minimum_degrees=minimum_sector_degrees,
    )
    return {
        "status": sectors["status"],
        "profile": "city_meshing_v1",
        "incident_sectors": sectors,
    }


class FidelityBudget:
    """Fixed occupied/open cores, prepared once from the original input.

    Candidate admission and reporting use the same numerical convention. Never
    construct a new budget from an intermediate repair: that would permit drift
    to accumulate beyond epsilon relative to the input.
    """

    def __init__(self, raw, epsilon):
        _validate_scale(epsilon, "epsilon", positive=False)
        occupied, self.interpretation = _interpret_input(raw)
        # GEOS overlay can lose thin differences at large map coordinates.
        # Keep one fixed frame for buffers, admission and reporting; never
        # choose a new origin from a candidate or an intermediate repair.
        self._origin = occupied.bounds[:2] if not occupied.is_empty else (0.0, 0.0)
        occupied = self._to_local(occupied)
        self.epsilon = epsilon
        quad_segs = 32
        self.uncertainty = (
            epsilon * (1 / math.cos(math.pi / (4 * quad_segs)) - 1) + 1e-7
        )
        self._local_bounds = [
            (
                occupied.buffer(-radius, quad_segs=quad_segs),
                occupied.buffer(radius, quad_segs=quad_segs),
            )
            for radius in (
                epsilon,
                max(0.0, epsilon - self.uncertainty),
                epsilon + self.uncertainty,
            )
        ]

    def _to_local(self, geometry):
        return translate(geometry, -self._origin[0], -self._origin[1])

    @property
    def _bounds(self):
        """World-coordinate views for research proposal geometry, not checking."""
        return [
            tuple(translate(g, *self._origin) for g in pair)
            for pair in self._local_bounds
        ]

    @staticmethod
    def _violations(bounds, candidate):
        protected, allowed = bounds
        return protected.difference(candidate).area, candidate.difference(allowed).area

    def _accepts_local_union(self, candidate):
        """Admit a union already expressed in this budget's fixed frame."""
        protected, allowed = self._local_bounds[1]
        if protected.difference(candidate).area > 1e-10:
            return False
        return candidate.difference(allowed).area <= 1e-10

    def accepts_union(self, candidate):
        """Admit an internally validated coverage union only on a definite pass."""
        return self._accepts_local_union(self._to_local(candidate))

    def can_drop_input_part(self, original):
        """Whether an original part contains no protected occupied interior.

        Use the fixed union's core: eroding this part alone would overlook
        protected space spanning shared walls. This tests complete disappearance
        only; it does not certify other changes made by an opening operation.
        """
        return (
            self._local_bounds[1][0].intersection(self._to_local(original)).area
            <= 1e-10
        )

    @property
    def protected_occupied_area(self):
        """Area that an accepted result must retain from the original union."""
        return float(self._local_bounds[1][0].area)

    def measure(self, cleaned):
        cleaned = list(cleaned)
        if any(
            not isinstance(g, Polygon)
            or g.is_empty
            or not g.is_valid
            or not np.isfinite(get_coordinates(g)).all()
            for g in cleaned
        ):
            return {"status": "not_checked", "reason": "invalid output"}
        candidate = unary_union([self._to_local(p) for p in cleaned])
        lost, added = self._violations(self._local_bounds[0], candidate)
        status = (
            "pass"
            if max(self._violations(self._local_bounds[1], candidate)) <= 1e-10
            else (
                "fail"
                if max(self._violations(self._local_bounds[2], candidate)) > 1e-10
                else "borderline"
            )
        )
        return {
            "status": status,
            "epsilon": self.epsilon,
            "input_interpretation": self.interpretation,
            "lost_protected_area": lost,
            "added_outside_budget_area": added,
            "buffer_uncertainty": self.uncertainty,
        }


def has_only_separation_defects(observation):
    """Classify an admissibility observation without changing its verdict."""
    return (
        observation.get("topology_ok") is True
        and observation.get("nonmanifold_boundary_vertices") == 0
        and observation.get("resolved") is False
        and observation.get("subscale_pairs", 0) > 0
    )


def fidelity(raw, cleaned, epsilon):
    """Measure the two-sided offset budget; numerical borderlines are inconclusive."""
    return FidelityBudget(raw, epsilon).measure(cleaned)


def check_cleaning_contract(raw, cleaned, *, delta, epsilon=None):
    """Return independent admissibility and fidelity observations.

    The initial fidelity budget is delta / 2; callers can choose epsilon
    independently. A borderline numerical result is never reported as a pass.
    Invalid raw inputs follow the documented GEOS interpretation, with non-area
    remnants counted. Malformed scales and unsupported inputs fail clearly.
    """
    _validate_scale(delta, "delta", positive=True)
    epsilon = delta / 2 if epsilon is None else epsilon
    _validate_scale(epsilon, "epsilon", positive=False)
    cleaned = list(cleaned)
    output = admissibility(cleaned, delta)
    drift = fidelity(raw, cleaned, epsilon)
    status = "fail" if not output["resolved"] else drift["status"]
    return {
        "status": status,
        "delta": delta,
        "epsilon": epsilon,
        "separation_tolerance": delta * _SEPARATION_RELATIVE_TOLERANCE,
        "admissibility": output,
        "fidelity": drift,
    }
