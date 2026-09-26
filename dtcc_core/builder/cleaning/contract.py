"""Independent checks of the footprint cleaning contract.

The contract (docs/design/footprint-cleaning-contract.md) follows the footprint
conditioning note. For the interpreted input occupied set P, an output
subdivision Q with occupied union U(Q), a resolution delta and a fidelity budget
epsilon, cleaning seeks

    Q in A_{delta,theta}  with  P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon}.

The checks here observe the three obligations separately:

* admissibility (A_delta): the union topology profile and a minimum feature
  size mfs(G_Q) >= delta of the essential boundary graph G_Q;
* fidelity: the two-sided erosion/dilation budget against the original occupied
  set, never against an intermediate repair;
* the mesher-input incident-sector profile (theta, 1 degree by default).

No repair operators or mesher heuristics participate. Distances are in the
input's planar coordinate units. A passing check is numerical evidence of
conformance, not an exact certificate.
"""

from __future__ import annotations

import math

import numpy as np
import shapely
from shapely import STRtree, distance, get_coordinates, linestrings, make_valid, points
from shapely.affinity import translate
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union

# Distances below delta * (1 - SEPARATION_RELATIVE_TOLERANCE) are violations.
SEPARATION_RELATIVE_TOLERANCE = 1e-6
# The topology profile implemented by `admissibility` (note, Section 2).
TOPOLOGY_PROFILE = "union"
_INCIDENT_SECTOR_DEGREE_TOLERANCE = 1e-6
DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES = 1.0
DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES = 3.0
# Area below which a fidelity violation counts as numerical noise.
_FIDELITY_AREA_TOLERANCE = 1e-10
_FIDELITY_QUAD_SEGMENTS = 32
# Bisection steps for the achieved budget: resolution (epsilon + u) / 2**7.
_ACHIEVED_BUDGET_STEPS = 7


def _validate_scale(value, name, *, positive):
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(
            f"{name} must be finite and {'positive' if positive else 'nonnegative'}"
        )


def interpret_input(raw):
    """Return the occupied set of raw footprints and how it was interpreted.

    Invalid geometries are interpreted with GEOS ``make_valid``; only polygonal
    parts are kept and lower-dimensional remnants are counted.

    Parameters
    ----------
    raw : iterable of Polygon or MultiPolygon
        Raw footprints with finite coordinates.

    Returns
    -------
    tuple[BaseGeometry, dict]
        The union of the interpreted polygons and a report with the method,
        the number of repaired inputs and the number of non-area remnants.

    Raises
    ------
    ValueError
        If an input is not polygonal or has nonfinite coordinates.
    """
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


def essential_boundary_graph(boundaries):
    """Return the essential boundary graph G_Q of region boundaries.

    Boundaries are noded so that shared pieces occur once; then every vertex of
    degree two whose edges continue along the same straight line in opposite
    directions is suppressed. The remaining vertices are exactly the bends and
    junctions. There is no snapping or cleanup tolerance, so a tiny genuine
    bend remains a feature; the collinearity decision uses floating-point
    arithmetic.

    Parameters
    ----------
    boundaries : list of BaseGeometry
        Region boundaries (lines or rings).

    Returns
    -------
    tuple[list, list]
        Sorted vertex coordinates and sorted edges as coordinate pairs.
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
        a, b = neighbors[p]
        if b < a:
            a, b = b, a
        ux, uy = a[0] - p[0], a[1] - p[1]
        vx, vy = b[0] - p[0], b[1] - p[1]
        if ux * vy - uy * vx == 0 and ux * vx + uy * vy < 0:
            neighbors[a].remove(p)
            neighbors[b].remove(p)
            neighbors[a].add(b)
            neighbors[b].add(a)
            del neighbors[p]
            pending.extend([a, b])
    vertices = sorted(neighbors)
    edges = sorted((a, b) for a in vertices for b in neighbors[a] if a < b)
    return vertices, edges


def _subscale_pairs(vertices, edges, delta, *, vertex_pairs):
    """Count feature pairs closer than delta and find the closest one.

    Nonincident vertex-edge pairs define mfs(G). Vertex-vertex pairs are
    counted only on request; by the note's Observation 1 they never change
    whether mfs(G) >= delta.
    """
    counts = [0, 0]  # vertex-vertex, vertex-edge
    minimum = delta
    witness = None
    if not vertices:
        return counts, minimum, witness
    tolerance = delta * SEPARATION_RELATIVE_TOLERANCE
    point_geometries = points(vertices)
    segments = linestrings(edges)
    vertex_ids = {xy: i for i, xy in enumerate(vertices)}
    endpoints = np.array([[vertex_ids[a], vertex_ids[b]] for a, b in edges])
    searches = [(1, STRtree(segments), segments)]
    if vertex_pairs:
        searches.insert(0, (0, STRtree(point_geometries), point_geometries))
    closest_key = (delta, len(vertices), 2)
    # Batch GEOS queries/distances instead of crossing Python for each pair.
    # Fixed-size batches bound query-result memory on dense input graphs.
    for start in range(0, len(vertices), 256):
        batch = point_geometries[start : start + 256]
        for kind, tree, targets in searches:
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
            counts[kind] += int(np.count_nonzero(distances < delta - tolerance))
            k = int(np.argmin(distances))
            i, j = int(source[k]), int(target[k])
            key = (float(distances[k]), i, kind)
            # Deterministic witness order: vertex order, then VV before VE,
            # then the tree's order within a vertex query.
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
    return counts, minimum, witness


def _invalid_polygons(polygons):
    """Whether any polygon is empty, invalid or has nonfinite coordinates."""
    return bool(
        polygons
        and (
            (shapely.is_empty(polygons) | ~shapely.is_valid(polygons)).any()
            or not np.isfinite(get_coordinates(polygons)).all()
        )
    )


def admissibility(polygons, delta, *, vertex_pairs=False, union_parts=False):
    """Check membership of a subdivision in A_delta (union topology profile).

    Topology: regions must be valid with pairwise disjoint interiors, and the
    boundary of their union must consist of disjoint simple closed curves, so
    two regions may share walls but not touch at a single point. Feature size:
    mfs(G_Q), the smallest distance between a vertex of the essential boundary
    graph and an edge not incident to it, must be at least ``delta``. Spatial
    queries only seek pairs closer than delta, so the reported minimum is
    capped at delta.

    Parameters
    ----------
    polygons : iterable of Polygon
        The output regions.
    delta : float
        The resolution, a positive distance.
    vertex_pairs : bool, optional
        Also count vertex-vertex pairs closer than delta, reported as
        ``subscale_vertex_pairs``. They do not change the verdict; the
        constructor uses the finer count to rank proposals.
    union_parts : bool, optional
        The polygons are the parts of one union, so their interiors are
        disjoint and their boundaries are the boundary of the union. The
        overlap and union-boundary computations are then skipped; the result
        is the same.

    Returns
    -------
    dict
        ``topology_profile``, ``topology_ok``, ``admissible``, the counts of
        interior overlaps and nonmanifold union-boundary vertices, the number
        of essential vertices, ``subscale_pairs`` (nonincident vertex-edge
        pairs closer than delta), ``mfs_capped_at_delta`` and a witness pair.
    """
    _validate_scale(delta, "delta", positive=True)
    polygons = list(polygons)
    failed = {
        "topology_profile": TOPOLOGY_PROFILE,
        "topology_ok": False,
        "admissible": False,
    }
    if any(not isinstance(p, Polygon) for p in polygons):
        return {**failed, "reason": "expected polygons"}
    if _invalid_polygons(polygons):
        return {**failed, "reason": "invalid polygon"}
    if len(polygons) == 1 or union_parts:
        vertices, edges = essential_boundary_graph([p.boundary for p in polygons])
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
        _, occupancy_edges = essential_boundary_graph(
            [occupied.boundary] if polygons else []
        )
        vertices, edges = essential_boundary_graph([p.boundary for p in polygons])
    degree = {}
    for a, b in occupancy_edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    # The occupied union's boundary must be disjoint simple cycles. Internal
    # parcel junctions may have higher valence in the full subdivision.
    nonmanifold = sum(d != 2 for d in degree.values())
    (vertex_vertex, vertex_edge), minimum, witness = _subscale_pairs(
        vertices, edges, delta, vertex_pairs=vertex_pairs
    )
    topology_ok = not overlap_pairs and nonmanifold == 0
    report = {
        "topology_profile": TOPOLOGY_PROFILE,
        "topology_ok": topology_ok,
        "admissible": topology_ok and vertex_edge + vertex_vertex == 0,
        "interior_overlap_pairs": overlap_pairs,
        "nonmanifold_boundary_vertices": nonmanifold,
        "essential_vertex_count": len(vertices),
        "subscale_pairs": vertex_edge,
        "mfs_capped_at_delta": minimum,
        "closest_pair": witness,
    }
    if vertex_pairs:
        report["subscale_vertex_pairs"] = vertex_vertex
    return report


def incident_sectors(
    polygons,
    *,
    minimum_degrees=DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES,
):
    """Measure incident sectors at the vertices of the essential boundary graph.

    Each polygon is a label-bearing region. Exact shared walls are represented
    once in the graph, while the point sampled inside each sector is tested
    against every region, so occupied/occupied source junctions are not
    mistaken for open ground. The smallest sector is alpha_min(Q) of the note;
    the angle is a property of consecutive incident graph rays, not an
    unsigned turn from one ring orientation.
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

    vertices, edges = essential_boundary_graph(
        [polygon.boundary for polygon in polygons]
    )
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
        occupied, self.interpretation = interpret_input(raw)
        # GEOS overlay can lose thin differences at large map coordinates.
        # Keep one fixed frame for buffers, admission and reporting; never
        # choose a new origin from a candidate or an intermediate repair.
        self._origin = occupied.bounds[:2] if not occupied.is_empty else (0.0, 0.0)
        occupied = self.to_local(occupied)
        self._occupied_local = occupied
        self.epsilon = epsilon
        self.uncertainty = (
            epsilon * (1 / math.cos(math.pi / (4 * _FIDELITY_QUAD_SEGMENTS)) - 1)
            + 1e-7
        )
        self._local_bounds = [
            (
                occupied.buffer(-radius, quad_segs=_FIDELITY_QUAD_SEGMENTS),
                occupied.buffer(radius, quad_segs=_FIDELITY_QUAD_SEGMENTS),
            )
            for radius in (
                epsilon,
                max(0.0, epsilon - self.uncertainty),
                epsilon + self.uncertainty,
            )
        ]

    @property
    def origin(self):
        """Origin of the fixed local frame used for all budget computations."""
        return self._origin

    def to_local(self, geometry):
        """Translate world-coordinate geometry into the budget's local frame."""
        return translate(geometry, -self._origin[0], -self._origin[1])

    def from_local(self, geometry):
        """Translate local-frame geometry back to world coordinates."""
        return translate(geometry, *self._origin)

    @property
    def local_core_and_envelope(self):
        """Protected core and allowed envelope used for admission, in the local frame.

        Both use the radius ``epsilon - buffer_uncertainty``, so admission is a
        definite pass.
        """
        return self._local_bounds[1]

    @staticmethod
    def _violations(bounds, candidate):
        protected, allowed = bounds
        return protected.difference(candidate).area, candidate.difference(allowed).area

    def accepts_local_union(self, candidate):
        """Admit a union already expressed in this budget's local frame."""
        protected, allowed = self._local_bounds[1]
        if protected.difference(candidate).area > _FIDELITY_AREA_TOLERANCE:
            return False
        return candidate.difference(allowed).area <= _FIDELITY_AREA_TOLERANCE

    def accepts_union(self, candidate):
        """Admit an internally validated coverage union only on a definite pass."""
        return self.accepts_local_union(self.to_local(candidate))

    @property
    def protected_occupied_area(self):
        """Area that an accepted result must retain from the original union."""
        return float(self._local_bounds[1][0].area)

    def achieved(self, candidate_local):
        """Smallest budget at which a local-frame union would pass, if at most ours.

        This is epsilon_Q of the note: the least epsilon with
        P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon}. It uses the same buffered test and
        area tolerance as :meth:`measure`, on the neighbourhood of the changed
        area only. The largest distance of a vertex of the changed area from
        the input boundary is a lower bound; the value is then bracketed to a
        resolution of ``epsilon / 128`` and rounded up, so it is an upper
        bound up to the buffer chord error.

        Returns
        -------
        tuple[float, float] or None
            The achieved budget and the resolution, or None when the union
            does not pass at ``epsilon``.
        """
        occupied = self._occupied_local
        added = candidate_local.difference(occupied)
        lost = occupied.difference(candidate_local)
        check_added = added.area > _FIDELITY_AREA_TOLERANCE
        check_lost = lost.area > _FIDELITY_AREA_TOLERANCE
        if not check_added and not check_lost:
            return 0.0, 0.0
        upper = self.epsilon
        resolution = upper / 2**_ACHIEVED_BUDGET_STEPS
        # A polygonal buffer boundary lies at least reach * cos(pi / 128) from
        # its source, so offsets by at most `upper` never see the clip edge.
        reach = upper * 1.001 + 1e-6
        near_added = (
            occupied.intersection(added.buffer(reach)) if check_added else None
        )
        near_lost = occupied.intersection(lost.buffer(reach)) if check_lost else None

        def passes(radius):
            if check_added and (
                added.difference(
                    near_added.buffer(radius, quad_segs=_FIDELITY_QUAD_SEGMENTS)
                ).area
                > _FIDELITY_AREA_TOLERANCE
            ):
                return False
            return not check_lost or (
                lost.intersection(
                    near_lost.buffer(-radius, quad_segs=_FIDELITY_QUAD_SEGMENTS)
                ).area
                <= _FIDELITY_AREA_TOLERANCE
            )

        # Added vertices lie outside P, lost vertices inside it; both distances
        # are to the input boundary, which is all the near sets can reach.
        lower = 0.0
        for region, near, check in (
            (added, near_added, check_added),
            (lost, near_lost, check_lost),
        ):
            if check and not near.is_empty:
                vertices = points(get_coordinates(region))
                lower = max(lower, float(distance(vertices, near.boundary).max()))
        low, high = min(lower, upper), upper
        if high - low > resolution and passes(low + resolution):
            high = low + resolution
        elif not passes(upper):
            return None
        else:
            low = min(low + resolution, high)
        while high - low > resolution:
            middle = (low + high) / 2
            if passes(middle):
                high = middle
            else:
                low = middle
        return high, resolution

    def measure(self, cleaned, *, achieved=False):
        """Report the fidelity of ``cleaned`` against this budget.

        Parameters
        ----------
        cleaned : iterable of Polygon
            Output regions in world coordinates.
        achieved : bool, optional
            Also measure the achieved budget epsilon_Q when the output passes
            (see :meth:`achieved`).

        Returns
        -------
        dict
            ``status`` (``pass``, ``borderline``, ``fail`` or ``not_checked``),
            the budget, the input interpretation, lost protected and added
            out-of-budget areas, the buffer uncertainty and, on request,
            ``achieved_epsilon`` and ``achieved_epsilon_resolution``.
        """
        cleaned = list(cleaned)
        if any(
            not isinstance(g, Polygon)
            or g.is_empty
            or not g.is_valid
            or not np.isfinite(get_coordinates(g)).all()
            for g in cleaned
        ):
            return {"status": "not_checked", "reason": "invalid output"}
        candidate = unary_union([self.to_local(p) for p in cleaned])
        lost, added = self._violations(self._local_bounds[0], candidate)
        status = (
            "pass"
            if max(self._violations(self._local_bounds[1], candidate))
            <= _FIDELITY_AREA_TOLERANCE
            else (
                "fail"
                if max(self._violations(self._local_bounds[2], candidate))
                > _FIDELITY_AREA_TOLERANCE
                else "borderline"
            )
        )
        report = {
            "status": status,
            "epsilon": self.epsilon,
            "input_interpretation": self.interpretation,
            "lost_protected_area": lost,
            "added_outside_budget_area": added,
            "buffer_uncertainty": self.uncertainty,
        }
        if achieved:
            value = self.achieved(candidate) if status == "pass" else None
            report["achieved_epsilon"] = None if value is None else value[0]
            report["achieved_epsilon_resolution"] = (
                None if value is None else value[1]
            )
        return report


def has_only_separation_defects(observation):
    """Return whether an admissibility report fails on feature size alone."""
    return (
        observation.get("topology_ok") is True
        and observation.get("nonmanifold_boundary_vertices") == 0
        and observation.get("admissible") is False
        and observation.get("subscale_pairs", 0) > 0
    )


def classify_candidate(contract, mesher_profile):
    """Classify a checked candidate by the outcomes of the note, Section 7.

    Parameters
    ----------
    contract : dict
        Report of :func:`check_cleaning_contract`.
    mesher_profile : dict
        Report of :func:`check_mesher_handoff_profile` for the same regions.

    Returns
    -------
    str
        ``"strict"`` when admissibility, fidelity and the profile pass;
        ``"warning"`` when only feature size fails (mfs < delta); ``"fail"``
        otherwise, including borderline and unchecked reports.
    """
    if not isinstance(contract, dict) or not isinstance(mesher_profile, dict):
        return "fail"
    if mesher_profile.get("status") != "pass":
        return "fail"
    if contract.get("status") == "pass":
        return "strict"
    if contract.get("fidelity", {}).get("status") == "pass" and (
        has_only_separation_defects(contract.get("admissibility", {}))
    ):
        return "warning"
    return "fail"


def fidelity(raw, cleaned, epsilon, *, achieved=False):
    """Measure the two-sided offset budget; numerical borderlines are inconclusive."""
    return FidelityBudget(raw, epsilon).measure(cleaned, achieved=achieved)


def check_cleaning_contract(raw, cleaned, *, delta, epsilon=None, achieved=False):
    """Return independent admissibility and fidelity observations.

    The default fidelity budget is delta / 2; callers can choose epsilon
    independently. A borderline numerical result is never reported as a pass.
    Invalid raw inputs follow the documented GEOS interpretation, with non-area
    remnants counted. Malformed scales and unsupported inputs fail clearly.

    Parameters
    ----------
    raw : iterable of Polygon or MultiPolygon
        The original input footprints.
    cleaned : iterable of Polygon
        The output regions.
    delta : float
        The resolution.
    epsilon : float or None, optional
        The fidelity budget; None uses ``delta / 2``.
    achieved : bool, optional
        Also measure the achieved budget epsilon_Q when fidelity passes.

    Returns
    -------
    dict
        ``status``, ``delta``, ``epsilon``, ``separation_tolerance`` and the
        ``admissibility`` and ``fidelity`` reports.
    """
    _validate_scale(delta, "delta", positive=True)
    epsilon = delta / 2 if epsilon is None else epsilon
    _validate_scale(epsilon, "epsilon", positive=False)
    cleaned = list(cleaned)
    output = admissibility(cleaned, delta)
    drift = fidelity(raw, cleaned, epsilon, achieved=achieved)
    status = "fail" if not output["admissible"] else drift["status"]
    return {
        "status": status,
        "delta": delta,
        "epsilon": epsilon,
        "separation_tolerance": delta * SEPARATION_RELATIVE_TOLERANCE,
        "admissibility": output,
        "fidelity": drift,
    }
