"""Bounded boundary-route research. No production callers or new dependencies.

Select directed noncrossing cycles from boundary-route proposals. The --focused
experiment adds whole-chain shortcuts through unresolved vertices to local
connections and preserves original edges outside the edit neighborhoods. The
default retains the broader local-chord family without whole-chain shortcuts.
Only selected routes are polygonized. The original contract certifies output;
limits and finite-model infeasibility return no accepted geometry.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix
from shapely import STRtree
from shapely.geometry import LineString, MultiPoint, Point, Polygon, mapping
from shapely.geometry.polygon import orient
from shapely.ops import polygonize, unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning.contract import (
    check_cleaning_contract,
    _SEPARATION_RELATIVE_TOLERANCE,
    _validate_scale,
)
from sandbox.cleaning_graph_prototype import canonical_faces
from sandbox.cleaning_pipeline_probe import (
    boundary_candidates,
    parts,
    topology_examples,
    saved_city_groups,
    plot_topology,
)


class WorkLimit(Exception):
    pass


class Rows:
    def __init__(self, variables, deadline):
        self.variables, self.deadline = variables, deadline
        self.i, self.j, self.values, self.low, self.high = [], [], [], [], []

    def add(self, terms, low=-np.inf, high=np.inf):
        if len(self.low) >= 100000:
            raise WorkLimit("constraint_limit")
        if time.perf_counter() >= self.deadline:
            raise WorkLimit("time_limit")
        for column, value in terms:
            if value:
                self.i.append(len(self.low))
                self.j.append(column)
                self.values.append(value)
        self.low.append(low)
        self.high.append(high)

    def constraint(self):
        matrix = coo_matrix(
            (self.values, (self.i, self.j)), shape=(len(self.low), self.variables)
        ).tocsc()
        return LinearConstraint(matrix, self.low, self.high)


def crossing_number(a, b, p):
    """Signed crossing of the rightward ray; half-open at segment endpoints."""
    cross = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
    return int(a[1] <= p[1] < b[1] and cross > 0) - int(
        b[1] <= p[1] < a[1] and cross < 0
    )


def route_graph(proposals):
    sites = proposals["sites"]
    index = {xy: i for i, xy in enumerate(sites)}
    preferred = {}
    for p in canonical_faces([orient(p, sign=1) for p in parts(proposals["occupied"])]):
        for ring in [p.exterior, *p.interiors]:
            xy = list(ring.coords)
            for a, b in zip(xy, xy[1:]):
                preferred[tuple(sorted((a, b)))] = int(a > b)
    edges, parents, costs = [], [], []
    node_parents = [set() for _ in sites]
    for parent, a, b, length in proposals["parent_routes"]:
        pair = (index[a], index[b])
        node_parents[index[a]].add(parent)
        node_parents[index[b]].add(parent)
        direction = preferred[proposals["parent_edges"][parent]]
        # Distinct parents can round to the same sampled segment. Keep their
        # identities; overlap constraints select at most one of these routes.
        edges.append(pair)
        parents.append(parent)
        costs.extend(([-length, length] if direction == 0 else [length, -length]))
    for chord in proposals["chords"]:
        a, b = map(tuple, chord.coords)
        pair = (index[a], index[b])
        edges.append(pair)
        parents.append(None)
        costs.extend([chord.length, chord.length])
    return sites, edges, parents, costs, node_parents


def add_turn_constraints(rows, sites, edges, parents, incoming, outgoing, v, active):
    """Require a corner unless the chosen outgoing route continues straight.

    Binary flow selects either no routes or one incoming/outgoing pair. Under
    that invariant this is equivalent to all incompatible-pair implications,
    with one row per incoming route. Same-edge reversal is excluded elsewhere.
    """
    for left in incoming:
        terms = [(left, 1), (active, -1)]
        for right in outgoing:
            if left // 2 == right // 2:
                continue
            p, q = parents[left // 2], parents[right // 2]
            a = edges[left // 2][left % 2]
            b = edges[right // 2][1 - right % 2]
            u, w = np.subtract(sites[a], sites[v]), np.subtract(sites[b], sites[v])
            straight = (p is not None and p == q) or (
                u[0] * w[1] - u[1] * w[0] == 0 and np.dot(u, w) < 0
            )
            if straight:
                terms.append((right, -1))
        rows.add(terms, high=0)


def build_model(proposals, delta, deadline):
    sites, edges, parents, costs, node_parents = route_graph(proposals)
    count = len(edges)
    incoming, outgoing = [[] for _ in sites], [[] for _ in sites]
    for k, (a, b) in enumerate(edges):
        outgoing[a].append(2 * k)
        incoming[b].append(2 * k)
        outgoing[b].append(2 * k + 1)
        incoming[a].append(2 * k + 1)
    rows = Rows(2 * count + len(sites), deadline)
    for k in range(count):
        rows.add([(2 * k, 1), (2 * k + 1, 1)], high=1)
        if parents[k] is not None and parents[k] not in proposals["editable_parents"]:
            # The negative-cost direction is the original occupied-side winding.
            direction = 2 * k + int(costs[2 * k + 1] < 0)
            rows.add([(direction, 1)], 1, 1)
    for v in range(len(sites)):
        rows.add([(i, 1) for i in incoming[v]] + [(i, -1) for i in outgoing[v]], 0, 0)
        rows.add([(i, 1) for i in outgoing[v]], high=1)
        active = 2 * count + v
        rows.add([(active, 1)] + [(i, -1) for i in outgoing[v]], high=0)
        add_turn_constraints(
            rows, sites, edges, parents, incoming[v], outgoing[v], v, active
        )
    lines = [LineString([sites[a], sites[b]]) for a, b in edges]
    tree = STRtree(lines)
    for k, line in enumerate(lines):
        for j in tree.query(line, predicate="intersects"):
            if j <= k:
                continue
            shared = set(edges[k]) & set(edges[j])
            endpoints = MultiPoint([sites[i] for i in shared])
            if not line.intersection(lines[j]).difference(endpoints).is_empty:
                rows.add(
                    [(2 * k, 1), (2 * k + 1, 1), (2 * j, 1), (2 * j + 1, 1)], high=1
                )
    points = [Point(xy) for xy in sites]
    point_tree = STRtree(points)
    threshold = delta - delta * _SEPARATION_RELATIVE_TOLERANCE
    for v, point in enumerate(points):
        active = 2 * count + v
        for w in point_tree.query(point, predicate="dwithin", distance=threshold):
            if w > v and point.distance(points[w]) < threshold:
                rows.add([(active, 1), (2 * count + int(w), 1)], high=1)
        for k in tree.query(point, predicate="dwithin", distance=threshold):
            if v in edges[k] or parents[k] in node_parents[v]:
                continue
            a, b = [sites[i] for i in edges[k]]
            # Collinear pieces might belong to the vertex's incident canonical
            # edge after selection. Leave these cases to the independent check.
            cross = (b[0] - a[0]) * (sites[v][1] - a[1]) - (b[1] - a[1]) * (
                sites[v][0] - a[0]
            )
            if cross and point.distance(lines[k]) < threshold:
                rows.add([(active, 1), (2 * int(k), 1), (2 * int(k) + 1, 1)], high=1)
    return sites, edges, parents, np.array(costs + [0.0] * len(sites)), rows


def winding_terms(sites, edges, point):
    p = tuple(point.coords[0])
    for k, (a, b) in enumerate(edges):
        value = crossing_number(sites[a], sites[b], p)
        if value:
            yield 2 * k, value
            yield 2 * k + 1, -value


def selected_regions(sites, edges, parents, selected):
    """Coalesce original-parent runs, then form regions from selected cycles."""
    next_edge, incoming = {}, {}
    for direction in selected:
        a, b = edges[direction // 2][:: 1 if direction % 2 == 0 else -1]
        if a in next_edge or b in incoming:
            return None, None
        next_edge[a] = (b, parents[direction // 2])
        incoming[b] = parents[direction // 2]
    if set(next_edge) != set(incoming):
        return None, None
    remaining = set(next_edge)
    rings = []
    while remaining:
        start = v = min(remaining)
        coords = []
        while True:
            remaining.remove(v)
            w, parent = next_edge[v]
            if parent is None or incoming[v] != parent:
                coords.append(sites[v])
            v = w
            if v == start:
                break
        if len(coords) < 3:
            return None, None
        ring = Polygon(coords)
        if not ring.is_valid or ring.is_empty:
            return None, None
        rings.append(ring)
    output = []
    for face in polygonize([r.exterior for r in rings]):
        point = face.representative_point()
        winding = sum(
            (1 if r.exterior.is_ccw else -1) for r in rings if r.contains(point)
        )
        if winding not in (0, 1):
            return None, point
        if winding == 1:
            output.append(face)
    return parts(unary_union(output)), None


def select_routes(
    raw,
    *,
    delta=0.5,
    epsilon=0.25,
    max_seconds=5.0,
    max_rounds=16,
    focus_unresolved=False,
):
    _validate_scale(max_seconds, "max_seconds", positive=True)
    if type(max_rounds) is not int or max_rounds < 1:
        raise ValueError("max_rounds must be a positive integer")
    raw = list(raw)
    started = time.perf_counter()
    proposals, report = boundary_candidates(
        raw, delta=delta, epsilon=epsilon, focus_unresolved=focus_unresolved
    )
    report.update(outcome="unresolved", rounds=0, focused=focus_unresolved)
    if proposals is None:
        return None, report
    if "chords" not in proposals:
        output = parts(proposals["occupied"])
        contract = check_cleaning_contract(raw, output, delta=delta, epsilon=epsilon)
        report["contract"] = contract
        if contract["status"] == "pass":
            return output, dict(report, outcome="conforming", reason="unchanged")
        return None, dict(report, reason="no_candidate_routes")
    deadline = started + max_seconds
    try:
        sites, edges, parents, costs, rows = build_model(proposals, delta, deadline)
        report.update(variables=len(costs), model_seconds=time.perf_counter() - started)
        protected, allowed = proposals["budget"]._bounds[1]
        outside = proposals["occupied"].convex_hull.difference(allowed)
        for geometry, label in ((protected, 1), (outside, 0)):
            for p in parts(geometry):
                rows.add(
                    winding_terms(sites, edges, p.representative_point()), label, label
                )
        for iteration in range(max_rounds):
            seconds = deadline - time.perf_counter()
            if seconds <= 0:
                raise WorkLimit("time_limit")
            result = milp(
                costs,
                integrality=np.ones(len(costs)),
                bounds=Bounds(0, 1),
                constraints=rows.constraint(),
                options={"time_limit": seconds},
            )
            report.update(
                rounds=iteration + 1,
                constraints=len(rows.low),
                solver_status=int(result.status),
            )
            if result.x is None or result.status not in (0, 1):
                reason = (
                    "model_infeasible" if result.status == 2 else "solver_unresolved"
                )
                return None, dict(report, reason=reason)
            if (
                not np.isfinite(result.x).all()
                or np.max(np.abs(result.x - np.rint(result.x))) > 1e-6
            ):
                return None, dict(report, reason="nonintegral_incumbent")
            selected = set(np.flatnonzero(result.x[: 2 * len(edges)] > 0.5))
            candidate, witness = selected_regions(sites, edges, parents, selected)
            if candidate is not None:
                contract = check_cleaning_contract(
                    raw, candidate, delta=delta, epsilon=epsilon
                )
                report["last_contract"] = contract
                if contract["status"] == "pass":
                    return candidate, dict(
                        report,
                        outcome="conforming",
                        reason="contract_pass",
                        selected_segments=len(selected),
                        contract=contract,
                    )
                union = unary_union(candidate)
                for bad, label in (
                    (protected.difference(union), 1),
                    (union.difference(allowed), 0),
                ):
                    if bad.area > 1e-10:
                        rows.add(
                            winding_terms(sites, edges, bad.representative_point()),
                            label,
                            label,
                        )
            elif witness is not None:
                rows.add(winding_terms(sites, edges, witness), 0, 1)
            # Reject this exact directed assignment. This guarantees finite
            # progress even where the sufficient linear constraints are weak.
            rows.add(
                ((i, 1 if i in selected else -1) for i in range(2 * len(edges))),
                high=len(selected) - 1,
            )
        return None, dict(report, reason="round_limit")
    except WorkLimit as error:
        return None, dict(report, reason=str(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--focused",
        action="store_true",
        help="Restrict proposals to unresolved input features (research comparison)",
    )
    parser.add_argument(
        "--run", type=Path, help="Evaluate the fixed saved raw city groups"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    reports, plotted = {}, {}
    for name, (raw, _) in topology_examples().items():
        start = time.perf_counter()
        output, report = select_routes(raw, focus_unresolved=args.focused)
        report["seconds"] = time.perf_counter() - start
        reports[name] = report
        plotted[name] = raw, output, []
        (args.output / "routes.json").write_text(json.dumps(reports, indent=2) + "\n")
        geometry_path = args.output / f"{name}.geojson"
        if output is not None:
            artifact = {
                "type": "FeatureCollection",
                "features": [
                    {"type": "Feature", "properties": {}, "geometry": mapping(p)}
                    for p in output
                ],
            }
            geometry_path.write_text(json.dumps(artifact) + "\n")
        else:
            # A repeated comparison must not leave an older successful output
            # beside a report that now says this example is unresolved.
            geometry_path.unlink(missing_ok=True)
        print(
            name,
            report["outcome"],
            report["reason"],
            report["rounds"],
            round(report["seconds"], 2),
            flush=True,
        )
    plot_topology(
        plotted,
        args.output / "route-examples.png",
        input_label="raw input",
        title="Automatic boundary-route search · δ = 0.5 m, ε = 0.25 m\n"
        "Only selected cycles are converted to regions; all outputs independently checked",
    )
    if args.run:
        rows = []
        for record, delta, epsilon, groups in saved_city_groups(args.run):
            metrics = []
            for group in groups:
                start = time.perf_counter()
                _, report = select_routes(
                    group, delta=delta, epsilon=epsilon, focus_unresolved=args.focused
                )
                report["seconds"] = time.perf_counter() - start
                metrics.append(report)
            row = {
                "case": record["case"],
                "delta": delta,
                "epsilon": epsilon,
                "groups": metrics,
                "seconds": sum(g["seconds"] for g in metrics),
            }
            rows.append(row)
            (args.output / "route-cities.json").write_text(
                json.dumps(rows, indent=2) + "\n"
            )
            print(
                record["case"]["city"],
                sum(g["outcome"] == "conforming" for g in metrics),
                "/",
                len(groups),
                "conforming groups",
                round(row["seconds"], 2),
                "seconds",
                flush=True,
            )


if __name__ == "__main__":
    main()
