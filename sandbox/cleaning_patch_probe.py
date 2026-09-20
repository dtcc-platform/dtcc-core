"""Constructive occupancy research; never called by production cleaning.

Local edits assign patch occupancy or separate conflicting features. Original
fidelity and strict integer defect descent admit edits. Internal source walls are outside
this experiment: merging permission is mandatory. No repair fallback.

  .venv/bin/python sandbox/cleaning_patch_probe.py --output /tmp/patch-probe
  .venv/bin/python sandbox/cleaning_patch_probe.py --output /tmp/patch-cities \
      --run benchmarks/runs/2026-09-17_100049_quick
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections import Counter
import sys
import time

import numpy as np
from scipy.optimize import LinearConstraint, minimize
from shapely import STRtree, get_coordinates, points as make_points
from shapely.errors import GEOSException
from shapely.geometry import LineString, MultiPoint, Point, Polygon, box, mapping
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning.contract import (
    FidelityBudget,
    _SEPARATION_RELATIVE_TOLERANCE,
    _canonical_graph,
    _interpret_input,
    _validate_scale,
    admissibility,
    check_cleaning_contract,
)
from sandbox.cleaning_graph_prototype import geometry_key, polygon_parts


def defect_tuple(state):
    return (
        state["nonmanifold_boundary_vertices"],
        state["subscale_pairs"],
        state["canonical_vertex_count"],
    )


def motion_reference(budget):
    """Original material boundaries and the allowed boundary-motion band."""
    core, envelope = budget._bounds[1]
    coordinates = np.unique(
        np.concatenate([get_coordinates(core), get_coordinates(envelope)]), axis=0
    )
    return coordinates, STRtree(make_points(coordinates)), envelope.difference(core)


def motion_interval(vertex, normal, band, radius):
    """Choose the nearest connected interval of the ORIGINAL band on this ray."""
    vertex = np.asarray(vertex)
    ray = LineString([vertex - radius * normal, vertex + radius * normal])
    pending = [band.intersection(ray)]
    intervals = []
    while pending:
        part = pending.pop()
        if part.is_empty:
            continue
        if hasattr(part, "geoms"):
            pending.extend(part.geoms)
            continue
        parameters = (np.asarray(part.coords)[:, :2] - vertex) @ normal
        lo, hi = max(-radius, parameters.min()), min(radius, parameters.max())
        if lo <= hi:
            intervals.append((lo, hi))
    return (
        min(intervals, key=lambda bounds: (max(bounds[0], -bounds[1], 0), bounds))
        if intervals
        else None
    )


def coupled_move(vertex, opposing, edges, incident, reference, delta, epsilon):
    """A bounded convex proposal in the original boundary-motion band.

    Move a conflicting pair and its graph neighbors along one normal. The
    degree-two restriction gives at most seven variables and fourteen affected
    edges. At most 512 inequalities and 128 solver iterations are allowed.
    Band sections and neighboring clearances constrain affected features.
    Numerical construction is not a certificate: the caller checks original
    fidelity and strict descent.
    """
    point = Point(vertex)
    feature = Point(opposing[0]) if len(opposing) == 1 else LineString(opposing)
    nearest = (
        feature if len(opposing) == 1 else feature.interpolate(feature.project(point))
    )
    normal = np.asarray(vertex) - nearest.coords[0]
    length = np.linalg.norm(normal)
    if length == 0 or length >= delta:
        return None, "no_direction"
    normal /= length
    seeds = {vertex, *opposing}
    if any(len(incident[v]) != 2 for v in seeds):
        return None, "junction"
    movable = sorted(seeds | {p for v in seeds for j in incident[v] for p in edges[j]})
    if any(len(incident[v]) != 2 for v in movable):
        return None, "junction"
    ids = {xy: i for i, xy in enumerate(movable)}
    affected = sorted({j for xy in movable for j in incident[xy]})
    rows, bounds = [], []
    for xy in opposing:
        row = np.zeros(len(ids))
        row[ids[vertex]], row[ids[xy]] = -1, 1
        rows.append(row)
        bounds.append(np.dot(np.asarray(vertex) - xy, normal) - delta)
    coordinates, tree, band = reference
    radius = 2 * epsilon + delta
    # Every section lies within radius of an affected edge. Include fixed
    # endpoints too, and round the query box outwards to retain boundary rays.
    # This restricts work, not the original fidelity band or allowed motion.
    support = np.asarray([p for j in affected for p in edges[j]])
    lower = np.nextafter(support.min(axis=0) - radius, -np.inf)
    upper = np.nextafter(support.max(axis=0) + radius, np.inf)
    band = band.intersection(box(*lower, *upper))
    movement_bounds = [motion_interval(xy, normal, band, radius) for xy in movable]
    if any(bound is None for bound in movement_bounds):
        return None, "no_motion_interval"
    tangent = np.array([-normal[1], normal[0]])
    for j in affected:
        aa, bb = edges[j]
        a, b = np.asarray(aa), np.asarray(bb)
        edge = b - a
        end = np.dot(edge, tangent)
        if end == 0:
            continue  # This segment moves along its line and sweeps no area.
        indices = tree.query(LineString(edges[j]), predicate="dwithin", distance=radius)
        offsets = coordinates[indices] - a
        along = offsets @ tangent
        indices = indices[(along >= min(0, end)) & (along <= max(0, end))]
        # Between projections of original boundary vertices, the band's
        # interval endpoints are linear. Constrain the moving edge at these
        # breakpoints and its endpoints, using the same connected corridor.
        fractions = sorted(
            set(float(np.dot(coordinates[k] - a, tangent) / end) for k in indices)
        )
        for fraction in fractions:
            location = a + fraction * edge
            interval = motion_interval(location, normal, band, radius)
            if interval is None:
                return None, "no_motion_interval"
            if len(rows) + 2 > 512:
                return None, "constraint_limit"
            row = np.zeros(len(ids))
            if aa in ids:
                row[ids[aa]] = 1 - fraction
            if bb in ids:
                row[ids[bb]] = fraction
            rows.extend([row, -row])
            bounds.extend([interval[1], -interval[0]])
    for row, bound in clearance_constraints(edges, ids, normal, radius, delta):
        if len(rows) >= 512:
            return None, "constraint_limit"
        rows.append(row)
        bounds.append(bound)
    result = minimize(
        lambda u: 0.5 * np.dot(u, u),
        np.zeros(len(ids)),
        jac=lambda u: u,
        constraints=LinearConstraint(np.asarray(rows), -np.inf, bounds),
        bounds=movement_bounds,
        method="SLSQP",
        options={"maxiter": 128, "ftol": 1e-12},
    )
    if not result.success or not np.isfinite(result.x).all():
        return None, (
            "iteration_limit"
            if getattr(result, "status", None) == 9
            else "solver_unresolved"
        )
    return {
        xy: tuple(np.asarray(xy) + result.x[i] * normal)
        for xy, i in ids.items()
        if result.x[i] != 0
    }, "proposed"


def clearance_constraints(edges, ids, normal, radius, delta):
    """Sufficient linear bounds preserving affected neighboring clearances.

    Each vertex and every point of an affected edge moves at most radius.
    Only pairs initially within delta + 2*radius can develop a new conflict.
    Fixed supporting directions preserve min(current distance, delta); the
    initiating pair is separately required to reach delta. These conservative
    constraints propose geometry; they do not replace the independent checker.
    """
    vertices = sorted({xy for edge in edges for xy in edge})
    points = [Point(xy) for xy in vertices]
    segments = [LineString(edge) for edge in edges]
    point_tree, edge_tree = STRtree(points), STRtree(segments)
    reach = np.nextafter(delta + 2 * radius, np.inf)

    def pairs():
        for xy in ids:
            point = Point(xy)
            for j in sorted(
                point_tree.query(point, predicate="dwithin", distance=reach)
            ):
                other = vertices[j]
                if other != xy and (other not in ids or xy < other):
                    yield xy, [other]
            for j in sorted(
                edge_tree.query(point, predicate="dwithin", distance=reach)
            ):
                if xy not in edges[j]:
                    yield xy, edges[j]
        # Fixed vertices can approach the interior of a moving edge too.
        for j, edge in enumerate(edges):
            if not any(xy in ids for xy in edge):
                continue
            for i in sorted(
                point_tree.query(segments[j], predicate="dwithin", distance=reach)
            ):
                xy = vertices[i]
                if xy not in ids and xy not in edge:
                    yield xy, edge

    # Parallel displacements make every inequality a bound on u_b - u_a.
    # Keep its strongest bound once: at most k*(k-1) + 2*k + 1 rows for k
    # movable vertices, including a possible inconsistent constant constraint.
    differences = {}
    for xy, opposing in pairs():
        point = Point(xy)
        feature = Point(opposing[0]) if len(opposing) == 1 else LineString(opposing)
        nearest = (
            feature
            if len(opposing) == 1
            else feature.interpolate(feature.project(point))
        )
        direction = np.asarray(xy) - nearest.coords[0]
        length = np.linalg.norm(direction)
        if length == 0:
            continue  # The existing zero clearance cannot decrease.
        direction /= length
        coefficient = np.dot(normal, direction)
        for other in opposing:
            bound = np.dot(np.asarray(xy) - other, direction) - min(length, delta)
            a, b = ids.get(xy), ids.get(other)
            if coefficient == 0 or a == b:
                if bound >= 0:
                    continue
                key = (None, None)
            else:
                key = (a, b) if coefficient > 0 else (b, a)
                bound /= abs(coefficient)
            differences[key] = min(differences.get(key, np.inf), bound)
    for (a, b), bound in differences.items():
        row = np.zeros(len(ids))
        if a is not None:
            row[a] -= 1
        if b is not None:
            row[b] += 1
        yield row, bound


def move_boundary(occupied, move, vertices):
    """Move canonical edges without pinning suppressed collinear ring points."""
    polygons = []
    for polygon in polygon_parts(occupied):
        rings = [
            [
                move.get(tuple(xy[:2]), tuple(xy[:2]))
                for xy in ring.coords
                if tuple(xy[:2]) in vertices
            ]
            for ring in [polygon.exterior, *polygon.interiors]
        ]
        moved = Polygon(rings[0], rings[1:])
        if not moved.is_valid or moved.is_empty:
            return None
        polygons.append(moved)
    return unary_union(polygons)


def patch_batches(occupied, delta, epsilon):
    """Occupancy and coordinate choices at each unresolved graph support."""
    vertices, edges = _canonical_graph([occupied.boundary])
    segments = [LineString(edge) for edge in edges]
    points = [Point(xy) for xy in vertices]
    point_tree, edge_tree = STRtree(points), STRtree(segments)
    incident = {xy: [] for xy in vertices}
    for i, edge in enumerate(edges):
        for xy in edge:
            incident[xy].append(i)
    threshold = delta * (1 - _SEPARATION_RELATIVE_TOLERANCE)
    seen = set()
    radii = [(2 * epsilon + delta) * k / 8 for k in range(1, 9)]
    # Canonical vertices and sorted edge indices fix the proposal order.
    # Prepare only supports reached by the consumer before its next edit.
    for i, (xy, point) in enumerate(zip(vertices, points)):
        near_points = [
            int(j)
            for j in point_tree.query(point, predicate="dwithin", distance=threshold)
            if j != i and point.distance(points[j]) < threshold
        ]
        near_edges = [
            int(j)
            for j in edge_tree.query(point, predicate="dwithin", distance=threshold)
            if xy not in edges[j] and point.distance(segments[j]) < threshold
        ]
        star = tuple(incident[xy])
        supports = set()
        if len(star) != 2 or near_points or near_edges:
            supports.add(star)
        for j in near_edges:
            for k in star:
                supports.add(tuple(sorted((j, k))))
        for indices in sorted(supports):
            found = {}
            for radius in [None, *radii]:
                coordinates = []
                for j in indices:
                    segment = segments[j]
                    if radius is None:
                        coordinates.extend(edges[j])
                    else:
                        distance = segment.project(point)
                        coordinates.extend(
                            segment.interpolate(d).coords[0]
                            for d in (
                                max(0, distance - radius),
                                min(segment.length, distance + radius),
                            )
                        )
                patch = MultiPoint(coordinates).convex_hull
                if patch.geom_type == "Polygon" and patch.area > 0:
                    key = geometry_key(patch)
                    if key not in seen:
                        found[key] = patch
                        seen.add(key)
            translations = []
            if indices == star:
                for j in near_points:
                    translations.append((xy, [vertices[j]], edges, incident))
                for j in near_edges:
                    translations.append((xy, edges[j], edges, incident))
            if found or translations:
                yield sorted(
                    found.values(), key=lambda p: (p.area, geometry_key(p))
                ), translations


def local_edits(patches, translations):
    # Both freedoms belong to one local choice; neither is a fallback stage.
    for patch in patches:
        yield "clear", patch
        yield "fill", patch
    for move in translations:
        yield "move", move


def construct(raw, *, delta=0.5, epsilon=None, merge_buildings=False, max_checks=2000):
    """Return only independently conforming output, or explicit unresolved.

    Each invocation is one interaction group. A cap bounds candidate admissions,
    not elapsed wall time. The original budget is never reset after an edit.
    """
    _validate_scale(delta, "delta", positive=True)
    epsilon = delta / 2 if epsilon is None else epsilon
    _validate_scale(epsilon, "epsilon", positive=False)
    if merge_buildings is not True:
        raise ValueError("Occupancy construction requires merge_buildings=True")
    if type(max_checks) is not int or max_checks <= 0:
        raise ValueError("max_checks must be a positive integer")
    raw = list(raw)
    occupied, interpretation = _interpret_input(raw)
    original = occupied
    budget = FidelityBudget(raw, epsilon)
    reference = None
    motion_outcomes = Counter()
    state = admissibility(polygon_parts(occupied), delta)
    history = [defect_tuple(state)]
    checks = rejected_fidelity = rejected_progress = numerical = 0
    reason = "local_minimum"
    while not state["resolved"]:
        best = None
        for batch, translations in patch_batches(occupied, delta, epsilon):
            for operation, geometry in local_edits(batch, translations):
                if checks >= max_checks:
                    reason = "work_limit"
                    break
                checks += 1
                try:
                    if operation == "move":
                        vertex, opposing, edges, incident = geometry
                        if reference is None:
                            reference = motion_reference(budget)
                        move, status = coupled_move(
                            vertex, opposing, edges, incident, reference, delta, epsilon
                        )
                        motion_outcomes[status] += 1
                        candidate = (
                            move_boundary(occupied, move, incident)
                            if move is not None
                            else None
                        )
                    elif operation == "fill":
                        candidate = occupied.union(geometry)
                    else:
                        candidate = occupied.difference(geometry)
                    if candidate is None or candidate.equals(occupied):
                        continue
                    if not budget.accepts_union(candidate):
                        rejected_fidelity += 1
                        continue
                    candidate_state = admissibility(polygon_parts(candidate), delta)
                    if "canonical_vertex_count" not in candidate_state:
                        numerical += 1
                        continue
                    if defect_tuple(candidate_state) >= history[-1]:
                        rejected_progress += 1
                        continue
                except GEOSException:
                    numerical += 1
                    continue
                rank = (
                    defect_tuple(candidate_state),
                    original.symmetric_difference(candidate).area,
                )
                if best is None or rank < best[0]:
                    best = rank, candidate, candidate_state
            # Finish one local support before committing, not a global search
            # over every defect in a potentially large interaction group.
            if best is not None or reason == "work_limit":
                break
        if best is None:
            break
        _, occupied, state = best
        history.append(defect_tuple(state))
        if reason == "work_limit":
            break
    output = polygon_parts(occupied)
    contract = check_cleaning_contract(raw, output, delta=delta, epsilon=epsilon)
    passed = contract["status"] == "pass"
    return output if passed else None, {
        "outcome": "conforming" if passed else "unresolved",
        "reason": "contract_pass" if passed else reason,
        "interpretation": interpretation,
        "checks": checks,
        "edits": len(history) - 1,
        "progress": history,
        "rejected_fidelity": rejected_fidelity,
        "rejected_progress": rejected_progress,
        "numerical_rejections": numerical,
        "motion_outcomes": dict(motion_outcomes),
        "contract": contract,
    }


def write_geometry(destination, output):
    if output is None:
        destination.unlink(missing_ok=True)
        return
    destination.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {"type": "Feature", "properties": {}, "geometry": mapping(p)}
                    for p in output
                ],
            }
        )
        + "\n"
    )


def main():
    from sandbox.cleaning_pipeline_probe import (
        plot_topology,
        saved_city_groups,
        topology_examples,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    examples, plotted = {}, {}
    for name, (raw, _) in topology_examples().items():
        started = time.perf_counter()
        output, report = construct(raw, merge_buildings=True)
        report["seconds"] = time.perf_counter() - started
        examples[name] = report
        plotted[name] = raw, output, []
        destination = args.output / f"{name}.geojson"
        write_geometry(destination, output)
        print(name, report["outcome"], report["reason"], report["checks"], flush=True)
    (args.output / "examples.json").write_text(json.dumps(examples, indent=2) + "\n")
    plot_topology(
        plotted,
        args.output / "examples.png",
        title="Local occupancy and coordinate construction",
        input_label="raw input",
    )
    if args.run:
        rows = []
        for record, delta, epsilon, groups in saved_city_groups(args.run):
            reports = []
            combined = []
            complete = True
            for i, raw in enumerate(groups):
                started = time.perf_counter()
                output, report = construct(
                    raw, delta=delta, epsilon=epsilon, merge_buildings=True
                )
                report["seconds"] = time.perf_counter() - started
                reports.append(report)
                complete &= output is not None
                if output is not None:
                    combined.extend(output)
                if report["seconds"] > 5:
                    print(
                        record["case"]["city"],
                        "group",
                        i,
                        report["reason"],
                        f"{report['seconds']:.2f}s",
                        flush=True,
                    )
            row = {
                "case": record["case"],
                "groups": reports,
                "seconds": sum(r["seconds"] for r in reports),
                "city_contract": (
                    check_cleaning_contract(
                        [p for group in groups for p in group],
                        combined,
                        delta=delta,
                        epsilon=epsilon,
                    )
                    if complete
                    else None
                ),
            }
            city = record["case"]["city"]
            if Path(city).name != city:
                raise ValueError("City name must be a filename component")
            write_geometry(
                args.output / f"{city}.geojson",
                (
                    combined
                    if complete and row["city_contract"]["status"] == "pass"
                    else None
                ),
            )
            rows.append(row)
            (args.output / "cities.json").write_text(json.dumps(rows, indent=2) + "\n")
            print(
                record["case"]["city"],
                sum(r["outcome"] == "conforming" for r in reports),
                "/",
                len(groups),
                f"{row['seconds']:.2f}s",
                flush=True,
            )


if __name__ == "__main__":
    main()
