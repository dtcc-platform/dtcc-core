"""Bounded research prototype; never called by production cleaning.

Only degree-two vertex removal is allowed. Shared faces change together; every
accepted edit stays inside the original fidelity band. Input normalization may
dissolve existing internal walls when merge_buildings permits it. Refinement has
no movement, merging of separated components, hole deletion or repair fallback.

  .venv/bin/python sandbox/cleaning_graph_prototype.py --output /tmp/graph-probe
  .venv/bin/python sandbox/cleaning_graph_prototype.py --output /tmp/graph-cities \
      --run benchmarks/runs/2026-09-17_100049_quick
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from shapely import STRtree, get_coordinates, make_valid
from shapely.geometry import LineString, MultiPoint, MultiPolygon, Point, Polygon
from shapely.ops import polygonize, unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning.contract import (
    FidelityBudget,
    _canonical_graph,
    _validate_scale,
    admissibility,
    check_cleaning_contract,
)


def polygon_parts(geometry):
    if geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    return [p for g in getattr(geometry, "geoms", ()) for p in polygon_parts(g)]


def geometry_key(polygon):
    return polygon.normalize().wkb


def interaction_groups(polygons, distance):
    """Connected components under the sufficient 2 epsilon + delta distance."""
    tree = STRtree(polygons)
    pending = set(range(len(polygons)))
    groups = []
    while pending:
        stack = [min(pending, key=lambda i: geometry_key(polygons[i]))]
        pending.remove(stack[0])
        group = []
        while stack:
            i = stack.pop()
            group.append(i)
            neighbors = (
                set(
                    map(
                        int,
                        tree.query(polygons[i], predicate="dwithin", distance=distance),
                    )
                )
                & pending
            )
            pending.difference_update(neighbors)
            stack.extend(sorted(neighbors))
        groups.append(sorted(group, key=lambda i: geometry_key(polygons[i])))
    return groups


def subdivision(polygons, sources, merge_buildings):
    """Node once, retaining each occupied face and its original source labels."""
    tree = STRtree(polygons)
    faces, labels = [], []
    if merge_buildings:
        faces = sorted(polygon_parts(unary_union(polygons)), key=geometry_key)
        labels = [
            sorted(
                {
                    s
                    for i in tree.query(face, predicate="intersects")
                    if face.intersection(polygons[int(i)]).area > 0
                    for s in sources[int(i)]
                }
            )
            for face in faces
        ]
        return canonical_faces(faces), labels
    linework = unary_union([p.boundary for p in polygons])
    for face in sorted(polygonize(linework), key=geometry_key):
        owners = tree.query(face.representative_point(), predicate="within")
        if len(owners):
            faces.append(face)
            labels.append(sorted({s for i in owners for s in sources[int(i)]}))
    return canonical_faces(faces), labels


def canonical_faces(faces):
    vertices, _ = _canonical_graph([p.boundary for p in faces])
    keep = set(vertices)
    result = []
    for p in faces:
        rings = [
            [tuple(c[:2]) for c in r.coords[:-1] if tuple(c[:2]) in keep]
            for r in [p.exterior, *p.interiors]
        ]
        result.append(Polygon(rings[0], rings[1:]))
    return result


def graph(faces, delta):
    """Discover editable vertices; the independent checker certifies output."""
    owners, neighbors, edge_set = {}, {}, set()
    for i, p in enumerate(faces):
        for ring in [p.exterior, *p.interiors]:
            coords = [tuple(c[:2]) for c in ring.coords[:-1]]
            for a, b in zip(coords, coords[1:] + coords[:1]):
                owners.setdefault(a, set()).add(i)
                neighbors.setdefault(a, set()).add(b)
                neighbors.setdefault(b, set()).add(a)
                edge_set.add(tuple(sorted((a, b))))
    vertices, edges = sorted(neighbors), sorted(edge_set)
    points = [Point(v) for v in vertices]
    segments = [LineString(e) for e in edges]
    point_tree, edge_tree = STRtree(points), STRtree(segments)
    affected = set()
    for i, (v, point) in enumerate(zip(vertices, points)):
        for j in point_tree.query(point, predicate="dwithin", distance=delta):
            if j > i and point.distance(points[j]) < delta:
                affected.update((v, vertices[j]))
        for j in edge_tree.query(point, predicate="dwithin", distance=delta):
            if v not in edges[j] and point.distance(segments[j]) < delta:
                affected.update((v, *edges[j]))
    candidates = sorted(v for v in affected if len(neighbors[v]) == 2)
    return owners, neighbors, edges, segments, point_tree, edge_tree, candidates


def shortcut(faces, vertex, context):
    owners, neighbors, edges, segments, point_tree, edge_tree, _ = context
    a, b = sorted(neighbors[vertex])
    chord = LineString([a, b])
    removed = {tuple(sorted((a, vertex))), tuple(sorted((vertex, b)))}
    endpoints = MultiPoint([a, b])
    for j in edge_tree.query(chord, predicate="intersects"):
        if edges[j] not in removed:
            intersection = chord.intersection(segments[j])
            if not intersection.difference(endpoints).is_empty:
                return None
    # A chord can sweep over a whole hole/component without crossing its edges.
    # Forbid swallowing any other graph vertex in the swept triangle.
    triangle = Polygon([a, vertex, b])
    if triangle.area and len(point_tree.query(triangle, predicate="contains")):
        return None
    candidate = list(faces)
    for i in owners[vertex]:
        p = faces[i]
        rings = [
            [tuple(c[:2]) for c in r.coords[:-1] if tuple(c[:2]) != vertex]
            for r in [p.exterior, *p.interiors]
        ]
        if any(len(r) < 3 for r in rings):
            return None
        changed = Polygon(rings[0], rings[1:])
        if changed.is_empty or not changed.is_valid:
            return None
        candidate[i] = changed
    return candidate


def simplify_group(original, sources, delta, epsilon, max_checks, merge_buildings):
    budget = FidelityBudget(original, epsilon)
    occupied = unary_union(original)
    faces, labels = subdivision(original, sources, merge_buildings)
    initial_vertices = admissibility(faces, delta)["canonical_vertex_count"]
    checks = accepted = 0
    reason = "no_admissible_shortcut"
    while True:
        state = admissibility(faces, delta)
        if not budget.accepts_union(unary_union(faces)):
            reason = "initial_fidelity_unresolved"
            break
        if state["resolved"]:
            reason = "resolved"
            break
        if not state["topology_ok"]:
            reason = "requires_topology_change"
            break
        context = graph(faces, delta)
        best = None
        for vertex in context[-1]:
            if checks >= max_checks:
                reason = "work_limit"
                break
            checks += 1
            candidate = shortcut(faces, vertex, context)
            if candidate is None:
                continue
            # Canonicalization is part of the proposal: admit exactly the
            # coordinates that would become the next state.
            candidate = canonical_faces(candidate)
            candidate_union = unary_union(candidate)
            if not budget.accepts_union(candidate_union):
                continue
            score = (occupied.symmetric_difference(candidate_union).area, vertex)
            if best is None or score < best[0]:
                best = score, candidate
        if reason == "work_limit" or best is None:
            break
        faces = best[1]
        accepted += 1
    final = admissibility(faces, delta)
    assert final["canonical_vertex_count"] <= initial_vertices - accepted
    return (
        faces,
        labels,
        {
            "reason": reason,
            "accepted_shortcuts": accepted,
            "candidate_checks": checks,
            "initial_vertices": initial_vertices,
            "final_vertices": final["canonical_vertex_count"],
            "remaining_subscale_pairs": final["subscale_pairs"],
        },
    )


def simplify_coverage(
    raw, *, delta=0.5, epsilon=None, merge_buildings=False, max_candidate_checks=2000
):
    """Return polygons, source maps and a conforming/unresolved research report.

    The work limit applies separately to each independent group. Hitting it is
    unresolved, never evidence of infeasibility. Source indices address raw.
    """
    _validate_scale(delta, "delta", positive=True)
    epsilon = delta / 2 if epsilon is None else epsilon
    _validate_scale(epsilon, "epsilon", positive=False)
    if type(max_candidate_checks) is not int or max_candidate_checks <= 0:
        raise ValueError("max_candidate_checks must be a positive integer")
    if type(merge_buildings) is not bool:
        raise ValueError("merge_buildings must be a boolean")
    raw = list(raw)
    atoms, sources = [], []
    for i, p in enumerate(raw):
        if not isinstance(p, (Polygon, MultiPolygon)):
            raise ValueError(
                "Raw footprints must be Polygon or MultiPolygon geometries"
            )
        if not np.isfinite(get_coordinates(p)).all():
            raise ValueError("Raw footprints must have finite coordinates")
        for part in polygon_parts(make_valid(p)):
            atoms.append(part)
            sources.append([i])
    output, labels, groups = [], [], []
    # Conservative margin belongs to grouping, not to physical fidelity.
    for indices in interaction_groups(atoms, 2 * epsilon + delta + delta * 1e-6):
        faces, face_sources, metrics = simplify_group(
            [atoms[i] for i in indices],
            [sources[i] for i in indices],
            delta,
            epsilon,
            max_candidate_checks,
            merge_buildings,
        )
        output.extend(faces)
        labels.extend(face_sources)
        groups.append(metrics)
    report = check_cleaning_contract(raw, output, delta=delta, epsilon=epsilon)
    report.update(
        outcome="conforming" if report["status"] == "pass" else "unresolved",
        groups=groups,
        merge_buildings=merge_buildings,
        accepted_shortcuts=sum(g["accepted_shortcuts"] for g in groups),
        candidate_checks=sum(g["candidate_checks"] for g in groups),
        max_candidate_checks_per_group=max_candidate_checks,
    )
    return output, labels, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    from sandbox.cleaning_contract_probe import examples

    cases = [(name, raw, True, 0.5) for name, raw in examples().items()]
    cases.append(
        (
            "shared_narrow_parts",
            [
                Polygon([(0, 0), (0.3, 0), (0.3, 10), (0, 10)]),
                Polygon([(0.3, 0), (0.6, 0), (0.6, 10), (0.3, 10)]),
            ],
            True,
            0.5,
        )
    )
    if args.run:
        from benchmarks.benchmark_phases import load_cleaning
        from dtcc_core.model import GeometryType

        cases = []
        results = json.loads((args.run / "results.json").read_text())["results"]
        for r in results:
            if r["dataset"] != "city_flat_mesh":
                continue
            artifact = (args.run / r["artifacts"]["cleaning_input"]["path"]).resolve()
            if not artifact.is_relative_to(args.run.resolve()):
                raise ValueError("Saved cleaning artifact must be inside its run")
            city, _ = load_cleaning(artifact, r)
            raw = [
                b.flatten_geometry(GeometryType.LOD0).to_polygon(simplify=0)
                for b in city.buildings
            ]
            cases.append(
                (
                    r["case"]["city"],
                    raw,
                    r["parameters"]["merge_buildings"],
                    r["parameters"]["min_building_detail"],
                )
            )
    rows = []
    for name, raw, merge_buildings, delta in cases:
        if Path(name).name != name:
            raise ValueError("Case name must be a filename component")
        started = time.perf_counter()
        polygons, sources, report = simplify_coverage(
            raw, delta=delta, merge_buildings=merge_buildings
        )
        row = {"case": name, "seconds": time.perf_counter() - started, "report": report}
        rows.append(row)
        from shapely.geometry import mapping

        artifact = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "properties": {"source_indices": ids},
                    "geometry": mapping(p),
                }
                for p, ids in zip(polygons, sources)
            ],
        }
        (args.output / f"{name}.geojson").write_text(json.dumps(artifact) + "\n")
        (args.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
        print(
            name,
            report["outcome"],
            f"{row['seconds']:.2f}s",
            "edits",
            report["accepted_shortcuts"],
            "checks",
            report["candidate_checks"],
            flush=True,
        )


if __name__ == "__main__":
    main()
