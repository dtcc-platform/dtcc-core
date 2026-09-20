"""Counterexamples and finite occupancy search, not a replacement cleaner.

Run from dtcc-core:
  .venv/bin/python sandbox/cleaning_pipeline_probe.py --output /tmp/pipeline-probe

Uses the independent contract checker. The Lund case replays an intermediate
stage snapshot; it does not claim to reproduce the entire raw-to-mesh workflow.
Explicit partitions demonstrate the finite search. Automatic local-chord
partitions are measured separately; required internal source walls are not
implemented. Add --run <saved benchmark run> for fixed raw city measurements.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from shapely import STRtree, get_coordinates
from shapely.geometry import LineString, Point, Polygon, box, shape
from shapely.ops import polygonize, unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning import footprints
from dtcc_core.builder.cleaning.contract import (
    FidelityBudget,
    _SEPARATION_RELATIVE_TOLERANCE,
    _canonical_graph,
    _interpret_input,
    _validate_scale,
    admissibility,
    check_cleaning_contract,
    fidelity,
)


def parts(geometry):
    if geometry.is_empty:
        return []
    return [geometry] if isinstance(geometry, Polygon) else list(geometry.geoms)


def counterexamples():
    square = box(0, 0, 10, 10)
    # Euclidean-style opening then closing. The smaller radius leaves numerical
    # headroom inside epsilon; it is not an alternative production parameter.
    rounded = square
    for radius in (-0.2, 0.2, 0.2, -0.2):
        rounded = rounded.buffer(radius, quad_segs=32, join_style="round")

    # Integer coordinates give separated vertices, but the third vertex is
    # only 1/sqrt(101) from the opposite edge.
    lattice_triangle = Polygon([(0, 0), (10, 1), (1, 0)])

    raw = [box(0, 0, 3, 3), box(3.6, 0, 6.6, 3)]
    occupied = unary_union(raw)
    merged = box(0, 0, 6.6, 3)
    projected = merged.intersection(occupied.buffer(0.2, quad_segs=32)).union(
        occupied.buffer(-0.2, quad_segs=32)
    )
    return {
        "offsets_do_not_resolve_polygon_graph": ([square], parts(rounded)),
        "grid_does_not_separate_vertex_from_edge": (
            [lattice_triangle],
            [lattice_triangle],
        ),
        "fidelity_projection_can_break_resolution": (raw, parts(projected)),
    }


def carrier_labels(raw, cells, *, delta=0.5, epsilon=0.25):
    """Classify a supplied partition using the oracle's original fidelity band."""
    raw, cells = list(raw), list(cells)
    _validate_scale(delta, "delta", positive=True)
    if any(
        not isinstance(p, Polygon)
        or p.is_empty
        or not p.is_valid
        or not np.isfinite(get_coordinates(p)).all()
        for p in cells
    ):
        raise ValueError("Carrier cells must be valid interior-disjoint polygons")
    tree = STRtree(cells)
    if any(
        i < j and cell.relate_pattern(cells[j], "T********")
        for i, cell in enumerate(cells)
        for j in tree.query(cell, predicate="intersects")
    ):
        raise ValueError("Carrier cells must be valid interior-disjoint polygons")
    cells.sort(key=lambda p: p.normalize().wkb)
    occupied, _ = _interpret_input(raw)
    budget = FidelityBudget(raw, epsilon)
    report = {"outcome": "unresolved", "cells": len(cells), "checked": 0}
    domain = unary_union(cells)
    if not budget.can_drop_input_part(occupied.difference(domain)):
        report["reason"] = "carrier_omits_protected_space"
    required, free, forbidden, conflicts = [], [], [], []
    for i, cell in enumerate(cells):
        must_keep = not budget.can_drop_input_part(cell)
        can_keep = budget.accepts_union(occupied.union(cell))
        if must_keep and not can_keep:
            conflicts.append(i)
        elif must_keep:
            required.append(i)
        elif can_keep:
            free.append(i)
        else:
            forbidden.append(i)
    report.update(
        required=len(required),
        forbidden=len(forbidden),
        free=len(free),
        conflicting=len(conflicts),
    )
    if conflicts:
        report["reason"] = "carrier_cell_conflict"
    return cells, occupied, required, free, report


def occupancy_search(raw, cells, *, delta=0.5, epsilon=0.25, max_free_cells=12):
    """Small finite-carrier oracle; all internal cell walls may dissolve.

    Search labels in increasing original symmetric-difference area until the
    independent contract passes. Unresolved returns no geometry. This studies
    occupancy only, not preservation of required internal source boundaries.
    """
    if type(max_free_cells) is not int or not 0 <= max_free_cells <= 12:
        raise ValueError("max_free_cells must be an integer between 0 and 12")
    raw = list(raw)
    cells, occupied, required, free, report = carrier_labels(
        raw, cells, delta=delta, epsilon=epsilon
    )
    if "reason" in report:
        return None, report
    if len(free) > max_free_cells:
        return None, dict(report, reason="work_limit")
    # Disjoint cells make occupied symmetric-difference area additive. The
    # omitted-domain contribution and fixed labels are constant across choices.
    costs = [
        (cells[i].intersection(occupied).area, cells[i].difference(occupied).area)
        for i in free
    ]
    choices = sorted(
        itertools.product((0, 1), repeat=len(free)),
        key=lambda bits: (sum(cost[b] for cost, b in zip(costs, bits)), bits),
    )
    report["assignments"] = len(choices)
    for bits in choices:
        selected = required + [i for i, b in zip(free, bits) if b]
        candidate = parts(unary_union([cells[i] for i in selected]))
        report["checked"] += 1
        contract = check_cleaning_contract(raw, candidate, delta=delta, epsilon=epsilon)
        if contract["status"] == "pass":
            return candidate, dict(
                report,
                outcome="conforming",
                reason="contract_pass",
                contract=contract,
                changed_area=occupied.symmetric_difference(unary_union(candidate)).area,
            )
    return None, dict(report, reason="carrier_search_exhausted")


def example_carrier(raw, cuts=()):
    """Node explicit example cuts with input boundaries inside their convex hull.

    These hand-specified cuts make the finite search falsifiable. They are not
    an automatic proposal rule for arbitrary inputs.
    """
    domain = unary_union(raw).convex_hull
    linework = unary_union(
        [domain.boundary, *[p.boundary for p in raw], *[p.boundary for p in cuts]]
    )
    return [p for p in polygonize(linework) if domain.covers(p.representative_point())]


def unresolved_vertices(vertices, edges, delta):
    """Proposal seeds only; acceptance remains the independent contract's job.

    Use every endpoint of a short vertex pair, every vertex near a nonincident
    edge, and every boundary junction. Do not chase only the worst witness.
    """
    points = [Point(xy) for xy in vertices]
    point_tree = STRtree(points)
    segments = [LineString(e) for e in edges]
    edge_tree = STRtree(segments)
    degree = dict.fromkeys(vertices, 0)
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1
    seeds = {xy for xy, count in degree.items() if count != 2}
    threshold = delta - delta * _SEPARATION_RELATIVE_TOLERANCE
    for i, (xy, point) in enumerate(zip(vertices, points)):
        for j in point_tree.query(point, predicate="dwithin", distance=threshold):
            if j > i and point.distance(points[j]) < threshold:
                seeds.update((xy, vertices[j]))
        for j in edge_tree.query(point, predicate="dwithin", distance=threshold):
            if xy not in edges[j] and point.distance(segments[j]) < threshold:
                seeds.add(xy)
    return sorted(seeds)


def boundary_chain_chords(vertices, edges, seeds, boundary, band, max_checks=2000):
    """Propose shortcuts across degree-two seed chains, with no length bound.

    Resolved vertices and junctions terminate a chain. Enumerate subchains,
    not pairs across unrelated boundaries; return no partial family on overflow.
    A chord is only a route proposal: joint selection must still preserve cores,
    topology and final feature separation. It never becomes an original parent.
    """
    neighbors = {v: [] for v in vertices}
    for a, b in edges:
        neighbors[a].append(b)
        neighbors[b].append(a)
    removable = {v for v in seeds if len(neighbors[v]) == 2}
    seen = set(edges)
    chords, checks = [], 0
    for a in vertices:
        for first in sorted(neighbors[a]):
            previous, v = a, first
            while v in removable:
                b = next(w for w in neighbors[v] if w != previous)
                if b == a:
                    break
                pair = tuple(sorted((a, b)))
                if pair not in seen:
                    if checks >= max_checks:
                        return None, checks
                    seen.add(pair)
                    checks += 1
                    chord = LineString(pair)
                    if not boundary.covers(chord) and band.covers(chord):
                        chords.append(chord)
                previous, v = v, b
    return chords, checks


def unsupported_short_chord(chord, coordinates, threshold):
    """A short new route needs a possible straight continuation at one end.

    No third collinear site means both endpoints must be active corners, which
    violates the model's vertex separation. Check both endpoint origins using
    the model's exact floating cross product. Any collinear third site suffices
    to retain the chord, even one between its endpoints or without an admitted
    connecting route: this is a conservative impossibility test, not acceptance.
    Original-parent segments must not be filtered by this rule.
    """
    a, b = np.asarray(chord.coords)
    if Point(a).distance(Point(b)) >= threshold:
        return False
    for v, w in ((a, b), (b, a)):
        direction = w - v
        offsets = coordinates - v
        cross = direction[0] * offsets[:, 1] - direction[1] * offsets[:, 0]
        others = np.any(offsets != 0, axis=1) & np.any(coordinates - w != 0, axis=1)
        if np.any((cross == 0) & others):
            return False
    return True


def boundary_candidates(
    raw,
    *,
    delta=0.5,
    epsilon=0.25,
    max_sites=2000,
    max_chords=2000,
    focus_unresolved=False,
    for_routes=True,
):
    """Bounded boundary-route proposals for one original interaction group.

    Sample nearby boundary edges at epsilon/2 from each vertex projection,
    within R = 2*epsilon + delta. Add all chords of length <= R contained in
    the strict original free band. Original boundaries remain in the partition.
    Caps return no partial proposal set. Keep original-edge provenance so a
    selected straight route need not acquire bends from interpolated samples.
    The route experiment focuses sampling on unresolved vertices and permits
    local chords only inside a common seed disk. Original edges outside every seed
    disk are fixed. The historical eager-carrier experiment retains all seeds.
    Route proposals exclude chords already represented by sampled original
    segments before applying the cap. The eager carrier uses raw boundaries,
    so it requests for_routes=False to retain its original linework semantics.
    Focused routes also include whole-chain shortcuts through degree-two seeds,
    with at most 2,000 chain checks and the same combined chord cap.
    New route chords with two unavoidable subscale corners are excluded before
    the cap. The site set and original-parent segments remain unchanged.
    """
    _validate_scale(delta, "delta", positive=True)
    _validate_scale(epsilon, "epsilon", positive=False)
    radius = 2 * epsilon + delta
    _validate_scale(radius, "interaction distance", positive=True)
    for name, value in (("max_sites", max_sites), ("max_chords", max_chords)):
        if type(value) is not int or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    raw = list(raw)
    occupied, _ = _interpret_input(raw)
    original = parts(occupied)
    state = admissibility(original, delta)
    report = {
        "initially_resolved": state["resolved"],
        "sites": 0,
        "chords": 0,
        "duplicate_chords": 0,
        "unsupported_chords": 0,
    }
    if occupied.is_empty:
        return {"occupied": occupied}, dict(report, reason="empty_input")
    if state["resolved"] or epsilon == 0:
        return {"occupied": occupied}, dict(report, reason="input_partition")
    spacing = epsilon / 2
    if spacing == 0 or not math.isfinite(radius / spacing):
        raise ValueError("epsilon is too small for representable sampling")
    if spacing <= np.spacing(max(1.0, np.abs(get_coordinates(occupied)).max())):
        return None, dict(report, reason="sampling_precision")
    boundary = occupied.boundary
    vertices, edges = _canonical_graph([boundary])
    segments = [LineString(e) for e in edges]
    parent_sites = [{a: 0.0, b: e.length} for (a, b), e in zip(edges, segments)]
    edge_tree = STRtree(segments)
    sites = set(vertices)
    if len(sites) > max_sites:
        return None, dict(report, reason="site_limit", sites=len(sites))
    seeds = (
        unresolved_vertices(vertices, edges, delta) if focus_unresolved else vertices
    )
    report["seed_vertices"] = len(seeds)
    editable_parents = set()
    steps = math.ceil(radius / spacing)
    for xy in seeds:
        point = Point(xy)
        for j in sorted(edge_tree.query(point, predicate="dwithin", distance=radius)):
            editable_parents.add(int(j))
            edge = segments[j]
            center = edge.project(point)
            first = max(-steps, math.ceil(-center / spacing))
            last = min(steps, math.floor((edge.length - center) / spacing))
            for k in range(first, last + 1):
                distance = center + k * spacing
                if not 0 < distance < edge.length:
                    continue
                sample = edge.interpolate(distance)
                if sample.distance(point) <= radius:
                    xy_sample = tuple(sample.coords[0])
                    sites.add(xy_sample)
                    parent_sites[j][xy_sample] = distance
                    if len(sites) > max_sites:
                        return None, dict(report, reason="site_limit", sites=len(sites))
    parent_routes = []
    for parent, samples in enumerate(parent_sites):
        ordered = sorted(samples, key=lambda xy: (samples[xy], xy))
        parent_routes.extend(
            (parent, a, b, samples[b] - samples[a])
            for a, b in zip(ordered, ordered[1:])
        )
    original_pairs = {(a, b) for _, a, b, _ in parent_routes}
    points = [Point(xy) for xy in sorted(sites)]
    coordinates = np.asarray(sorted(sites))
    threshold = delta - delta * _SEPARATION_RELATIVE_TOLERANCE
    report["sites"] = len(points)
    point_tree = STRtree(points)
    neighborhoods = None
    if focus_unresolved:
        seed_tree = STRtree([Point(xy) for xy in seeds])
        neighborhoods = [
            set(seed_tree.query(p, predicate="dwithin", distance=radius))
            for p in points
        ]
    budget = FidelityBudget(raw, epsilon)
    # Reuse the checker's declared strict band; never create another tolerance.
    protected, allowed = budget._bounds[1]
    band = allowed.difference(protected)
    chords = []
    chord_pairs = set()
    if focus_unresolved and for_routes:
        chords, report["chain_checks"] = boundary_chain_chords(
            vertices, edges, seeds, boundary, band
        )
        if chords is None:
            return None, dict(report, reason="chain_limit")
        chord_pairs = {tuple(c.coords) for c in chords}
        retained = [
            c for c in chords if not unsupported_short_chord(c, coordinates, threshold)
        ]
        report["unsupported_chords"] += len(chords) - len(retained)
        chords = retained
        report["chain_chords"] = len(chords)
        if len(chords) > max_chords:
            return None, dict(report, reason="chord_limit", chords=len(chords))
    for i, point in enumerate(points):
        for j in sorted(point_tree.query(point, predicate="dwithin", distance=radius)):
            if j <= i:
                continue
            if neighborhoods is not None and not neighborhoods[i] & neighborhoods[j]:
                continue
            chord = LineString([point.coords[0], points[j].coords[0]])
            if boundary.covers(chord) or not band.covers(chord):
                continue
            if for_routes and tuple(chord.coords) in original_pairs:
                report["duplicate_chords"] += 1
                continue
            if tuple(chord.coords) in chord_pairs:
                continue
            if for_routes and unsupported_short_chord(chord, coordinates, threshold):
                report["unsupported_chords"] += 1
                continue
            chords.append(chord)
            if len(chords) > max_chords:
                return None, dict(report, reason="chord_limit", chords=len(chords))
    return {
        "occupied": occupied,
        "budget": budget,
        "sites": sorted(sites),
        "parent_edges": edges,
        "parent_routes": parent_routes,
        "chords": chords,
        "editable_parents": (
            editable_parents if focus_unresolved else set(range(len(edges)))
        ),
    }, dict(report, reason="proposed", chords=len(chords))


def automatic_carrier(raw, *, delta=0.5, epsilon=0.25, max_sites=2000, max_chords=2000):
    """Historical eager arrangement experiment, capped at 20,000 cells."""
    proposals, report = boundary_candidates(
        raw,
        delta=delta,
        epsilon=epsilon,
        max_sites=max_sites,
        max_chords=max_chords,
        for_routes=False,
    )
    if proposals is None:
        return None, report
    occupied = proposals["occupied"]
    if occupied.is_empty:
        return [], dict(report, cells=0)
    if "chords" not in proposals:
        cells = example_carrier(parts(occupied))
        return cells, dict(report, cells=len(cells))
    boundary, chords = occupied.boundary, proposals["chords"]
    domain = occupied.convex_hull
    linework = unary_union([boundary, domain.boundary, *chords])
    cells = []
    for p in polygonize(linework):
        if domain.covers(p.representative_point()):
            cells.append(p)
            if len(cells) > 20000:
                return None, dict(
                    report, reason="cell_limit", chords=len(chords), cells=len(cells)
                )
    cells.sort(key=lambda p: p.normalize().wkb)
    return cells, dict(
        report, reason="constructed", chords=len(chords), cells=len(cells)
    )


def topology_examples():
    from sandbox.cleaning_contract_probe import examples

    simple = examples()
    touching = [box(0, 0, 2, 2), box(2, 2, 4, 4)]
    corner_cuts = [
        Polygon([(2, 2), (1.6, 2), (2, 1.6)]),
        Polygon([(2, 2), (2.4, 2), (2, 2.4)]),
    ]
    cases = {
        "point_contact": (touching, corner_cuts),
        "near_buildings": (simple["near_buildings"], []),
        "courtyard_passage": (simple["courtyard_passage"], [box(5.9, 7, 6.1, 10)]),
        "tiny_hole": (simple["tiny_hole"], []),
    }
    return {
        name: (raw, example_carrier(raw, cuts)) for name, (raw, cuts) in cases.items()
    }


def topology_probe():
    reports, plotted = {}, {}
    examples = topology_examples()
    for name, (raw, cells) in examples.items():
        output, report = occupancy_search(raw, cells)
        reports[name] = report
        plotted[name] = raw, output, cells
    touching = examples["point_contact"][0]
    _, reports["point_contact_without_candidate_cuts"] = occupancy_search(
        touching, example_carrier(touching)
    )
    # The same raw input is feasible on the refined carrier. Exhaustion of the
    # coarser carrier therefore cannot be a proof of geometric infeasibility.
    return reports, plotted


def automatic_probe():
    """Measure automatic partitions; witnesses are evidence, never search input."""
    reports = {}
    for name, (raw, manual_cells) in topology_examples().items():
        started = time.perf_counter()
        cells, report = automatic_carrier(raw)
        if cells is not None:
            report["labels"] = carrier_labels(raw, cells)[-1]
            target, _ = occupancy_search(raw, manual_cells)
            if name == "point_contact":
                # An independent analytic witness at a sampled length. This is
                # not passed to automatic_carrier or used to generate its cuts.
                t = 0.375
                target = parts(
                    unary_union(raw).difference(
                        unary_union(
                            [
                                Polygon([(2, 2), (2 - t, 2), (2, 2 - t)]),
                                Polygon([(2, 2), (2 + t, 2), (2, 2 + t)]),
                            ]
                        )
                    )
                )
            for label, reference in (
                ("original", unary_union(raw)),
                ("witness", unary_union(target)),
            ):
                rebuilt = parts(
                    unary_union(
                        [c for c in cells if reference.covers(c.representative_point())]
                    )
                )
                report[label] = {
                    "reconstruction_difference_area": reference.symmetric_difference(
                        unary_union(rebuilt)
                    ).area,
                    "contract": check_cleaning_contract(raw, rebuilt, delta=0.5),
                }
        report["seconds"] = time.perf_counter() - started
        reports[name] = report
    return reports


def saved_city_groups(run):
    """Read the same saved raw input and partition for both bounded experiments."""
    from benchmarks.benchmark_phases import load_cleaning
    from dtcc_core.model import GeometryType
    from sandbox.cleaning_graph_prototype import interaction_groups

    for record in json.loads((run / "results.json").read_text())["results"]:
        if record["dataset"] != "city_flat_mesh":
            continue
        if record["parameters"]["merge_buildings"] is not True:
            raise ValueError(
                "Occupancy-only carrier survey requires merging permission"
            )
        artifact = (run / record["artifacts"]["cleaning_input"]["path"]).resolve()
        if not artifact.is_relative_to(run.resolve()):
            raise ValueError("Saved cleaning artifact must be inside its run")
        city, _ = load_cleaning(artifact, record)
        raw = [
            b.flatten_geometry(GeometryType.LOD0).to_polygon(simplify=0)
            for b in city.buildings
        ]
        atoms = [p for original in raw for p in parts(_interpret_input([original])[0])]
        delta = record["parameters"]["min_building_detail"]
        epsilon = delta / 2
        groups = interaction_groups(atoms, 2 * epsilon + delta + delta * 1e-6)
        yield record, delta, epsilon, [
            [atoms[i] for i in indices] for indices in groups
        ]


def carrier_survey(run, destination):
    """Measure fixed saved raw groups; do not search or produce cleaned output."""
    rows = []
    for record, delta, epsilon, groups in saved_city_groups(run):
        metrics = []
        for group in groups:
            started = time.perf_counter()
            cells, report = automatic_carrier(group, delta=delta, epsilon=epsilon)
            if cells is not None and not report["initially_resolved"]:
                try:
                    report["labels"] = carrier_labels(
                        group, cells, delta=delta, epsilon=epsilon
                    )[-1]
                except ValueError as error:
                    # A generated invalid partition is research evidence, not
                    # permission to repair it or weaken the boundary validator.
                    report.update(
                        reason="invalid_carrier",
                        error=str(error),
                        invalid_cells=sum(not c.is_valid for c in cells),
                    )
            report["seconds"] = time.perf_counter() - started
            metrics.append(report)
        row = {
            "case": record["case"],
            "delta": delta,
            "epsilon": epsilon,
            "groups": metrics,
            "seconds": sum(g["seconds"] for g in metrics),
        }
        rows.append(row)
        destination.write_text(json.dumps(rows, indent=2) + "\n")
        print(
            record["case"]["city"],
            len(groups),
            "groups",
            round(row["seconds"], 2),
            "seconds",
            flush=True,
        )
    return rows


def plot_topology(
    cases, destination, *, title=None, input_label="input + supplied partition"
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from dtcc_core.plotting.style import plot_polygon_geometries, set_axes_extent

    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 12), layout="constrained")
    for (name, (raw, output, cells)), pair in zip(cases.items(), axes):
        for ax, polygons in zip(pair, (raw, output or [])):
            plot_polygon_geometries(
                ax, polygons, palette=["#66a8bd"], edgecolor="#234653"
            )
            ax.set_aspect("equal")
        for cell in cells:
            for ring in [cell.exterior, *cell.interiors]:
                x, y = ring.xy
                pair[0].plot(x, y, color="#777777", lw=0.7, ls="--")
        pair[0].set_title(name.replace("_", " ") + " — " + input_label)
        pair[1].set_title(
            "Conforming occupancy" if output is not None else "Unresolved"
        )
        set_axes_extent(pair, raw + (output or []))
    fig.suptitle(
        title
        or (
            "One occupancy rule on four explicit partitions · δ = 0.5 m, ε = 0.25 m\n"
            "Finite search demonstration; candidate partition generation is not solved"
        )
    )
    fig.savefig(destination, dpi=150)
    plt.close(fig)


def contact_stage(polygons, sources, metadata):
    diagnostics = footprints._empty_diagnostics(len(polygons))
    diagnostics.update(collect_stage_metrics=False, enable_logging=False)
    output, source_map = footprints._regularize_coverage_contacts(
        polygons,
        sources,
        min_segment_length=metadata["delta"],
        grid=metadata["grid"],
        min_area=0,
        min_hole_area=metadata["min_hole_area"],
        diagnostics=diagnostics,
    )
    return output, source_map


def lund_case():
    fixture = ROOT / "tests/data/cleaning/lund-contact-coupling.geojson"
    data = json.loads(fixture.read_text())
    metadata = data["metadata"]
    raw_features = [f for f in data["features"] if f["properties"]["kind"] == "raw"]
    before_features = [
        f for f in data["features"] if f["properties"]["kind"] == "before"
    ]
    raw = [shape(f["geometry"]) for f in raw_features]
    before = [shape(f["geometry"]) for f in before_features]
    sources = [f["properties"]["source_indices"] for f in before_features]
    after, after_sources = contact_stage(before, sources, metadata)
    alone, _ = contact_stage(before[-1:], sources[-1:], metadata)
    distant_ids = set(sources[-1])
    distant_raw = [
        shape(f["geometry"])
        for f in raw_features
        if distant_ids.intersection(f["properties"]["source_indices"])
    ]
    trigger_raw = [
        shape(f["geometry"])
        for f in raw_features
        if not distant_ids.intersection(f["properties"]["source_indices"])
    ]
    distant_after = [
        p
        for p, indices in zip(after, after_sources)
        if distant_ids.intersection(indices)
    ]
    delta, epsilon = metadata["delta"], metadata["epsilon"]
    report = {
        "fixture": str(fixture.relative_to(ROOT)),
        "distance_from_trigger_to_distant_footprint": min(
            p.distance(before[-1]) for p in before[:-1]
        ),
        "distance_between_original_source_groups": unary_union(trigger_raw).distance(
            unary_union(distant_raw)
        ),
        "independent_above_distance": 2 * epsilon + delta,
        "distant_footprint_unchanged_when_alone": unary_union(alone).equals(before[-1]),
        "distant_raw_fidelity_before": fidelity(distant_raw, before[-1:], epsilon),
        "distant_raw_fidelity_after": fidelity(distant_raw, distant_after, epsilon),
        "additional_loss_relative_to_stage_input": before[-1]
        .buffer(-epsilon, quad_segs=32)
        .difference(unary_union(after))
        .area,
        "before": check_cleaning_contract(raw, before, delta=delta, epsilon=epsilon),
        "after": check_cleaning_contract(raw, after, delta=delta, epsilon=epsilon),
    }
    return report, before, after


def plot_lund(before, after, destination):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from dtcc_core.plotting.style import plot_polygon_geometries, set_axes_extent

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
    for ax, polygons, title in zip(
        axes, (before, after), ("Before contact repair", "After contact repair")
    ):
        plot_polygon_geometries(ax, polygons, palette=["#66a8bd"], edgecolor="#234653")
        ax.set_title(title)
        ax.set_aspect("equal")
    lost = before[-1].buffer(-0.25, quad_segs=32).difference(unary_union(after))
    plot_polygon_geometries(
        axes[1], parts(lost), palette=["#d95c4f"], edgecolor="#b63e32"
    )
    set_axes_extent(axes, before + after)
    point = before[-1].representative_point()
    for ax in axes:
        ax.annotate("Independent footprint", (point.x, point.y), ha="center")
    fig.suptitle(
        "Lund: a close pair triggers changes to a separate footprint\n"
        "Red: protected interior lost from the incoming footprint (epsilon = 0.25 m)"
    )
    fig.savefig(destination, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--run", type=Path, help="Measure automatic partitions on saved raw city groups"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    reports = {
        name: check_cleaning_contract(raw, candidate, delta=0.5, epsilon=0.25)
        for name, (raw, candidate) in counterexamples().items()
    }
    # These are mathematical counterexamples, not expected production failures.
    assert all(r["fidelity"]["status"] == "pass" for r in reports.values())
    assert all(not r["admissibility"]["resolved"] for r in reports.values())
    lund, before, after = lund_case()
    topology, plotted = topology_probe()
    report = {
        "counterexamples": reports,
        "lund": lund,
        "topology": topology,
        "automatic_carriers": automatic_probe(),
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_lund(before, after, args.output / "lund-contact-coupling.png")
    plot_topology(plotted, args.output / "cleaning-topology-search.png")
    print(json.dumps(report, indent=2))
    if args.run:
        carrier_survey(args.run, args.output / "carrier-cities.json")


if __name__ == "__main__":
    main()
