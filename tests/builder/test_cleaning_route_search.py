"""Boundary-route outcomes; no production branches or city-specific goldens."""

from types import SimpleNamespace
from itertools import product
import time

import numpy as np
import pytest
from shapely.geometry import LineString, Point, Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

from sandbox import cleaning_route_probe as route
from sandbox.cleaning_pipeline_probe import (
    boundary_chain_chords,
    unresolved_vertices,
    unsupported_short_chord,
)
from sandbox.cleaning_graph_prototype import simplify_coverage
from dtcc_core.builder.cleaning.contract import (
    _SEPARATION_RELATIVE_TOLERANCE,
    _canonical_graph,
    check_cleaning_contract,
)


def test_unsupported_short_chords_are_impossible_in_the_existing_model():
    a, b, c = (0, 0), (0.25, 0), (0, 2)
    chord = LineString([a, b])
    sites = sorted([a, b, c])
    threshold = 0.5 - 0.5 * _SEPARATION_RELATIVE_TOLERANCE
    assert unsupported_short_chord(chord, np.array(sites), threshold)
    parents = [tuple(sorted(pair)) for pair in ((a, c), (b, c))]
    proposals = {
        "occupied": Polygon([a, b, c]),
        "sites": sites,
        "parent_edges": parents,
        "parent_routes": [
            (i, u, v, LineString([u, v]).length) for i, (u, v) in enumerate(parents)
        ],
        "chords": [chord],
        "editable_parents": {0, 1},
    }
    _, edges, labels, costs, rows = route.build_model(
        proposals, 0.5, time.perf_counter() + 5
    )
    constraint = rows.constraint()
    assignments = np.array(list(product((0, 1), repeat=len(costs))))
    values = constraint.A @ assignments.T
    feasible = np.all(values >= constraint.lb[:, None], axis=0) & np.all(
        values <= constraint.ub[:, None], axis=0
    )
    assert feasible.any()
    k = labels.index(None)
    assert not assignments[feasible, 2 * k : 2 * k + 2].any()
    # Support is deliberately conservative: a collinear site can be on either
    # side, even between the endpoints, and need not have a connecting route.
    for third in ((0.5, 0), (0.125, 0)):
        assert not unsupported_short_chord(chord, np.array([a, b, third]), threshold)
    assert unsupported_short_chord(chord, np.array([a, b, (0.5, 1e-12)]), threshold)
    boundary_chord = LineString([a, (threshold, 0)])
    assert not unsupported_short_chord(
        boundary_chord, np.array(boundary_chord.coords), threshold
    )


def test_compact_turn_constraints_preserve_discrete_route_choices():
    # Horizontal continuation, a turn, rounded samples sharing an original
    # parent, and a near-collinear chord without that parent guarantee.
    sites = [(0, 0), (-1, 0), (1, 0), (0, 1), (1, 1e-12), (-2, 0)]
    edges = [(1, 0), (0, 2), (0, 3), (0, 4), (5, 0)]
    parents = [0, None, None, 0, None]
    incoming, outgoing = [0, 3, 5, 7, 8], [1, 2, 4, 6, 9]
    active = 10
    rows = route.Rows(11, time.perf_counter() + 5)
    route.add_turn_constraints(
        rows, sites, edges, parents, incoming, outgoing, 0, active
    )
    constraint = rows.constraint()
    straight_pairs = {(0, 1), (0, 3), (1, 4)}
    for left, right, corner in product(range(5), range(5), (0, 1)):
        if left == right:  # The existing edge-direction row forbids reversal.
            continue
        x = np.zeros(11)
        x[incoming[left]] = x[outgoing[right]] = 1
        x[active] = corner
        admitted = np.all(constraint.A @ x <= constraint.ub)
        assert admitted == (
            bool(corner) or tuple(sorted((left, right))) in straight_pairs
        )
    # At an unused vertex, existing flow/activation rows force all these bits to zero.
    assert np.all(constraint.A @ np.zeros(11) <= constraint.ub)


def test_rounded_parallel_routes_keep_distinct_parent_identities_and_costs():
    # A valid thin triangle in metric world coordinates. Sampling its two long
    # edges produces a coincident segment although the original edges differ.
    raw = [Polygon([(1e6, 6e6), (1e6 + 2, 6e6 + 3), (np.nextafter(1e6, np.inf), 6e6)])]
    proposals, _ = route.boundary_candidates(raw, focus_unresolved=True)
    sites, edges, parents, costs, _ = route.route_graph(proposals)
    coincident = [p for p in set(edges) if edges.count(p) > 1]
    assert coincident
    for pair in coincident:
        identities = [parent for edge, parent in zip(edges, parents) if edge == pair]
        assert len(identities) == len(set(identities))
        assert None not in identities
    assert len(costs) == 2 * len(edges)
    # Chords aliasing those original routes were already excluded by the graph.
    # They must not consume the proposal cap before that graph can be built.
    limited, report = route.boundary_candidates(
        raw, focus_unresolved=True, max_chords=len(proposals["chords"])
    )
    assert report["duplicate_chords"] > 0
    assert limited is not None
    assert route.route_graph(limited) == route.route_graph(proposals)
    _, _, _, objective, rows = route.build_model(
        proposals, 0.5, time.perf_counter() + 5
    )
    assert rows.constraint().A.shape[1] == len(objective)


def test_automatic_point_contact_search_satisfies_contract_without_mutating_input():
    # One real solve protects the workflow. The four-case timed survey records
    # search coverage separately: finite time limits cannot guarantee discovery.
    # All four geometric witnesses are checked against the actual model below.
    raw = route.topology_examples()["point_contact"][0]
    original = [p.wkb for p in raw]
    output, report = route.select_routes(raw)
    assert report["outcome"] == "conforming", report
    assert report["contract"]["status"] == "pass"
    assert output is not None
    assert [p.wkb for p in raw] == original


def test_resolved_geometry_is_unchanged_and_zero_budget_is_not_relaxed():
    raw = [
        Polygon([(0, 0), (30, -2), (30, 2)]),
        Polygon(
            box(40, 0, 50, 10).exterior.coords, [box(43, 3, 47, 7).exterior.coords]
        ),
    ]
    output, report = route.select_routes(raw)
    assert report["outcome"] == "conforming"
    assert unary_union(output).equals(unary_union(raw))
    touching = route.topology_examples()["point_contact"][0]
    output, report = route.select_routes(touching, epsilon=0)
    assert output is None and report["outcome"] == "unresolved"
    with pytest.raises(ValueError, match="epsilon"):
        route.select_routes(touching, epsilon=-1)


def test_input_permutation_does_not_change_point_contact_solution():
    raw = route.topology_examples()["point_contact"][0]
    forward, first = route.select_routes(raw)
    reverse, second = route.select_routes(raw[::-1])
    assert first["outcome"] == second["outcome"] == "conforming"
    assert unary_union(forward).equals(unary_union(reverse))


def test_no_solver_incumbent_is_unresolved_and_time_budget_is_validated(monkeypatch):
    raw = route.topology_examples()["point_contact"][0]
    monkeypatch.setattr(
        route, "milp", lambda *a, **k: SimpleNamespace(status=1, x=None)
    )
    output, report = route.select_routes(raw)
    assert output is None and report["outcome"] == "unresolved"
    # Even a purported integral incumbent must pass the independent contract.
    # Empty output here would erase both protected occupied cores.
    monkeypatch.setattr(
        route,
        "milp",
        lambda costs, **k: SimpleNamespace(status=0, x=np.zeros(len(costs))),
    )
    output, report = route.select_routes(raw, max_rounds=1)
    assert output is None and report["outcome"] == "unresolved"
    assert report["last_contract"]["fidelity"]["status"] == "fail"
    with pytest.raises(ValueError, match="max_seconds"):
        route.select_routes(raw, max_seconds=0)


def test_selected_parent_run_does_not_turn_sampling_error_into_a_feature():
    # The middle sample belongs to the original bottom edge. Its rounded XY
    # value must not replace that original straight edge with two new bends.
    sites = [(0, 0), (1, 1e-12), (2, 0), (2, 2), (0, 2)]
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    output, witness = route.selected_regions(
        sites, edges, [0, 0, 1, 2, 3], {0, 2, 4, 6, 8}
    )
    assert witness is None
    assert unary_union(output).equals(box(0, 0, 2, 2))


def test_collinear_replacement_must_preserve_original_parent_transitions():
    # A is a rounded sample of original parent O--C. AB equals AC + CB as
    # linework, but replacing AB with original-parent AC removes A on emission.
    o, c = (1e6, 6e6), (1e6 + 0.75, 6e6 + 1)
    a = tuple(LineString([o, c]).interpolate(1).coords[0])
    b = (2 * c[0] - a[0], 2 * c[1] - a[1])
    d, e = (b[0], b[1] + 4), (o[0], b[1] + 4)
    sites = [o, a, c, b, d, e]
    edges = [(0, 1), (1, 3), (3, 4), (4, 5), (5, 0), (1, 2), (2, 3)]
    parents = [0, None, 1, 2, 3, 0, 4]
    raw = [Polygon([o, c, b, d, e])]
    assert LineString([sites[1], sites[3]]).equals(
        LineString([sites[1], sites[2], sites[3]])
    )
    before, _ = route.selected_regions(sites, edges, parents, {0, 2, 4, 6, 8})
    after, _ = route.selected_regions(sites, edges, parents, {0, 10, 12, 4, 6, 8})
    assert not unary_union(before).equals(unary_union(after))
    assert check_cleaning_contract(raw, before, delta=0.5)["status"] == "pass"
    replacement = check_cleaning_contract(raw, after, delta=0.5)
    assert replacement["fidelity"]["status"] == "pass"
    assert not replacement["admissibility"]["resolved"]


def test_focused_family_represents_valid_long_chain_shortcut():
    # The interaction distance bounds independent groups, not replacement-edge
    # length. A small bend on a long wall can be removed inside the original band.
    raw = [Polygon([(0, 0), (0.15, 0.04), (5, 0), (5, 5), (0, 5)])]
    output, _, report = simplify_coverage(raw, merge_buildings=True)
    assert report["status"] == "pass"
    assert unary_union(output).equals(box(0, 0, 5, 5))
    proposals, proposal_report = route.boundary_candidates(raw, focus_unresolved=True)
    sites, edges, _, _, node_parents = route.route_graph(proposals)
    shortcut = LineString([(0, 0), (5, 0)])
    a, b = map(tuple, shortcut.coords)
    assert shortcut.length > 2 * report["epsilon"] + report["delta"]
    assert [xy for xy in sites if shortcut.covers(Point(xy))] == [a, b]
    i, j = sites.index(a), sites.index(b)
    assert (i, j) in edges
    assert not node_parents[i] & node_parents[j]
    assert proposal_report["chain_chords"] > 0
    # It must remain a new chord, not inherit original-parent coalescing.
    assert proposals["parent_edges"].count((a, b)) == 0
    assert any(c.equals(shortcut) for c in proposals["chords"])
    assert_model_represents(raw, unary_union(output), proposals)


def test_chain_candidates_stop_at_resolved_vertices_and_report_work_limit():
    raw = Polygon([(0, 0), (0.1, 0.02), (0.2, 0.03), (5, 0), (5, 5), (0, 5)])
    vertices, edges = _canonical_graph([raw.boundary])
    seeds = [(0.1, 0.02), (0.2, 0.03)]
    # A generous band isolates chain discovery from geometric admission.
    band = raw.buffer(10)
    chords, checks = boundary_chain_chords(vertices, edges, seeds, raw.boundary, band)
    assert LineString([(0, 0), (5, 0)]) in chords
    assert all(set(c.coords) <= {(0, 0), *seeds, (5, 0)} for c in chords)
    partial, used = boundary_chain_chords(
        vertices, edges, seeds, raw.boundary, band, max_checks=checks - 1
    )
    assert partial is None and used == checks - 1
    junction = (1, 0)
    tips = [(0, 0), (2, 0), (1, 1)]
    candidates, _ = boundary_chain_chords(
        sorted([junction, *tips]),
        sorted(tuple(sorted((junction, p))) for p in tips),
        [junction],
        raw.boundary,
        band,
    )
    assert candidates == []


def assert_model_represents(raw, target, proposals):
    occupied = unary_union(raw)
    sites, edges, parents, costs, rows = route.build_model(
        proposals, 0.5, time.perf_counter() + 30
    )
    target_parts = route.parts(target)
    target_vertices, _ = _canonical_graph([target.boundary])
    directed = []
    for p in target_parts:
        p = orient(p, sign=1)
        for ring in [p.exterior, *p.interiors]:
            directed.extend(zip(ring.coords, list(ring.coords)[1:]))
    assignment = np.zeros(len(costs))
    selected = set()
    for k, (a, b) in enumerate(edges):
        line = LineString([sites[a], sites[b]])
        if not target.boundary.covers(line):
            continue
        if parents[k] is None and line.intersection(occupied.boundary).length > 0:
            continue
        u, v = next((u, v) for u, v in directed if LineString([u, v]).covers(line))
        reverse = np.dot(np.subtract(sites[b], sites[a]), np.subtract(v, u)) < 0
        index = 2 * k + int(reverse)
        selected.add(index)
        assignment[index] = 1
    for v, xy in enumerate(sites):
        assignment[2 * len(edges) + v] = xy in target_vertices
    protected, allowed = proposals["budget"]._bounds[1]
    for geometry, label in (
        (protected, 1),
        (occupied.convex_hull.difference(allowed), 0),
    ):
        for p in route.parts(geometry):
            rows.add(
                route.winding_terms(sites, edges, p.representative_point()),
                label,
                label,
            )
    constraint = rows.constraint()
    values = constraint.A @ assignment
    assert np.all(values >= constraint.lb)
    assert np.all(values <= constraint.ub)
    output, witness = route.selected_regions(sites, edges, parents, selected)
    assert witness is None and unary_union(output).equals(target)
    assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"
    # An unselected original edge must have been admitted to an edit region.
    for k, parent in enumerate(parents):
        if parent is not None and parent not in proposals["editable_parents"]:
            assert 2 * k in selected or 2 * k + 1 in selected


def test_focused_model_retains_analytic_solutions_and_fixes_remote_edges():
    # Independent witnesses prove representation, not the solver's ability to
    # discover an incumbent within its time limit. None is supplied to search.
    for name, (raw, _) in route.topology_examples().items():
        occupied = unary_union(raw)
        if name == "point_contact":
            target = occupied.difference(
                unary_union(
                    [
                        Polygon([(2, 2), (1.625, 2), (2, 1.625)]),
                        Polygon([(2, 2), (2.375, 2), (2, 2.375)]),
                    ]
                )
            )
        elif name == "near_buildings":
            target = box(0, 0, 20.2, 6)
        elif name == "courtyard_passage":
            target = box(0, 0, 12, 10).difference(box(4, 3, 8, 7))
        else:
            target = box(0, 0, 10, 10)
        proposals, _ = route.boundary_candidates(raw, focus_unresolved=True)
        assert_model_represents(raw, target, proposals)


def test_focused_seeds_include_vertex_edge_conflicts_without_short_vertex_pairs():
    # Every vertex pair is >= 1 m apart, but (1, 0) is close to the opposite
    # edge. A vertex-pair-only discovery rule would miss the unresolved feature.
    vertices, edges = _canonical_graph([Polygon([(0, 0), (10, 1), (1, 0)]).boundary])
    assert unresolved_vertices(vertices, edges, 0.5) == [(1.0, 0.0)]


def test_unresolved_cli_run_removes_an_older_success_artifact(tmp_path, monkeypatch):
    stale = tmp_path / "point_contact.geojson"
    stale.write_text("older successful output")
    raw = route.topology_examples()["point_contact"][0]
    monkeypatch.setattr(
        route, "topology_examples", lambda: {"point_contact": (raw, [])}
    )
    monkeypatch.setattr(route, "plot_topology", lambda *a, **k: None)
    monkeypatch.setattr(
        route,
        "select_routes",
        lambda *a, **k: (
            None,
            {"outcome": "unresolved", "reason": "solver_unresolved", "rounds": 1},
        ),
    )
    monkeypatch.setattr(
        route.sys, "argv", ["probe", "--focused", "--output", str(tmp_path)]
    )
    route.main()
    assert not stale.exists()
    assert '"unresolved"' in (tmp_path / "routes.json").read_text()
