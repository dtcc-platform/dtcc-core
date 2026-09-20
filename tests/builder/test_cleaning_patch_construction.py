"""Construction invariants, independent of production repair branches."""

import numpy as np
import pytest
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

from dtcc_core.builder.cleaning.contract import (
    FidelityBudget,
    _canonical_graph,
    check_cleaning_contract,
)
from sandbox.cleaning_patch_probe import (
    clearance_constraints,
    construct,
    coupled_move,
    motion_interval,
    motion_reference,
    move_boundary,
    polygon_parts,
)
from sandbox.cleaning_pipeline_probe import topology_examples


def test_local_construction_resolves_topology_with_fixed_fidelity_and_strict_progress():
    for name, (raw, _) in topology_examples().items():
        original = [p.wkb for p in raw]
        output, report = construct(raw, merge_buildings=True)
        assert output is not None, (name, report)
        assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"
        assert [p.wkb for p in raw] == original
        assert all(b < a for a, b in zip(report["progress"], report["progress"][1:]))
        if name == "point_contact":
            assert report["progress"][-1][2] > report["progress"][0][2]
        if name == "courtyard_passage":
            # Closing or widening the entrance are both permitted. The
            # protected courtyard itself must remain open in either case.
            assert not unary_union(output).covers(Point(6, 5))


def test_long_shallow_wall_is_not_excluded_by_a_local_edge_length_cap():
    raw = [Polygon([(0, 0), (0.15, 0.04), (5, 0), (5, 5), (0, 5)])]
    output, report = construct(raw, merge_buildings=True)
    assert output is not None and report["outcome"] == "conforming"
    assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"


def test_clearance_can_be_resolved_by_moving_a_corner_without_changing_topology():
    # Patch-only edits stall here. A 0.1 m corner displacement fits the original
    # budget and separates the vertex from the long opposing wall.
    raw = [box(0, 0, 10, 20), Polygon([(10.4, 5), (15, 5), (15, 15), (11, 15)])]
    output, report = construct(raw, merge_buildings=True)
    assert output is not None and len(output) == 2
    assert report["progress"] == [(0, 1, 8), (0, 0, 8)]
    assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"

    # An extra collinear ring point must not pin the middle of a moved wall.
    densified = Polygon([(10, 0), (10, 10), (10, 20), (0, 20), (0, 0)])
    repeated, _ = construct([densified, raw[1]], merge_buildings=True)
    assert repeated is not None
    assert unary_union(output).symmetric_difference(unary_union(repeated)).area < 1e-9


def test_coupled_edit_can_restore_an_adjacent_corner_within_original_fidelity():
    # A three-rectangle witness: an earlier corner move can obstruct completion
    # even though the remaining vertex/edge conflict has an obvious direction.
    raw = [box(0, 0, 2, 10), box(-3, 1, -0.1, 3), box(-3, 6, -0.02, 8)]
    stalled = [
        box(0.2, 0, 2, 10),
        box(-3, 1, -0.3, 3),
        Polygon([(-3, 6), (-0.3, 6), (-0.02, 8), (-3, 8)]),
    ]
    state = check_cleaning_contract(raw, stalled, delta=0.5)
    assert state["fidelity"]["status"] == "pass"
    assert state["admissibility"]["subscale_pairs"] == 1

    # Resolving the last gap by moving only its corner consumes protected core.
    isolated = [*stalled[:2], box(-3, 6, -0.3, 8)]
    rejected = check_cleaning_contract(raw, isolated, delta=0.5)
    assert rejected["admissibility"]["resolved"]
    assert rejected["fidelity"]["status"] == "fail"

    # Restore the lower corner by 0.04 m while moving the upper corner left
    # and the opposing wall right. The whole edit uses the ORIGINAL budget.
    coupled = [box(0.24, 0, 2, 10), stalled[1], box(-3, 6, -0.26, 8)]
    accepted = check_cleaning_contract(raw, coupled, delta=0.5)
    assert accepted["status"] == "pass"
    assert accepted["admissibility"]["canonical_vertex_count"] == 12

    # Generate a correction automatically from the already-stalled state,
    # retaining the ORIGINAL reference rather than cleaning that state as raw.
    occupied = unary_union(stalled)
    vertices, edges = _canonical_graph([occupied.boundary])
    incident = {v: [] for v in vertices}
    for i, edge in enumerate(edges):
        for v in edge:
            incident[v].append(i)
    move, _ = coupled_move(
        (-0.02, 8),
        [(0.2, 0), (0.2, 10)],
        edges,
        incident,
        motion_reference(FidelityBudget(raw, 0.25)),
        0.5,
        0.25,
    )
    assert move is not None
    generated = polygon_parts(move_boundary(occupied, move, incident))
    assert check_cleaning_contract(raw, generated, delta=0.5)["status"] == "pass"

    output, report = construct(raw, merge_buildings=True)
    assert output is not None and report["outcome"] == "conforming"
    assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"


def test_near_contact_motion_preserves_original_material_at_large_coordinates():
    # The wall is one floating-point step inside the strict protected core.
    # Its tiny area error is admitted, but preserving the corner's measured
    # side would let the wall move a further 5 cm into protected material.
    x, y = 374000, 6163700
    raw = [box(x, y, x + 10, y + 2), box(x + 10.1, y + 0.6, x + 13, y + 1.4)]
    budget = FidelityBudget(raw, 0.25)
    wall = np.nextafter(budget._bounds[1][0].intersection(raw[0]).bounds[2], -np.inf)
    occupied = unary_union([box(x, y, wall, y + 2), raw[1]])
    assert budget.accepts_union(occupied)
    before = check_cleaning_contract(raw, polygon_parts(occupied), delta=0.5)
    vertices, edges = _canonical_graph([occupied.boundary])
    incident = {v: [] for v in vertices}
    for i, edge in enumerate(edges):
        for v in edge:
            incident[v].append(i)
    move, _ = coupled_move(
        (x + 10.1, y + 0.6),
        [(wall, y), (wall, y + 2)],
        edges,
        incident,
        motion_reference(budget),
        0.5,
        0.25,
    )
    assert move is not None
    candidate = move_boundary(occupied, move, incident)
    assert budget.accepts_union(candidate)
    after = check_cleaning_contract(raw, polygon_parts(candidate), delta=0.5)
    assert after["fidelity"]["status"] == "pass"
    assert (
        after["admissibility"]["subscale_pairs"]
        < before["admissibility"]["subscale_pairs"]
    )


def test_motion_interval_cannot_jump_across_original_protected_material():
    budget = FidelityBudget([box(0, 0, 2, 2)], 0.25)
    _, _, band = motion_reference(budget)
    normal = np.array([1.0, 0.0])
    lo, hi = motion_interval((0, 1), normal, band, 3)
    # This ray also crosses a second allowed interval at the opposite wall.
    # Only the connected interval at the original corner may be used.
    assert -0.25 <= lo < 0 < hi <= 0.25
    assert motion_interval((10, 10), normal, band, 1) is None


def test_motion_resolves_gap_without_shortening_neighboring_edge():
    # Minimum movement for the initiating gap alone shortens the bevel below
    # delta. A feasible motion exists within the same variable support/budget.
    raw = [
        box(-3, -3, -0.44, 3),
        Polygon([(0, 0), (0.45, -0.25), (1.5, -5), (5, -5), (5, 5), (0.8, 5)]),
    ]
    output, report = construct(raw, merge_buildings=True)
    assert output is not None
    assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"
    assert report["progress"] == [(0, 1, 10), (0, 0, 10)]


def test_motion_constraints_protect_fixed_vertices_opposite_moving_edge():
    raw = [box(0, -3, 4, 0), box(1.9, 0.6, 2.1, 2.6)]
    _, edges = _canonical_graph([unary_union(raw).boundary])
    ids = {(0.0, 0.0): 0, (4.0, 0.0): 1}
    constraints = list(clearance_constraints(edges, ids, np.array([0, 1]), 1, 0.5))
    # Parallel motion needs one strongest bound per directed variable pair
    # (including the fixed zero), rather than one row per geometric neighbor.
    assert len(constraints) <= 2 * 1 + 2 * 2 + 1
    rows, bounds = map(np.asarray, zip(*constraints))
    # Both moving endpoints remain far from the other polygon, but its fixed
    # bottom corners would be only 0.4 m from the raised edge's interior.
    assert np.all(rows @ np.zeros(2) <= bounds + 1e-12)
    assert np.any(rows @ np.array([0.2, 0.2]) > bounds + 1e-12)


def test_local_band_retains_complete_motion_sections(monkeypatch):
    from sandbox import cleaning_patch_probe as probe

    raw = [box(0, 0, 2, 2), box(2.25, 0, 4, 2), box(100, 100, 110, 110)]
    _, _, original_band = motion_reference(FidelityBudget(raw, 0.25))
    queried = []

    def compare_sections(vertex, normal, local_band, radius):
        point = np.asarray(vertex)
        ray = LineString([point - radius * normal, point + radius * normal])
        assert original_band.intersection(ray).equals(local_band.intersection(ray))
        queried.append(local_band.area < original_band.area)
        return motion_interval(vertex, normal, local_band, radius)

    monkeypatch.setattr(probe, "motion_interval", compare_sections)
    output, _ = construct(raw, merge_buildings=True)
    assert queried and all(queried)
    assert output is not None
    assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"


def test_resolved_input_is_unchanged_and_order_does_not_choose_a_different_repair():
    courtyard = Polygon(
        box(0, 0, 10, 10).exterior.coords, [box(3, 3, 7, 7).exterior.coords]
    )
    output, report = construct([courtyard], merge_buildings=True)
    assert report["edits"] == report["checks"] == 0
    assert unary_union(output).equals(courtyard)
    raw = topology_examples()["point_contact"][0]
    forward, _ = construct(raw, merge_buildings=True)
    reverse, _ = construct(raw[::-1], merge_buildings=True)
    assert unary_union(forward).equals(unary_union(reverse))


def test_unresolved_is_not_delivered_as_cleaned_and_budget_is_not_relaxed():
    raw = topology_examples()["point_contact"][0]
    output, report = construct(raw, epsilon=0, merge_buildings=True)
    assert output is None and report["outcome"] == "unresolved"
    assert report["edits"] == 0
    assert report["contract"]["fidelity"]["status"] == "pass"
    output, report = construct(raw, max_checks=1, merge_buildings=True)
    assert output is None and report["reason"] == "work_limit"


def test_unsupported_source_policy_and_invalid_parameters_fail_clearly():
    with pytest.raises(ValueError, match="merge_buildings"):
        construct([box(0, 0, 2, 2)])
    with pytest.raises(ValueError, match="epsilon"):
        construct([], epsilon=-1, merge_buildings=True)
    with pytest.raises(ValueError, match="max_checks"):
        construct([], max_checks=True, merge_buildings=True)
