"""Geometric obligations of the bounded research method, not legacy branches."""

import pytest
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from sandbox.cleaning_graph_prototype import simplify_coverage


def noisy_rectangle():
    return Polygon([(0, 0), (0.15, 0.04), (5, 0), (5, 5), (0, 5)])


def test_shortcut_changes_shared_wall_coherently():
    left = Polygon([(0, 0), (2, 0), (2.05, 0.1), (2, 4), (0, 4)])
    right = Polygon([(2, 0), (4, 0), (4, 4), (2, 4), (2.05, 0.1)])
    original = [left.wkb, right.wkb]
    polygons, sources, report = simplify_coverage([left, right])
    by_source = {tuple(ids): p for p, ids in zip(polygons, sources)}
    assert report["outcome"] == "conforming"
    assert by_source[(0,)].equals(box(0, 0, 2, 4))
    assert by_source[(1,)].equals(box(2, 0, 4, 4))
    assert unary_union(polygons).equals(unary_union([left, right]))
    assert [left.wkb, right.wkb] == original
    assert all(
        g["final_vertices"] <= g["initial_vertices"] - g["accepted_shortcuts"]
        for g in report["groups"]
    )


def test_already_resolved_corner_and_courtyard_are_unchanged():
    sharp = Polygon([(0, 0), (30, -2), (30, 2)])
    courtyard = Polygon(
        box(40, 0, 50, 10).exterior.coords, [box(43, 3, 47, 7).exterior.coords]
    )
    polygons, sources, report = simplify_coverage([sharp, courtyard])
    assert report["outcome"] == "conforming"
    assert report["accepted_shortcuts"] == 0
    assert all(
        p.equals([sharp, courtyard][ids[0]]) for p, ids in zip(polygons, sources)
    )


def test_original_band_is_not_relaxed_to_resolve_geometry():
    raw = [noisy_rectangle()]
    polygons, _, strict = simplify_coverage(raw, epsilon=0)
    assert strict["outcome"] == "unresolved"
    assert strict["fidelity"]["status"] == "pass"
    assert polygons[0].equals(raw[0])
    _, _, allowed = simplify_coverage(raw, epsilon=0.25)
    assert allowed["outcome"] == "conforming"
    assert allowed["fidelity"]["lost_protected_area"] == 0
    assert allowed["fidelity"]["added_outside_budget_area"] == 0


def test_topology_changes_remain_explicitly_unresolved():
    courtyard = Polygon(
        box(0, 0, 10, 10).exterior.coords, [box(4.9, 4.9, 5.1, 5.1).exterior.coords]
    )
    polygons, _, report = simplify_coverage([courtyard])
    assert report["outcome"] == "unresolved"
    assert len(polygons) == 1 and len(polygons[0].interiors) == 1
    assert report["fidelity"]["status"] == "pass"
    touching = [box(0, 0, 2, 2), box(2, 2, 4, 4)]
    polygons, _, report = simplify_coverage(touching)
    assert report["outcome"] == "unresolved"
    assert unary_union(polygons).equals(unary_union(touching))


def test_permutation_and_distant_input_do_not_change_local_result():
    near, distant = noisy_rectangle(), box(20, 0, 25, 5)
    alone, _, _ = simplify_coverage([near])
    for raw in ([near, distant], [distant, near]):
        polygons, sources, report = simplify_coverage(raw)
        assert report["outcome"] == "conforming"
        by_source = {tuple(ids): p for p, ids in zip(polygons, sources)}
        assert by_source[(raw.index(near),)].equals(alone[0])
        assert by_source[(raw.index(distant),)].equals(distant)


def test_input_subdivision_and_merge_permission_preserve_sources():
    raw = [box(0, 0, 3, 3), box(2, 0, 5, 3)]
    polygons, sources, report = simplify_coverage(raw)
    assert report["outcome"] == "conforming"
    assert sorted(sources) == [[0], [0, 1], [1]]
    assert sum(p.area for p in polygons) == unary_union(raw).area
    polygons, sources, report = simplify_coverage(raw, merge_buildings=True)
    assert report["outcome"] == "conforming"
    assert sources == [[0, 1]]
    assert polygons[0].equals(unary_union(raw))


def test_work_limit_is_unresolved_and_invalid_options_fail():
    polygons, _, report = simplify_coverage([noisy_rectangle()], max_candidate_checks=1)
    assert report["outcome"] == "unresolved"
    assert report["groups"][0]["reason"] == "work_limit"
    assert report["accepted_shortcuts"] == 0
    assert polygons[0].equals(noisy_rectangle())
    with pytest.raises(ValueError, match="epsilon"):
        simplify_coverage([], epsilon=-1)
    with pytest.raises(ValueError, match="merge_buildings"):
        simplify_coverage([], merge_buildings="false")


def test_shortcut_cannot_swallow_a_separate_region():
    # Filling the narrow notch would respect the occupancy band, but would
    # swallow an island without the new chord crossing any island boundary.
    notched = Polygon(
        [(0, 0), (10, 0), (10, 10), (5.1, 10), (5, 8), (4.9, 10), (0, 10)]
    )
    island = box(4.99, 9.4, 5.01, 9.5)
    polygons, sources, report = simplify_coverage([notched, island])
    assert report["admissibility"]["topology_ok"]
    assert sorted(sources) == [[0], [1]]
    assert polygons[0].disjoint(polygons[1])
    assert report["fidelity"]["status"] == "pass"
