"""Regression tests for the bounded footprint constructor."""

import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union

from dtcc_core.builder.cleaning import construction
from dtcc_core.builder.cleaning.contract import FidelityBudget, check_cleaning_contract
from tests.builder.cleaning_fixtures import survey_merge_group, topology_examples


def test_required_construction_witnesses():
    cases = topology_examples()
    cases.update(
        coupled=[box(0, 0, 2, 10), box(-3, 1, -0.1, 3), box(-3, 6, -0.02, 8)],
        long_wall=[Polygon([(0, 0), (0.15, 0.04), (5, 0), (5, 5), (0, 5)])],
        sampled_square=[
            Polygon([(x, 0) for x in np.linspace(0, 5, 101)] + [(5, 5), (0, 5)])
        ],
        opposing_wall=[
            box(-3, -3, -0.44, 3),
            Polygon([(0, 0), (0.45, -0.25), (1.5, -5), (5, -5), (5, 5), (0.8, 5)]),
        ],
    )
    for name, raw in cases.items():
        original = [p.wkb for p in raw]
        output, report = construction.construct(iter(raw))
        assert output is not None, (name, report)
        assert check_cleaning_contract(raw, output, delta=0.5)["status"] == "pass"
        assert [p.wkb for p in raw] == original
        if name == "courtyard_passage":
            assert not unary_union(output).covers(Point(6, 5))


@pytest.mark.parametrize(
    "case_id,group_id",
    [
        ("city_grid:helsingborg:015", 2),
        ("city_grid:norrkoping:039", 23),
    ],
)
def test_full_survey_conformance_regressions_use_general_structural_snap(
    case_id, group_id
):
    raw = survey_merge_group(case_id, group_id)
    for candidate in (raw, [translate(polygon, 17.125, -9.875) for polygon in raw]):
        output, report = construction.construct(candidate)
        assert output is not None, report
        assert check_cleaning_contract(candidate, output, delta=0.5)["status"] == "pass"
        assert any(
            attempt.get("accepted_for_ranking")
            for attempt in report["fallback_attempts"]
        )


def test_dense_structural_proposal_repairs_one_cap_case_and_explains_the_residual():
    repaired, report = construction.construct(
        survey_merge_group("city_grid:helsingborg:007", 42)
    )
    assert repaired is not None, report
    assert report["evaluations"] < 1000

    unresolved, report = construction.construct(
        survey_merge_group("city_grid:linkoping:066", 45)
    )
    assert unresolved is None
    assert report["contract"]["admissibility"]["subscale_pairs"] == 2
    assert report["selected_attempt_work"]["evaluations"] < 500
    assert report["terminal_rejection"]["reason"] in {
        "fidelity",
        "global_fidelity",
        "global_progress",
        "no_candidate_helped",
    }
    assert report["fallback_attempts"][1]["reason"] == "residual_defects"
    assert report["total_work"]["evaluations"] > report["selected_attempt_work"][
        "evaluations"
    ]


def test_search_stops_at_the_first_candidate_with_preferred_angles():
    raw = survey_merge_group("city_grid:norrkoping:003", 0)
    output, report = construction.construct(raw)
    assert output is not None, report
    assert report["fallback_attempts"][0]["accepted_for_ranking"] is True
    assert report["mesher_profile"]["preferred"]["status"] == "pass"
    assert report["reason"] == "preferred_candidate"
    assert [row["reason"] for row in report["fallback_attempts"][1:]] == [
        "preferred_candidate_found",
        "preferred_candidate_found",
    ]
    assert report["evaluations"] < 100


def test_a_conforming_result_without_preferred_angles_does_not_stop_the_search(
    monkeypatch,
):
    raw = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]
    variants = [[box(0, 0, 10, 5 + 0.05 * index)] for index in range(3)]
    calls = []

    def attempt(*args, **kwargs):
        repaired = variants[len(calls)]
        calls.append(len(calls))
        contract = check_cleaning_contract(raw, repaired, delta=0.5)
        return repaired, dict(
            contract=contract, outcome="conforming", stages={},
            final=construction.defect_tuple(contract["admissibility"]),
            edits=0, sites=0, evaluations=1, admissions=0, global_evaluations=0,
        )

    def finish(raw, output, **kwargs):
        preferred = {"status": "pass" if len(calls) == 3 else "fail"}
        return output, {"status": "pass", "preferred": preferred}

    monkeypatch.setattr(construction, "attempt", attempt)
    monkeypatch.setattr(construction, "_finish_mesher_profile", finish)
    monkeypatch.setattr(
        construction, "_fallback_preprocessing_is_equivalent", lambda *args: False
    )
    output, report = construction.construct(raw)
    assert calls == [0, 1, 2]
    assert output == variants[2]
    assert report["reason"] == "preferred_candidate"


def test_local_improvement_cannot_hide_new_distant_conflicts():
    raw = [Polygon([(0, 0), (0.1, 0.01), (100, 1), (100, 10), (0, 10)])]
    raw += [
        box(x, x * 0.01 - 0.49997 - 2, x + 2, x * 0.01 - 0.49997)
        for x in (10, 20, 30, 40, 50)
    ]
    occupied = unary_union(raw)
    candidate = construction.rebuild(occupied, {(0.1, 0.01): None})
    judge = construction.LocalJudge(FidelityBudget(raw, 0.25), 0.5)
    coordinates = np.array([(0, 0), (0.1, 0.01)])
    clip = judge.clip(judge.window(coordinates))
    assert construction.separation_guard(
        judge.state(occupied, clip), judge.state(candidate, clip)
    )
    assert judge.admits(candidate)
    accepted, reason = construction.repair_site(
        occupied,
        judge,
        coordinates,
        0.5,
        construction.separation_guard,
        lambda *args: iter([("remove_vertex", candidate)]),
    )
    assert accepted is None and reason == "global_progress"


def test_equal_occupancy_is_recognised_across_representations():
    square = box(0, 0, 10, 10)
    resampled = Polygon([(0, 0), (5, 0), (10, 0), (10, 10), (0, 10)])
    moved = box(0, 0, 10, 10.001)
    assert construction._same_occupancy(resampled, square, square.area)
    assert not construction._same_occupancy(moved, square, square.area)
    # An equal area alone does not make two sets equal.
    shifted = box(0.001, 0, 10.001, 10)
    assert not construction._same_occupancy(shifted, square, square.area)


def test_local_fidelity_clip_cache_is_invalidated_between_sites():
    raw = [box(0, 0, 10, 10)]
    judge = construction.LocalJudge(FidelityBudget(raw, 0.25), 0.5)
    candidate = box(0.4, 0, 10, 10)
    unaffected = box(8, -1, 11, 11)
    affected = box(-1, -1, 2, 11)

    assert judge.keeps_fidelity(candidate, unaffected)
    assert not judge.keeps_fidelity(candidate, affected)
    assert judge.keeps_fidelity(candidate, unaffected)
    # The same holds for material added beyond the envelope.
    grown = box(0, 0, 10.4, 10)
    assert judge.keeps_fidelity(grown, affected)
    assert not judge.keeps_fidelity(grown, unaffected)
    assert judge.keeps_fidelity(grown, affected)


def test_rebuilder_cache_is_bounded_to_its_occupied_geometry():
    occupied = unary_union([box(0, 0, 2, 2), box(5, 0, 7, 2)])
    rebuilder = construction._Rebuilder(occupied)
    replacement = {(2.0, 0.0): (2.0, 0.25)}

    first = rebuilder.rebuild(replacement)
    assert rebuilder.rebuild(replacement) is first

    changed = construction._Rebuilder(first).rebuild({(7.0, 0.0): (7.0, 0.25)})
    fresh = construction.rebuild(first, {(7.0, 0.0): (7.0, 0.25)})
    assert changed.equals_exact(fresh, 0)


def test_collinear_points_do_not_pin_a_wall():
    other = Polygon([(10.4, 5), (15, 5), (15, 15), (11, 15)])
    simple, _ = construction.construct([box(0, 0, 10, 20), other])
    sampled, _ = construction.construct(
        [Polygon([(10, 0), (10, 10), (10, 20), (0, 20), (0, 0)]), other]
    )
    assert simple is not None and sampled is not None
    assert unary_union(simple).symmetric_difference(unary_union(sampled)).area < 1e-9


def test_zero_budget_is_not_silently_relaxed():
    raw = topology_examples()["point_contact"]
    output, report = construction.construct(raw, epsilon=0)
    assert output is None and report["outcome"] == "unresolved"
    assert len(report["sampling_attempts"]) == 1
    assert [row["reason"] for row in report["fallback_attempts"][1:]] == [
        "canonical_preprocessing_equivalence",
        "canonical_preprocessing_equivalence",
    ]


def test_retry_work_is_counted(monkeypatch):
    def attempt(*args, **kwargs):
        return None, dict(
            outcome="unresolved",
            stages={},
            edits=1,
            sites=2,
            evaluations=3,
            admissions=4,
            global_evaluations=5,
        )

    monkeypatch.setattr(construction, "attempt", attempt)
    monkeypatch.setattr(
        construction, "_fallback_preprocessing_is_equivalent", lambda *args: False
    )
    _, report = construction.construct(topology_examples()["point_contact"], epsilon=0)
    assert (
        report["edits"],
        report["sites"],
        report["evaluations"],
        report["admissions"],
        report["global_evaluations"],
    ) == (3, 6, 9, 12, 15)


def test_warning_candidate_never_outranks_a_conforming_fallback(monkeypatch):
    raw = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]
    repaired = [box(0, 0, 10, 5)]
    calls = []

    def attempt(*args, **kwargs):
        calls.append(kwargs["use_bulk_simplify"])
        output = raw if len(calls) == 1 else repaired
        contract = check_cleaning_contract(raw, output, delta=0.5)
        return output, dict(
            contract=contract, outcome="unresolved" if len(calls) == 1 else "conforming",
            stages={}, final=construction.defect_tuple(contract["admissibility"]),
            edits=0, sites=0, evaluations=1, admissions=0, global_evaluations=0,
        )

    monkeypatch.setattr(construction, "attempt", attempt)
    output, report = construction.construct(raw, allow_residual_separation=True)
    assert output == repaired
    assert report["outcome"] == "conforming"
    assert calls[:2] == [True, False]


def test_failed_bulk_attempt_retries_without_bulk_simplification(monkeypatch):
    raw = topology_examples()["point_contact"]
    rescue, _ = construction.construct(raw)
    assert rescue is not None
    calls = []

    def attempt(*args, **kwargs):
        calls.append(kwargs["use_bulk_simplify"])
        output = None if kwargs["use_bulk_simplify"] else rescue
        return output, dict(
            outcome="unresolved" if output is None else "conforming",
            reason="residual_defects" if output is None else "contract_pass",
            final=(1, 0, 8) if output is None else (0, 0, 8),
            stages={"dense_sampling": {"tolerance": None}},
            edits=0,
            sites=0,
            evaluations=1,
            admissions=0,
            global_evaluations=0,
            contract=check_cleaning_contract(raw, output or raw, delta=0.5),
        )

    monkeypatch.setattr(construction, "attempt", attempt)
    monkeypatch.setattr(
        construction, "_fallback_preprocessing_is_equivalent", lambda *args: True
    )
    output, report = construction.construct(raw)
    assert output is not None, report
    assert calls == [True, False]
    assert report["fallback_attempts"][1]["accepted_for_ranking"] is True


@pytest.mark.parametrize(
    "delta,epsilon", [(0, 0.25), (float("nan"), 0.25), (0.5, -1), (0.5, float("inf"))]
)
def test_invalid_parameters(delta, epsilon):
    with pytest.raises(ValueError):
        construction.construct([], delta=delta, epsilon=epsilon)


def _budget_retry_fixture(monkeypatch, outcomes):
    """Two merge groups: a pair 0.1 m apart and a separate square."""
    pair = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]
    single = box(20, 0, 30, 10)
    calls = []

    def construct(raw, *, delta, epsilon, allow_residual_separation, **kwargs):
        calls.append(epsilon)
        if len(raw) == 1:
            return list(raw), {"outcome": "unchanged"}
        outcome = outcomes[epsilon]
        output = {
            "conforming": [box(0, 0, 10, 5.3)],  # adds points 0.3 m out
            "warning": list(raw),
            "unresolved": None,
        }[outcome]
        return output, {"outcome": outcome}

    monkeypatch.setattr(construction, "construct", construct)
    return pair + [single], calls


def test_group_that_fails_is_retried_at_the_larger_budget(monkeypatch):
    raw, calls = _budget_retry_fixture(
        monkeypatch, {0.25: "unresolved", 0.375: "conforming"}
    )
    polygons, sources, report = construction.construct_coverage(
        raw, [[0], [1], [2]], delta=0.5, epsilon=0.25, retry_epsilon=0.375
    )
    assert sorted(calls) == [0.25, 0.25, 0.375]
    assert report["outcome"] == "conforming"
    assert report["fidelity_budget_used"] == 0.375
    contract = report["before_selection_contract"]
    assert contract["epsilon"] == 0.375 and contract["status"] == "pass"
    # The farthest added point is above the middle of the 0.1 m gap.
    assert contract["fidelity"]["achieved_epsilon"] == pytest.approx(
        (0.3**2 + 0.05**2) ** 0.5, abs=0.004
    )
    (retried,) = [row for row in report["group_reports"] if "budget_attempts" in row]
    assert retried["fidelity_budget"] == 0.375
    assert [row["outcome"] for row in retried["budget_attempts"]] == [
        "unresolved",
        "conforming",
    ]
    assert sorted(sources) == [[0, 1], [2]]


def test_a_retry_that_does_not_improve_keeps_the_declared_budget(monkeypatch):
    raw, calls = _budget_retry_fixture(
        monkeypatch, {0.25: "warning", 0.375: "unresolved"}
    )
    polygons, _, report = construction.construct_coverage(
        raw,
        [[0], [1], [2]],
        delta=0.5,
        epsilon=0.25,
        retry_epsilon=0.375,
        allow_residual_separation=True,
    )
    assert report["outcome"] == "warning"
    assert report["fidelity_budget_used"] == 0.25
    assert report["before_selection_contract"]["epsilon"] == 0.25

    calls.clear()
    _, _, report = construction.construct_coverage(
        raw, [[0], [1], [2]], delta=0.5, epsilon=0.25,
        allow_residual_separation=True,
    )
    assert 0.375 not in calls and report["fidelity_budget_used"] == 0.25


def test_early_stop_is_revisited_when_groups_end_up_too_close(monkeypatch):
    # Two merge groups 0.6 m apart. The first preferred candidate of the left
    # group grows towards the right one; ranking every attempt keeps it apart.
    left = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]
    right = box(10.6, 0, 20, 5)
    calls = []

    def construct(raw, *, delta, epsilon, allow_residual_separation, stop_at_preferred):
        calls.append(stop_at_preferred)
        if len(raw) == 1:
            return list(raw), {"outcome": "unchanged"}
        grown = box(0, 0, 10.2, 5) if stop_at_preferred else box(0, 0, 10, 5)
        skipped = [{"reason": "preferred_candidate_found"}] if stop_at_preferred else []
        return [grown], {"outcome": "conforming", "fallback_attempts": skipped}

    monkeypatch.setattr(construction, "construct", construct)
    polygons, _, report = construction.construct_coverage(
        left + [right], [[0], [1], [2]], delta=0.5, epsilon=0.25,
        allow_residual_separation=True,
    )
    assert report["outcome"] == "conforming"
    assert report["cross_group_rebuild"] == {"groups": [0], "adopted": True}
    assert box(0, 0, 10, 5) in polygons
    assert calls.count(False) == 1
