"""The staged research construction must earn its own acceptance evidence."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from shapely import from_wkb
from shapely.affinity import translate
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union

from dtcc_core.builder.cleaning.contract import FidelityBudget, check_cleaning_contract
from sandbox import cleaning_staged_probe as staged
from sandbox.cleaning_delta_budget_probe import measured_costs
from sandbox.cleaning_mesher_precondition_probe import configurations
from sandbox.cleaning_pipeline_probe import occupancy_search, topology_examples


FOCUSED_FIXTURE = (
    Path(__file__).parents[1]
    / "data/cleaning/full-survey-focused-groups.json"
)


def focused_group(case_id, group_id):
    payload = json.loads(FOCUSED_FIXTURE.read_text())
    case = next(row for row in payload["cases"] if row["case_id"] == case_id)
    group = next(row for row in case["groups"] if row["group"] == group_id)
    return [from_wkb(bytes.fromhex(value)) for value in group["polygon_wkb_hex"]]


def test_required_construction_witnesses():
    cases = {name: raw for name, (raw, _) in topology_examples().items()}
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
        output, report = staged.construct(iter(raw))
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
    raw = focused_group(case_id, group_id)
    for candidate in (raw, [translate(polygon, 17.125, -9.875) for polygon in raw]):
        output, report = staged.construct(candidate)
        assert output is not None, report
        assert check_cleaning_contract(candidate, output, delta=0.5)["status"] == "pass"
        assert any(
            attempt.get("accepted_for_ranking")
            for attempt in report["fallback_attempts"]
        )


def test_dense_structural_proposal_repairs_one_cap_case_and_explains_the_residual():
    repaired, report = staged.construct(
        focused_group("city_grid:helsingborg:007", 42)
    )
    assert repaired is not None, report
    assert report["evaluations"] < 1000

    unresolved, report = staged.construct(
        focused_group("city_grid:linkoping:066", 45)
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


def test_identical_fixed_fallback_starts_are_evaluated_once():
    raw = focused_group("city_grid:norrkoping:003", 0)
    output, report = staged.construct(raw)
    assert output is not None, report
    assert report["fallback_attempts"][0]["accepted_for_ranking"] is True
    assert [row["reason"] for row in report["fallback_attempts"][1:]] == [
        "canonical_preprocessing_equivalence",
        "canonical_preprocessing_equivalence",
    ]
    assert report["evaluations"] < 100


def test_local_improvement_cannot_hide_new_distant_conflicts():
    raw = [Polygon([(0, 0), (0.1, 0.01), (100, 1), (100, 10), (0, 10)])]
    raw += [
        box(x, x * 0.01 - 0.49997 - 2, x + 2, x * 0.01 - 0.49997)
        for x in (10, 20, 30, 40, 50)
    ]
    occupied = unary_union(raw)
    candidate = staged.rebuild(occupied, {(0.1, 0.01): None})
    judge = staged.LocalJudge(FidelityBudget(raw, 0.25), 0.5)
    coordinates = np.array([(0, 0), (0.1, 0.01)])
    clip = judge.clip(judge.window(coordinates))
    assert staged.separation_guard(
        judge.state(occupied, clip), judge.state(candidate, clip)
    )
    assert judge.admits(candidate)
    accepted, reason = staged.repair_site(
        occupied,
        judge,
        coordinates,
        0.5,
        staged.separation_guard,
        lambda *args: iter([("remove_vertex", candidate)]),
    )
    assert accepted is None and reason == "global_progress"


def test_local_fidelity_clip_cache_is_invalidated_between_sites():
    raw = [box(0, 0, 10, 10)]
    judge = staged.LocalJudge(FidelityBudget(raw, 0.25), 0.5)
    candidate = box(0.4, 0, 10, 10)
    unaffected = box(8, -1, 11, 11)
    affected = box(-1, -1, 2, 11)

    assert judge.keeps_fidelity(candidate, unaffected)
    assert not judge.keeps_fidelity(candidate, affected)
    assert judge.keeps_fidelity(candidate, unaffected)


def test_rebuilder_cache_is_bounded_to_its_occupied_geometry():
    from dtcc_core.builder.cleaning import construction

    occupied = unary_union([box(0, 0, 2, 2), box(5, 0, 7, 2)])
    rebuilder = construction._Rebuilder(occupied)
    replacement = {(2.0, 0.0): (2.0, 0.25)}

    first = rebuilder.rebuild(replacement)
    assert rebuilder.rebuild(replacement) is first

    changed = construction._Rebuilder(first).rebuild({(7.0, 0.0): (7.0, 0.25)})
    fresh = staged.rebuild(first, {(7.0, 0.0): (7.0, 0.25)})
    assert changed.equals_exact(fresh, 0)


def test_collinear_points_do_not_pin_a_wall():
    other = Polygon([(10.4, 5), (15, 5), (15, 15), (11, 15)])
    simple, _ = staged.construct([box(0, 0, 10, 20), other])
    sampled, _ = staged.construct(
        [Polygon([(10, 0), (10, 10), (10, 20), (0, 20), (0, 0)]), other]
    )
    assert simple is not None and sampled is not None
    assert unary_union(simple).symmetric_difference(unary_union(sampled)).area < 1e-9


def test_zero_budget_is_not_silently_relaxed():
    raw = topology_examples()["point_contact"][0]
    output, report = staged.construct(raw, epsilon=0)
    assert output is None and report["outcome"] == "unresolved"
    assert len(report["sampling_attempts"]) == 1
    assert [row["reason"] for row in report["fallback_attempts"][1:]] == [
        "canonical_preprocessing_equivalence",
        "canonical_preprocessing_equivalence",
    ]


def test_retry_work_is_counted(monkeypatch):
    from dtcc_core.builder.cleaning import construction

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
    _, report = staged.construct(topology_examples()["point_contact"][0], epsilon=0)
    assert (
        report["edits"],
        report["sites"],
        report["evaluations"],
        report["admissions"],
        report["global_evaluations"],
    ) == (3, 6, 9, 12, 15)


def test_warning_candidate_never_outranks_a_conforming_fallback(monkeypatch):
    from dtcc_core.builder.cleaning import construction

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
    from dtcc_core.builder.cleaning import construction

    raw, cells = topology_examples()["point_contact"]
    rescue, _ = occupancy_search(raw, cells)
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
    output, report = staged.construct(raw)
    assert output is not None, report
    assert calls == [True, False]
    assert report["fallback_attempts"][1]["accepted_for_ranking"] is True


def test_reused_output_fails_before_work(tmp_path, monkeypatch, capsys):
    sentinel = tmp_path / "groups-delta0.5.jsonl"
    sentinel.write_text('{"epsilon":0.25}\n')
    monkeypatch.setattr(
        sys, "argv", ["probe", "--output", str(tmp_path), "--epsilon", "0"]
    )
    with pytest.raises(SystemExit) as error:
        staged.main()
    assert error.value.code == 2
    assert "fresh directory" in capsys.readouterr().err
    assert json.loads(sentinel.read_text())["epsilon"] == 0.25


@pytest.mark.parametrize(
    "delta,epsilon", [(0, 0.25), (float("nan"), 0.25), (0.5, -1), (0.5, float("inf"))]
)
def test_invalid_parameters(delta, epsilon):
    with pytest.raises(ValueError):
        staged.construct([], delta=delta, epsilon=epsilon)


def test_thin_passage_has_the_declared_width():
    for width in (1.0, 0.001):
        polygon = configurations(width)["thin_passage"][0]
        passage = box(3, 0, 10, 10).difference(polygon)
        assert passage.area == pytest.approx(7 * width)
        assert passage.bounds[3] - passage.bounds[1] == pytest.approx(width)


def test_cost_report_uses_actual_faces_without_cross_family_calibration():
    def row(family, width, **attempt):
        return dict(
            family=family,
            separation=width,
            attempts=[dict(requested_spacing=10, **attempt)],
        )

    rows = measured_costs(
        dict(
            rows=[
                row("near_walls", 1, outcome="meshed", faces=100),
                row("small_hole", 1, outcome="meshed", faces=30),
                row("near_walls", 0.1, outcome="meshed", faces=900),
                row("small_hole", 0.1, outcome="meshed", faces=32),
                row("thin_passage", 0.1, outcome="raised"),
            ]
        )
    )
    assert [r["extra_faces"] for r in rows] == [0, 0, 800, 2, None]
