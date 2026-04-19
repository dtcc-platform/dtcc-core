from pathlib import Path

from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.runner import Runner, BuildingResult


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_runner_produces_one_result_per_sample():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    assert len(results) == 1
    assert isinstance(results[0], BuildingResult)


def test_runner_result_has_expected_fields_for_flat_roof_fixture():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    r = results[0]
    assert r.building_id == "b001"
    assert r.stage_outcome in ("success", "validation_failed")
    assert r.roof_type_predicted is not None
    assert r.roof_type_truth == "FLAT"
    assert r.timings is not None
    assert "filter" in r.timings
    assert r.plane_count >= 0


def test_runner_records_total_wall_clock():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    r = results[0]
    assert r.total_ms >= 0.0
