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


def test_runner_exposes_built_buildings_by_id():
    runner = Runner()
    runner.run(EvalDataset(FIXTURE_ROOT))
    assert "b001" in runner.buildings
    b = runner.buildings["b001"]
    assert b.id == "b001"
    assert b.lod0 is not None
    assert b.lod1 is not None


def test_runner_buildings_reset_between_runs():
    runner = Runner()
    runner.run(EvalDataset(FIXTURE_ROOT))
    first_ids = set(runner.buildings.keys())
    runner.run(EvalDataset(FIXTURE_ROOT))
    second_ids = set(runner.buildings.keys())
    assert first_ids == second_ids


def test_runner_default_wall_height_is_used_when_gt_omits_eave():
    import numpy as np
    from dtcc_core.model.geometry.pointcloud import PointCloud
    from dtcc_core.model.geometry.surface import Surface
    from dtcc_core.builder.evaluation.dataset import EvalSample, GroundTruth

    fp = Surface(vertices=np.array(
        [[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]], dtype=float
    ))
    pc = PointCloud(points=np.array([[1, 1, 5]], dtype=float))
    gt = GroundTruth(roof_type="FLAT")  # no eave_height

    sample = EvalSample(
        building_id="bx", footprint=fp, point_cloud=pc, ground_truth=gt,
    )

    runner = Runner(default_wall_height=7.5)
    building = runner._sample_to_building(sample)
    assert building.attributes["height"] == 7.5
