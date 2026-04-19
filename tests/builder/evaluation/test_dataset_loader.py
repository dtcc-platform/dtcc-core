from pathlib import Path

from dtcc_core.builder.evaluation.dataset import EvalDataset, EvalSample


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_eval_dataset_discovers_buildings():
    ds = EvalDataset(FIXTURE_ROOT)
    ids = [sample.building_id for sample in ds]
    assert ids == ["b001"]


def test_eval_dataset_len_matches_iter():
    ds = EvalDataset(FIXTURE_ROOT)
    assert len(ds) == 1


def test_eval_dataset_loads_sample_fields():
    ds = EvalDataset(FIXTURE_ROOT)
    sample = next(iter(ds))
    assert isinstance(sample, EvalSample)
    assert sample.building_id == "b001"
    assert sample.footprint.vertices.shape[1] == 3
    assert sample.point_cloud.points.shape == (200, 3)
    assert sample.ground_truth.roof_type == "FLAT"
    assert sample.ground_truth.eave_height == 5.0
    assert sample.ground_truth.expected_plane_count == 1
