from pathlib import Path

from dtcc_core.builder.evaluation.runner import Runner
from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.report import FailureThresholds, write_failure_gallery


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_gallery_exports_vtk_for_threshold_crossing(tmp_path: Path):
    # max_plane_count_delta=-1 forces any non-None delta to trip (abs always >= 0)
    thresholds = FailureThresholds(
        max_plane_count_delta=-1,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=False,
    )
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))

    gallery_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        dataset_root=FIXTURE_ROOT,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )
    assert len(paths) >= 1
    assert any(gallery_dir.iterdir())


def test_gallery_skips_passing_buildings(tmp_path: Path):
    # permissive thresholds — fixture should pass, gallery stays empty
    thresholds = FailureThresholds(
        max_plane_count_delta=100,
        min_coverage=0.0,
        min_semantic_iou=0.0,
        require_watertight=False,
    )
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    # Skip if fixture happens to fail by stage_outcome (it shouldn't on a flat roof)
    if any(r.stage_outcome != "success" for r in results):
        import pytest
        pytest.skip("fixture did not reach success; test assumption violated")

    gallery_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        dataset_root=FIXTURE_ROOT,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )
    assert paths == []
