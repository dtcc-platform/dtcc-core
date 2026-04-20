import csv
from pathlib import Path

from dtcc_core.builder.evaluation.runner import BuildingResult
from dtcc_core.builder.evaluation.report import write_per_building_csv


def _fake_result(bid="b001", outcome="success"):
    return BuildingResult(
        building_id=bid,
        stage_outcome=outcome,
        roof_type_predicted="FLAT",
        roof_type_truth="FLAT",
        roof_type_correct=True,
        confidence=0.9,
        classifier_source=None,
        fallback_reason_attr=None,
        plane_count=1,
        plane_count_delta=0,
        coverage_ratio=0.95,
        top2_area_share=1.0,
        symmetry_error_deg=None,
        watertight=True,
        ridge_height_error=None,
        eave_height_error=None,
        semantic_iou={"ROOF": 0.98, "GROUND": 0.99},
        timings={
            "normals": 1.0,
            "filter": 2.0,
            "ransac": 8.0,
            "region_growing": 3.0,
            "merge": 0.5,
            "classify": 1.5,
            "geometry": 2.0,
            "validate": 0.5,
        },
        total_ms=20.5,
    )


def test_csv_one_row_per_result(tmp_path: Path):
    out = tmp_path / "per_building.csv"
    write_per_building_csv(
        [_fake_result(), _fake_result(bid="b002", outcome="validation_failed")], out
    )
    with out.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 2
    assert rows[0]["building_id"] == "b001"
    assert rows[0]["stage_outcome"] == "success"
    assert rows[1]["stage_outcome"] == "validation_failed"


def test_csv_flattens_semantic_iou_columns(tmp_path: Path):
    out = tmp_path / "per_building.csv"
    write_per_building_csv([_fake_result()], out)
    with out.open() as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
    assert "iou_ROOF" in headers
    assert "iou_GROUND" in headers


def test_csv_flattens_stage_timings(tmp_path: Path):
    out = tmp_path / "per_building.csv"
    write_per_building_csv([_fake_result()], out)
    with out.open() as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
    assert "timing_filter_ms" in headers
    assert "timing_ransac_ms" in headers
    assert "timing_region_growing_ms" in headers
    assert "timing_normals_ms" in headers
