from pathlib import Path

import numpy as np

from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface, Surface
from dtcc_core.builder.evaluation.runner import BuildingResult
from dtcc_core.builder.evaluation.report import FailureThresholds, write_failure_gallery


def _square(z: float = 0.0) -> Surface:
    return Surface(vertices=np.array([[0,0,z],[1,0,z],[1,1,z],[0,1,z]], dtype=float))


def _building_with_lod2(bid: str) -> Building:
    b = Building(id=bid)
    ms = MultiSurface()
    ms.surfaces.append(_square())
    b.add_geometry(ms, GeometryType.LOD2)
    b.attributes["roof_type"] = "FLAT"
    b.attributes["roof_confidence"] = 0.9
    return b


def _fake_failing_result(bid: str) -> BuildingResult:
    return BuildingResult(
        building_id=bid,
        stage_outcome="validation_failed",
        roof_type_predicted="FLAT",
        roof_type_truth="FLAT",
        roof_type_correct=True,
        confidence=0.9,
        classifier_source=None,
        fallback_reason_attr="validation_failed",
        plane_count=1,
        plane_count_delta=0,
        coverage_ratio=0.9,
        top2_area_share=1.0,
        symmetry_error_deg=None,
        watertight=False,
        ridge_height_error=None,
        eave_height_error=None,
        semantic_iou={},
        timings={"filter": 1.0},
        total_ms=5.0,
    )


def _fake_passing_result(bid: str) -> BuildingResult:
    return BuildingResult(
        building_id=bid,
        stage_outcome="success",
        roof_type_predicted="FLAT",
        roof_type_truth="FLAT",
        roof_type_correct=True,
        confidence=0.95,
        classifier_source=None,
        fallback_reason_attr=None,
        plane_count=1,
        plane_count_delta=0,
        coverage_ratio=0.98,
        top2_area_share=1.0,
        symmetry_error_deg=None,
        watertight=True,
        ridge_height_error=None,
        eave_height_error=None,
        semantic_iou={},
        timings={"filter": 1.0},
        total_ms=5.0,
    )


def test_gallery_writes_vtk_for_each_failing_building(tmp_path: Path):
    buildings = {"b1": _building_with_lod2("b1"), "b2": _building_with_lod2("b2")}
    results = [_fake_failing_result("b1"), _fake_failing_result("b2")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    assert len(paths) == 2
    assert (out_dir / "b1.vtk").exists()
    assert (out_dir / "b2.vtk").exists()


def test_gallery_skips_passing_buildings(tmp_path: Path):
    buildings = {"b1": _building_with_lod2("b1")}
    results = [_fake_passing_result("b1")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    assert paths == []


def test_gallery_skips_buildings_not_in_map(tmp_path: Path):
    buildings = {"b1": _building_with_lod2("b1")}
    results = [_fake_failing_result("b1"), _fake_failing_result("ghost")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    assert len(paths) == 1
    assert (out_dir / "b1.vtk").exists()


def test_gallery_skips_building_with_no_lod2(tmp_path: Path):
    # Failing result references an id whose Building has no LoD2 geometry —
    # should be silently skipped (with a logged warning) and not written.
    b_no_lod2 = Building(id="empty")
    # deliberately no add_geometry(..., LOD2)
    buildings = {"empty": b_no_lod2, "b1": _building_with_lod2("b1")}
    results = [_fake_failing_result("empty"), _fake_failing_result("b1")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    # Only b1 should be exported; 'empty' is dropped.
    assert len(paths) == 1
    assert (out_dir / "b1.vtk").exists()
    assert not (out_dir / "empty.vtk").exists()
