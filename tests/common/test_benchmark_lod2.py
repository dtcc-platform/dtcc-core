"""Offline LOD2 acceptance tests through the benchmark's actual task runner."""

from copy import deepcopy

import numpy as np
import pytest
from shapely.geometry import box, mapping

import dtcc_core
from benchmarks import benchmark_datasets, benchmark_phases
from benchmarks.benchmark_catalog import build_tasks, task_to_json
from dtcc_core.builder.geometry_builders.buildings import (
    compute_building_heights,
    extract_roof_points,
)
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.io.model import load_model, save_model
from dtcc_core.model import City, GeometryType, PointCloud, Surface
from tests.builder.test_roof_catalogue import building, city_of, model


def test_lod_selection_keeps_defaults_and_distinguishes_methods():
    default = build_tasks("quick")
    assert len(default) == 20
    assert all(t.parameters.get("lod", "lod1") == "lod1" for t in default)
    assert all(t.id.endswith(":baseline") for t in default)
    ids = set()
    for method in ("catalogue", "ransac"):
        tasks = build_tasks(
            "quick", datasets=["city_surface_mesh"], lod="lod2", lod2_method=method
        )
        assert len(tasks) == 10
        assert len({t.case.city for t in tasks}) == 10
        assert all(t.case.id.endswith(":056") for t in tasks)
        assert all(t.parameters["lod2_method"] == method for t in tasks)
        assert not ids.intersection(t.id for t in tasks)
        ids.update(t.id for t in tasks)
    with pytest.raises(ValueError, match="requires --lod lod2"):
        build_tasks("quick", lod2_method="ransac")
    with pytest.raises(ValueError, match="requires a city surface or volume"):
        build_tasks("quick", datasets=["city_flat_mesh"], lod="lod2")


def test_unmerged_cleaned_outline_does_not_change_roof_observations():
    original = box(0, 0, 10, 10)
    source = building(original, model("flat", [20.0], original))
    source.attributes["source_footprint"] = mapping(original)
    cleaned = Surface()
    cleaned.from_polygon(box(2, 0, 12, 10), 0.0)
    source.add_geometry(cleaned, GeometryType.LOD0)
    points = np.array([[0.25, 5, 20], [5, 5, 20], [10.75, 5, 50]], dtype=float)
    extract_roof_points(
        [source], PointCloud(points=points), statistical_outlier_remover=False
    )
    np.testing.assert_array_equal(source.point_cloud.points, points[:2])
    terrain = city_of([], bounds=(-5, -5, 20, 15), slope=10.0).terrain.raster
    assert terrain.get_value(5.0, 5.0) != terrain.get_value(7.0, 5.0)
    compute_building_heights([source], terrain, overwrite=True)
    assert source.attributes["ground_height"] == terrain.get_value(5.0, 5.0)
    assert source.lod0.zmax == pytest.approx(20.0)


@pytest.fixture
def roof_inputs(monkeypatch):
    """A hip/gable merge and an isolated shed; only data acquisition is replaced."""
    footprints = [box(-10, -5, 10, 5), box(10.2, -5, 20, 5), box(25, -5, 40, 5)]
    buildings = [
        building(p, model(shape, params, p, key=f"source-{i}"))
        for i, (p, shape, params) in enumerate(
            zip(
                footprints,
                ["hip", "gable", "shed"],
                [[0.0, 8.0, 0.6, 0.6], [0.0, 11.0, 0.4], [9.0, 0.15, -0.25]],
            )
        )
    ]
    city = city_of(buildings)
    roofs = np.vstack([b.point_cloud.points for b in buildings])
    x, y = np.meshgrid(np.linspace(-25, 45, 36), np.linspace(-20, 25, 24))
    ground = np.c_[x.ravel(), y.ravel(), np.zeros(x.size)]
    cloud = PointCloud(
        points=np.vstack([ground, roofs]),
        classification=np.r_[np.full(len(ground), 2), np.full(len(roofs), 6)].astype(
            np.uint8
        ),
    )
    footprints_only = deepcopy(buildings)
    for source in footprints_only:
        source.remove_geometry(GeometryType.POINT_CLOUD)
    monkeypatch.setattr(
        dtcc_core.io.data,
        "download_footprints",
        lambda **kwargs: deepcopy(footprints_only),
    )
    monkeypatch.setattr(
        dtcc_core.io.data, "download_pointcloud", lambda **kwargs: deepcopy(cloud)
    )
    return city


def test_failed_meshing_preserves_audit_and_inputs(roof_inputs, tmp_path, monkeypatch):
    task = task_to_json(
        build_tasks("quick", city="lund", datasets=["city_surface_mesh"], lod="lod2")[0]
    )
    task["case"]["bounds"] = [-25.0, -20.0, 45.0, 25.0]
    task["task_dir"] = str(tmp_path)

    def fail(city, *, stage_audit, **kwargs):
        stage_audit["roof_catalogue"] = {"rejection": "test failure"}
        raise ValueError("test failure")

    monkeypatch.setattr(benchmark_phases.meshes, "build_city_surface_mesh", fail)
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] == "failed"
    assert result["metrics"]["meshing"]["stage_audit"]["roof_catalogue"] == {
        "rejection": "test failure"
    }
    assert result["metrics"]["mesh_input"]["roof_evidence"]["sources_with_points"] == 3
    assert (tmp_path / "meshing" / "input.json").exists()


@pytest.mark.parametrize("dataset", ["city_surface_mesh", "city_volume_mesh"])
def test_quick_lod2_runner_and_exact_evidence_replay(
    dataset, roof_inputs, tmp_path, monkeypatch
):
    if dataset == "city_volume_mesh" and not is_tetgen_available():
        pytest.skip("TetGen is not installed")
    task = task_to_json(
        build_tasks("quick", city="lund", datasets=[dataset], lod="lod2")[0]
    )
    task["case"]["bounds"] = [-25.0, -20.0, 45.0, 25.0]
    task["parameters"].update(max_mesh_size=8.0, domain_height=40.0)
    task["task_dir"] = str(tmp_path / "catalogue")
    first = benchmark_datasets.run_dataset(task)
    assert first["status"] in ("success", "warning"), first["error"]
    assert first["metrics"]["meshing"]["lod2"]["retained_nonflat_fits"] > 0
    evidence = first["metrics"]["mesh_input"]["roof_evidence"]
    assert evidence["source_count"] == evidence["sources_with_points"] == 3
    assert len({s["points_sha256"] for s in evidence["sources"]}) == 3
    saved = first["artifacts"]["prepared_city"]["path"]
    prepared = load_model(saved, expected_type=City)
    assert len(prepared.buildings) == 2
    assert sorted(len(b.building_parts) for b in prepared.buildings) == [0, 2]
    assert benchmark_phases.roof_evidence(prepared) == evidence

    def forbidden(*args, **kwargs):
        raise AssertionError("Saved replay must not acquire data or clean footprints")

    monkeypatch.setattr(benchmark_phases, "prepare_city_from_bounds", forbidden)
    monkeypatch.setattr(
        benchmark_phases, "prepare_footprint_city_from_bounds", forbidden
    )
    monkeypatch.setattr(City, "clean_footprints", forbidden)
    task.update(
        phase="meshing",
        task_dir=str(tmp_path / "ransac"),
        prepared_city=saved,
        cleaning_input=first["artifacts"]["cleaning_input"]["path"],
    )
    task["parameters"]["lod2_method"] = "ransac"
    second = benchmark_datasets.run_dataset(task)
    assert second["status"] in ("success", "warning"), second["error"]
    assert second["metrics"]["mesh_input"]["source"] == "saved"
    assert second["metrics"]["mesh_input"]["roof_evidence"] == evidence
    assert (
        second["metrics"]["cleaning"]["contract"]
        == first["metrics"]["cleaning"]["contract"]
    )
    for result in (first, second):
        mesh = result["metrics"]["meshing"]
        assert mesh["degenerate_cell_count"] == 0
        assert mesh["num_cells" if dataset == "city_volume_mesh" else "num_faces"] > 0
        assert mesh["lod2"]["retained_nonflat_fits"] > 0
        assert mesh["lod2"]["block_fallbacks"] == []

    # A parent cloud cannot replace an individual source, and changed saved
    # observations cannot silently invalidate a catalogue/RANSAC comparison.
    block = next(b for b in prepared.buildings if b.building_parts)
    part = block.building_parts[0]
    part.point_cloud.points[0, 2] += 1
    altered = tmp_path / "altered.dtcc"
    save_model(prepared, altered)
    task.update(prepared_city=str(altered), task_dir=str(tmp_path / "altered"))
    rejected = benchmark_datasets.run_dataset(task)
    assert rejected["status"] == "failed"
    assert "matching source roof evidence" in rejected["error"]["message"]
    block.add_geometry(part.point_cloud, GeometryType.POINT_CLOUD)
    part.remove_geometry(GeometryType.POINT_CLOUD)
    with pytest.raises(ValueError, match="Missing source roof evidence"):
        benchmark_phases.roof_evidence(prepared)


def test_raw_saved_input_replay_recleans_and_preserves_source_evidence(
    roof_inputs, tmp_path, monkeypatch
):
    task = task_to_json(
        build_tasks("quick", city="lund", datasets=["city_surface_mesh"], lod="lod2")[0]
    )
    task["case"]["bounds"] = [-25.0, -20.0, 45.0, 25.0]
    task["task_dir"] = str(tmp_path / "initial")
    first = benchmark_datasets.run_dataset(task)
    assert first["status"] in ("success", "warning"), first["error"]

    def forbidden(*args, **kwargs):
        raise AssertionError("Raw replay must not acquire provider data")

    monkeypatch.setattr(benchmark_phases, "prepare_city_from_bounds", forbidden)
    monkeypatch.setattr(
        benchmark_phases, "prepare_footprint_city_from_bounds", forbidden
    )
    original = City.clean_footprints
    calls = []

    def cleaned(city, *args, **kwargs):
        calls.append(len(city.buildings))
        return original(city, *args, **kwargs)

    monkeypatch.setattr(City, "clean_footprints", cleaned)
    task.update(
        phase="replay",
        task_dir=str(tmp_path / "replay"),
        cleaning_input=first["artifacts"]["cleaning_input"]["path"],
        prepared_city=first["artifacts"]["prepared_city"]["path"],
    )
    second = benchmark_datasets.run_dataset(task)
    assert second["status"] in ("success", "warning"), second["error"]
    assert calls == [3]
    assert (
        second["metrics"]["mesh_input"]["roof_evidence"]
        == first["metrics"]["mesh_input"]["roof_evidence"]
    )
    assert second["metrics"]["input"]["source"] == "saved raw"
