"""Acceptance checks for the benchmark commands and cleaning/meshing boundary."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import box

from benchmarks import benchmark_datasets, benchmark_phases
from benchmarks.benchmark_catalog import build_tasks, task_to_json
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.model import Building, City, GeometryType, Surface

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "benchmarks" / "bench"


def cli(*args):
    return subprocess.run(
        [sys.executable, str(BENCH), *map(str, args)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )


@pytest.fixture
def task(tmp_path):
    value = task_to_json(build_tasks("quick", city="lund", phase="cleaning")[0])
    value.update(phase="cleaning", task_dir=str(tmp_path / "tasks" / "clean"))
    return value


@pytest.fixture
def raw_city(monkeypatch):
    city = City()
    for i, polygon in enumerate(
        [box(386300, 6174680, 386310, 6174700), box(386310.2, 6174680, 386320, 6174700)]
    ):
        surface = Surface()
        surface.from_polygon(polygon, 10)
        building = Building(id=f"raw-{i}")
        building.add_geometry(surface, GeometryType.LOD0)
        building.estimated_height = 10
        city.add_building(building)
    monkeypatch.setattr(
        benchmark_phases, "prepare_footprint_city_from_bounds", lambda bounds: city
    )
    return city


@pytest.fixture
def cleaning_run(task, raw_city, tmp_path):
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] == "success", result["error"]
    benchmark_phases.write_json(
        tmp_path / "results.json",
        {
            "manifest": {"suite": "quick", "phase": "cleaning", "tasks": [task]},
            "results": [result],
        },
    )
    return tmp_path, task, result


def test_simple_cli_and_spatial_scopes():
    help_result = cli("--help")
    assert help_result.returncode == 0
    assert "_run-case" not in help_result.stdout
    assert "SUPPRESS" not in help_result.stdout
    assert cli("run", "quick").returncode == 2
    quick = json.loads(cli("quick", "--dry-run").stdout)
    assert quick["task_count"] == 20
    survey = json.loads(cli("survey", "--dry-run").stdout)
    assert survey["task_count"] == 2000
    assert {t["scenario"]["id"] for t in survey["tasks"]} == {"baseline"}
    cleaning = json.loads(
        cli("survey", "--city", "lund", "--phase", "cleaning", "--dry-run").stdout
    )
    assert cleaning["task_count"] == 100
    assert {t["dataset"] for t in cleaning["tasks"]} == {"city_footprints"}


def test_sweeps_are_relevant_to_phase_and_respect_volume_envelope():
    cleaning = build_tasks("sweep", city="lund", phase="cleaning")
    assert all(
        not set(t.scenario.parameters) & {"max_mesh_size", "raster_cell_size"}
        for t in cleaning
    )
    volume = build_tasks("sweep", city="lund", datasets=["city_volume_mesh"])
    assert min(t.parameters["max_mesh_size"] for t in volume) == 5
    selected = cli(
        "sweep",
        "--city",
        "lund",
        "--dataset",
        "city_surface_mesh",
        "--scenario",
        "max_mesh_size_2",
        "--dry-run",
    )
    assert selected.returncode == 0, selected.stderr
    assert json.loads(selected.stdout)["task_count"] == 1


def test_invalid_selections_fail_before_work(tmp_path):
    missing = cli("quick", "--phase", "meshing", "--dry-run")
    assert missing.returncode == 2
    assert "requires --input" in missing.stderr
    assert cli("quick", "--input", tmp_path, "--dry-run").returncode == 2
    assert (
        cli(
            "quick",
            "--phase",
            "cleaning",
            "--dataset",
            "city_surface_mesh",
            "--dry-run",
        ).returncode
        == 2
    )


def test_cleaning_is_independent_and_measures_fidelity(cleaning_run, raw_city):
    directory, task, result = cleaning_run
    assert len(raw_city.buildings) == 2
    metrics = result["metrics"]["cleaning"]
    assert metrics["output_count"] == 1
    assert metrics["added_area"] > 0
    assert metrics["invalid_output_count"] == 0
    assert metrics["unrepresented_source_count"] == 0
    assert "mesh_input" not in result["metrics"]
    saved = Path(result["artifacts"]["cleaning_input"]["path"])
    city, conditioned = benchmark_phases.load_cleaning(saved, task)
    assert [b.id for b in city.buildings] == ["raw-0", "raw-1"]
    assert conditioned.source_map == [[0, 1]]
    assert (saved.parent / "metrics.json").exists()


@pytest.mark.parametrize("dataset", ["city_flat_mesh", "city_surface_mesh"])
def test_meshing_replays_cleaned_geometry_without_cleaner(
    cleaning_run, monkeypatch, dataset
):
    directory, clean_task, clean_result = cleaning_run
    task = task_to_json(build_tasks("quick", city="lund", datasets=[dataset])[0])
    task.update(
        phase="meshing",
        task_dir=str(directory / dataset),
        cleaning_input=clean_result["artifacts"]["cleaning_input"]["path"],
    )
    original = json.loads(Path(task["cleaning_input"]).read_text())["features"]

    def forbidden(*args, **kwargs):
        raise AssertionError("Cleaner or footprint download called during meshing")

    monkeypatch.setattr(meshes, "condition_polygon_coverage", forbidden)
    monkeypatch.setattr(meshes, "build_conditioned_footprints", forbidden)
    monkeypatch.setattr(
        benchmark_phases, "prepare_footprint_city_from_bounds", forbidden
    )

    def terrain(bounds, *, buildings, **kwargs):
        city = City()
        city.add_buildings(buildings)
        from dtcc_core.model import Raster

        raster = Raster()
        raster.data = np.zeros((250, 250))
        raster.set_bounds(bounds)
        city.add_terrain(raster)
        return city

    monkeypatch.setattr(benchmark_phases, "prepare_city_from_bounds", terrain)
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] in {"success", "warning"}, result["error"]
    assert result["metrics"]["cleaning"]["status"] == "reused"
    assert result["metrics"]["meshing"]["num_faces"] > 0
    assert result["metrics"]["meshing"]["element_quality_p01"] > 0
    assert (
        json.loads(Path(result["artifacts"]["cleaning_input"]["path"]).read_text())[
            "features"
        ]
        == original
    )
    # Prepared surface input is reusable without new terrain/height acquisition.
    task["prepared_city"] = result["artifacts"]["prepared_city"]["path"]
    task["task_dir"] = str(directory / (dataset + "-replayed"))
    monkeypatch.setattr(benchmark_phases, "prepare_city_from_bounds", forbidden)
    replay = benchmark_datasets.run_dataset(task)
    assert replay["status"] in {"success", "warning"}, replay["error"]
    assert (
        replay["metrics"]["meshing"]["num_faces"]
        == result["metrics"]["meshing"]["num_faces"]
    )


def test_cli_saved_run_replay_report_compare_and_no_overwrite(cleaning_run):
    directory, task, result = cleaning_run
    output = directory / "replay"
    args = (
        "quick",
        "--city",
        "lund",
        "--phase",
        "meshing",
        "--dataset",
        "city_flat_mesh",
        "--input",
        directory,
        "--output",
        output,
    )
    completed = cli(*args)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = json.loads((output / "results.json").read_text())
    assert payload["results"][0]["metrics"]["cleaning"]["status"] == "reused"
    assert payload["manifest"]["elapsed_seconds"] > 0
    assert (output / "summary.md").exists()
    assert cli(*args).returncode == 2
    assert "meshing=success" in cli("report", output).stdout
    comparison = cli("compare", output, output)
    assert comparison.returncode == 0
    assert "meshing.quality.element_quality.min" in comparison.stdout
    rerun = cli(
        "rerun", output, "--task", payload["results"][0]["task_id"], "--dry-run"
    )
    assert rerun.returncode == 0, rerun.stderr
    assert json.loads(rerun.stdout)["task_count"] == 1


@pytest.mark.parametrize("corruption", ["crs", "sources", "bounds"])
def test_corrupt_cleaning_artifact_fails_clearly(cleaning_run, corruption):
    _, task, result = cleaning_run
    path = Path(result["artifacts"]["cleaning_input"]["path"])
    payload = json.loads(path.read_text())
    if corruption == "crs":
        payload["metadata"]["crs"] = "EPSG:4326"
    elif corruption == "sources":
        payload["features"][0]["properties"]["source_indices"] = [999]
    else:
        payload["metadata"]["bounds"][0] += 1
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="Invalid cleaning artifact"):
        benchmark_phases.load_cleaning(path, task)


def test_mesh_failure_retains_completed_cleaning(cleaning_run, monkeypatch):
    directory, task, _ = cleaning_run
    task = {
        **task,
        "dataset": "city_flat_mesh",
        "phase": "both",
        "task_dir": str(directory / "failure"),
    }

    def fail(*args, **kwargs):
        raise RuntimeError("mesh failure")

    monkeypatch.setattr(meshes, "build_city_flat_mesh", fail)
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] == "failed"
    assert result["error"]["phase"] == "meshing"
    assert result["metrics"]["cleaning"]["status"] == "success"
    assert Path(result["artifacts"]["cleaning_input"]["path"]).exists()


def test_quality_failure_is_not_a_success(cleaning_run, monkeypatch):
    directory, task, _ = cleaning_run
    task = {
        **task,
        "dataset": "city_flat_mesh",
        "phase": "both",
        "task_dir": str(directory / "quality"),
    }
    monkeypatch.setattr(
        benchmark_datasets,
        "_mesh_metrics",
        lambda mesh: (_ for _ in ()).throw(ValueError("invalid mesh quality")),
    )
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] == "failed"
    assert result["metrics"]["cleaning"]["status"] == "success"
    assert result["error"]["phase"] == "meshing"


def test_optional_export_failure_does_not_erase_mesh_success(cleaning_run, monkeypatch):
    directory, task, _ = cleaning_run
    task = {
        **task,
        "dataset": "city_flat_mesh",
        "phase": "both",
        "task_dir": str(directory / "export"),
        "artifact_dir": str(directory / "artifacts"),
    }
    monkeypatch.setattr(
        benchmark_datasets,
        "_save_result_artifacts",
        lambda *args: (_ for _ in ()).throw(OSError("cannot export")),
    )
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] == "warning"
    assert result["error"]["failure_class"] == "artifact_export"
    assert result["metrics"]["meshing"]["status"] == "success"


def test_contract_warning_classification_retains_quality_and_terrain_only():
    warnings = [
        {"stage": "ground_mesh", "status": "warning", "message": "bad quality"},
        {
            "stage": "conditioned_footprints",
            "status": "warning",
            "message": benchmark_datasets.TERRAIN_ONLY_CONDITIONED_FOOTPRINTS_WARNING,
        },
    ]
    status, informational = benchmark_datasets._split_contract_warnings(warnings)
    assert len(status) == len(informational) == 1
    assert (
        benchmark_datasets._warning_error_payload(status)["failure_class"]
        == "mesh_quality_warning"
    )


def test_cleaning_warning_does_not_change_mesh_status(task, monkeypatch):
    def execute(task, metrics, artifacts):
        metrics.update(
            cleaning={
                "status": "success",
                "contract": {"status": "warn", "warnings": ["Residual clearance"]},
            },
            meshing={"status": "success"},
        )

    monkeypatch.setattr(benchmark_phases, "execute_phases", execute)
    result = benchmark_datasets.run_dataset(task)
    assert result["status"] == "warning"
    assert result["metrics"]["cleaning"]["status"] == "warning"
    assert result["metrics"]["meshing"]["status"] == "success"


def test_benchmark_failure_classification_separates_data_from_geometry() -> None:
    assert (
        benchmark_datasets.classify_failure(
            "ValueError",
            "ground_only=True requires one integer LAS classification per point",
        )
        == "lidar_input"
    )
    assert (
        benchmark_datasets.classify_failure(
            "RuntimeError",
            "Footprint download failed for bounds (1, 2, 3, 4).",
        )
        == "footprint_coverage"
    )
    assert (
        benchmark_datasets.classify_failure(
            "NoFootprintTilesError",
            "No footprint tiles intersect the requested bounding box.",
        )
        == "footprint_coverage"
    )
    assert (
        benchmark_datasets.classify_failure(
            "FootprintDownloadError",
            "Footprint tile download did not produce all expected files: tile.gpkg",
        )
        == "footprint_cache"
    )
    assert (
        benchmark_datasets.classify_failure(
            "FootprintDownloadError",
            "Footprint tile lookup failed for bounds (1, 2, 3, 4): connection refused",
        )
        == "footprint_download"
    )
    assert (
        benchmark_datasets.classify_failure(
            "FileNotFoundError",
            "File /Users/example/Library/Caches/dtcc-data/downloaded-gpkg/tile.gpkg not found",
        )
        == "footprint_cache"
    )
    assert (
        benchmark_datasets.classify_failure(
            "ConnectionError",
            "HTTPConnectionPool(host='compute.dtcc.chalmers.se', port=8000): "
            "Max retries exceeded with url: /get_lidar",
        )
        == "lidar_download"
    )
    assert (
        benchmark_datasets.classify_failure(
            "RuntimeError",
            'Request failed with status 404: {"detail":"No lidar tiles intersect the requested bounding box."}',
        )
        == "lidar_coverage"
    )
    assert (
        benchmark_datasets.classify_failure(
            "LazrsError",
            "IoError: failed to fill whole buffer",
        )
        == "lidar_cache"
    )
    assert (
        benchmark_datasets.classify_failure(
            "RuntimeError",
            "Conditioned footprints contract failed: short_edge_count=2",
        )
        == "conditioned_footprint_contract"
    )


def test_data_coverage_and_cache_failures_are_warning_statuses() -> None:
    assert (
        benchmark_datasets.result_status_for_failure("footprint_coverage") == "warning"
    )
    assert benchmark_datasets.result_status_for_failure("footprint_cache") == "warning"
    assert (
        benchmark_datasets.result_status_for_failure("footprint_download") == "failed"
    )
    assert benchmark_datasets.result_status_for_failure("lidar_coverage") == "warning"
    assert benchmark_datasets.result_status_for_failure("lidar_cache") == "warning"
    assert benchmark_datasets.result_status_for_failure("lidar_download") == "failed"
    assert (
        benchmark_datasets.result_status_for_failure("conditioned_footprint_warning")
        == "warning"
    )
    assert (
        benchmark_datasets.result_status_for_failure("mesh_quality_warning")
        == "warning"
    )
    assert (
        benchmark_datasets.result_status_for_failure("stage_contract_warning")
        == "warning"
    )
    assert benchmark_datasets.result_status_for_failure("pipeline") == "failed"
