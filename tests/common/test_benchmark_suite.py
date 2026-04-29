from __future__ import annotations

import json
import os
import runpy
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

from shapely.geometry import Polygon

from benchmarks import benchmark_datasets
from benchmarks.benchmark_datasets import dataset_parameters
from benchmarks.benchmark_catalog import DATASET_NAMES, CITIES, build_tasks


def test_benchmark_entrypoint_is_executable() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    assert bench.is_file()
    assert os.access(bench, os.X_OK)


def test_benchmark_entrypoint_exposes_show_output_flag() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    completed = subprocess.run(
        [sys.executable, str(bench), "run", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "--show-output" in completed.stdout


def test_benchmark_entrypoint_exposes_save_artifacts_flag() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    completed = subprocess.run(
        [sys.executable, str(bench), "run", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "--save-artifacts" in completed.stdout


def test_benchmark_entrypoint_uses_dataset_flag() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    completed = subprocess.run(
        [sys.executable, str(bench), "run", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "--dataset" in completed.stdout


def test_benchmark_entrypoint_uses_scenario_filter() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    completed = subprocess.run(
        [
            sys.executable,
            str(bench),
            "run",
            "stress",
            "--city",
            "lund",
            "--dataset",
            "city_surface_mesh",
            "--scenario",
            "max_mesh_size_2",
            "--dry-run",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    manifest = json.loads(completed.stdout)
    assert manifest["scenarios"] == ["max_mesh_size_2"]
    assert manifest["task_count"] == 1
    assert manifest["tasks"][0]["dataset"] == "city_surface_mesh"
    assert manifest["tasks"][0]["scenario"]["id"] == "max_mesh_size_2"


def test_benchmark_list_scenarios_marks_default_parameters() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    completed = subprocess.run(
        [sys.executable, str(bench), "list", "scenarios"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Default parameters:" in completed.stdout
    assert "raster_cell_size" in completed.stdout
    assert "max_mesh_size" in completed.stdout
    assert "baseline" in completed.stdout
    assert "default" in completed.stdout
    assert "raster_cell_size_0.5" in completed.stdout
    assert "overrides: raster_cell_size=0.5" in completed.stdout


def test_benchmark_datasets_use_dataset_names() -> None:
    assert DATASET_NAMES == (
        "city_footprints",
        "terrain_surface_mesh",
        "city_flat_mesh",
        "city_surface_mesh",
        "city_volume_mesh",
    )


def test_city_catalog_has_one_automatic_grid_model() -> None:
    assert "lund" in CITIES
    assert len(build_tasks("grid", city="malmo")) == 500


def test_grid_suite_runs_100_spatial_cases_across_all_datasets() -> None:
    tasks = build_tasks("grid", city="stockholm")
    assert len({task.case.id for task in tasks}) == 100
    assert len({task.dataset for task in tasks}) == len(DATASET_NAMES)
    assert len({task.scenario.id for task in tasks}) == 1
    assert len(tasks) == 100 * len(DATASET_NAMES)


def test_grid_suite_requires_explicit_city() -> None:
    try:
        build_tasks("grid")
    except ValueError as exc:
        assert "--city" in str(exc)
    else:
        raise AssertionError("grid suite should require an explicit city")


def test_benchmark_tasks_do_not_define_mesher_dimension() -> None:
    tasks = build_tasks("smoke")
    assert tasks

    for task in tasks:
        assert task.dataset in DATASET_NAMES
        assert "mesher" not in task.parameters
        assert "triangle" not in task.id
        assert "spade" not in task.id


def test_dataset_parameters_map_normalized_sweeps_to_dataset_arguments() -> None:
    parameters = {
        "raster_cell_size": 0.5,
        "max_mesh_size": 2.0,
        "outlier_threshold": 7.0,
        "mesher": "triangle",
    }

    terrain = dataset_parameters("terrain_surface_mesh", parameters)
    assert terrain["raster_resolution"] == 0.5
    assert terrain["mesh_resolution"] == 2.0
    assert terrain["remove_outlier_threshold"] == 7.0
    assert "mesher" not in terrain

    flat = dataset_parameters("city_flat_mesh", parameters)
    assert flat["raster_cell_size"] == 0.5
    assert flat["max_mesh_size"] == 2.0
    assert flat["outlier_threshold"] == 7.0
    assert "mesher" not in flat


def test_sweep_suite_keeps_parameter_sweeps_on_surface_mesh_dataset() -> None:
    tasks = build_tasks("sweep", city="lund")
    assert tasks
    assert {task.dataset for task in tasks} == {"city_surface_mesh"}
    assert {task.scenario.id for task in tasks} >= {
        "raster_cell_size_0.5",
        "max_mesh_size_1",
        "bbox_size_m_500",
    }


def test_survey_suite_covers_spatial_baseline_and_parameter_envelope() -> None:
    tasks = build_tasks("survey")
    assert len(tasks) == 2400
    assert {task.dataset for task in tasks} == {
        "city_flat_mesh",
        "city_surface_mesh",
    }

    baseline_grid_tasks = [
        task
        for task in tasks
        if task.case.kind == "city_grid" and task.scenario.id == "baseline"
    ]
    assert len(baseline_grid_tasks) == 10 * 100 * 2

    envelope_tasks = [
        task
        for task in tasks
        if task.case.kind == "city_grid" and task.scenario.id != "baseline"
    ]
    assert len(envelope_tasks) == 10 * 15 * 2
    assert {task.scenario.id for task in envelope_tasks} >= {
        "raster_cell_size_0.5",
        "max_mesh_size_1",
        "max_mesh_size_2",
        "min_building_detail_0.25",
        "min_building_area_1",
    }

    bbox_tasks = [task for task in tasks if task.case.kind == "city_center"]
    assert len(bbox_tasks) == 10 * 5 * 2
    assert {task.scenario.id for task in bbox_tasks} == {
        "bbox_size_m_50",
        "bbox_size_m_100",
        "bbox_size_m_200",
        "bbox_size_m_350",
        "bbox_size_m_500",
    }


def test_survey_suite_can_be_reduced_to_one_city() -> None:
    tasks = build_tasks("survey", city="lund")
    assert len(tasks) == 240
    assert {task.case.city for task in tasks} == {"lund"}


def test_stress_suite_uses_dataset_specific_mesh_size_limits() -> None:
    tasks = build_tasks("stress", city="lund")
    assert tasks

    scenarios_by_dataset = {}
    for task in tasks:
        scenarios_by_dataset.setdefault(task.dataset, set()).add(task.scenario.id)

    assert scenarios_by_dataset == {
        "city_surface_mesh": {
            "raster_cell_size_0.5",
            "max_mesh_size_1",
            "max_mesh_size_2",
        },
        "city_volume_mesh": {
            "raster_cell_size_0.5",
            "max_mesh_size_5",
        },
    }


def test_stress_suite_keeps_volume_mesh_size_at_or_above_five_meters() -> None:
    tasks = build_tasks("stress", city="lund", datasets=["city_volume_mesh"])
    assert tasks
    assert {task.dataset for task in tasks} == {"city_volume_mesh"}
    for task in tasks:
        max_mesh_size = task.parameters["max_mesh_size"]
        assert task.scenario.id == "raster_cell_size_0.5" or max_mesh_size >= 5.0


def test_benchmark_failure_classification_separates_data_from_geometry() -> None:
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
    assert benchmark_datasets.result_status_for_failure("footprint_coverage") == "warning"
    assert benchmark_datasets.result_status_for_failure("footprint_cache") == "warning"
    assert benchmark_datasets.result_status_for_failure("footprint_download") == "failed"
    assert benchmark_datasets.result_status_for_failure("lidar_coverage") == "warning"
    assert benchmark_datasets.result_status_for_failure("lidar_cache") == "warning"
    assert benchmark_datasets.result_status_for_failure("lidar_download") == "failed"
    assert benchmark_datasets.result_status_for_failure("pipeline") == "failed"


def test_run_dataset_promotes_stage_contract_warnings(monkeypatch) -> None:
    class FakeDataset:
        class ArgsModel:
            model_fields = {}

        def __call__(self, *, bounds):
            return SimpleNamespace(
                vertices=[],
                faces=[],
                cells=[],
                markers=[],
                stage_audit={
                    "selected_attempt_index": 0,
                    "attempts": [
                        {
                            "stages": {
                                "ground_mesh": {
                                    "contract": {
                                        "status": "warn",
                                        "warnings": ["Flat mesh has short edge tail."],
                                    }
                                }
                            }
                        }
                    ],
                },
            )

    monkeypatch.setattr(
        benchmark_datasets.dtcc.datasets,
        "city_flat_mesh",
        FakeDataset(),
    )

    result = benchmark_datasets.run_dataset(
        {
            "id": "task",
            "dataset": "city_flat_mesh",
            "case": {"bounds": [0, 0, 1, 1]},
            "scenario": {"id": "baseline"},
            "parameters": {},
        }
    )

    assert result["status"] == "warning"
    assert result["error"]["failure_class"] == "stage_contract_warning"
    assert result["error"]["severity"] == "warning"
    assert result["error"]["warnings"][0]["stage"] == "ground_mesh"


def test_summary_markdown_reports_status_counts_and_quality_metrics() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    bench_module = runpy.run_path(str(bench))
    summary_markdown = bench_module["_summary_markdown"]

    summary = summary_markdown(
        [
            {
                "task_id": "smoke:city_footprints:city_center:lund:500m:baseline",
                "dataset": "city_footprints",
                "case": {"id": "city_center:lund:500m"},
                "scenario": {"id": "baseline"},
                "status": "success",
                "elapsed_seconds": 1.25,
                "metrics": {
                    "footprint_count": 12,
                    "polygon_count": 12,
                    "source_map_max_size": 2,
                    "contract": {
                        "status": "passed",
                        "metrics": {
                            "min_clearance": 0.5,
                            "pair_issue_count": 0,
                        },
                    },
                },
            },
            {
                "task_id": "smoke:city_surface_mesh:city_center:malmo:500m:baseline",
                "dataset": "city_surface_mesh",
                "case": {"id": "city_center:malmo:500m"},
                "scenario": {"id": "baseline"},
                "status": "failed",
                "elapsed_seconds": 2.5,
                "metrics": {},
                "error": {
                    "type": "RuntimeError",
                    "message": "mesh failed",
                    "failure_class": "pipeline",
                },
                "stdout_log": "tasks/example/stdout.log",
                "stderr_log": "tasks/example/stderr.log",
            },
            {
                "task_id": "smoke:city_flat_mesh:city_center:helsingborg:500m:baseline",
                "dataset": "city_flat_mesh",
                "case": {"id": "city_center:helsingborg:500m"},
                "scenario": {"id": "baseline"},
                "status": "warning",
                "elapsed_seconds": 2.75,
                "metrics": {},
                "error": {
                    "type": "RuntimeError",
                    "message": "No lidar tiles intersect the requested bounding box.",
                    "failure_class": "lidar_coverage",
                    "severity": "warning",
                },
                "stdout_log": "tasks/warning/stdout.log",
                "stderr_log": "tasks/warning/stderr.log",
            },
            {
                "task_id": "smoke:terrain_surface_mesh:city_center:uppsala:500m:baseline",
                "dataset": "terrain_surface_mesh",
                "case": {"id": "city_center:uppsala:500m"},
                "scenario": {"id": "baseline"},
                "status": "success",
                "elapsed_seconds": 3.0,
                "metrics": {
                    "num_vertices": 1234,
                    "num_faces": 2000,
                    "num_cells": 0,
                    "building_faces": 0,
                    "quality": {
                        "element_quality": {"min": 0.5, "mean": 0.75},
                        "aspect_ratio": {"max": 2.0},
                        "skewness": {"max": 0.1},
                    },
                },
                "artifacts": {
                    "mesh": {
                        "format": "vtu",
                        "path": "tasks/example/artifacts/mesh.vtu",
                    }
                },
            },
        ],
        {
            "suite": "smoke",
            "plan": "4 cases x 1 dataset x 1 scenario = 4 tasks",
        },
    )

    assert "| Spatial cases | 2 | 1 | 1 | 4 |" in summary
    assert "| Execution tasks | 2 | 1 | 1 | 4 |" in summary
    assert "| Scope | ✓ Success | ⚠ Warning | ✗ Fail | Total |" in summary
    assert "## Failure Classes" in summary
    assert "| ✗ failed | pipeline | 1 |" in summary
    assert "| ⚠ warning | lidar_coverage | 1 |" in summary
    assert summary.index("## Results") < summary.index("## Status Summary")
    assert "✓ success" in summary
    assert "⚠ warning" in summary
    assert "✗ failed" in summary
    assert "pipeline / RuntimeError: mesh failed" in summary
    assert "lidar_coverage / RuntimeError: No lidar tiles intersect" in summary
    assert "footprints=12, polygons=12" in summary
    assert "contract=passed" in summary
    assert "V=1,234, F=2,000, C=0, building_faces=0" in summary
    assert "q_min=0.5, q_mean=0.75, aspect_max=2, skew_max=0.1" in summary
    assert "Artifacts" in summary
    assert "mesh: tasks/example/artifacts/mesh.vtu" in summary


def test_summary_markdown_reports_total_run_time() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    bench_module = runpy.run_path(str(bench))
    summary_markdown = bench_module["_summary_markdown"]

    summary = summary_markdown(
        [],
        {
            "suite": "survey",
            "plan": "0 cases x 0 datasets x 0 scenarios = 0 tasks",
            "started_at": "2026-04-29T10:00:00Z",
            "finished_at": "2026-04-29T11:01:05Z",
            "elapsed_seconds": 3665.432,
        },
    )

    assert "Started: 2026-04-29T10:00:00Z" in summary
    assert "Finished: 2026-04-29T11:01:05Z" in summary
    assert "Total time: 1h 01m 05.432s" in summary


def test_benchmark_run_persists_total_run_time(tmp_path, monkeypatch) -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    bench_module = runpy.run_path(str(bench))
    command_run = bench_module["command_run"]

    def fake_run_task_subprocess(
        task,
        run_dir,
        *,
        show_output=False,
        save_artifacts=False,
    ):
        return {
            "task_id": task["id"],
            "dataset": task["dataset"],
            "case": task["case"],
            "scenario": task["scenario"],
            "status": "success",
            "elapsed_seconds": 0.001,
            "bounds": task["case"]["bounds"],
            "parameters": task["parameters"],
            "metrics": {},
            "artifacts": {},
            "error": None,
        }

    monkeypatch.setitem(bench_module, "_run_task_subprocess", fake_run_task_subprocess)
    monkeypatch.setitem(bench_module, "_load_table_helpers", lambda: None)
    monkeypatch.setitem(bench_module, "_print_run_tables", lambda results: None)

    rc = command_run(
        SimpleNamespace(
            suite="smoke",
            city="lund",
            dataset=["city_footprints"],
            scenario=None,
            run_id="timed",
            output_dir=tmp_path,
            dry_run=False,
            show_output=False,
            save_artifacts=False,
        )
    )

    run_dir = tmp_path / "timed"
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    results_payload = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
    summary = (run_dir / "summary.md").read_text(encoding="utf-8")

    assert rc == 0
    assert manifest["started_at"].endswith("Z")
    assert manifest["finished_at"].endswith("Z")
    assert manifest["elapsed_seconds"] >= 0
    assert results_payload["manifest"]["elapsed_seconds"] == manifest["elapsed_seconds"]
    assert "Total time:" in summary


def test_live_status_labels_include_checkmarks_and_crosses() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    bench_module = runpy.run_path(str(bench))
    terminal_status_label = bench_module["_terminal_status_label"]

    assert "✓ success" in terminal_status_label("success")
    assert "⚠ warning" in terminal_status_label("warning")
    assert "✗ failed" in terminal_status_label("failed")


def test_warning_status_is_not_a_hard_benchmark_failure() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    bench_module = runpy.run_path(str(bench))
    is_hard_failure = bench_module["_is_hard_failure"]

    assert not is_hard_failure("success")
    assert not is_hard_failure("warning")
    assert is_hard_failure("failed")


def test_save_result_artifacts_writes_mesh_artifact(tmp_path, monkeypatch) -> None:
    saved_paths = []

    def fake_save_mesh(_mesh, path):
        saved_paths.append(Path(path))
        Path(path).write_text("mesh", encoding="utf-8")

    monkeypatch.setattr(
        benchmark_datasets.dtcc,
        "io",
        SimpleNamespace(save_mesh=fake_save_mesh),
        raising=False,
    )

    artifacts = benchmark_datasets._save_result_artifacts(
        "city_surface_mesh",
        object(),
        tmp_path,
    )

    assert saved_paths == [tmp_path / "mesh.vtu"]
    assert artifacts == {
        "mesh": {
            "format": "vtu",
            "path": str(tmp_path / "mesh.vtu"),
        }
    }


def test_save_result_artifacts_writes_footprint_geojson(tmp_path) -> None:
    footprints = SimpleNamespace(
        polygons=[Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])],
        source_map=[[3, 4]],
        subdomain_resolution=[2.5],
        diagnostics={"output_grid": 0.5},
        contract={"status": "passed"},
    )

    artifacts = benchmark_datasets._save_result_artifacts(
        "city_footprints",
        footprints,
        tmp_path,
    )

    artifact_path = tmp_path / "footprints.geojson"
    payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert artifacts["footprints"]["path"] == str(artifact_path)
    assert payload["type"] == "FeatureCollection"
    assert payload["features"][0]["properties"]["source_indices"] == [3, 4]
    assert payload["features"][0]["properties"]["subdomain_resolution"] == 2.5
