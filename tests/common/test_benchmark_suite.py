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
                "error": {"type": "RuntimeError", "message": "mesh failed"},
                "stdout_log": "tasks/example/stdout.log",
                "stderr_log": "tasks/example/stderr.log",
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
            "plan": "3 cases x 1 dataset x 1 scenario = 3 tasks",
        },
    )

    assert "| Spatial cases | 2 | 1 | 3 |" in summary
    assert "| Execution tasks | 2 | 1 | 3 |" in summary
    assert "✓ success" in summary
    assert "✗ failed" in summary
    assert "footprints=12, polygons=12" in summary
    assert "contract=passed" in summary
    assert "V=1,234, F=2,000, C=0, building_faces=0" in summary
    assert "q_min=0.5, q_mean=0.75, aspect_max=2, skew_max=0.1" in summary
    assert "Artifacts" in summary
    assert "mesh: tasks/example/artifacts/mesh.vtu" in summary


def test_live_status_labels_include_checkmarks_and_crosses() -> None:
    bench = Path(__file__).resolve().parents[2] / "benchmarks" / "bench"
    bench_module = runpy.run_path(str(bench))
    terminal_status_label = bench_module["_terminal_status_label"]

    assert "✓ success" in terminal_status_label("success")
    assert "✗ failed" in terminal_status_label("failed")


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
