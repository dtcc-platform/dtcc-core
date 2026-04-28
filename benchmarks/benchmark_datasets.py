from __future__ import annotations

import json
import time
import traceback
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

import dtcc_core as dtcc

from benchmarks.benchmark_catalog import DATASET_NAMES


BLOCKED_PARAMETER_NAMES = {"bounds", "strict_live", "mesher"}

DATASET_PARAMETER_ALIASES: dict[str, dict[str, str]] = {
    "terrain_surface_mesh": {
        "raster_cell_size": "raster_resolution",
        "max_mesh_size": "mesh_resolution",
        "outlier_threshold": "remove_outlier_threshold",
    },
}

ARTIFACT_FORMATS: dict[str, str] = {
    "city_footprints": "geojson",
    "terrain_surface_mesh": "vtu",
    "city_flat_mesh": "vtu",
    "city_surface_mesh": "vtu",
    "city_volume_mesh": "xdmf",
}


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def dataset_parameters(dataset_name: str, parameters: dict[str, Any]) -> dict[str, Any]:
    dataset_fn = getattr(dtcc.datasets, dataset_name)
    fields = set(dataset_fn.ArgsModel.model_fields)
    aliases = DATASET_PARAMETER_ALIASES.get(dataset_name, {})
    kwargs: dict[str, Any] = {}
    for key, value in parameters.items():
        dataset_key = aliases.get(key, key)
        if dataset_key in fields and dataset_key not in BLOCKED_PARAMETER_NAMES:
            kwargs[dataset_key] = value
    return kwargs


def _marker_histogram(markers: Any) -> dict[str, int]:
    if markers is None:
        return {}
    array = np.asarray(markers)
    if array.size == 0:
        return {}
    return {str(key): int(value) for key, value in Counter(array.astype(int).tolist()).items()}


def _building_face_count(markers: Any) -> int | None:
    if markers is None:
        return None
    array = np.asarray(markers)
    if array.size == 0:
        return 0
    return int(np.count_nonzero(array >= 0))


def _mesh_metrics(mesh: Any) -> dict[str, Any]:
    metrics: dict[str, Any] = {
        "num_vertices": int(getattr(mesh, "num_vertices", 0) or len(getattr(mesh, "vertices", []))),
        "num_faces": int(getattr(mesh, "num_faces", 0) or len(getattr(mesh, "faces", []))),
        "num_cells": int(getattr(mesh, "num_cells", 0) or len(getattr(mesh, "cells", []))),
    }
    markers = getattr(mesh, "markers", None)
    boundary_markers = getattr(mesh, "boundary_markers", None)
    metrics["marker_histogram"] = _marker_histogram(markers)
    metrics["boundary_marker_histogram"] = _marker_histogram(boundary_markers)
    metrics["building_faces"] = _building_face_count(markers)
    if hasattr(mesh, "quality"):
        try:
            metrics["quality"] = json_ready(mesh.quality())
        except Exception as exc:
            metrics["quality_error"] = {
                "type": type(exc).__name__,
                "message": str(exc),
            }
    stage_audit = getattr(mesh, "stage_audit", None)
    if stage_audit is not None:
        metrics["stage_audit"] = json_ready(stage_audit)
    return metrics


def _footprint_metrics(footprints: Any) -> dict[str, Any]:
    source_map = getattr(footprints, "source_map", []) or []
    source_sizes = [len(indices) for indices in source_map]
    return {
        "footprint_count": int(len(getattr(footprints, "footprints", []) or [])),
        "polygon_count": int(len(getattr(footprints, "polygons", []) or [])),
        "source_map_count": int(len(source_map)),
        "source_map_max_size": int(max(source_sizes) if source_sizes else 0),
        "diagnostics": json_ready(getattr(footprints, "diagnostics", {}) or {}),
        "contract": json_ready(getattr(footprints, "contract", {}) or {}),
    }


def _write_footprint_geojson(footprints: Any, path: Path) -> None:
    from shapely.geometry import mapping

    source_map = getattr(footprints, "source_map", []) or []
    subdomain_resolution = getattr(footprints, "subdomain_resolution", []) or []
    features = []
    for index, polygon in enumerate(getattr(footprints, "polygons", []) or []):
        source_indices = source_map[index] if index < len(source_map) else []
        feature = {
            "type": "Feature",
            "geometry": mapping(polygon),
            "properties": {
                "id": index,
                "source_indices": list(source_indices),
                "source_count": len(source_indices),
            },
        }
        if index < len(subdomain_resolution):
            feature["properties"]["subdomain_resolution"] = subdomain_resolution[index]
        features.append(feature)

    payload = {
        "type": "FeatureCollection",
        "features": features,
        "metadata": {
            "diagnostics": json_ready(getattr(footprints, "diagnostics", {}) or {}),
            "contract": json_ready(getattr(footprints, "contract", {}) or {}),
        },
    }
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _save_result_artifacts(dataset_name: str, result: Any, artifact_dir: Path) -> dict[str, dict[str, str]]:
    artifact_dir.mkdir(parents=True, exist_ok=True)
    artifact_format = ARTIFACT_FORMATS[dataset_name]
    artifacts: dict[str, dict[str, str]] = {}

    if dataset_name == "city_footprints":
        path = artifact_dir / f"footprints.{artifact_format}"
        _write_footprint_geojson(result, path)
        artifacts["footprints"] = {"format": artifact_format, "path": str(path)}
        return artifacts

    stem = "volume_mesh" if dataset_name == "city_volume_mesh" else "mesh"
    path = artifact_dir / f"{stem}.{artifact_format}"
    dtcc.io.save_mesh(result, path)
    artifacts[stem] = {"format": artifact_format, "path": str(path)}

    if artifact_format == "xdmf":
        h5_path = path.with_suffix(".h5")
        if h5_path.exists():
            artifacts[f"{stem}_h5"] = {"format": "h5", "path": str(h5_path)}
    return artifacts


def classify_failure(exc_type: str, message: str) -> str:
    """Return a coarse benchmark failure class for triage summaries."""
    text = f"{exc_type}: {message}".lower()
    if (
        exc_type == "FootprintDownloadError"
        or "footprint download failed" in text
        or "footprint tile lookup failed" in text
        or "footprint tile download did not produce" in text
    ):
        return "footprint_download"
    if "downloaded-gpkg" in text or "gpkg" in text and "not found" in text:
        return "footprint_cache"
    if "lidar" in text and ("404" in text or "not found" in text):
        return "lidar_coverage"
    if (
        "lidar" in text
        and (
            "connectionerror" in text
            or "max retries exceeded" in text
            or "failed to establish" in text
            or "get_lidar" in text
        )
    ):
        return "lidar_download"
    if exc_type == "TimeoutError" or "timeouterror" in text:
        return "timeout"
    if "conditioned footprint" in text and "contract" in text:
        return "conditioned_footprint_contract"
    if "point outside" in text and "domain" in text:
        return "point_outside_domain"
    if exc_type == "IndexError" and "terrain" in text:
        return "terrain_lookup"
    return "pipeline"


def run_dataset(task: dict[str, Any]) -> dict[str, Any]:
    dataset_name = task["dataset"]
    if dataset_name not in DATASET_NAMES:
        raise ValueError(f"unknown benchmark dataset: {dataset_name}")

    bounds = task["case"]["bounds"]
    parameters = dict(task.get("parameters", {}))
    kwargs = dataset_parameters(dataset_name, parameters)

    started = time.perf_counter()
    try:
        dataset_fn = getattr(dtcc.datasets, dataset_name)
        result = dataset_fn(bounds=bounds, **kwargs)
        elapsed = time.perf_counter() - started

        if dataset_name == "city_footprints":
            metrics = _footprint_metrics(result)
        else:
            metrics = _mesh_metrics(result)
        artifacts = {}
        artifact_dir = task.get("artifact_dir")
        if artifact_dir:
            artifacts = _save_result_artifacts(dataset_name, result, Path(artifact_dir))

        return {
            "task_id": task["id"],
            "dataset": dataset_name,
            "case": task["case"],
            "scenario": task["scenario"],
            "status": "success",
            "elapsed_seconds": round(elapsed, 3),
            "bounds": bounds,
            "parameters": kwargs,
            "metrics": json_ready(metrics),
            "artifacts": artifacts,
            "error": None,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - started
        exc_type = type(exc).__name__
        message = str(exc)
        return {
            "task_id": task["id"],
            "dataset": dataset_name,
            "case": task["case"],
            "scenario": task["scenario"],
            "status": "failed",
            "elapsed_seconds": round(elapsed, 3),
            "bounds": bounds,
            "parameters": kwargs,
            "metrics": {},
            "artifacts": {},
            "error": {
                "type": exc_type,
                "message": message,
                "failure_class": classify_failure(exc_type, message),
                "traceback": traceback.format_exc(),
            },
        }


def result_to_line(result: dict[str, Any]) -> str:
    return "BENCHMARK_RESULT " + json.dumps(json_ready(result), sort_keys=True)
