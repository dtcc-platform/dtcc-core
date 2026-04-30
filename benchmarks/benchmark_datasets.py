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

WARNING_FAILURE_CLASSES = {
    "footprint_coverage",
    "footprint_cache",
    "lidar_coverage",
    "lidar_cache",
    "conditioned_footprint_warning",
    "mesh_quality_warning",
    "stage_contract_warning",
}

TERRAIN_ONLY_CONDITIONED_FOOTPRINTS_WARNING = (
    "No conditioned building footprints remain; downstream meshing will run terrain-only."
)
CONDITIONED_FOOTPRINT_WARNING_STAGES = {"contract", "conditioned_footprints"}
MESH_QUALITY_WARNING_STAGES = {
    "ground_mesh",
    "surface_shell",
    "plc",
    "volume_mesh",
}
WARNING_CLASS_PRIORITY = (
    "mesh_quality_warning",
    "conditioned_footprint_warning",
    "stage_contract_warning",
)


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


def _contract_warning_messages(
    *,
    label: str,
    contract: Any,
) -> list[dict[str, Any]]:
    if not isinstance(contract, dict):
        return []

    warnings = [str(message) for message in contract.get("warnings", [])]
    status = str(contract.get("status", "") or "")
    warning_statuses = {"warn", "warning"}
    if not warnings and status.lower() not in warning_statuses:
        return []
    if not warnings:
        warnings = [f"{label} contract reported warning status."]

    return [
        {
            "stage": label,
            "status": status or "warn",
            "message": message,
        }
        for message in warnings
    ]


def _metric_contract_warnings(metrics: dict[str, Any]) -> list[dict[str, Any]]:
    warnings: list[dict[str, Any]] = []
    warnings.extend(
        _contract_warning_messages(
            label="contract",
            contract=metrics.get("contract"),
        )
    )

    stage_audit = metrics.get("stage_audit")
    if not isinstance(stage_audit, dict):
        return warnings

    attempts = stage_audit.get("attempts")
    if not isinstance(attempts, list):
        return warnings

    selected_attempt_index = stage_audit.get("selected_attempt_index")
    if (
        isinstance(selected_attempt_index, int)
        and 0 <= selected_attempt_index < len(attempts)
    ):
        attempts_to_check = [attempts[selected_attempt_index]]
    else:
        attempts_to_check = attempts

    for attempt in attempts_to_check:
        if not isinstance(attempt, dict):
            continue
        stages = attempt.get("stages")
        if not isinstance(stages, dict):
            continue
        for stage_name, stage in stages.items():
            if not isinstance(stage, dict):
                continue
            warnings.extend(
                _contract_warning_messages(
                    label=str(stage_name),
                    contract=stage.get("contract"),
                )
            )

    return warnings


def _is_informational_stage_warning(warning: dict[str, Any]) -> bool:
    stage = str(warning.get("stage", "") or "")
    message = str(warning.get("message", "") or "")
    return (
        stage in CONDITIONED_FOOTPRINT_WARNING_STAGES
        and message == TERRAIN_ONLY_CONDITIONED_FOOTPRINTS_WARNING
    )


def _split_contract_warnings(
    warnings: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    status_warnings: list[dict[str, Any]] = []
    informational_warnings: list[dict[str, Any]] = []
    for warning in warnings:
        if _is_informational_stage_warning(warning):
            informational = dict(warning)
            informational["classification"] = "terrain_only"
            informational_warnings.append(informational)
        else:
            status_warnings.append(warning)
    return status_warnings, informational_warnings


def _contract_warning_class(warning: dict[str, Any]) -> str:
    stage = str(warning.get("stage", "") or "")
    message = str(warning.get("message", "") or "").lower()
    if stage in MESH_QUALITY_WARNING_STAGES:
        return "mesh_quality_warning"
    if (
        stage in CONDITIONED_FOOTPRINT_WARNING_STAGES
        or "conditioned footprint" in message
    ):
        return "conditioned_footprint_warning"
    return "stage_contract_warning"


def _warning_classes(warnings: list[dict[str, Any]]) -> list[str]:
    classes = {_contract_warning_class(warning) for warning in warnings}
    return sorted(
        classes,
        key=lambda name: (
            WARNING_CLASS_PRIORITY.index(name)
            if name in WARNING_CLASS_PRIORITY
            else len(WARNING_CLASS_PRIORITY)
        ),
    )


def _warning_error_payload(
    warnings: list[dict[str, Any]],
) -> dict[str, Any] | None:
    if not warnings:
        return None

    first = warnings[0]
    warning_classes = _warning_classes(warnings)
    extra_count = len(warnings) - 1
    suffix = f" (+{extra_count} more)" if extra_count else ""
    return {
        "type": "BenchmarkWarning",
        "message": f"{first['stage']}: {first['message']}{suffix}",
        "failure_class": warning_classes[0],
        "severity": "warning",
        "warning_classes": warning_classes,
        "warnings": warnings,
    }


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
    if exc_type == "NoFootprintTilesError" or "no footprint tiles intersect" in text:
        return "footprint_coverage"
    if "footprint download failed for bounds" in text:
        return "footprint_coverage"
    if (
        "downloaded-gpkg" in text
        or ("gpkg" in text and "not found" in text)
        or "footprint tile download did not produce" in text
    ):
        return "footprint_cache"
    if (
        exc_type == "FootprintDownloadError"
        or "footprint tile lookup failed" in text
        or "failed to download footprint tile" in text
    ):
        return "footprint_download"
    if "lidar" in text and (
        "404" in text
        or "not found" in text
        or "no lidar data" in text
        or "no lidar tiles intersect" in text
    ):
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
    if (
        exc_type == "LazrsError"
        or "lazrserror" in text
        or "failed to fill whole buffer" in text
    ):
        return "lidar_cache"
    if exc_type == "TimeoutError" or "timeouterror" in text:
        return "timeout"
    if "conditioned footprint" in text and "contract" in text:
        return "conditioned_footprint_contract"
    if "point outside" in text and "domain" in text:
        return "point_outside_domain"
    if exc_type == "IndexError" and "terrain" in text:
        return "terrain_lookup"
    return "pipeline"


def result_status_for_failure(failure_class: str) -> str:
    """Return the benchmark result status for a classified failure."""
    if failure_class in WARNING_FAILURE_CLASSES:
        return "warning"
    return "failed"


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
        contract_warnings = _metric_contract_warnings(metrics)
        status_warnings, informational_warnings = _split_contract_warnings(
            contract_warnings
        )
        if status_warnings or informational_warnings:
            metrics = dict(metrics)
        if status_warnings:
            metrics["stage_contract_warnings"] = status_warnings
        if informational_warnings:
            metrics["informational_stage_warnings"] = informational_warnings
        warning_error = _warning_error_payload(status_warnings)
        artifacts = {}
        artifact_dir = task.get("artifact_dir")
        if artifact_dir:
            artifacts = _save_result_artifacts(dataset_name, result, Path(artifact_dir))

        return {
            "task_id": task["id"],
            "dataset": dataset_name,
            "case": task["case"],
            "scenario": task["scenario"],
            "status": "warning" if warning_error is not None else "success",
            "elapsed_seconds": round(elapsed, 3),
            "bounds": bounds,
            "parameters": kwargs,
            "metrics": json_ready(metrics),
            "artifacts": artifacts,
            "error": warning_error,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - started
        exc_type = type(exc).__name__
        message = str(exc)
        failure_class = classify_failure(exc_type, message)
        status = result_status_for_failure(failure_class)
        return {
            "task_id": task["id"],
            "dataset": dataset_name,
            "case": task["case"],
            "scenario": task["scenario"],
            "status": status,
            "elapsed_seconds": round(elapsed, 3),
            "bounds": bounds,
            "parameters": kwargs,
            "metrics": {},
            "artifacts": {},
            "error": {
                "type": exc_type,
                "message": message,
                "failure_class": failure_class,
                "severity": "warning" if status == "warning" else "error",
                "traceback": traceback.format_exc(),
            },
        }


def result_to_line(result: dict[str, Any]) -> str:
    return "BENCHMARK_RESULT " + json.dumps(json_ready(result), sort_keys=True)
