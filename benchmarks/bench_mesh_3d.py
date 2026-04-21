"""
Benchmark TetGen-backed 3D volume-mesh quality across central Stockholm tiles.

This script reuses the same city-preparation pipeline as the 2D survey, then
builds a tetrahedral volume mesh with dtcc-core's TetGen path for each tile.
It caches per-case results, saves per-case XDMF/HDF5 meshes, and writes
overview heatmaps for the key 3D metrics.

Typical usage:
    python benchmarks/bench_mesh_3d.py --cases 55
    python benchmarks/bench_mesh_3d.py
    python benchmarks/bench_mesh_3d.py --cases 55 --max-mesh-size 8
    python benchmarks/bench_mesh_3d.py --cases 55 --quality-ratio 1.4
    python benchmarks/bench_mesh_3d.py --cases 55 --size-only --preserve-surface
    python benchmarks/bench_mesh_3d.py --cases 55 --preserve-surface --max-added-points 10000
    python benchmarks/bench_mesh_3d.py --no-merge-buildings

Coordinate system: SWEREF99 TM (EPSG:3006)
"""

from __future__ import annotations

import argparse
import json
import re
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

import dtcc_core
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.model import GeometryType, VolumeMesh
from dtcc_core.model.mixins.mesh.quality import tet_element_quality

try:
    from _benchmark_conditioning import CONDITIONING_MODES, benchmark_conditioning_override
    from _stockholm_common import (
        DEFAULT_DELAY_BETWEEN_CASES,
        DEFAULT_MAX_MESH_SIZE,
        MIN_BUILDING_AREA,
        MIN_BUILDING_DETAIL,
        MIN_MESH_ANGLE,
        MERGE_BUILDINGS,
        NX,
        NY,
        add_cases_argument,
        annotate_heatmap,
        bounds_to_dict,
        case_to_grid,
        format_console_table,
        json_ready,
        load_plot_modules,
        load_results,
        load_results_bundle,
        make_bounds,
        mesh_size_label,
        positive_float,
        positive_int,
        prepare_city,
        quantile,
        resolve_case_numbers,
        save_results,
        slug_token,
        stockholm_output_dir,
    )
except ImportError:
    from benchmarks._benchmark_conditioning import (
        CONDITIONING_MODES,
        benchmark_conditioning_override,
    )
    from benchmarks._stockholm_common import (
        DEFAULT_DELAY_BETWEEN_CASES,
        DEFAULT_MAX_MESH_SIZE,
        MIN_BUILDING_AREA,
        MIN_BUILDING_DETAIL,
        MIN_MESH_ANGLE,
        MERGE_BUILDINGS,
        NX,
        NY,
        add_cases_argument,
        annotate_heatmap,
        bounds_to_dict,
        case_to_grid,
        format_console_table,
        json_ready,
        load_plot_modules,
        load_results,
        load_results_bundle,
        make_bounds,
        mesh_size_label,
        positive_float,
        positive_int,
        prepare_city,
        quantile,
        resolve_case_numbers,
        save_results,
        slug_token,
        stockholm_output_dir,
    )

# Configuration ---------------------------------------------------------------

DEFAULT_DOMAIN_HEIGHT = 100.0
DEFAULT_QUALITY_RATIO = 1.6
DEFAULT_LOD = GeometryType.LOD0
DEFAULT_MESHER = "auto"
DEFAULT_MAX_ADDED_POINTS = 50_000

SHORT_EDGE_THRESHOLDS = (0.5, 1.0)
LOW_QUALITY_THRESHOLDS = (0.05, 0.10)

OVERVIEW_METRICS = (
    ("element_quality_mean", "ElemQ mean", "RdYlGn", ".3f"),
    ("element_quality_min", "ElemQ min", "RdYlGn", ".3f"),
    ("aspect_ratio_max", "AR max", "RdYlGn_r", ".2f"),
    ("edge_length_p01", "Edge p01", "RdYlGn", ".2f"),
    ("volume_p01", "Vol p01", "RdYlGn", ".2f"),
)

REQUIRED_METRIC_KEYS = (
    "num_cells",
    "element_quality_min",
    "element_quality_mean",
    "aspect_ratio_max",
    "edge_ratio_max",
    "edge_length_min",
    "edge_length_p01",
    "edge_length_p05",
    "volume_min",
    "volume_p01",
    "volume_p05",
    "short_edges_lt_0_5_count",
    "short_edges_lt_1_0_count",
    "low_quality_lt_0_05_count",
    "low_quality_lt_0_10_count",
)


# Helpers --------------------------------------------------------------------


def parse_lod(value: str) -> GeometryType:
    normalized = value.strip().upper()
    try:
        return GeometryType[normalized]
    except KeyError as exc:
        valid = ", ".join(member.name for member in GeometryType)
        raise argparse.ArgumentTypeError(f"lod must be one of: {valid}") from exc


def regular_tet_volume(max_edge_length: float) -> float:
    # Volume of a regular tetrahedron with edge length h.
    return (np.sqrt(2.0) / 12.0) * float(max_edge_length) ** 3


def config_slug(
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> str:
    return ".".join(
        [
            f"maxh-{slug_token(max_mesh_size)}",
            f"mina-{slug_token(min_mesh_angle)}",
            f"height-{slug_token(domain_height)}",
            (
                f"q-{slug_token(quality_ratio)}"
                if quality_enabled
                else "size-only"
            ),
            "preserve" if preserve_surface else "split-surface",
            (
                f"S-{max_added_points}"
                if max_added_points is not None
                else "S-unlimited"
            ),
            lod.name.lower(),
            "merged" if merge_buildings else "split",
            mesher,
            f"pipeline-{pipeline_mode}",
            "audit" if stage_audit_enabled else "noaudit",
        ]
    )


def results_file_path(
    output_dir: Path,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> Path:
    del (
        max_mesh_size,
        min_mesh_angle,
        domain_height,
        quality_ratio,
        quality_enabled,
        preserve_surface,
        max_added_points,
        lod,
        merge_buildings,
        mesher,
        pipeline_mode,
        stage_audit_enabled,
    )
    return output_dir / "results.json"


def overview_plot_path(
    output_dir: Path,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> Path:
    del (
        max_mesh_size,
        min_mesh_angle,
        domain_height,
        quality_ratio,
        quality_enabled,
        preserve_surface,
        max_added_points,
        lod,
        merge_buildings,
        mesher,
        pipeline_mode,
        stage_audit_enabled,
    )
    return output_dir / "overview.png"


def summary_text_path(
    output_dir: Path,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> Path:
    del (
        max_mesh_size,
        min_mesh_angle,
        domain_height,
        quality_ratio,
        quality_enabled,
        preserve_surface,
        max_added_points,
        lod,
        merge_buildings,
        mesher,
        pipeline_mode,
        stage_audit_enabled,
    )
    return output_dir / "summary.txt"


def case_output_dir(
    output_dir: Path,
    number: int,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> Path:
    del (
        max_mesh_size,
        min_mesh_angle,
        domain_height,
        quality_ratio,
        quality_enabled,
        preserve_surface,
        max_added_points,
        lod,
        merge_buildings,
        mesher,
        pipeline_mode,
        stage_audit_enabled,
    )
    return output_dir / f"{number:03d}"


def case_mesh_path(
    output_dir: Path,
    number: int,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> Path:
    return (
        case_output_dir(
            output_dir,
            number,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            domain_height=domain_height,
            quality_ratio=quality_ratio,
            quality_enabled=quality_enabled,
            preserve_surface=preserve_surface,
            max_added_points=max_added_points,
            lod=lod,
            merge_buildings=merge_buildings,
            mesher=mesher,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        / "volume_mesh.xdmf"
    )


def case_tetgen_input_stem(
    number: int,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> str:
    del (
        number,
        max_mesh_size,
        min_mesh_angle,
        domain_height,
        quality_ratio,
        quality_enabled,
        preserve_surface,
        max_added_points,
        lod,
        merge_buildings,
        mesher,
        pipeline_mode,
        stage_audit_enabled,
    )
    return "tetgen_input"


def case_tetgen_input_paths(
    output_dir: Path,
    number: int,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> dict[str, Path]:
    case_dir = case_output_dir(
        output_dir,
        number,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=domain_height,
        quality_ratio=quality_ratio,
        quality_enabled=quality_enabled,
        preserve_surface=preserve_surface,
        max_added_points=max_added_points,
        lod=lod,
        merge_buildings=merge_buildings,
        mesher=mesher,
        pipeline_mode=pipeline_mode,
        stage_audit_enabled=stage_audit_enabled,
    )
    return {
        "ground": case_dir / "tetgen_input_ground.xdmf",
        "shell": case_dir / "tetgen_input_shell.xdmf",
        "plc": case_dir / "tetgen_input_plc.xdmf",
    }


def case_tetgen_quality_failure_report_path(
    output_dir: Path,
    number: int,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
    stage_audit_enabled: bool = False,
) -> Path:
    return (
        case_output_dir(
            output_dir,
            number,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            domain_height=domain_height,
            quality_ratio=quality_ratio,
            quality_enabled=quality_enabled,
            preserve_surface=preserve_surface,
            max_added_points=max_added_points,
            lod=lod,
            merge_buildings=merge_buildings,
            mesher=mesher,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        / "tetgen_quality_failure.json"
    )


def selected_stage_audit_attempt(stage_audit: dict[str, Any] | None) -> dict[str, Any] | None:
    if not stage_audit:
        return None

    attempts = stage_audit.get("attempts", [])
    if not isinstance(attempts, list) or not attempts:
        return None

    selected_index = stage_audit.get("selected_attempt_index")
    if isinstance(selected_index, int) and 0 <= selected_index < len(attempts):
        attempt = attempts[selected_index]
        return attempt if isinstance(attempt, dict) else None

    selected_label = stage_audit.get("selected_attempt_label")
    if selected_label is not None:
        for attempt in attempts:
            if isinstance(attempt, dict) and attempt.get("label") == selected_label:
                return attempt

    attempt = attempts[-1]
    return attempt if isinstance(attempt, dict) else None


def selected_stage_contracts(stage_audit: dict[str, Any] | None) -> dict[str, str]:
    attempt = selected_stage_audit_attempt(stage_audit)
    if attempt is None:
        return {}

    statuses: dict[str, str] = {}
    for stage_name, stage_data in attempt.get("stages", {}).items():
        if not isinstance(stage_data, dict):
            continue
        contract = stage_data.get("contract", {})
        if not isinstance(contract, dict):
            continue
        status = str(contract.get("status", "")).strip().lower()
        if status:
            statuses[str(stage_name)] = status
    return statuses


def format_stage_contracts(stage_audit: dict[str, Any] | None) -> str:
    statuses = selected_stage_contracts(stage_audit)
    if not statuses:
        return "-"

    preferred_order = ("conditioned_footprints", "ground_mesh", "surface_shell", "plc", "volume_mesh")
    ordered_stage_names = [
        *[name for name in preferred_order if name in statuses],
        *sorted(name for name in statuses if name not in preferred_order),
    ]
    return ", ".join(f"{name}={statuses[name]}" for name in ordered_stage_names)


def tetgen_switches(
    max_mesh_size: float,
    min_mesh_angle: float,
    quality_ratio: float,
    *,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
) -> dict[str, Any]:
    switches: dict[str, Any] = {
        "max_volume": regular_tet_volume(max_mesh_size),
        "quiet": True,
    }
    if quality_enabled:
        switches["quality"] = (quality_ratio, min_mesh_angle)
    if preserve_surface:
        switches["preserve_surface"] = True
    if max_added_points is not None:
        switches["max_added_points"] = int(max_added_points)
    return switches


def boundary_marker_histogram(mesh: VolumeMesh) -> dict[str, int]:
    markers = getattr(mesh, "boundary_markers", None)
    if markers is None or len(markers) == 0:
        return {}

    values, counts = np.unique(np.asarray(markers, dtype=np.int64), return_counts=True)
    return {str(int(value)): int(count) for value, count in zip(values, counts)}


def grading_metrics(mesh: VolumeMesh) -> dict[str, float | int]:
    cells = np.asarray(mesh.cells, dtype=np.int64)
    if cells.size == 0:
        return {
            "unique_edge_count": 0,
            "edge_length_min": float("nan"),
            "edge_length_p01": float("nan"),
            "edge_length_p05": float("nan"),
            "edge_length_mean": float("nan"),
            "volume_min": float("nan"),
            "volume_p01": float("nan"),
            "volume_p05": float("nan"),
            "volume_mean": float("nan"),
            "short_edges_lt_0_5_count": 0,
            "short_edges_lt_1_0_count": 0,
            "low_quality_lt_0_05_count": 0,
            "low_quality_lt_0_10_count": 0,
            "boundary_face_count": 0,
            "boundary_marker_count": 0,
        }

    vertices = np.asarray(mesh.vertices[:, :3], dtype=np.float64)
    tets = vertices[cells]
    v0 = tets[:, 0]
    v1 = tets[:, 1]
    v2 = tets[:, 2]
    v3 = tets[:, 3]

    triple_products = np.einsum("ij,ij->i", v1 - v0, np.cross(v2 - v0, v3 - v0))
    volumes = np.abs(triple_products) / 6.0

    edges = np.concatenate(
        [
            cells[:, [0, 1]],
            cells[:, [0, 2]],
            cells[:, [0, 3]],
            cells[:, [1, 2]],
            cells[:, [1, 3]],
            cells[:, [2, 3]],
        ],
        axis=0,
    )
    edges = np.sort(edges, axis=1)
    unique_edges = np.unique(edges, axis=0)
    edge_vectors = vertices[unique_edges[:, 1]] - vertices[unique_edges[:, 0]]
    edge_lengths = np.linalg.norm(edge_vectors, axis=1)

    qualities = tet_element_quality(vertices, cells)

    boundary_faces = getattr(mesh, "boundary_faces", None)
    boundary_markers = getattr(mesh, "boundary_markers", None)
    boundary_face_count = 0 if boundary_faces is None else int(len(boundary_faces))
    boundary_marker_count = 0 if boundary_markers is None else int(len(np.unique(boundary_markers)))

    return {
        "unique_edge_count": int(len(unique_edges)),
        "edge_length_min": float(edge_lengths.min()),
        "edge_length_p01": quantile(edge_lengths, 1.0),
        "edge_length_p05": quantile(edge_lengths, 5.0),
        "edge_length_mean": float(edge_lengths.mean()),
        "volume_min": float(volumes.min()),
        "volume_p01": quantile(volumes, 1.0),
        "volume_p05": quantile(volumes, 5.0),
        "volume_mean": float(volumes.mean()),
        "short_edges_lt_0_5_count": int(np.count_nonzero(edge_lengths < SHORT_EDGE_THRESHOLDS[0])),
        "short_edges_lt_1_0_count": int(np.count_nonzero(edge_lengths < SHORT_EDGE_THRESHOLDS[1])),
        "low_quality_lt_0_05_count": int(np.count_nonzero(qualities < LOW_QUALITY_THRESHOLDS[0])),
        "low_quality_lt_0_10_count": int(np.count_nonzero(qualities < LOW_QUALITY_THRESHOLDS[1])),
        "boundary_face_count": boundary_face_count,
        "boundary_marker_count": boundary_marker_count,
    }


def mesh_metrics(mesh: VolumeMesh, quality: dict[str, Any]) -> dict[str, float | int]:
    metrics = {
        "num_cells": int(quality["num_cells"]),
        "element_quality_min": float(quality["element_quality"]["min"]),
        "element_quality_mean": float(quality["element_quality"]["mean"]),
        "aspect_ratio_max": float(quality["aspect_ratio"]["max"]),
        "aspect_ratio_mean": float(quality["aspect_ratio"]["mean"]),
        "edge_ratio_max": float(quality["edge_ratio"]["max"]),
        "edge_ratio_mean": float(quality["edge_ratio"]["mean"]),
        "skewness_max": float(quality["skewness"]["max"]),
        "skewness_mean": float(quality["skewness"]["mean"]),
    }
    metrics.update(grading_metrics(mesh))
    return metrics


def error_details(exc: Exception) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def format_case_row(result: dict[str, Any]) -> str:
    if result.get("status") != "success":
        error = result.get("error", {})
        return (
            f"{'FAILED':<8}  "
            f"{error.get('type', 'Error')}: {error.get('message', 'unknown error')}"
        )

    metrics = result["metrics"]
    return (
        f"{'ok':<8}  "
        f"{metrics['num_cells']:>8}  "
        f"{metrics['element_quality_min']:>9.4f}  "
        f"{metrics['element_quality_mean']:>10.4f}  "
        f"{metrics['aspect_ratio_max']:>8.4f}  "
        f"{metrics['edge_ratio_max']:>8.4f}  "
        f"{metrics['skewness_max']:>8.4f}  "
        f"{result['time']:>7.2f}s"
    )


def format_case_scale_row(result: dict[str, Any]) -> str:
    if result.get("status") != "success":
        return f"{'FAILED':<8}"

    metrics = result["metrics"]
    return (
        f"{'ok':<8}  "
        f"{metrics['edge_length_min']:>9.4f}  "
        f"{metrics['edge_length_p01']:>9.4f}  "
        f"{metrics['edge_length_p05']:>9.4f}  "
        f"{metrics['volume_min']:>10.4g}  "
        f"{metrics['volume_p01']:>10.4f}  "
        f"{metrics['volume_p05']:>10.4f}  "
        f"{metrics['low_quality_lt_0_05_count']:>11}  "
        f"{metrics['low_quality_lt_0_10_count']:>11}"
    )


def print_case_summary(case_record: dict[str, Any]) -> None:
    number = case_record["number"]
    bounds = case_record["bounds"]
    result = case_record["result"]

    if result.get("status") == "success":
        metrics = result["metrics"]
        quality_rows = [[
            "ok",
            metrics["num_cells"],
            f"{metrics['element_quality_min']:.4f}",
            f"{metrics['element_quality_mean']:.4f}",
            f"{metrics['aspect_ratio_max']:.2f}",
            f"{metrics['edge_ratio_max']:.2f}",
            f"{metrics['skewness_max']:.2f}",
            f"{result['time']:.2f}s",
        ]]
        scale_rows = [[
            "ok",
            f"{metrics['edge_length_min']:.4f}",
            f"{metrics['edge_length_p01']:.4f}",
            f"{metrics['edge_length_p05']:.4f}",
            f"{metrics['volume_min']:.4g}",
            f"{metrics['volume_p01']:.4f}",
            f"{metrics['volume_p05']:.4f}",
            metrics["low_quality_lt_0_05_count"],
            metrics["low_quality_lt_0_10_count"],
        ]]
    else:
        quality_rows = [["failed", "-", "-", "-", "-", "-", "-", "-"]]
        scale_rows = [["failed", "-", "-", "-", "-", "-", "-", "-", "-"]]

    title = (
        f"Case {number}  ({bounds['xmin']:.0f}, {bounds['ymin']:.0f}) -> "
        f"({bounds['xmax']:.0f}, {bounds['ymax']:.0f})"
    )
    print()
    print(
        format_console_table(
            ["Status", "Cells", "EQ min", "EQ mean", "AR max", "ER max", "Skew max", "Time"],
            quality_rows,
            title=title,
        )
    )
    print()
    print(
        format_console_table(
            [
                "Status",
                "Edge min",
                "Edge p01",
                "Edge p05",
                "Vol min",
                "Vol p01",
                "Vol p05",
                "Q<0.05",
                "Q<0.10",
            ],
            scale_rows,
            title="Scale / grading",
        )
    )
    print(f"Preparation time: {case_record['prepare_time']:.2f}s")
    if result.get("file"):
        print(f"Volume mesh: {result['file']}")
    if result.get("tetgen_input"):
        tetgen_input = result["tetgen_input"]
        print(
            "TetGen input: "
            f"ground={tetgen_input.get('ground', '-')}, "
            f"shell={tetgen_input.get('shell', '-')}, "
            f"plc={tetgen_input.get('plc', '-')}"
        )
    if result.get("tetgen_switches"):
        print(f"TetGen switches: {result['tetgen_switches']}")
    if result.get("boundary_markers"):
        print(f"Boundary markers: {result['boundary_markers']}")
    if result.get("stage_audit"):
        stage_audit = result["stage_audit"]
        attempts = stage_audit.get("attempts", [])
        selected_label = stage_audit.get("selected_attempt_label", "-")
        print(f"Stage audit: {len(attempts)} attempt(s), selected={selected_label}")
        print(f"Selected contracts: {format_stage_contracts(stage_audit)}")


def print_survey_summary(results: dict[int, dict[str, Any]]) -> None:
    successes = [
        record["result"]
        for record in results.values()
        if record.get("result", {}).get("status") == "success"
    ]
    if not successes:
        print()
        print(
            format_console_table(
                [
                    "Success",
                    "EQ mean",
                    "EQ min",
                    "AR max",
                    "ER max",
                    "Edge p01",
                    "Vol p01",
                    "Q<0.05",
                    "Time mean",
                ],
                [[0, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a", "n/a", "n/a"]],
                title="TetGen summary",
            )
        )
        return

    eq_means = [item["metrics"]["element_quality_mean"] for item in successes]
    eq_mins = [item["metrics"]["element_quality_min"] for item in successes]
    ar_maxs = [item["metrics"]["aspect_ratio_max"] for item in successes]
    er_maxs = [item["metrics"]["edge_ratio_max"] for item in successes]
    edge_p01s = [item["metrics"]["edge_length_p01"] for item in successes]
    volume_p01s = [item["metrics"]["volume_p01"] for item in successes]
    low_quality = [item["metrics"]["low_quality_lt_0_05_count"] for item in successes]
    times = [item["time"] for item in successes]

    print()
    print(
        format_console_table(
            [
                "Success",
                "EQ mean",
                "EQ min",
                "AR max",
                "ER max",
                "Edge p01",
                "Vol p01",
                "Q<0.05",
                "Time mean",
            ],
            [[
                len(successes),
                f"{np.mean(eq_means):.4f}",
                f"{np.min(eq_mins):.4f}",
                f"{np.max(ar_maxs):.4f}",
                f"{np.max(er_maxs):.4f}",
                f"{np.mean(edge_p01s):.4f}",
                f"{np.mean(volume_p01s):.4f}",
                np.max(low_quality),
                f"{np.mean(times):.2f}s",
            ]],
            title="TetGen summary",
        )
    )


def case_complete(record: dict[str, Any]) -> bool:
    result = record.get("result", {})
    if result.get("status") != "success":
        return False
    metrics = result.get("metrics", {})
    return all(key in metrics for key in REQUIRED_METRIC_KEYS)


def build_grid(results: dict[int, dict[str, Any]], metric_key: str) -> np.ndarray:
    grid = np.full((NY, NX), np.nan)

    for number, record in results.items():
        result = record.get("result", {})
        if result.get("status") != "success":
            continue

        value = float(result["metrics"][metric_key])
        ix, iy = case_to_grid(number)
        grid[iy, ix] = value

    return grid


def plot_overview(
    results: dict[int, dict[str, Any]],
    output_path: Path,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    pipeline_mode: str,
) -> None:
    plt, _, _, _ = load_plot_modules()

    fig, axes = plt.subplots(
        1,
        len(OVERVIEW_METRICS),
        figsize=(5.4 * len(OVERVIEW_METRICS), 4.4),
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)

    for ax, (metric_key, title, cmap, fmt) in zip(axes, OVERVIEW_METRICS):
        grid = build_grid(results, metric_key)
        image = ax.imshow(grid, origin="lower", cmap=cmap, aspect="equal")
        annotate_heatmap(ax, grid, fmt)
        ax.set_title(title)
        ax.set_xticks(range(NX))
        ax.set_yticks(range(NY))
        ax.set_xlabel("Grid X")
        ax.set_ylabel("Grid Y")
        fig.colorbar(image, ax=ax, shrink=0.82, pad=0.03)

    fig.suptitle(
        "TetGen 3D mesh quality overview\n"
        f"{mesh_size_label(max_mesh_size)}, min angle={min_mesh_angle:g}, "
        f"height={domain_height:g}, "
        f"{'q=' + f'{quality_ratio:g}' if quality_enabled else 'size-only'}, "
        f"{'preserve' if preserve_surface else 'split-surface'}, "
        f"{'S=' + str(max_added_points) if max_added_points is not None else 'S=unlimited'}, "
        f"{lod.name.lower()}, mesher={mesher}, "
        f"merge={'on' if merge_buildings else 'off'}, "
        f"pipeline={pipeline_mode}"
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def compact_case_status(record: dict[str, Any]) -> str:
    result = record["result"]
    if result.get("status") != "success":
        error = result.get("error", {})
        failure = f"FAIL {error.get('type', 'Error')}: {error.get('message', 'unknown error')}"
        stage_contracts = format_stage_contracts(result.get("stage_audit"))
        if stage_contracts != "-":
            failure += f" [{stage_contracts}]"
        return failure

    metrics = result["metrics"]
    summary = (
        f"ARmax {metrics['aspect_ratio_max']:.2f}, "
        f"edge p01 {metrics['edge_length_p01']:.2f}m, "
        f"vol p01 {metrics['volume_p01']:.2f}, "
        f"Q<0.05 {metrics['low_quality_lt_0_05_count']}"
    )
    stage_contracts = format_stage_contracts(result.get("stage_audit"))
    if stage_contracts != "-":
        summary += f", contracts [{stage_contracts}]"
    return summary


def status_emoji(status: str | None) -> str:
    return "✅" if status == "success" else "❌"


def maybe_rename_case_artifact(source: Path, target: Path) -> Path:
    if not source.exists():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()
    source.rename(target)
    sibling_source = source.with_suffix(".h5")
    sibling_target = target.with_suffix(".h5")
    if sibling_source.exists():
        if sibling_target.exists():
            sibling_target.unlink()
        sibling_source.rename(sibling_target)
    return target


def maybe_rename_json_artifact(source: Path, target: Path) -> Path:
    if not source.exists():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()
    source.rename(target)
    return target


def _simplified_retry_label(label: str) -> str:
    return label.replace("-", "_")


def simplify_case_artifact_name(name: str) -> str:
    tetgen_input_match = re.fullmatch(
        r"tetgen_input(?:\.(retry-[^.]+))?\.tetgen-input-(ground|shell|plc)\.xdmf",
        name,
    )
    if tetgen_input_match:
        retry_label, kind = tetgen_input_match.groups()
        if retry_label is None:
            return f"tetgen_input_{kind}.xdmf"
        return f"tetgen_input_{_simplified_retry_label(retry_label)}_{kind}.xdmf"

    quality_match = re.fullmatch(
        r"tetgen_quality_failure(?:\.(retry-[^.]+))?\.tetgen-quality-failure\.json",
        name,
    )
    if quality_match:
        retry_label = quality_match.group(1)
        if retry_label is None:
            return "tetgen_quality_failure.json"
        return f"tetgen_quality_failure_{_simplified_retry_label(retry_label)}.json"

    return name


def simplify_case_artifacts(case_dir: Path) -> list[str]:
    renamed: list[str] = []
    for path in sorted(case_dir.glob("*.xdmf")):
        target_name = simplify_case_artifact_name(path.name)
        if target_name == path.name:
            continue
        renamed_path = maybe_rename_case_artifact(path, case_dir / target_name)
        renamed.append(renamed_path.name)
    for path in sorted(case_dir.glob("*.json")):
        if path.name == "summary.json":
            continue
        target_name = simplify_case_artifact_name(path.name)
        if target_name == path.name:
            continue
        renamed_path = maybe_rename_json_artifact(path, case_dir / target_name)
        renamed.append(renamed_path.name)
    return renamed


def aggregate_status_emoji(successes: int, failures: int) -> str:
    if failures == 0:
        return "✅"
    if successes == 0:
        return "❌"
    return "⚠️"


def build_summary_report(
    results: dict[int, dict[str, Any]],
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    conditioning_mode: str,
    pipeline_mode: str,
    stage_audit_enabled: bool,
    elapsed_seconds: float,
) -> str:
    detailed_rows: list[list[Any]] = []
    successes = 0
    failures = 0
    for number in sorted(results):
        result = results[number].get("result", {})
        if result.get("status") == "success":
            successes += 1
            metrics = result["metrics"]
            detailed_rows.append(
                [
                    f"{number:03d}",
                    f"{status_emoji('success')} ok",
                    metrics["num_cells"],
                    f"{metrics['element_quality_min']:.4f}",
                    f"{metrics['aspect_ratio_max']:.2f}",
                    f"{metrics['edge_length_p01']:.2f}",
                    f"{metrics['volume_p01']:.2f}",
                    metrics["low_quality_lt_0_05_count"],
                    f"{result['time']:.2f}s",
                ]
            )
        else:
            failures += 1
            error = result.get("error", {})
            detailed_rows.append(
                [
                    f"{number:03d}",
                    f"{status_emoji('failed')} fail",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    f"{error.get('type', '-')}: {error.get('message', '-')}",
                ]
            )

    successful_results = [
        record["result"]
        for record in results.values()
        if record.get("result", {}).get("status") == "success"
    ]
    if successful_results:
        summary_rows = [[
            aggregate_status_emoji(successes, failures),
            successes,
            failures,
            f"{min(r['metrics']['element_quality_min'] for r in successful_results):.4f}",
            f"{max(r['metrics']['aspect_ratio_max'] for r in successful_results):.2f}",
            f"{float(np.mean([r['metrics']['edge_length_p01'] for r in successful_results])):.2f}",
            f"{float(np.mean([r['metrics']['volume_p01'] for r in successful_results])):.2f}",
            int(max(r['metrics']['low_quality_lt_0_05_count'] for r in successful_results)),
        ]]
    else:
        summary_rows = [[aggregate_status_emoji(0, failures), 0, failures, "-", "-", "-", "-", "-"]]

    config_line = (
        f"{mesh_size_label(max_mesh_size)}, min angle={min_mesh_angle:g}, "
        f"height={domain_height:g}, "
        f"{'q=' + f'{quality_ratio:g}' if quality_enabled else 'size-only'}, "
        f"{'preserve' if preserve_surface else 'split-surface'}, "
        f"{'S=' + str(max_added_points) if max_added_points is not None else 'S=unlimited'}, "
        f"{lod.name.lower()}, mesher={mesher}, merge={'on' if merge_buildings else 'off'}, "
        f"conditioning={conditioning_mode}, "
        f"pipeline={pipeline_mode}, stage-audit={'on' if stage_audit_enabled else 'off'}"
    )

    return "\n".join(
        [
            "bench_mesh_3d",
            f"Config: {config_line}",
            f"Elapsed time: {elapsed_seconds:.1f}s",
            "",
            format_console_table(
                ["Status", "Success", "Fail", "Worst EQ", "Worst AR", "Mean edge p01", "Mean vol p01", "Max Q<0.05"],
                summary_rows,
                title="Summary",
            ),
            "",
            format_console_table(
                ["Case", "Status", "Cells", "EQ min", "AR max", "Edge p01", "Vol p01", "Q<0.05", "Time/Error"],
                detailed_rows,
                title="Detailed results",
            ),
        ]
    )


def run_case(
    number: int,
    output_dir: Path,
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    domain_height: float,
    quality_ratio: float,
    quality_enabled: bool,
    preserve_surface: bool,
    max_added_points: int | None,
    lod: GeometryType,
    merge_buildings: bool,
    mesher: str,
    conditioning_mode: str,
    pipeline_mode: str,
    save_tetgen_input: bool,
    stage_audit_enabled: bool,
) -> dict[str, Any]:
    ix, iy = case_to_grid(number)
    bounds = make_bounds(ix, iy)
    case_dir = case_output_dir(
        output_dir,
        number,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=domain_height,
        quality_ratio=quality_ratio,
        quality_enabled=quality_enabled,
        preserve_surface=preserve_surface,
        max_added_points=max_added_points,
        lod=lod,
        merge_buildings=merge_buildings,
        mesher=mesher,
        pipeline_mode=pipeline_mode,
        stage_audit_enabled=stage_audit_enabled,
    )
    case_dir.mkdir(parents=True, exist_ok=True)
    switches_params = tetgen_switches(
        max_mesh_size,
        min_mesh_angle,
        quality_ratio,
        quality_enabled=quality_enabled,
        preserve_surface=preserve_surface,
        max_added_points=max_added_points,
    )
    tetgen_input_paths = case_tetgen_input_paths(
        output_dir,
        number,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=domain_height,
        quality_ratio=quality_ratio,
        quality_enabled=quality_enabled,
        preserve_surface=preserve_surface,
        max_added_points=max_added_points,
        lod=lod,
        merge_buildings=merge_buildings,
        mesher=mesher,
        pipeline_mode=pipeline_mode,
        stage_audit_enabled=stage_audit_enabled,
    )
    tetgen_quality_failure_report = case_tetgen_quality_failure_report_path(
        output_dir,
        number,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=domain_height,
        quality_ratio=quality_ratio,
        quality_enabled=quality_enabled,
        preserve_surface=preserve_surface,
        max_added_points=max_added_points,
        lod=lod,
        merge_buildings=merge_buildings,
        mesher=mesher,
        pipeline_mode=pipeline_mode,
        stage_audit_enabled=stage_audit_enabled,
    )

    case_record: dict[str, Any] = {
        "number": number,
        "bounds": bounds_to_dict(bounds),
        "config": {
            "max_mesh_size": max_mesh_size,
            "min_mesh_angle": min_mesh_angle,
            "domain_height": domain_height,
            "quality_ratio": quality_ratio,
            "quality_enabled": quality_enabled,
            "preserve_surface": preserve_surface,
            "max_added_points": max_added_points,
            "lod": lod.name,
            "merge_buildings": merge_buildings,
            "mesher": mesher,
            "conditioning_mode": conditioning_mode,
            "pipeline_mode": pipeline_mode,
            "save_tetgen_input": save_tetgen_input,
            "stage_audit_enabled": stage_audit_enabled,
        },
    }

    start = time.perf_counter()
    try:
        city, prepare_time = prepare_city(bounds)
        case_record["prepare_time"] = round(prepare_time, 2)
    except Exception as exc:
        case_record["prepare_time"] = None
        case_record["result"] = {
            "status": "failed",
            "time": round(time.perf_counter() - start, 2),
            "tetgen_switches": json_ready(switches_params),
            "error": error_details(exc),
        }
        with (case_dir / "summary.json").open("w", encoding="utf-8") as handle:
            json.dump(json_ready(case_record), handle, indent=2)
        return json_ready(case_record)

    stage_audit: dict[str, Any] | None = {} if stage_audit_enabled else None
    try:
        with benchmark_conditioning_override(conditioning_mode):
            volume_mesh = dtcc_core.builder.build_city_volume_mesh(
                city,
                lod=lod,
                domain_height=domain_height,
                max_mesh_size=max_mesh_size,
                min_mesh_angle=min_mesh_angle,
                merge_buildings=merge_buildings,
                min_building_detail=MIN_BUILDING_DETAIL,
                min_building_area=MIN_BUILDING_AREA,
                smoothing=0,
                boundary_face_markers=True,
                tetgen_switches=switches_params,
                report_mesh_quality=False,
                mesher=mesher,
                tetgen_debug_output_dir=case_dir if save_tetgen_input else None,
                tetgen_debug_output_stem=(
                    case_tetgen_input_stem(
                        number,
                        max_mesh_size=max_mesh_size,
                        min_mesh_angle=min_mesh_angle,
                        domain_height=domain_height,
                        quality_ratio=quality_ratio,
                        quality_enabled=quality_enabled,
                        preserve_surface=preserve_surface,
                        max_added_points=max_added_points,
                        lod=lod,
                        merge_buildings=merge_buildings,
                        mesher=mesher,
                        pipeline_mode=pipeline_mode,
                        stage_audit_enabled=stage_audit_enabled,
                    )
                    if save_tetgen_input
                    else None
                ),
                tetgen_quality_failure_output_dir=case_dir,
                tetgen_quality_failure_output_stem="tetgen_quality_failure",
                pipeline_mode=pipeline_mode,
                stage_audit=stage_audit,
            )
        quality = json_ready(volume_mesh.quality())
        metrics = mesh_metrics(volume_mesh, quality)
        mesh_path = case_mesh_path(
            output_dir,
            number,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            domain_height=domain_height,
            quality_ratio=quality_ratio,
            quality_enabled=quality_enabled,
            preserve_surface=preserve_surface,
            max_added_points=max_added_points,
            lod=lod,
            merge_buildings=merge_buildings,
            mesher=mesher,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        volume_mesh.save(str(mesh_path))
        simplify_case_artifacts(case_dir)

        case_record["result"] = {
            "status": "success",
            "time": round(time.perf_counter() - start, 2),
            "file": mesh_path.name,
            "conditioning_mode": conditioning_mode,
            "tetgen_switches": json_ready(switches_params),
            "quality": quality,
            "metrics": metrics,
            "boundary_markers": boundary_marker_histogram(volume_mesh),
        }
        if stage_audit_enabled and stage_audit:
            case_record["result"]["stage_audit"] = json_ready(stage_audit)
        tetgen_inputs = {
            key: path.name for key, path in tetgen_input_paths.items() if path.exists()
        }
        if tetgen_inputs:
            case_record["result"]["tetgen_input"] = tetgen_inputs
        if tetgen_quality_failure_report.exists():
            case_record["result"]["tetgen_quality_failure_report"] = (
                tetgen_quality_failure_report.name
            )
    except Exception as exc:
        simplify_case_artifacts(case_dir)
        case_record["result"] = {
            "status": "failed",
            "time": round(time.perf_counter() - start, 2),
            "conditioning_mode": conditioning_mode,
            "tetgen_switches": json_ready(switches_params),
            "error": error_details(exc),
        }
        if stage_audit_enabled and stage_audit:
            case_record["result"]["stage_audit"] = json_ready(stage_audit)
        tetgen_inputs = {
            key: path.name for key, path in tetgen_input_paths.items() if path.exists()
        }
        if tetgen_inputs:
            case_record["result"]["tetgen_input"] = tetgen_inputs
        if tetgen_quality_failure_report.exists():
            case_record["result"]["tetgen_quality_failure_report"] = (
                tetgen_quality_failure_report.name
            )

    with (case_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(json_ready(case_record), handle, indent=2)

    return json_ready(case_record)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Survey TetGen 3D volume-mesh quality across Stockholm tiles."
    )
    add_cases_argument(parser)
    parser.add_argument(
        "--max-mesh-size",
        type=positive_float,
        default=DEFAULT_MAX_MESH_SIZE,
        help="Maximum target edge length in meters for the surface mesh and TetGen volume cap.",
    )
    parser.add_argument(
        "--min-mesh-angle",
        type=positive_float,
        default=MIN_MESH_ANGLE,
        help="Minimum angle / dihedral constraint passed through the 3D meshing pipeline.",
    )
    parser.add_argument(
        "--domain-height",
        type=positive_float,
        default=DEFAULT_DOMAIN_HEIGHT,
        help="Height of the volume domain above terrain in meters.",
    )
    parser.add_argument(
        "--quality-ratio",
        type=positive_float,
        default=DEFAULT_QUALITY_RATIO,
        help="TetGen radius-edge quality target used in the q switch.",
    )
    parser.add_argument(
        "--size-only",
        action="store_true",
        help="Disable TetGen's q-switch refinement and keep only the volume cap.",
    )
    parser.add_argument(
        "--preserve-surface",
        action="store_true",
        help="Pass TetGen's Y switch to preserve the input shell surface during refinement.",
    )
    parser.add_argument(
        "--max-added-points",
        type=positive_int,
        default=DEFAULT_MAX_ADDED_POINTS,
        help="TetGen S limit for added Steiner points.",
    )
    parser.add_argument(
        "--lod",
        type=parse_lod,
        default=DEFAULT_LOD,
        help="Building geometry level of detail. Example: LOD0 or LOD1.",
    )
    parser.add_argument(
        "--mesher",
        choices=("auto", "dtcc_mesher", "triangle", "spade"),
        default=DEFAULT_MESHER,
        help="2D backend for the intermediate ground/surface mesh stages.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=DEFAULT_DELAY_BETWEEN_CASES,
        help="Delay in seconds between cases during multi-case runs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=stockholm_output_dir("output_mesh_3d"),
        help="Directory for JSON, XDMF/HDF5, and PNG outputs.",
    )
    parser.add_argument(
        "--conditioning-mode",
        choices=CONDITIONING_MODES,
        default="new",
        help="Footprint conditioning benchmark mode used before volume meshing.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip overview PNG generation.",
    )
    parser.add_argument(
        "--save-tetgen-input",
        dest="save_tetgen_input",
        action="store_true",
        help="Save the flat ground mesh, surface shell, and PLC shell handed to TetGen.",
    )
    parser.add_argument(
        "--no-save-tetgen-input",
        dest="save_tetgen_input",
        action="store_false",
        help="Do not save intermediate TetGen input surface meshes.",
    )
    parser.add_argument(
        "--stage-audit",
        action="store_true",
        help="Record per-attempt stage metrics for conditioned footprints, meshes, PLC, and volume mesh.",
    )
    parser.add_argument(
        "--merge-buildings",
        dest="merge_buildings",
        action="store_true",
        help="Merge nearby/overlapping footprints before volume meshing.",
    )
    parser.add_argument(
        "--no-merge-buildings",
        dest="merge_buildings",
        action="store_false",
        help="Disable footprint merging before volume meshing.",
    )
    parser.set_defaults(merge_buildings=MERGE_BUILDINGS)
    parser.set_defaults(save_tetgen_input=None)
    args = parser.parse_args()
    cases_explicit = args.cases is not None
    args.cases = resolve_case_numbers(args.cases)
    args.cases_explicit = cases_explicit
    args.pipeline_mode = "strict"
    args.stage_audit = True

    if args.save_tetgen_input is None:
        args.save_tetgen_input = len(args.cases) == 1

    return args


def main() -> None:
    if not is_tetgen_available():
        raise RuntimeError(
            "TetGen is not available in the active environment. "
            "Install dtcc-tetgen-wrapper in the same venv first."
        )

    args = parse_args()
    benchmark_start = time.perf_counter()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "benchmark": "bench_mesh_3d",
        "config": {
            "max_mesh_size": args.max_mesh_size,
            "min_mesh_angle": args.min_mesh_angle,
            "domain_height": args.domain_height,
            "quality_ratio": None if args.size_only else args.quality_ratio,
            "quality_enabled": not args.size_only,
            "preserve_surface": args.preserve_surface,
            "max_added_points": args.max_added_points,
            "lod": args.lod.name,
            "merge_buildings": args.merge_buildings,
            "mesher": args.mesher,
            "conditioning_mode": args.conditioning_mode,
            "pipeline_mode": args.pipeline_mode,
            "save_tetgen_input": args.save_tetgen_input,
            "stage_audit_enabled": args.stage_audit,
        },
    }

    results_path = results_file_path(
        output_dir,
        max_mesh_size=args.max_mesh_size,
        min_mesh_angle=args.min_mesh_angle,
        domain_height=args.domain_height,
        quality_ratio=args.quality_ratio,
        quality_enabled=not args.size_only,
        preserve_surface=args.preserve_surface,
        max_added_points=args.max_added_points,
        lod=args.lod,
        merge_buildings=args.merge_buildings,
        mesher=args.mesher,
        pipeline_mode=args.pipeline_mode,
        stage_audit_enabled=args.stage_audit,
    )
    saved_metadata, results = load_results_bundle(results_path)
    if (
        saved_metadata is not None
        and (
            saved_metadata.get("benchmark") != metadata["benchmark"]
            or saved_metadata.get("config") != json_ready(metadata["config"])
        )
    ):
        print(
            f"Existing results at {results_path} use a different configuration; "
            "starting with a fresh cache."
        )
        print()
        results = {}
    total = NX * NY

    loaded_count = len(results)
    if loaded_count > 0:
        print(f"Loaded {loaded_count} previous case record(s) from {results_path}")
        print("  Delete that file to recompute everything from scratch.")
        print()

    print(
        f"Configuration: {mesh_size_label(args.max_mesh_size)}, "
        f"min angle={args.min_mesh_angle:g}, height={args.domain_height:g}, "
        f"{'q=' + f'{args.quality_ratio:g}' if not args.size_only else 'size-only'}, "
        f"{'preserve' if args.preserve_surface else 'split-surface'}, "
        f"{'S=' + str(args.max_added_points) if args.max_added_points is not None else 'S=unlimited'}, "
        f"{args.lod.name.lower()}, mesher={args.mesher}, "
        f"merge={'on' if args.merge_buildings else 'off'}, "
        f"conditioning={args.conditioning_mode}, "
        f"pipeline={args.pipeline_mode}, "
        f"tetgen-input={'on' if args.save_tetgen_input else 'off'}, "
        f"stage-audit={'on' if args.stage_audit else 'off'}"
    )
    print()

    if args.cases_explicit:
        cases_to_run = args.cases
        if len(cases_to_run) == 1 and cases_to_run[0] in results:
            print(f"Case {cases_to_run[0]} exists in cache and will be recomputed.")
            print()
    else:
        cases_to_run = [
            number for number in args.cases if not case_complete(results.get(number, {}))
        ]
        if len(args.cases) == total and not cases_to_run:
            print(f"All {total} cases are already complete.")
            print()
        else:
            print(
                f"{len(cases_to_run)} case(s) to compute: "
                f"{', '.join(str(case) for case in cases_to_run)}"
            )
            print()

    success_count = 0
    failure_count = 0
    start_all = time.perf_counter()

    for index, number in enumerate(cases_to_run):
        if index > 0 and args.delay > 0:
            time.sleep(args.delay)

        bounds = make_bounds(*case_to_grid(number))
        label = (
            f"[{index + 1}/{len(cases_to_run)}] case {number:3d}  "
            f"({bounds.xmin:.0f}, {bounds.ymin:.0f}) -> "
            f"({bounds.xmax:.0f}, {bounds.ymax:.0f})"
        )
        print(f"{label}  generating...", end="", flush=True)

        case_record = run_case(
            number,
            output_dir=output_dir,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            domain_height=args.domain_height,
            quality_ratio=args.quality_ratio,
            quality_enabled=not args.size_only,
            preserve_surface=args.preserve_surface,
            max_added_points=args.max_added_points,
            lod=args.lod,
            merge_buildings=args.merge_buildings,
            mesher=args.mesher,
            conditioning_mode=args.conditioning_mode,
            pipeline_mode=args.pipeline_mode,
            save_tetgen_input=args.save_tetgen_input,
            stage_audit_enabled=args.stage_audit,
        )
        results[number] = case_record
        save_results(results_path, results, metadata=metadata)

        if case_record["result"].get("status") == "success":
            success_count += 1
        else:
            failure_count += 1

        print(f"  {compact_case_status(case_record)}")

    elapsed = time.perf_counter() - start_all

    summary_results = {number: results[number] for number in args.cases if number in results}

    if not args.no_plots:
        overview_path = overview_plot_path(
            output_dir,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            domain_height=args.domain_height,
            quality_ratio=args.quality_ratio,
            quality_enabled=not args.size_only,
            preserve_surface=args.preserve_surface,
            max_added_points=args.max_added_points,
            lod=args.lod,
            merge_buildings=args.merge_buildings,
            mesher=args.mesher,
            pipeline_mode=args.pipeline_mode,
            stage_audit_enabled=args.stage_audit,
        )
        plot_overview(
            summary_results,
            overview_path,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            domain_height=args.domain_height,
            quality_ratio=args.quality_ratio,
            quality_enabled=not args.size_only,
            preserve_surface=args.preserve_surface,
            max_added_points=args.max_added_points,
            lod=args.lod,
            merge_buildings=args.merge_buildings,
            mesher=args.mesher,
            pipeline_mode=args.pipeline_mode,
        )

    total_elapsed = time.perf_counter() - benchmark_start
    save_results(results_path, results, metadata=metadata)
    report = build_summary_report(
        summary_results,
        max_mesh_size=args.max_mesh_size,
        min_mesh_angle=args.min_mesh_angle,
        domain_height=args.domain_height,
        quality_ratio=args.quality_ratio,
        quality_enabled=not args.size_only,
        preserve_surface=args.preserve_surface,
        max_added_points=args.max_added_points,
        lod=args.lod,
        merge_buildings=args.merge_buildings,
        mesher=args.mesher,
        conditioning_mode=args.conditioning_mode,
        pipeline_mode=args.pipeline_mode,
        stage_audit_enabled=args.stage_audit,
        elapsed_seconds=total_elapsed,
    )
    summary_path = summary_text_path(
        output_dir,
        max_mesh_size=args.max_mesh_size,
        min_mesh_angle=args.min_mesh_angle,
        domain_height=args.domain_height,
        quality_ratio=args.quality_ratio,
        quality_enabled=not args.size_only,
        preserve_surface=args.preserve_surface,
        max_added_points=args.max_added_points,
        lod=args.lod,
        merge_buildings=args.merge_buildings,
        mesher=args.mesher,
        pipeline_mode=args.pipeline_mode,
        stage_audit_enabled=args.stage_audit,
    )
    summary_path.write_text(report + "\n", encoding="utf-8")
    print()
    print(report)
    print()
    print(f"Results file: {results_path}")
    if not args.no_plots:
        print(f"Overview plot: {overview_path}")
    print(f"Summary text: {summary_path}")


if __name__ == "__main__":
    main()
