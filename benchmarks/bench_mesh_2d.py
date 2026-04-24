"""
Benchmark flat-mesh quality across a benchmark area for one or more 2D meshers.

This script builds the same prepared City input once per tile, then runs one or
more mesh generators against that identical input. It caches per-case results,
saves per-mesher VTU files, and writes visual comparison plots.

Typical usage:
    python benchmarks/bench_mesh_2d.py --cases 55
    python benchmarks/bench_mesh_2d.py
    python benchmarks/bench_mesh_2d.py --cases 55 --meshers triangle dtcc_mesher spade
    python benchmarks/bench_mesh_2d.py --cases 55 --max-mesh-size none
    python benchmarks/bench_mesh_2d.py --per-case-plots

With no explicit meshers, the script compares all available meshers from:
    dtcc_mesher, triangle, spade

Coordinate system: SWEREF99 TM (EPSG:3006)
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

import dtcc_core
from dtcc_core.model import Bounds, City, Mesh
try:
    from _benchmark_conditioning import CONDITIONING_MODES, benchmark_conditioning_override
    from _stockholm_common import (
        DEFAULT_DELAY_BETWEEN_CASES,
        DEFAULT_MAX_MESH_SIZE,
        MERGE_BUILDINGS,
        MIN_BUILDING_AREA,
        MIN_BUILDING_DETAIL,
        MIN_MESH_ANGLE,
        active_area_label,
        add_area_argument,
        add_cases_argument,
        annotate_heatmap,
        benchmark_output_dir,
        bounds_to_dict,
        box_size,
        case_to_grid,
        format_console_table,
        grid_shape,
        json_ready,
        load_plot_modules,
        load_results,
        load_results_bundle,
        make_bounds,
        mesh_size_label,
        normalize_max_mesh_size,
        parse_max_mesh_size_argument,
        prepare_city,
        quantile,
        resolve_case_numbers,
        save_results,
        set_active_area,
        slug_token,
    )
except ImportError:
    from benchmarks._benchmark_conditioning import (
        CONDITIONING_MODES,
        benchmark_conditioning_override,
    )
    from benchmarks._stockholm_common import (
        DEFAULT_DELAY_BETWEEN_CASES,
        DEFAULT_MAX_MESH_SIZE,
        MERGE_BUILDINGS,
        MIN_BUILDING_AREA,
        MIN_BUILDING_DETAIL,
        MIN_MESH_ANGLE,
        active_area_label,
        add_area_argument,
        add_cases_argument,
        annotate_heatmap,
        benchmark_output_dir,
        bounds_to_dict,
        box_size,
        case_to_grid,
        format_console_table,
        grid_shape,
        json_ready,
        load_plot_modules,
        load_results,
        load_results_bundle,
        make_bounds,
        mesh_size_label,
        normalize_max_mesh_size,
        parse_max_mesh_size_argument,
        prepare_city,
        quantile,
        resolve_case_numbers,
        save_results,
        set_active_area,
        slug_token,
    )

# Configuration ---------------------------------------------------------------

DEFAULT_MESHER_ORDER = ("dtcc_mesher", "triangle", "spade")
OVERVIEW_METRICS = (
    ("element_quality_mean", "ElemQ mean", "RdYlGn", ".3f"),
    ("aspect_ratio_max", "AR max", "RdYlGn_r", ".2f"),
    ("edge_ratio_max", "ER max", "RdYlGn_r", ".2f"),
    ("edge_length_p01", "Edge p01", "RdYlGn", ".2f"),
    ("area_p01", "Area p01", "RdYlGn", ".2f"),
)

MARKER_CATEGORY_COLORS = ("#eceff4", "#f2cc8f", "#5b8fd1")
EDGE_FACE_LIMIT = 40_000
SHORT_EDGE_THRESHOLDS = (0.5, 1.0)
REQUIRED_METRIC_KEYS = (
    "num_cells",
    "element_quality_min",
    "element_quality_mean",
    "aspect_ratio_max",
    "edge_ratio_max",
    "edge_length_min",
    "edge_length_p01",
    "edge_length_p05",
    "area_min",
    "area_p01",
    "area_p05",
    "short_edges_lt_0_5_count",
    "short_edges_lt_1_0_count",
)


# Helpers --------------------------------------------------------------------


def meshers_slug(meshers: list[str]) -> str:
    ordered = [mesher for mesher in DEFAULT_MESHER_ORDER if mesher in meshers]
    extras = sorted(mesher for mesher in meshers if mesher not in DEFAULT_MESHER_ORDER)
    return "-".join([*ordered, *extras])


def config_slug(max_mesh_size: float | None, min_mesh_angle: float) -> str:
    normalized = normalize_max_mesh_size(max_mesh_size)
    size_slug = "maxh-unrestricted" if normalized is None else f"maxh-{slug_token(normalized)}"
    return f"{size_slug}.mina-{slug_token(min_mesh_angle)}"


def run_slug(
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> str:
    return ".".join(
        [
            meshers_slug(meshers),
            config_slug(max_mesh_size, min_mesh_angle),
            f"pipeline-{pipeline_mode}",
            "audit" if stage_audit_enabled else "noaudit",
        ]
    )


def results_file_path(
    output_dir: Path,
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> Path:
    del meshers, max_mesh_size, min_mesh_angle, pipeline_mode, stage_audit_enabled
    return output_dir / "results.json"


def overview_plot_path(
    output_dir: Path,
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> Path:
    del meshers, max_mesh_size, min_mesh_angle, pipeline_mode, stage_audit_enabled
    return output_dir / "overview.png"


def summary_text_path(
    output_dir: Path,
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> Path:
    del meshers, max_mesh_size, min_mesh_angle, pipeline_mode, stage_audit_enabled
    return output_dir / "summary.txt"


def case_output_dir(
    output_dir: Path,
    number: int,
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> Path:
    del meshers, max_mesh_size, min_mesh_angle, pipeline_mode, stage_audit_enabled
    return output_dir / f"{number:03d}"


def case_plot_path(
    output_dir: Path,
    number: int,
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> Path:
    return (
        case_output_dir(
            output_dir,
            number,
            meshers,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        / "comparison.png"
    )


def case_mesh_path(
    output_dir: Path,
    number: int,
    mesher: str,
    *,
    meshers: list[str],
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> Path:
    return (
        case_output_dir(
            output_dir,
            number,
            meshers,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        / f"mesh_{mesher}.vtu"
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


def grading_metrics(mesh: Mesh) -> dict[str, float | int]:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if faces.size == 0:
        return {
            "unique_edge_count": 0,
            "edge_length_min": float("nan"),
            "edge_length_p01": float("nan"),
            "edge_length_p05": float("nan"),
            "edge_length_mean": float("nan"),
            "area_min": float("nan"),
            "area_p01": float("nan"),
            "area_p05": float("nan"),
            "area_mean": float("nan"),
            "short_edges_lt_0_5_count": 0,
            "short_edges_lt_1_0_count": 0,
        }

    vertices = np.asarray(mesh.vertices[:, :2], dtype=np.float64)
    triangles = vertices[faces]

    doubled_area = np.abs(
        (triangles[:, 1, 0] - triangles[:, 0, 0]) * (triangles[:, 2, 1] - triangles[:, 0, 1])
        - (triangles[:, 1, 1] - triangles[:, 0, 1]) * (triangles[:, 2, 0] - triangles[:, 0, 0])
    )
    areas = 0.5 * doubled_area

    edges = np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]],
        axis=0,
    )
    edges = np.sort(edges, axis=1)
    unique_edges = np.unique(edges, axis=0)
    edge_vectors = vertices[unique_edges[:, 1]] - vertices[unique_edges[:, 0]]
    edge_lengths = np.linalg.norm(edge_vectors, axis=1)

    return {
        "unique_edge_count": int(len(unique_edges)),
        "edge_length_min": float(edge_lengths.min()),
        "edge_length_p01": quantile(edge_lengths, 1.0),
        "edge_length_p05": quantile(edge_lengths, 5.0),
        "edge_length_mean": float(edge_lengths.mean()),
        "area_min": float(areas.min()),
        "area_p01": quantile(areas, 1.0),
        "area_p05": quantile(areas, 5.0),
        "area_mean": float(areas.mean()),
        "short_edges_lt_0_5_count": int(np.count_nonzero(edge_lengths < SHORT_EDGE_THRESHOLDS[0])),
        "short_edges_lt_1_0_count": int(np.count_nonzero(edge_lengths < SHORT_EDGE_THRESHOLDS[1])),
    }


def mesh_metrics(mesh: Mesh, quality: dict[str, Any]) -> dict[str, float | int]:
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


def format_case_row(mesher: str, result: dict[str, Any]) -> str:
    if result.get("status") != "success":
        error = result.get("error", {})
        return (
            f"{mesher:<12}  {'FAILED':<8}  "
            f"{error.get('type', 'Error')}: {error.get('message', 'unknown error')}"
        )

    metrics = result["metrics"]
    return (
        f"{mesher:<12}  {'ok':<8}  "
        f"{metrics['num_cells']:>8}  "
        f"{metrics['element_quality_min']:>9.4f}  "
        f"{metrics['element_quality_mean']:>10.4f}  "
        f"{metrics['aspect_ratio_max']:>8.4f}  "
        f"{metrics['aspect_ratio_mean']:>8.4f}  "
        f"{metrics['edge_ratio_max']:>8.4f}  "
        f"{metrics['edge_ratio_mean']:>8.4f}  "
        f"{metrics['skewness_max']:>8.4f}  "
        f"{result['time']:>7.2f}s"
    )


def format_case_scale_row(mesher: str, result: dict[str, Any]) -> str:
    if result.get("status") != "success":
        return f"{mesher:<12}  {'FAILED':<8}"

    metrics = result["metrics"]
    return (
        f"{mesher:<12}  {'ok':<8}  "
        f"{metrics['edge_length_min']:>9.4f}  "
        f"{metrics['edge_length_p01']:>9.4f}  "
        f"{metrics['edge_length_p05']:>9.4f}  "
        f"{metrics['area_min']:>10.4g}  "
        f"{metrics['area_p01']:>10.4f}  "
        f"{metrics['area_p05']:>10.4f}  "
        f"{metrics['short_edges_lt_0_5_count']:>8}  "
        f"{metrics['short_edges_lt_1_0_count']:>8}"
    )


def print_case_summary(case_record: dict[str, Any], meshers: list[str]) -> None:
    number = case_record["number"]
    bounds = case_record["bounds"]

    quality_rows: list[list[Any]] = []
    scale_rows: list[list[Any]] = []
    for mesher in meshers:
        result = case_record["meshers"].get(mesher, {})
        if result.get("status") != "success":
            quality_rows.append([mesher, "failed", "-", "-", "-", "-", "-", "-", "-", "-"])
            scale_rows.append([mesher, "failed", "-", "-", "-", "-", "-", "-", "-", "-"])
            continue

        metrics = result["metrics"]
        quality_rows.append(
            [
                mesher,
                "ok",
                metrics["num_cells"],
                f"{metrics['element_quality_min']:.4f}",
                f"{metrics['element_quality_mean']:.4f}",
                f"{metrics['aspect_ratio_max']:.2f}",
                f"{metrics['aspect_ratio_mean']:.2f}",
                f"{metrics['edge_ratio_max']:.2f}",
                f"{metrics['skewness_max']:.2f}",
                f"{result['time']:.2f}s",
            ]
        )
        scale_rows.append(
            [
                mesher,
                "ok",
                f"{metrics['edge_length_min']:.4f}",
                f"{metrics['edge_length_p01']:.4f}",
                f"{metrics['edge_length_p05']:.4f}",
                f"{metrics['area_min']:.4g}",
                f"{metrics['area_p01']:.4f}",
                f"{metrics['area_p05']:.4f}",
                metrics["short_edges_lt_0_5_count"],
                metrics["short_edges_lt_1_0_count"],
            ]
        )

    title = (
        f"Case {number}  ({bounds['xmin']:.0f}, {bounds['ymin']:.0f}) -> "
        f"({bounds['xmax']:.0f}, {bounds['ymax']:.0f})"
    )
    print()
    print(
        format_console_table(
            [
                "Mesher",
                "Status",
                "Cells",
                "EQ min",
                "EQ mean",
                "AR max",
                "AR mean",
                "ER max",
                "Skew max",
                "Time",
            ],
            quality_rows,
            title=title,
        )
    )
    print()
    print(
        format_console_table(
            [
                "Mesher",
                "Status",
                "Edge min",
                "Edge p01",
                "Edge p05",
                "Area min",
                "Area p01",
                "Area p05",
                "E<0.5m",
                "E<1.0m",
            ],
            scale_rows,
            title="Scale / grading",
        )
    )
    for mesher in meshers:
        result = case_record["meshers"].get(mesher, {})
        stage_contracts = format_stage_contracts(result.get("stage_audit"))
        if stage_contracts != "-":
            print(f"{mesher} contracts: {stage_contracts}")
    print(f"Preparation time: {case_record['prepare_time']:.2f}s")
    if case_record.get("plot_file"):
        print(f"Comparison plot: {case_record['plot_file']}")


def print_mesher_summary(results: dict[int, dict[str, Any]], meshers: list[str]) -> None:
    quality_rows: list[list[Any]] = []
    scale_rows: list[list[Any]] = []
    for mesher in meshers:
        successes = []
        for record in results.values():
            result = record.get("meshers", {}).get(mesher)
            if result and result.get("status") == "success":
                successes.append(result)

        if not successes:
            quality_rows.append([mesher, 0, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a", "n/a"])
            scale_rows.append([mesher, 0, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a"])
            continue

        eq_means = [item["metrics"]["element_quality_mean"] for item in successes]
        eq_mins = [item["metrics"]["element_quality_min"] for item in successes]
        ar_means = [item["metrics"]["aspect_ratio_mean"] for item in successes]
        ar_maxs = [item["metrics"]["aspect_ratio_max"] for item in successes]
        er_maxs = [item["metrics"]["edge_ratio_max"] for item in successes]
        sk_maxs = [item["metrics"]["skewness_max"] for item in successes]
        times = [item["time"] for item in successes]

        quality_rows.append(
            [
                mesher,
                len(successes),
                f"{np.mean(eq_means):.4f}",
                f"{np.min(eq_mins):.4f}",
                f"{np.mean(ar_means):.4f}",
                f"{np.max(ar_maxs):.4f}",
                f"{np.max(er_maxs):.4f}",
                f"{np.max(sk_maxs):.4f}",
                f"{np.mean(times):.2f}s",
            ]
        )

        edge_mins = [item["metrics"]["edge_length_min"] for item in successes]
        edge_p01s = [item["metrics"]["edge_length_p01"] for item in successes]
        area_mins = [item["metrics"]["area_min"] for item in successes]
        area_p01s = [item["metrics"]["area_p01"] for item in successes]
        short_half = [item["metrics"]["short_edges_lt_0_5_count"] for item in successes]
        short_one = [item["metrics"]["short_edges_lt_1_0_count"] for item in successes]

        scale_rows.append(
            [
                mesher,
                len(successes),
                f"{np.min(edge_mins):.4f}",
                f"{np.mean(edge_p01s):.4f}",
                f"{np.min(area_mins):.4g}",
                f"{np.mean(area_p01s):.4f}",
                np.max(short_half),
                np.max(short_one),
            ]
        )

    print()
    print(
        format_console_table(
            [
                "Mesher",
                "Success",
                "EQ mean",
                "EQ min",
                "AR mean",
                "AR max",
                "ER max",
                "Skew max",
                "Time mean",
            ],
            quality_rows,
            title="Mesher summary",
        )
    )
    print()
    print(
        format_console_table(
            [
                "Mesher",
                "Success",
                "Edge min",
                "Edge p01",
                "Area min",
                "Area p01",
                "E<0.5m max",
                "E<1.0m max",
            ],
            scale_rows,
            title="Scale / grading summary",
        )
    )


def case_complete(record: dict[str, Any], meshers: list[str]) -> bool:
    mesher_results = record.get("meshers", {})
    for mesher in meshers:
        result = mesher_results.get(mesher, {})
        if result.get("status") != "success":
            return False
        metrics = result.get("metrics", {})
        if not all(key in metrics for key in REQUIRED_METRIC_KEYS):
            return False
    return True


def resolve_requested_meshers(requested: list[str] | None) -> list[str]:
    available = set(dtcc_core.builder.available_2d_meshers())
    if requested is None:
        selected = [mesher for mesher in DEFAULT_MESHER_ORDER if mesher in available]
    else:
        selected = []
        for mesher in requested:
            normalized = mesher.strip().lower()
            if normalized not in selected:
                selected.append(normalized)

    if not selected:
        raise RuntimeError("No requested 2D meshers are available.")

    missing = [mesher for mesher in selected if mesher not in available]
    if missing:
        raise RuntimeError(
            f"Requested meshers are not available: {', '.join(missing)}. "
            f"Available meshers: {', '.join(sorted(available)) or 'none'}."
        )

    return selected


def run_mesher_for_case(
    number: int,
    city: City,
    bounds: Bounds,
    mesher: str,
    output_dir: Path,
    *,
    meshers: list[str],
    max_mesh_size: float | None,
    min_mesh_angle: float,
    conditioning_mode: str,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> tuple[dict[str, Any], Mesh | None]:
    start = time.perf_counter()

    try:
        with benchmark_conditioning_override(conditioning_mode):
            mesh = dtcc_core.datasets.city_flat_mesh.build_from_city(
                city,
                bounds=bounds.tuple,
                max_mesh_size=max_mesh_size,
                min_mesh_angle=min_mesh_angle,
                merge_buildings=MERGE_BUILDINGS,
                min_building_detail=MIN_BUILDING_DETAIL,
                min_building_area=MIN_BUILDING_AREA,
                report_mesh_quality=False,
                show_footprints=False,
                footprint_cleaning_plot_block=True,
                mesher=mesher,
                pipeline_mode=pipeline_mode,
                stage_audit_enabled=stage_audit_enabled,
            )
        quality = json_ready(mesh.quality())
        metrics = mesh_metrics(mesh, quality)
        mesh_path = case_mesh_path(
            output_dir,
            number,
            mesher,
            meshers=meshers,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        mesh.save(str(mesh_path))

        result: dict[str, Any] = {
            "status": "success",
            "time": round(time.perf_counter() - start, 2),
            "file": mesh_path.name,
            "conditioning_mode": conditioning_mode,
            "quality": quality,
            "metrics": metrics,
        }
        if stage_audit_enabled and getattr(mesh, "stage_audit", None) is not None:
            result["stage_audit"] = json_ready(mesh.stage_audit)
        return (
            result,
            mesh,
        )
    except Exception as exc:
        result = {
            "status": "failed",
            "time": round(time.perf_counter() - start, 2),
            "conditioning_mode": conditioning_mode,
            "error": error_details(exc),
        }
        return (
            result,
            None,
        )


def marker_categories(markers: np.ndarray | None, num_faces: int) -> np.ndarray:
    if markers is None or len(markers) != num_faces:
        return np.zeros(num_faces, dtype=np.int32)

    markers_array = np.asarray(markers, dtype=np.int32)
    categories = np.zeros(num_faces, dtype=np.int32)
    categories[markers_array == -1] = 1
    categories[markers_array >= 0] = 2
    return categories


def draw_mesh_panel(ax, mesh: Mesh, line_collection_cls, poly_collection_cls) -> None:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    vertices = np.asarray(mesh.vertices[:, :2], dtype=np.float64)
    polygons = vertices[faces]

    categories = marker_categories(mesh.markers, len(faces))
    facecolors = [MARKER_CATEGORY_COLORS[index] for index in categories]

    poly_collection = poly_collection_cls(
        polygons,
        facecolors=facecolors,
        edgecolors="none",
        alpha=0.92,
        rasterized=len(faces) > 5_000,
    )
    ax.add_collection(poly_collection)

    if len(faces) <= EDGE_FACE_LIMIT:
        edges = np.concatenate(
            [polygons[:, [0, 1]], polygons[:, [1, 2]], polygons[:, [2, 0]]],
            axis=0,
        )
        edge_collection = line_collection_cls(
            edges,
            colors="black",
            linewidths=0.14,
            alpha=0.45,
            rasterized=len(faces) > 5_000,
        )
        ax.add_collection(edge_collection)
    else:
        ax.text(
            0.98,
            0.02,
            "edges omitted",
            ha="right",
            va="bottom",
            fontsize=8,
            transform=ax.transAxes,
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "alpha": 0.85,
                "edgecolor": "none",
            },
        )

    ax.set_facecolor("#fcfcfc")


def mesher_panel_text(result: dict[str, Any]) -> str:
    if result.get("status") != "success":
        error = result.get("error", {})
        return "\n".join(
            [
                "FAILED",
                error.get("type", "Error"),
                error.get("message", "unknown error")[:180],
            ]
        )

    metrics = result["metrics"]
    return "\n".join(
        [
            f"cells {metrics['num_cells']}",
            f"EQ {metrics['element_quality_min']:.3f} / {metrics['element_quality_mean']:.3f}",
            f"AR {metrics['aspect_ratio_max']:.3f} / {metrics['aspect_ratio_mean']:.3f}",
            f"ER {metrics['edge_ratio_max']:.3f} / {metrics['edge_ratio_mean']:.3f}",
            f"Sk {metrics['skewness_max']:.3f} / {metrics['skewness_mean']:.3f}",
            f"edge {metrics['edge_length_min']:.3f} / {metrics['edge_length_p01']:.3f}",
            f"area {metrics['area_min']:.3g} / {metrics['area_p01']:.3f}",
            f"short <0.5 {metrics['short_edges_lt_0_5_count']}  <1.0 {metrics['short_edges_lt_1_0_count']}",
            f"time {result['time']:.2f}s",
        ]
    )


def set_panel_extent(axes, bounds: Bounds) -> None:
    padding = 0.05 * max(bounds.xmax - bounds.xmin, bounds.ymax - bounds.ymin)
    xmin = bounds.xmin - padding
    xmax = bounds.xmax + padding
    ymin = bounds.ymin - padding
    ymax = bounds.ymax + padding
    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])


def plot_case_comparison(
    output_path: Path,
    number: int,
    bounds: Bounds,
    meshers: list[str],
    mesh_objects: dict[str, Mesh],
    mesher_results: dict[str, dict[str, Any]],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    show_plot: bool = False,
) -> None:
    plt, line_collection_cls, poly_collection_cls, patch_cls = load_plot_modules()

    fig, axes = plt.subplots(
        1,
        len(meshers),
        figsize=(6.2 * len(meshers), 7.2),
        constrained_layout=True,
    )
    if len(meshers) == 1:
        axes = [axes]

    for ax, mesher in zip(axes, meshers):
        result = mesher_results.get(mesher, {})
        mesh = mesh_objects.get(mesher)

        if mesh is not None:
            draw_mesh_panel(ax, mesh, line_collection_cls, poly_collection_cls)
        else:
            ax.set_facecolor("#fcfcfc")

        ax.set_title(mesher)
        ax.text(
            0.02,
            0.02,
            mesher_panel_text(result),
            ha="left",
            va="bottom",
            fontsize=9,
            family="monospace",
            transform=ax.transAxes,
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "white",
                "alpha": 0.9,
                "edgecolor": "black",
                "linewidth": 0.4,
            },
        )

    set_panel_extent(axes, bounds)

    legend_handles = [
        patch_cls(facecolor=MARKER_CATEGORY_COLORS[0], edgecolor="black", label="ground"),
        patch_cls(facecolor=MARKER_CATEGORY_COLORS[1], edgecolor="black", label="halo"),
        patch_cls(facecolor=MARKER_CATEGORY_COLORS[2], edgecolor="black", label="building"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        f"Case {number} mesh comparison  ({bounds.xmin:.0f}, {bounds.ymin:.0f}) -> "
        f"({bounds.xmax:.0f}, {bounds.ymax:.0f})\n"
        f"{mesh_size_label(max_mesh_size)}, min angle={min_mesh_angle:g}°"
    )
    fig.savefig(output_path, dpi=200)
    if show_plot:
        plt.show(block=True)
    plt.close(fig)


def build_grid(
    results: dict[int, dict[str, Any]],
    mesher: str,
    metric_key: str,
) -> np.ndarray:
    nx, ny = grid_shape()
    grid = np.full((ny, nx), np.nan)

    for number, record in results.items():
        mesher_result = record.get("meshers", {}).get(mesher)
        if not mesher_result or mesher_result.get("status") != "success":
            continue

        value = float(mesher_result["metrics"][metric_key])
        ix, iy = case_to_grid(number)
        grid[iy, ix] = value

    return grid


def plot_overview(
    results: dict[int, dict[str, Any]],
    meshers: list[str],
    output_path: Path,
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    pipeline_mode: str,
) -> None:
    plt, _, _, _ = load_plot_modules()

    fig, axes = plt.subplots(
        len(meshers),
        len(OVERVIEW_METRICS),
        figsize=(5.4 * len(OVERVIEW_METRICS), 3.9 * len(meshers)),
        constrained_layout=True,
    )
    if len(meshers) == 1:
        axes = np.asarray([axes])
    nx, ny = grid_shape()

    for row, mesher in enumerate(meshers):
        for col, (metric_key, title, cmap, fmt) in enumerate(OVERVIEW_METRICS):
            ax = axes[row, col]
            grid = build_grid(results, mesher, metric_key)
            image = ax.imshow(grid, origin="lower", cmap=cmap, aspect="equal")
            annotate_heatmap(ax, grid, fmt)
            ax.set_title(f"{mesher} - {title}")
            ax.set_xticks(range(nx))
            ax.set_yticks(range(ny))
            ax.set_xlabel("Grid X")
            ax.set_ylabel("Grid Y")
            fig.colorbar(image, ax=ax, shrink=0.82, pad=0.03)

    fig.suptitle(
        f"{active_area_label()} mesh quality overview by mesher\n"
        f"{mesh_size_label(max_mesh_size)}, min angle={min_mesh_angle:g}°, "
        f"pipeline={pipeline_mode}"
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def compact_case_status(record: dict[str, Any], meshers: list[str]) -> str:
    parts = []
    for mesher in meshers:
        result = record["meshers"].get(mesher, {})
        if result.get("status") != "success":
            error = result.get("error", {})
            failure = f"{mesher}: FAIL"
            stage_contracts = format_stage_contracts(result.get("stage_audit"))
            if stage_contracts != "-":
                failure += f" ({stage_contracts})"
            elif error.get("type"):
                failure += f" ({error.get('type')})"
            parts.append(failure)
            continue
        metrics = result["metrics"]
        summary = (
            f"{mesher}: ARmax {metrics['aspect_ratio_max']:.2f}, "
            f"edge p01 {metrics['edge_length_p01']:.2f}m"
        )
        stage_contracts = format_stage_contracts(result.get("stage_audit"))
        if stage_contracts != "-":
            summary += f", contracts [{stage_contracts}]"
        parts.append(summary)
    return "  ".join(parts)


def status_emoji(status: str | None) -> str:
    return "✅" if status == "success" else "❌"


def aggregate_status_emoji(successes: int, failures: int) -> str:
    if failures == 0:
        return "✅"
    if successes == 0:
        return "❌"
    return "⚠️"


def build_summary_report(
    results: dict[int, dict[str, Any]],
    meshers: list[str],
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    conditioning_mode: str,
    pipeline_mode: str,
    stage_audit_enabled: bool,
    elapsed_seconds: float,
) -> str:
    case_rows: list[list[Any]] = []
    for number in sorted(results):
        record = results[number]
        for mesher in meshers:
            result = record.get("meshers", {}).get(mesher, {})
            if result.get("status") == "success":
                metrics = result["metrics"]
                case_rows.append(
                    [
                        f"{number:03d}",
                        mesher,
                        f"{status_emoji('success')} ok",
                        metrics["num_cells"],
                        f"{metrics['element_quality_min']:.4f}",
                        f"{metrics['aspect_ratio_max']:.2f}",
                        f"{metrics['edge_length_p01']:.2f}",
                        metrics["short_edges_lt_0_5_count"],
                        f"{result['time']:.2f}s",
                    ]
                )
            else:
                error = result.get("error", {})
                case_rows.append(
                    [
                        f"{number:03d}",
                        mesher,
                        f"{status_emoji('failed')} fail",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        f"{error.get('type', '-')}: {error.get('message', '-')}",
                    ]
                )

    summary_rows: list[list[Any]] = []
    for mesher in meshers:
        mesher_results = [
            record.get("meshers", {}).get(mesher, {})
            for record in results.values()
        ]
        successes = [result for result in mesher_results if result.get("status") == "success"]
        failures = len(mesher_results) - len(successes)
        if successes:
            eq_min = min(result["metrics"]["element_quality_min"] for result in successes)
            ar_max = max(result["metrics"]["aspect_ratio_max"] for result in successes)
            edge_p01 = float(np.mean([result["metrics"]["edge_length_p01"] for result in successes]))
            mean_time = float(np.mean([result["time"] for result in successes]))
            summary_rows.append(
                [
                    aggregate_status_emoji(len(successes), failures),
                    mesher,
                    f"{len(successes)}/{len(mesher_results)}",
                    failures,
                    f"{eq_min:.4f}",
                    f"{ar_max:.2f}",
                    f"{edge_p01:.2f}",
                    f"{mean_time:.2f}s",
                ]
            )
        else:
            summary_rows.append(
                [
                    aggregate_status_emoji(0, failures),
                    mesher,
                    f"0/{len(mesher_results)}",
                    failures,
                    "-",
                    "-",
                    "-",
                    "-",
                ]
            )

    lines = [
        "bench_mesh_2d",
        (
            f"Config: {mesh_size_label(max_mesh_size)}, min angle={min_mesh_angle:g}°, "
            f"conditioning={conditioning_mode}, "
            f"pipeline={pipeline_mode}, stage-audit={'on' if stage_audit_enabled else 'off'}"
        ),
        f"Elapsed time: {elapsed_seconds:.1f}s",
        "",
        format_console_table(
            ["Status", "Mesher", "Success", "Fail", "Worst EQ", "Worst AR", "Mean edge p01", "Mean time"],
            summary_rows,
            title="Summary",
        ),
        "",
        format_console_table(
            ["Case", "Mesher", "Status", "Cells", "EQ min", "AR max", "Edge p01", "E<0.5m", "Time/Error"],
            case_rows,
            title="Detailed results",
        ),
    ]
    return "\n".join(lines)


def run_case(
    number: int,
    meshers: list[str],
    output_dir: Path,
    create_case_plot: bool,
    show_plot: bool,
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    conditioning_mode: str,
    pipeline_mode: str,
    stage_audit_enabled: bool,
) -> dict[str, Any]:
    ix, iy = case_to_grid(number)
    bounds = make_bounds(ix, iy)
    case_dir = case_output_dir(
        output_dir,
        number,
        meshers,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        pipeline_mode=pipeline_mode,
        stage_audit_enabled=stage_audit_enabled,
    )
    case_dir.mkdir(parents=True, exist_ok=True)

    city, prepare_time = prepare_city(bounds)

    case_record: dict[str, Any] = {
        "number": number,
        "bounds": bounds_to_dict(bounds),
        "prepare_time": round(prepare_time, 2),
        "config": {
            "max_mesh_size": max_mesh_size,
            "min_mesh_angle": min_mesh_angle,
            "conditioning_mode": conditioning_mode,
            "pipeline_mode": pipeline_mode,
            "stage_audit_enabled": stage_audit_enabled,
        },
        "meshers": {},
    }
    mesh_objects: dict[str, Mesh] = {}

    for mesher in meshers:
        result, mesh = run_mesher_for_case(
            number,
            city,
            bounds,
            mesher,
            output_dir,
            meshers=meshers,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            conditioning_mode=conditioning_mode,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        case_record["meshers"][mesher] = json_ready(result)
        if mesh is not None:
            mesh_objects[mesher] = mesh

    if create_case_plot:
        plot_path = case_plot_path(
            output_dir,
            number,
            meshers,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            pipeline_mode=pipeline_mode,
            stage_audit_enabled=stage_audit_enabled,
        )
        plot_case_comparison(
            plot_path,
            number=number,
            bounds=bounds,
            meshers=meshers,
            mesh_objects=mesh_objects,
            mesher_results=case_record["meshers"],
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            show_plot=show_plot,
        )
        case_record["plot_file"] = plot_path.name

    summary_path = case_dir / "summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(json_ready(case_record), handle, indent=2)

    return case_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare flat-mesh quality across a benchmark area for multiple meshers."
    )
    add_area_argument(parser)
    add_cases_argument(parser)
    parser.add_argument(
        "--meshers",
        nargs="+",
        default=None,
        help="Meshers to compare. Defaults to all available from: dtcc_mesher triangle spade.",
    )
    parser.add_argument(
        "--max-mesh-size",
        type=parse_max_mesh_size_argument,
        default=DEFAULT_MAX_MESH_SIZE,
        help=(
            "Maximum target triangle edge length in meters. "
            "Use 'none' for unrestricted size."
        ),
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
        default=None,
        help="Directory for JSON, VTU, and PNG outputs.",
    )
    parser.add_argument(
        "--conditioning-mode",
        choices=CONDITIONING_MODES,
        default="new",
        help="Footprint conditioning benchmark mode used before meshing.",
    )
    parser.add_argument(
        "--per-case-plots",
        action="store_true",
        help="When running many cases, also save a side-by-side PNG for every case.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip all PNG plot generation.",
    )
    parser.add_argument(
        "--show-plot",
        action="store_true",
        help="Display the side-by-side comparison plot for a single case.",
    )
    parser.add_argument(
        "--stage-audit",
        action="store_true",
        help="Record per-stage contract/audit data for conditioned footprints and ground mesh.",
    )
    args = parser.parse_args()
    args.pipeline_mode = "strict"
    args.stage_audit = True

    return args


def main() -> None:
    args = parse_args()
    set_active_area(args.area)
    benchmark_start = time.perf_counter()
    output_dir = args.output_dir or benchmark_output_dir("output_mesh_2d")
    output_dir.mkdir(parents=True, exist_ok=True)
    cases_explicit = args.cases is not None
    selected_cases = resolve_case_numbers(args.cases)

    meshers = resolve_requested_meshers(args.meshers)
    metadata = {
        "benchmark": "bench_mesh_2d",
        "config": {
            "meshers": meshers,
            "area": args.area,
            "max_mesh_size": args.max_mesh_size,
            "min_mesh_angle": MIN_MESH_ANGLE,
            "conditioning_mode": args.conditioning_mode,
            "pipeline_mode": args.pipeline_mode,
            "stage_audit_enabled": args.stage_audit,
        },
    }
    results_path = results_file_path(
        output_dir,
        meshers,
        max_mesh_size=args.max_mesh_size,
        min_mesh_angle=MIN_MESH_ANGLE,
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
    nx, ny = grid_shape()
    total = nx * ny

    loaded_count = len(results)
    if loaded_count > 0:
        print(f"Loaded {loaded_count} previous case record(s) from {results_path}")
        print("  Delete that file to recompute everything from scratch.")
        print()

    print(
        f"Area: {active_area_label()} | "
        f"Configuration: {mesh_size_label(args.max_mesh_size)}, "
        f"min angle={MIN_MESH_ANGLE:g}°, "
        f"conditioning={args.conditioning_mode}, "
        f"pipeline={args.pipeline_mode}, "
        f"stage-audit={'on' if args.stage_audit else 'off'}"
    )
    print()

    if cases_explicit:
        cases_to_run = selected_cases
        if len(cases_to_run) == 1 and cases_to_run[0] in results:
            print(f"Case {cases_to_run[0]} exists in cache and will be recomputed.")
            print()
    else:
        cases_to_run = [
            number
            for number in selected_cases
            if not case_complete(results.get(number, {}), meshers)
        ]
        if len(selected_cases) == total and not cases_to_run:
            print(f"All {total} cases are already complete for {', '.join(meshers)}.")
            print()
        else:
            print(
                f"{len(cases_to_run)} case(s) to compute for {', '.join(meshers)}: "
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

        create_case_plot = (
            not args.no_plots
            and (cases_explicit or args.per_case_plots)
        )
        case_record = run_case(
            number,
            meshers=meshers,
            output_dir=output_dir,
            create_case_plot=create_case_plot,
            show_plot=bool(args.show_plot and len(selected_cases) == 1),
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=MIN_MESH_ANGLE,
            conditioning_mode=args.conditioning_mode,
            pipeline_mode=args.pipeline_mode,
            stage_audit_enabled=args.stage_audit,
        )
        results[number] = case_record
        save_results(results_path, results, metadata=metadata)

        mesher_successes = sum(
            1
            for mesher in meshers
            if case_record["meshers"].get(mesher, {}).get("status") == "success"
        )
        if mesher_successes == len(meshers):
            success_count += 1
        else:
            failure_count += 1

        print(f"  {compact_case_status(case_record, meshers)}")

    elapsed = time.perf_counter() - start_all

    if cases_to_run:
        print()
        print(
            f"Finished in {elapsed:.0f}s  "
            f"({success_count} fully successful case(s), {failure_count} incomplete case(s), "
            f"{loaded_count} previously cached)"
        )

    summary_results = {number: results[number] for number in selected_cases if number in results}

    if not args.no_plots:
        overview_path = overview_plot_path(
            output_dir,
            meshers,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=MIN_MESH_ANGLE,
            pipeline_mode=args.pipeline_mode,
            stage_audit_enabled=args.stage_audit,
        )
        plot_overview(
            summary_results,
            meshers,
            overview_path,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=MIN_MESH_ANGLE,
            pipeline_mode=args.pipeline_mode,
        )
    total_elapsed = time.perf_counter() - benchmark_start
    save_results(results_path, results, metadata=metadata)
    report = build_summary_report(
        summary_results,
        meshers,
        max_mesh_size=args.max_mesh_size,
        min_mesh_angle=MIN_MESH_ANGLE,
        conditioning_mode=args.conditioning_mode,
        pipeline_mode=args.pipeline_mode,
        stage_audit_enabled=args.stage_audit,
        elapsed_seconds=total_elapsed,
    )
    summary_path = summary_text_path(
        output_dir,
        meshers,
        max_mesh_size=args.max_mesh_size,
        min_mesh_angle=MIN_MESH_ANGLE,
        pipeline_mode=args.pipeline_mode,
        stage_audit_enabled=args.stage_audit,
    )
    summary_path.write_text(report + "\n")
    print()
    print(report)
    print()
    print(f"Results file: {results_path}")
    if not args.no_plots:
        print(f"Overview plot: {overview_path}")
    print(f"Summary text: {summary_path}")


if __name__ == "__main__":
    main()
