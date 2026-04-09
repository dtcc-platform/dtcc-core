"""
Survey TetGen-backed 3D volume-mesh quality across central Stockholm tiles.

This script reuses the same city-preparation pipeline as the 2D survey, then
builds a tetrahedral volume mesh with dtcc-core's TetGen path for each tile.
It caches per-case results, saves per-case XDMF/HDF5 meshes, and writes
overview heatmaps for the key 3D metrics.

Typical usage:
    python mesh_quality_survey_3d.py 55
    python mesh_quality_survey_3d.py
    python mesh_quality_survey_3d.py 55 --max-mesh-size 8
    python mesh_quality_survey_3d.py 55 --quality-ratio 1.4
    python mesh_quality_survey_3d.py 55 --size-only --preserve-surface
    python mesh_quality_survey_3d.py 55 --preserve-surface --max-added-points 10000
    python mesh_quality_survey_3d.py --no-merge-buildings

Coordinate system: SWEREF99 TM (EPSG:3006)
"""

from __future__ import annotations

import argparse
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

import dtcc_core
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.model import GeometryType, VolumeMesh
from dtcc_core.model.mixins.mesh.quality import tet_element_quality

from mesh_quality_survey_2d import (
    DEFAULT_DELAY_BETWEEN_CASES,
    DEFAULT_MAX_MESH_SIZE,
    MIN_BUILDING_AREA,
    MIN_BUILDING_DETAIL,
    MERGE_BUILDINGS,
    NX,
    NY,
    _load_plot_modules,
    _quantile,
    _slug_token,
    annotate_heatmap,
    bounds_to_dict,
    case_to_grid,
    json_ready,
    load_results,
    make_bounds,
    mesh_size_label,
    prepare_city,
    save_results,
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


def positive_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be a positive number") from exc
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be a positive integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


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
) -> str:
    return ".".join(
        [
            f"maxh-{_slug_token(max_mesh_size)}",
            f"mina-{_slug_token(min_mesh_angle)}",
            f"height-{_slug_token(domain_height)}",
            (
                f"q-{_slug_token(quality_ratio)}"
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
) -> Path:
    return output_dir / (
        f"mesh_quality_survey_3d."
        f"{config_slug(max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, domain_height=domain_height, quality_ratio=quality_ratio, quality_enabled=quality_enabled, preserve_surface=preserve_surface, max_added_points=max_added_points, lod=lod, merge_buildings=merge_buildings, mesher=mesher)}."
        f"json"
    )


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
) -> Path:
    return output_dir / (
        f"mesh_quality_survey_3d."
        f"{config_slug(max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, domain_height=domain_height, quality_ratio=quality_ratio, quality_enabled=quality_enabled, preserve_surface=preserve_surface, max_added_points=max_added_points, lod=lod, merge_buildings=merge_buildings, mesher=mesher)}."
        f"png"
    )


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
) -> Path:
    return output_dir / (
        f"{number:03d}.tetgen."
        f"{config_slug(max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, domain_height=domain_height, quality_ratio=quality_ratio, quality_enabled=quality_enabled, preserve_surface=preserve_surface, max_added_points=max_added_points, lod=lod, merge_buildings=merge_buildings, mesher=mesher)}."
        f"xdmf"
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
) -> str:
    return (
        f"{number:03d}.tetgen-input."
        f"{config_slug(max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, domain_height=domain_height, quality_ratio=quality_ratio, quality_enabled=quality_enabled, preserve_surface=preserve_surface, max_added_points=max_added_points, lod=lod, merge_buildings=merge_buildings, mesher=mesher)}"
    )


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
) -> dict[str, Path]:
    stem = case_tetgen_input_stem(
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
    )
    return {
        "ground": output_dir / f"{stem}.tetgen-input-ground.xdmf",
        "shell": output_dir / f"{stem}.tetgen-input-shell.xdmf",
        "plc": output_dir / f"{stem}.tetgen-input-plc.xdmf",
    }


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
        "edge_length_p01": _quantile(edge_lengths, 1.0),
        "edge_length_p05": _quantile(edge_lengths, 5.0),
        "edge_length_mean": float(edge_lengths.mean()),
        "volume_min": float(volumes.min()),
        "volume_p01": _quantile(volumes, 1.0),
        "volume_p05": _quantile(volumes, 5.0),
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

    print()
    print(
        f"Case {number}  ({bounds['xmin']:.0f}, {bounds['ymin']:.0f}) -> "
        f"({bounds['xmax']:.0f}, {bounds['ymax']:.0f})"
    )
    print(
        f"{'Status':<8}  {'Cells':>8}  {'ElemQ min':>9}  {'ElemQ mean':>10}  "
        f"{'AR max':>8}  {'ER max':>8}  {'Skew max':>8}  {'Time':>8}"
    )
    print("-" * 82)
    print(format_case_row(result))
    print("-" * 82)
    print("Scale / grading metrics")
    print(
        f"{'Status':<8}  {'Edge min':>9}  {'Edge p01':>9}  {'Edge p05':>9}  "
        f"{'Vol min':>10}  {'Vol p01':>10}  {'Vol p05':>10}  {'Q<0.05':>11}  {'Q<0.10':>11}"
    )
    print("-" * 110)
    print(format_case_scale_row(result))
    print("-" * 110)
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


def print_survey_summary(results: dict[int, dict[str, Any]]) -> None:
    print()
    print("TetGen summary across cached cases")
    print(
        f"{'Success':>7}  {'EQ mean':>9}  {'EQ min':>9}  {'AR max':>9}  "
        f"{'ER max':>9}  {'Edge p01':>10}  {'Vol p01':>10}  {'Q<0.05':>11}  {'Time mean':>10}"
    )
    print("-" * 100)

    successes = [
        record["result"]
        for record in results.values()
        if record.get("result", {}).get("status") == "success"
    ]
    if not successes:
        print(
            f"{0:>7}  {'n/a':>9}  {'n/a':>9}  {'n/a':>9}  "
            f"{'n/a':>9}  {'n/a':>10}  {'n/a':>10}  {'n/a':>11}  {'n/a':>10}"
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

    print(
        f"{len(successes):>7}  "
        f"{np.mean(eq_means):>9.4f}  "
        f"{np.min(eq_mins):>9.4f}  "
        f"{np.max(ar_maxs):>9.4f}  "
        f"{np.max(er_maxs):>9.4f}  "
        f"{np.mean(edge_p01s):>10.4f}  "
        f"{np.mean(volume_p01s):>10.4f}  "
        f"{np.max(low_quality):>11}  "
        f"{np.mean(times):>10.2f}s"
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
) -> None:
    plt, _, _, _ = _load_plot_modules()

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
        f"merge={'on' if merge_buildings else 'off'}"
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def compact_case_status(record: dict[str, Any]) -> str:
    result = record["result"]
    if result.get("status") != "success":
        error = result.get("error", {})
        return f"FAIL {error.get('type', 'Error')}: {error.get('message', 'unknown error')}"

    metrics = result["metrics"]
    return (
        f"ARmax {metrics['aspect_ratio_max']:.2f}, "
        f"edge p01 {metrics['edge_length_p01']:.2f}m, "
        f"vol p01 {metrics['volume_p01']:.2f}, "
        f"Q<0.05 {metrics['low_quality_lt_0_05_count']}"
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
    save_tetgen_input: bool,
) -> dict[str, Any]:
    ix, iy = case_to_grid(number)
    bounds = make_bounds(ix, iy)
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
            "save_tetgen_input": save_tetgen_input,
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
        return case_record

    try:
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
            tetgen_debug_output_dir=output_dir if save_tetgen_input else None,
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
                )
                if save_tetgen_input
                else None
            ),
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
        )
        volume_mesh.save(str(mesh_path))

        case_record["result"] = {
            "status": "success",
            "time": round(time.perf_counter() - start, 2),
            "file": mesh_path.name,
            "tetgen_switches": json_ready(switches_params),
            "quality": quality,
            "metrics": metrics,
            "boundary_markers": boundary_marker_histogram(volume_mesh),
        }
        if save_tetgen_input:
            case_record["result"]["tetgen_input"] = {
                key: path.name for key, path in tetgen_input_paths.items() if path.exists()
            }
    except Exception as exc:
        case_record["result"] = {
            "status": "failed",
            "time": round(time.perf_counter() - start, 2),
            "tetgen_switches": json_ready(switches_params),
            "error": error_details(exc),
        }
        if save_tetgen_input:
            case_record["result"]["tetgen_input"] = {
                key: path.name for key, path in tetgen_input_paths.items() if path.exists()
            }

    return json_ready(case_record)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Survey TetGen 3D volume-mesh quality across Stockholm tiles."
    )
    parser.add_argument(
        "case_number",
        nargs="?",
        type=int,
        help=f"Single case to recompute (1-{NX * NY}). Omit to process all missing cases.",
    )
    parser.add_argument(
        "--max-mesh-size",
        type=positive_float,
        default=DEFAULT_MAX_MESH_SIZE,
        help="Maximum target edge length in meters for the surface mesh and TetGen volume cap.",
    )
    parser.add_argument(
        "--min-mesh-angle",
        type=positive_float,
        default=25.0,
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
        default=Path(__file__).resolve().parent / "output_3d",
        help="Directory for JSON, XDMF/HDF5, and PNG outputs.",
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

    total = NX * NY
    if args.case_number is not None and not (1 <= args.case_number <= total):
        parser.error(f"case_number must be between 1 and {total}")
    if args.save_tetgen_input is None:
        args.save_tetgen_input = args.case_number is not None

    return args


def main() -> None:
    if not is_tetgen_available():
        raise RuntimeError(
            "TetGen is not available in the active environment. "
            "Install dtcc-tetgen-wrapper in the same venv first."
        )

    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

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
    )
    results = load_results(results_path)
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
        f"tetgen-input={'on' if args.save_tetgen_input else 'off'}"
    )
    print()

    if args.case_number is not None:
        cases_to_run = [args.case_number]
        if args.case_number in results:
            print(f"Case {args.case_number} exists in cache and will be recomputed.")
            print()
    else:
        cases_to_run = [
            number for number in range(1, total + 1) if not case_complete(results.get(number, {}))
        ]
        if not cases_to_run:
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
            save_tetgen_input=args.save_tetgen_input,
        )
        results[number] = case_record
        save_results(results_path, results)

        if case_record["result"].get("status") == "success":
            success_count += 1
        else:
            failure_count += 1

        print(f"  {compact_case_status(case_record)}")

    elapsed = time.perf_counter() - start_all

    if cases_to_run:
        print()
        print(
            f"Finished in {elapsed:.0f}s  "
            f"({success_count} successful case(s), {failure_count} failed case(s), "
            f"{loaded_count} previously cached)"
        )

    print()
    print(f"Results file: {results_path}")

    if args.case_number is not None:
        print_case_summary(results[args.case_number])
    else:
        incomplete_cases = [
            number for number, record in results.items() if not case_complete(record)
        ]
        print_survey_summary(results)
        if incomplete_cases:
            print()
            print(
                "Incomplete cases: "
                + ", ".join(str(number) for number in sorted(incomplete_cases))
            )

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
        )
        plot_overview(
            results,
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
        )
        print(f"Overview plot: {overview_path}")


if __name__ == "__main__":
    main()
