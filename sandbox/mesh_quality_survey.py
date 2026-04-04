"""
Compare flat-mesh quality across central Stockholm for one or more 2D meshers.

This script builds the same prepared City input once per tile, then runs one or
more mesh generators against that identical input. It caches per-case results,
saves per-mesher VTU files, and writes visual comparison plots.

Typical usage:
    python mesh_quality_survey.py 55
    python mesh_quality_survey.py
    python mesh_quality_survey.py 55 --meshers triangle dtcc_mesher spade
    python mesh_quality_survey.py --per-case-plots

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
from dtcc_core.model import Bounds, City, GeometryType, Mesh

# Configuration ---------------------------------------------------------------

X_MIN = 673_000
Y_MIN = 6_578_500
NX, NY = 10, 10
BOX_SIZE = 500

MAX_MESH_SIZE = 10.0
MIN_MESH_ANGLE = 25.0
MIN_BUILDING_DETAIL = 0.5
MIN_BUILDING_AREA = 15.0
MERGE_BUILDINGS = True

RASTER_CELL_SIZE = 2.0
RASTER_RADIUS = 3.0
OUTLIER_THRESHOLD = 3.0

DEFAULT_DELAY_BETWEEN_CASES = 8.0
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


def case_to_grid(number: int) -> tuple[int, int]:
    iy, ix = divmod(number - 1, NX)
    return ix, iy


def make_bounds(ix: int, iy: int) -> Bounds:
    x0 = X_MIN + ix * BOX_SIZE
    y0 = Y_MIN + iy * BOX_SIZE
    return Bounds(x0, y0, x0 + BOX_SIZE, y0 + BOX_SIZE)


def bounds_to_dict(bounds: Bounds) -> dict[str, float]:
    return {
        "xmin": float(bounds.xmin),
        "ymin": float(bounds.ymin),
        "xmax": float(bounds.xmax),
        "ymax": float(bounds.ymax),
    }


def meshers_slug(meshers: list[str]) -> str:
    ordered = [mesher for mesher in DEFAULT_MESHER_ORDER if mesher in meshers]
    extras = sorted(mesher for mesher in meshers if mesher not in DEFAULT_MESHER_ORDER)
    return "-".join([*ordered, *extras])


def results_file_path(output_dir: Path, meshers: list[str]) -> Path:
    return output_dir / f"mesh_quality_survey.{meshers_slug(meshers)}.json"


def overview_plot_path(output_dir: Path, meshers: list[str]) -> Path:
    return output_dir / f"mesh_quality_survey.{meshers_slug(meshers)}.png"


def case_plot_path(output_dir: Path, number: int, meshers: list[str]) -> Path:
    return output_dir / f"{number:03d}.compare.{meshers_slug(meshers)}.png"


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_results(path: Path) -> dict[int, dict[str, Any]]:
    if not path.exists():
        return {}
    with path.open() as handle:
        data = json.load(handle)
    return {int(key): value for key, value in data.items()}


def save_results(path: Path, results: dict[int, dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = {str(key): results[key] for key in sorted(results)}
    with path.open("w") as handle:
        json.dump(ordered, handle, indent=2)


def _quantile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values, q))


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
        "edge_length_p01": _quantile(edge_lengths, 1.0),
        "edge_length_p05": _quantile(edge_lengths, 5.0),
        "edge_length_mean": float(edge_lengths.mean()),
        "area_min": float(areas.min()),
        "area_p01": _quantile(areas, 1.0),
        "area_p05": _quantile(areas, 5.0),
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

    print()
    print(
        f"Case {number}  ({bounds['xmin']:.0f}, {bounds['ymin']:.0f}) -> "
        f"({bounds['xmax']:.0f}, {bounds['ymax']:.0f})"
    )
    print(
        f"{'Mesher':<12}  {'Status':<8}  {'Cells':>8}  "
        f"{'ElemQ min':>9}  {'ElemQ mean':>10}  "
        f"{'AR max':>8}  {'AR mean':>8}  "
        f"{'ER max':>8}  {'ER mean':>8}  "
        f"{'Skew max':>8}  {'Time':>8}"
    )
    print("-" * 110)
    for mesher in meshers:
        print(format_case_row(mesher, case_record["meshers"].get(mesher, {})))
    print("-" * 110)
    print("Scale / grading metrics")
    print(
        f"{'Mesher':<12}  {'Status':<8}  {'Edge min':>9}  {'Edge p01':>9}  {'Edge p05':>9}  "
        f"{'Area min':>10}  {'Area p01':>10}  {'Area p05':>10}  {'E<0.5m':>8}  {'E<1.0m':>8}"
    )
    print("-" * 110)
    for mesher in meshers:
        print(format_case_scale_row(mesher, case_record["meshers"].get(mesher, {})))
    print("-" * 110)
    print(f"Preparation time: {case_record['prepare_time']:.2f}s")
    if case_record.get("plot_file"):
        print(f"Comparison plot: {case_record['plot_file']}")


def print_mesher_summary(results: dict[int, dict[str, Any]], meshers: list[str]) -> None:
    print()
    print("Mesher summary across cached cases")
    print(
        f"{'Mesher':<12}  {'Success':>7}  {'EQ mean':>9}  {'EQ min':>9}  "
        f"{'AR mean':>9}  {'AR max':>9}  {'ER max':>9}  {'Skew max':>10}  {'Time mean':>10}"
    )
    print("-" * 100)

    for mesher in meshers:
        successes = []
        for record in results.values():
            result = record.get("meshers", {}).get(mesher)
            if result and result.get("status") == "success":
                successes.append(result)

        if not successes:
            print(f"{mesher:<12}  {0:>7}  {'n/a':>9}  {'n/a':>9}  {'n/a':>9}  {'n/a':>9}  {'n/a':>9}  {'n/a':>10}  {'n/a':>10}")
            continue

        eq_means = [item["metrics"]["element_quality_mean"] for item in successes]
        eq_mins = [item["metrics"]["element_quality_min"] for item in successes]
        ar_means = [item["metrics"]["aspect_ratio_mean"] for item in successes]
        ar_maxs = [item["metrics"]["aspect_ratio_max"] for item in successes]
        er_maxs = [item["metrics"]["edge_ratio_max"] for item in successes]
        sk_maxs = [item["metrics"]["skewness_max"] for item in successes]
        times = [item["time"] for item in successes]

        print(
            f"{mesher:<12}  {len(successes):>7}  "
            f"{np.mean(eq_means):>9.4f}  {np.min(eq_mins):>9.4f}  "
            f"{np.mean(ar_means):>9.4f}  {np.max(ar_maxs):>9.4f}  "
            f"{np.max(er_maxs):>9.4f}  {np.max(sk_maxs):>10.4f}  "
            f"{np.mean(times):>10.2f}s"
        )

    print()
    print("Scale / grading summary across cached cases")
    print(
        f"{'Mesher':<12}  {'Success':>7}  {'Edge min':>9}  {'Edge p01':>9}  {'Area min':>10}  "
        f"{'Area p01':>10}  {'E<0.5m max':>11}  {'E<1.0m max':>11}"
    )
    print("-" * 90)

    for mesher in meshers:
        successes = []
        for record in results.values():
            result = record.get("meshers", {}).get(mesher)
            if result and result.get("status") == "success":
                successes.append(result)

        if not successes:
            print(f"{mesher:<12}  {0:>7}  {'n/a':>9}  {'n/a':>9}  {'n/a':>10}  {'n/a':>10}  {'n/a':>11}  {'n/a':>11}")
            continue

        edge_mins = [item["metrics"]["edge_length_min"] for item in successes]
        edge_p01s = [item["metrics"]["edge_length_p01"] for item in successes]
        area_mins = [item["metrics"]["area_min"] for item in successes]
        area_p01s = [item["metrics"]["area_p01"] for item in successes]
        short_half = [item["metrics"]["short_edges_lt_0_5_count"] for item in successes]
        short_one = [item["metrics"]["short_edges_lt_1_0_count"] for item in successes]

        print(
            f"{mesher:<12}  {len(successes):>7}  "
            f"{np.min(edge_mins):>9.4f}  {np.mean(edge_p01s):>9.4f}  "
            f"{np.min(area_mins):>10.4g}  {np.mean(area_p01s):>10.4f}  "
            f"{np.max(short_half):>11}  {np.max(short_one):>11}"
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


def prepare_city(bounds: Bounds) -> tuple[City, float]:
    start = time.perf_counter()

    pointcloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
    buildings = dtcc_core.io.data.download_footprints(bounds=bounds)
    pointcloud = pointcloud.remove_global_outliers(OUTLIER_THRESHOLD)

    raster = dtcc_core.builder.build_terrain_raster(
        pointcloud,
        cell_size=RASTER_CELL_SIZE,
        radius=RASTER_RADIUS,
        ground_only=True,
    )
    buildings = dtcc_core.builder.extract_roof_points(buildings, pointcloud)
    buildings = dtcc_core.builder.compute_building_heights(buildings, raster, overwrite=True)

    city = City()
    city.add_terrain(raster)
    city.add_buildings(buildings, remove_outside_terrain=True)

    return city, time.perf_counter() - start


def run_mesher_for_case(
    number: int,
    city: City,
    mesher: str,
    output_dir: Path,
) -> tuple[dict[str, Any], Mesh | None]:
    start = time.perf_counter()

    try:
        mesh = dtcc_core.builder.build_city_flat_mesh(
            city,
            lod=GeometryType.LOD0,
            max_mesh_size=MAX_MESH_SIZE,
            min_mesh_angle=MIN_MESH_ANGLE,
            merge_buildings=MERGE_BUILDINGS,
            min_building_detail=MIN_BUILDING_DETAIL,
            min_building_area=MIN_BUILDING_AREA,
            report_mesh_quality=False,
            mesher=mesher,
        )
        quality = json_ready(mesh.quality())
        metrics = mesh_metrics(mesh, quality)
        filename = f"{number:03d}.{mesher}.vtu"
        mesh.save(str(output_dir / filename))

        return (
            {
                "status": "success",
                "time": round(time.perf_counter() - start, 2),
                "file": filename,
                "quality": quality,
                "metrics": metrics,
            },
            mesh,
        )
    except Exception as exc:
        return (
            {
                "status": "failed",
                "time": round(time.perf_counter() - start, 2),
                "error": error_details(exc),
            },
            None,
        )


def _load_plot_modules():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection, PolyCollection
        from matplotlib.patches import Patch
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting") from exc

    return plt, LineCollection, PolyCollection, Patch


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
    padding = 0.05 * BOX_SIZE
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
    show_plot: bool = False,
) -> None:
    plt, line_collection_cls, poly_collection_cls, patch_cls = _load_plot_modules()

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
        f"({bounds.xmax:.0f}, {bounds.ymax:.0f})"
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
    grid = np.full((NY, NX), np.nan)

    for number, record in results.items():
        mesher_result = record.get("meshers", {}).get(mesher)
        if not mesher_result or mesher_result.get("status") != "success":
            continue

        value = float(mesher_result["metrics"][metric_key])
        ix, iy = case_to_grid(number)
        grid[iy, ix] = value

    return grid


def annotate_heatmap(ax, grid: np.ndarray, fmt: str) -> None:
    valid = grid[~np.isnan(grid)]
    if valid.size == 0:
        return
    midpoint = 0.5 * (float(valid.min()) + float(valid.max()))
    for iy in range(NY):
        for ix in range(NX):
            value = grid[iy, ix]
            if np.isnan(value):
                continue
            ax.text(
                ix,
                iy,
                format(value, fmt),
                ha="center",
                va="center",
                fontsize=6,
                color="black" if value >= midpoint else "white",
            )


def plot_overview(results: dict[int, dict[str, Any]], meshers: list[str], output_path: Path) -> None:
    plt, _, _, _ = _load_plot_modules()

    fig, axes = plt.subplots(
        len(meshers),
        len(OVERVIEW_METRICS),
        figsize=(5.4 * len(OVERVIEW_METRICS), 3.9 * len(meshers)),
        constrained_layout=True,
    )
    if len(meshers) == 1:
        axes = np.asarray([axes])

    for row, mesher in enumerate(meshers):
        for col, (metric_key, title, cmap, fmt) in enumerate(OVERVIEW_METRICS):
            ax = axes[row, col]
            grid = build_grid(results, mesher, metric_key)
            image = ax.imshow(grid, origin="lower", cmap=cmap, aspect="equal")
            annotate_heatmap(ax, grid, fmt)
            ax.set_title(f"{mesher} - {title}")
            ax.set_xticks(range(NX))
            ax.set_yticks(range(NY))
            ax.set_xlabel("Grid X")
            ax.set_ylabel("Grid Y")
            fig.colorbar(image, ax=ax, shrink=0.82, pad=0.03)

    fig.suptitle("Mesh quality overview by mesher")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def compact_case_status(record: dict[str, Any], meshers: list[str]) -> str:
    parts = []
    for mesher in meshers:
        result = record["meshers"].get(mesher, {})
        if result.get("status") != "success":
            parts.append(f"{mesher}: FAIL")
            continue
        metrics = result["metrics"]
        parts.append(
            f"{mesher}: ARmax {metrics['aspect_ratio_max']:.2f}, "
            f"edge p01 {metrics['edge_length_p01']:.2f}m"
        )
    return "  ".join(parts)


def run_case(
    number: int,
    meshers: list[str],
    output_dir: Path,
    create_case_plot: bool,
    show_plot: bool,
) -> dict[str, Any]:
    ix, iy = case_to_grid(number)
    bounds = make_bounds(ix, iy)

    city, prepare_time = prepare_city(bounds)

    case_record: dict[str, Any] = {
        "number": number,
        "bounds": bounds_to_dict(bounds),
        "prepare_time": round(prepare_time, 2),
        "meshers": {},
    }
    mesh_objects: dict[str, Mesh] = {}

    for mesher in meshers:
        result, mesh = run_mesher_for_case(number, city, mesher, output_dir)
        case_record["meshers"][mesher] = json_ready(result)
        if mesh is not None:
            mesh_objects[mesher] = mesh

    if create_case_plot:
        plot_path = case_plot_path(output_dir, number, meshers)
        plot_case_comparison(
            plot_path,
            number=number,
            bounds=bounds,
            meshers=meshers,
            mesh_objects=mesh_objects,
            mesher_results=case_record["meshers"],
            show_plot=show_plot,
        )
        case_record["plot_file"] = plot_path.name

    return case_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare flat-mesh quality across Stockholm tiles for multiple meshers."
    )
    parser.add_argument(
        "case_number",
        nargs="?",
        type=int,
        help=f"Single case to recompute (1-{NX * NY}). Omit to process all missing cases.",
    )
    parser.add_argument(
        "--meshers",
        nargs="+",
        default=None,
        help="Meshers to compare. Defaults to all available from: dtcc_mesher triangle spade.",
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
        default=Path(__file__).resolve().parent / "output",
        help="Directory for JSON, VTU, and PNG outputs.",
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
    args = parser.parse_args()

    total = NX * NY
    if args.case_number is not None and not (1 <= args.case_number <= total):
        parser.error(f"case_number must be between 1 and {total}")

    return args


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    meshers = resolve_requested_meshers(args.meshers)
    results_path = results_file_path(output_dir, meshers)
    results = load_results(results_path)
    total = NX * NY

    loaded_count = len(results)
    if loaded_count > 0:
        print(f"Loaded {loaded_count} previous case record(s) from {results_path}")
        print("  Delete that file to recompute everything from scratch.")
        print()

    if args.case_number is not None:
        cases_to_run = [args.case_number]
        if args.case_number in results:
            print(f"Case {args.case_number} exists in cache and will be recomputed.")
            print()
    else:
        cases_to_run = [
            number
            for number in range(1, total + 1)
            if not case_complete(results.get(number, {}), meshers)
        ]
        if not cases_to_run:
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
            and (args.case_number is not None or args.per_case_plots)
        )
        case_record = run_case(
            number,
            meshers=meshers,
            output_dir=output_dir,
            create_case_plot=create_case_plot,
            show_plot=bool(args.show_plot and args.case_number is not None),
        )
        results[number] = case_record
        save_results(results_path, results)

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

    print()
    print(f"Results file: {results_path}")

    if args.case_number is not None:
        print_case_summary(results[args.case_number], meshers)
    else:
        incomplete_cases = [
            number for number, record in results.items() if not case_complete(record, meshers)
        ]
        print_mesher_summary(results, meshers)
        if incomplete_cases:
            print()
            print(
                "Incomplete cases: "
                + ", ".join(str(number) for number in sorted(incomplete_cases))
            )

    if not args.no_plots:
        overview_path = overview_plot_path(output_dir, meshers)
        plot_overview(results, meshers, overview_path)
        print(f"Overview plot: {overview_path}")


if __name__ == "__main__":
    main()
