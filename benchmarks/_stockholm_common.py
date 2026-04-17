from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

import dtcc_core
from dtcc_core.model import Bounds, City

X_MIN = 673_000
Y_MIN = 6_578_500
NX = 10
NY = 10
BOX_SIZE = 500

DEFAULT_MAX_MESH_SIZE = 10.0
MIN_MESH_ANGLE = 25.0
MIN_BUILDING_DETAIL = 0.5
MIN_BUILDING_AREA = 15.0
MERGE_BUILDINGS = True

RASTER_CELL_SIZE = 2.0
RASTER_RADIUS = 3.0
OUTLIER_THRESHOLD = 3.0

DEFAULT_DELAY_BETWEEN_CASES = 8.0


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def benchmark_root() -> Path:
    return repo_root() / "benchmarks"


def stockholm_output_dir(name: str) -> Path:
    return benchmark_root() / name


def total_case_count() -> int:
    return NX * NY


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


def validate_case_number(parser: argparse.ArgumentParser, number: int | None) -> None:
    if number is None:
        return
    total = total_case_count()
    if not (1 <= number <= total):
        parser.error(f"case_number must be between 1 and {total}")


def validate_case_numbers(numbers: list[int]) -> None:
    total = total_case_count()
    invalid = [number for number in numbers if number < 1 or number > total]
    if invalid:
        joined = ", ".join(str(number) for number in invalid)
        raise ValueError(f"Invalid case numbers: {joined}. Expected values in 1..{total}.")


def add_cases_argument(
    parser: argparse.ArgumentParser,
    *,
    default: list[int] | None = None,
) -> None:
    parser.add_argument(
        "--cases",
        nargs="+",
        type=int,
        default=default,
        help="One or more Stockholm grid case numbers (1-100). Omit to use the script default.",
    )


def resolve_case_numbers(cases: list[int] | None) -> list[int]:
    if not cases:
        return list(range(1, total_case_count() + 1))
    ordered: list[int] = []
    for number in cases:
        if number not in ordered:
            ordered.append(int(number))
    validate_case_numbers(ordered)
    return ordered


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


def parse_max_mesh_size_argument(raw: str) -> float | None:
    value = raw.strip().lower()
    if value in {"none", "inf", "unrestricted"}:
        return None
    return positive_float(raw)


def normalize_max_mesh_size(max_mesh_size: float | None) -> float | None:
    if max_mesh_size is None:
        return None
    value = float(max_mesh_size)
    if value <= 0.0:
        return None
    return value


def mesh_size_label(max_mesh_size: float | None) -> str:
    normalized = normalize_max_mesh_size(max_mesh_size)
    if normalized is None:
        return "unrestricted"
    return f"maxh={normalized:g}m"


def slug_token(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


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
    if isinstance(data, dict) and isinstance(data.get("cases"), dict):
        data = data["cases"]
    return {int(key): value for key, value in data.items()}


def save_results(
    path: Path,
    results: dict[int, dict[str, Any]],
    *,
    metadata: dict[str, Any] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = {str(key): results[key] for key in sorted(results)}
    payload: dict[str, Any] | dict[str, dict[str, Any]]
    if metadata is None:
        payload = ordered
    else:
        payload = {**json_ready(metadata), "cases": ordered}
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2)


def quantile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values, q))


def load_plot_modules():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection, PolyCollection
        from matplotlib.patches import Patch
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting") from exc

    return plt, LineCollection, PolyCollection, Patch


def format_console_table(
    headers: list[str],
    rows: list[list[Any]],
    *,
    title: str | None = None,
) -> str:
    formatted_rows = [[str(item) for item in row] for row in rows]
    widths = [len(header) for header in headers]
    for row in formatted_rows:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))

    def render_row(values: list[str]) -> str:
        return " | ".join(value.ljust(widths[index]) for index, value in enumerate(values))

    separator = "-+-".join("-" * width for width in widths)
    lines: list[str] = []
    if title:
        lines.append(title)
    lines.append(render_row(headers))
    lines.append(separator)
    for row in formatted_rows:
        lines.append(render_row(row))
    return "\n".join(lines)


def annotate_heatmap(ax, grid: np.ndarray, fmt: str) -> None:
    finite_values = grid[np.isfinite(grid)]
    if finite_values.size == 0:
        threshold = 0.0
    else:
        threshold = float(np.nanmedian(finite_values))

    for row in range(grid.shape[0]):
        for col in range(grid.shape[1]):
            value = grid[row, col]
            if not np.isfinite(value):
                continue
            color = "white" if value >= threshold else "black"
            ax.text(col, row, format(value, fmt), ha="center", va="center", color=color)


def prepare_city(
    bounds: Bounds,
    *,
    raster_cell_size: float = RASTER_CELL_SIZE,
    raster_radius: float = RASTER_RADIUS,
    outlier_threshold: float = OUTLIER_THRESHOLD,
) -> tuple[City, float]:
    start = time.perf_counter()

    pointcloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
    buildings = dtcc_core.io.data.download_footprints(bounds=bounds)
    pointcloud = pointcloud.remove_global_outliers(outlier_threshold)

    raster = dtcc_core.builder.build_terrain_raster(
        pointcloud,
        cell_size=raster_cell_size,
        radius=raster_radius,
        ground_only=True,
    )
    buildings = dtcc_core.builder.extract_roof_points(buildings, pointcloud)
    buildings = dtcc_core.builder.compute_building_heights(buildings, raster, overwrite=True)

    city = City()
    city.add_terrain(raster)
    city.add_buildings(buildings, remove_outside_terrain=True)

    return city, time.perf_counter() - start
