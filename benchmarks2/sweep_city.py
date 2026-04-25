"""Swedish city parameter sweep harness.

Runs exactly one parameter-combination case per invocation against a supported
city center bbox, printing a machine-parseable SWEEP_RESULT sentinel on the last
line of stdout. Designed to be called by sweep_city.sh, which wraps invocations
in a watchdog.
"""

from __future__ import annotations

import signal
import sys
import time
import traceback
from typing import Iterable

import numpy as np

import dtcc_core as dtcc

# Explicitly install Python's default SIGINT handler. When this script is run
# as a backgrounded subprocess of a non-interactive shell, POSIX job control
# sets SIGINT = SIG_IGN on the async command. Python inherits SIG_IGN and will
# not override it during normal startup, so reinstall the default handler.
signal.signal(signal.SIGINT, signal.default_int_handler)

# City centers in EPSG:3006 / SWEREF99TM.
CITY_CENTERS: dict[str, tuple[str, float, float]] = {
    "lund": ("Lund", 386325.0, 6174697.0),
    "stockholm": ("Stockholm", 674571.9, 6580743.0),
    "gothenburg": ("Gothenburg", 319758.0, 6400326.0),
    "malmo": ("Malmo", 374243.8, 6163926.6),
    "uppsala": ("Uppsala", 647793.5, 6638608.1),
    "linkoping": ("Linkoping", 536308.5, 6474615.0),
    "orebro": ("Orebro", 512162.3, 6570727.0),
    "vasteras": ("Vasteras", 587172.6, 6608981.7),
    "helsingborg": ("Helsingborg", 356398.3, 6213652.0),
    "norrkoping": ("Norrkoping", 569321.0, 6494759.3),
}

# Fixed per-call arguments - identical in every case.
FIXED_MIN_MESH_ANGLE = 25.0

# Baseline values for axes not currently being swept.
BASELINE_RASTER_CELL_SIZE = 2.0
BASELINE_MIN_BUILDING_DETAIL = 0.5
BASELINE_MIN_BUILDING_AREA = 15.0
BASELINE_BBOX_SIZE_M = 500.0
BASELINE_MAX_MESH_SIZE = 10.0

RASTER_CELL_SIZE_SWEEP = [0.5, 1.0, 2.0, 5.0, 10.0]
MIN_BUILDING_DETAIL_SWEEP = [0.25, 0.5, 1.0, 2.0]
MIN_BUILDING_AREA_SWEEP = [1.0, 5.0, 15.0, 25.0, 100.0]
BBOX_SIZE_M_SWEEP = [50, 100, 200, 350, 500]
MAX_MESH_SIZE_SWEEP = [1, 2, 5, 10, 20]

GROUPS = ("raster", "detail", "area", "bbox", "max")


def bounds_for(city: str, bbox_size_m: float) -> "dtcc.Bounds":
    _, x0, y0 = CITY_CENTERS[city]
    half = bbox_size_m / 2.0
    return dtcc.Bounds(x0 - half, y0 - half, x0 + half, y0 + half)


def _build_case_registry() -> dict[str, tuple[str, float]]:
    registry: dict[str, tuple[str, float]] = {}
    for v in RASTER_CELL_SIZE_SWEEP:
        registry[f"raster_{v}"] = ("raster", v)
    for v in MIN_BUILDING_DETAIL_SWEEP:
        registry[f"detail_{v}"] = ("detail", v)
    for v in MIN_BUILDING_AREA_SWEEP:
        registry[f"area_{v}"] = ("area", v)
    for v in BBOX_SIZE_M_SWEEP:
        registry[f"bbox_{v}"] = ("bbox", v)
    for v in MAX_MESH_SIZE_SWEEP:
        registry[f"max_{v}"] = ("max", v)
    return registry


CASES: dict[str, tuple[str, float]] = _build_case_registry()


def _params_for(axis: str, value: float) -> dict[str, float]:
    params = {
        "raster_cell_size": BASELINE_RASTER_CELL_SIZE,
        "min_building_detail": BASELINE_MIN_BUILDING_DETAIL,
        "min_building_area": BASELINE_MIN_BUILDING_AREA,
        "bbox_size_m": BASELINE_BBOX_SIZE_M,
        "max_mesh_size": BASELINE_MAX_MESH_SIZE,
    }
    if axis == "raster":
        params["raster_cell_size"] = value
    elif axis == "detail":
        params["min_building_detail"] = value
    elif axis == "area":
        params["min_building_area"] = value
    elif axis == "bbox":
        params["bbox_size_m"] = value
    elif axis == "max":
        params["max_mesh_size"] = value
    else:
        raise ValueError(f"unknown axis: {axis}")
    return params


def _count_building_faces(mesh) -> int:
    markers = np.asarray(mesh.markers, dtype=np.int64)
    if markers.size == 0:
        return 0
    return int(np.sum(markers >= 0))


def run_case(city: str, case_id: str) -> int:
    if city not in CITY_CENTERS:
        valid = ", ".join(CITY_CENTERS)
        print(f"unknown city: {city}", file=sys.stderr)
        print(f"valid cities: {valid}", file=sys.stderr)
        print(
            f"SWEEP_RESULT city={city} case={case_id} status=error elapsed=0.00",
            flush=True,
        )
        return 2

    if case_id not in CASES:
        print(f"unknown case: {case_id}", file=sys.stderr)
        print(
            f"SWEEP_RESULT city={city} case={case_id} status=error elapsed=0.00",
            flush=True,
        )
        return 2

    axis, value = CASES[case_id]
    params = _params_for(axis, value)
    bounds = bounds_for(city, params["bbox_size_m"])

    start = time.perf_counter()
    try:
        mesh = dtcc.datasets.city_surface_mesh(
            bounds=bounds,
            max_mesh_size=params["max_mesh_size"],
            min_mesh_angle=FIXED_MIN_MESH_ANGLE,
            raster_cell_size=params["raster_cell_size"],
            min_building_detail=params["min_building_detail"],
            min_building_area=params["min_building_area"],
        )
    except BaseException:
        elapsed = time.perf_counter() - start
        traceback.print_exc()
        print(
            f"SWEEP_RESULT city={city} case={case_id} "
            f"status=error elapsed={elapsed:.2f}",
            flush=True,
        )
        return 1

    elapsed = time.perf_counter() - start
    vertices = int(mesh.num_vertices)
    faces = int(mesh.num_faces)
    building_faces = _count_building_faces(mesh)
    status = "ok" if building_faces > 0 else "ok_terrain_only"
    print(
        f"SWEEP_RESULT city={city} case={case_id} status={status} "
        f"elapsed={elapsed:.2f} vertices={vertices} faces={faces} "
        f"building_faces={building_faces}",
        flush=True,
    )
    return 0


def _cases_in_group(selector: str) -> Iterable[str]:
    if selector == "all":
        return list(CASES)
    if selector in GROUPS:
        return [cid for cid, (axis, _) in CASES.items() if axis == selector]
    return []


def main() -> int:
    if len(sys.argv) < 3:
        print(
            f"usage: {sys.argv[0]} <city> <case-id|{'|'.join(GROUPS)}|all>",
            file=sys.stderr,
        )
        print(f"cities: {', '.join(CITY_CENTERS)}", file=sys.stderr)
        print(f"cases: {', '.join(CASES)}", file=sys.stderr)
        return 2

    city = sys.argv[1]
    selector = sys.argv[2]

    if selector in CASES:
        return run_case(city, selector)

    group_cases = list(_cases_in_group(selector))
    if not group_cases:
        valid = ", ".join(list(CASES) + list(GROUPS) + ["all"])
        print(f"unknown selector: {selector}", file=sys.stderr)
        print(f"valid: {valid}", file=sys.stderr)
        return 2

    overall = 0
    for cid in group_cases:
        rc = run_case(city, cid)
        if rc != 0 and overall == 0:
            overall = rc
    return overall


if __name__ == "__main__":
    sys.exit(main())
