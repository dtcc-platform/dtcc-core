"""
Survey of flat mesh quality across central Stockholm.

Overview
--------
This script generates flat 2D meshes (at z = 0) for a 10x10 grid of 500 m
bounding boxes covering central Stockholm, evaluates mesh quality for each,
and prints a summary table with global statistics.

The grid spans 5 km x 5 km in SWEREF99 TM (EPSG:3006):
    x: 673 000 -> 678 000   (Kungsholmen -> Ostermalm)
    y: 6 578 500 -> 6 583 500   (Sodermalm -> Vasastan)

How it works
------------
1.  On each run the script loads previous results from
    output/mesh_quality_survey.json (if the file exists).

2.  It determines which cases still need to be computed:
    - Full run (no arguments):  all cases NOT already in the results file.
    - Single case (e.g. ``python mesh_quality_survey.py 55``):  only that
      case, even if it is already in the results file (useful for debugging).

3.  For each case to compute it calls ``dtcc_core.datasets.city_flat_mesh()``
    to download data and build the mesh, then evaluates ``mesh.quality()``.

4.  Each successfully built mesh is saved to output/<NNN>.vtu so it can be
    opened in ParaView for visual inspection.

5.  Successfully computed quality metrics are merged into the results file
    and saved back.  Failed cases are NOT saved, so the next full run
    will automatically retry them.

6.  Finally the script presents ALL results (loaded + newly computed) in a
    table with per-mesh metrics and global statistics.

Rerunning
---------
- To retry only failed cases:  just run the script again with no arguments.
- To rerun a specific case:    ``python mesh_quality_survey.py <n>``
- To rerun everything:         delete output/mesh_quality_survey.json first.

Coordinate system: SWEREF99 TM (EPSG:3006)
"""

import json
import os
import sys
import time
import traceback

import numpy as np

import dtcc_core
from dtcc_core.model import Bounds
from dtcc_core.model.mixins.mesh.quality import format_quality

# ── Configuration ────────────────────────────────────────────────────────

# Central Stockholm extent (SWEREF99 TM)
X_MIN = 673_000
Y_MIN = 6_578_500
NX, NY = 10, 10
BOX_SIZE = 500  # metres

# Meshing parameters
MAX_MESH_SIZE = 10.0
MIN_MESH_ANGLE = 25.0

# Delay between meshes (seconds) to avoid rate-limiting by dtcc-data.
# The server defaults to 5 requests per 30 s per IP.  Each mesh needs
# ~2 requests (pointcloud + footprints), so 8 s is a safe minimum.
# Set to 0 if rate-limiting is disabled (ENABLE_RATE_LIMIT=false).
DELAY_BETWEEN_MESHES = 8

# Output paths
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
RESULTS_FILE = os.path.join(OUTPUT_DIR, "mesh_quality_survey.json")

# ── Helpers ──────────────────────────────────────────────────────────────


def case_to_grid(number):
    """Convert 1-based case number to (ix, iy) grid indices."""
    iy, ix = divmod(number - 1, NX)
    return ix, iy


def make_bounds(ix, iy):
    """Return a Bounds for grid cell (ix, iy)."""
    x0 = X_MIN + ix * BOX_SIZE
    y0 = Y_MIN + iy * BOX_SIZE
    return Bounds(x0, y0, x0 + BOX_SIZE, y0 + BOX_SIZE)


def bounds_to_dict(b):
    """Serialize a Bounds to a plain dict."""
    return {"xmin": b.xmin, "ymin": b.ymin, "xmax": b.xmax, "ymax": b.ymax}


def dict_to_bounds(d):
    """Deserialize a Bounds from a plain dict."""
    return Bounds(d["xmin"], d["ymin"], d["xmax"], d["ymax"])


def load_results():
    """Load previous results from the JSON file.  Returns dict keyed by
    case number (int)."""
    if not os.path.exists(RESULTS_FILE):
        return {}
    with open(RESULTS_FILE) as f:
        data = json.load(f)
    # JSON keys are strings; convert to int
    return {int(k): v for k, v in data.items()}


def save_results(results_dict):
    """Save results dict to JSON (keys are case numbers)."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    # Sort by case number for readability
    ordered = dict(sorted(results_dict.items()))
    with open(RESULTS_FILE, "w") as f:
        json.dump(ordered, f, indent=2)


def run_case(number):
    """Build mesh for a single case number.  Returns result dict on success."""
    ix, iy = case_to_grid(number)
    bounds = make_bounds(ix, iy)

    t1 = time.time()
    mesh = dtcc_core.datasets.city_flat_mesh(
        bounds=bounds,
        max_mesh_size=MAX_MESH_SIZE,
        min_mesh_angle=MIN_MESH_ANGLE,
    )
    q = mesh.quality()
    dt = time.time() - t1

    # Save mesh file
    filename = f"{number:03d}.vtu"
    mesh.save(os.path.join(OUTPUT_DIR, filename))

    return {
        "number": number,
        "bounds": bounds_to_dict(bounds),
        "quality": q,
        "time": round(dt, 2),
        "file": filename,
    }


def print_summary_table(results_dict):
    """Print a combined table of all mesh quality results."""

    # Sort by case number
    items = sorted(results_dict.items())
    if not items:
        print("No results to display.")
        return

    sep = "-" * 115

    print()
    print(
        f"{'#':>4}  {'Bounds (xmin, ymin)':>22}  {'Cells':>8}  "
        f"{'ElemQ min':>9}  {'ElemQ mean':>10}  "
        f"{'AR max':>8}  {'AR mean':>8}  "
        f"{'ER max':>8}  {'ER mean':>8}  "
        f"{'Skew max':>8}"
    )
    print(sep)

    eq_mins, eq_means = [], []
    ar_maxs, ar_means = [], []
    er_maxs, er_means = [], []
    sk_maxs = []

    for num, r in items:
        q = r["quality"]
        eq = q["element_quality"]
        ar = q["aspect_ratio"]
        er = q["edge_ratio"]
        sk = q["skewness"]

        eq_mins.append(eq["min"])
        eq_means.append(eq["mean"])
        ar_maxs.append(ar["max"])
        ar_means.append(ar["mean"])
        er_maxs.append(er["max"])
        er_means.append(er["mean"])
        sk_maxs.append(sk["max"])

        b = r["bounds"]
        print(
            f"{num:>4}  ({b['xmin']:10.0f}, {b['ymin']:10.0f})  "
            f"{q['num_cells']:>8}  "
            f"{eq['min']:>9.4f}  {eq['mean']:>10.4f}  "
            f"{ar['max']:>8.4f}  {ar['mean']:>8.4f}  "
            f"{er['max']:>8.4f}  {er['mean']:>8.4f}  "
            f"{sk['max']:>8.4f}"
        )

    print(sep)

    # Global statistics
    total_cells = sum(r["quality"]["num_cells"] for _, r in items)
    print()
    print("Global statistics across all meshes")
    print(sep)
    print(f"  Total meshes           : {len(items)}")
    print(f"  Total cells            : {total_cells}")
    print()
    print(f"  Element quality  min   : {np.min(eq_mins):.4f}  (worst single element)")
    print(
        f"  Element quality  mean  : {np.mean(eq_means):.4f}  (avg of per-mesh means)"
    )
    print(f"  Aspect ratio     max   : {np.max(ar_maxs):.4f}  (worst single element)")
    print(
        f"  Aspect ratio     mean  : {np.mean(ar_means):.4f}  (avg of per-mesh means)"
    )
    print(f"  Edge ratio       max   : {np.max(er_maxs):.4f}  (worst single element)")
    print(
        f"  Edge ratio       mean  : {np.mean(er_means):.4f}  (avg of per-mesh means)"
    )
    print(f"  Skewness         max   : {np.max(sk_maxs):.4f}  (worst single element)")
    print(sep)


# ── Main ─────────────────────────────────────────────────────────────────


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    total = NX * NY

    # ── Parse command-line arguments ─────────────────────────────────
    single_case = None
    if len(sys.argv) > 1:
        try:
            single_case = int(sys.argv[1])
        except ValueError:
            print(f"Usage: python {os.path.basename(__file__)} [<case_number>]")
            print(f"  case_number: 1-{total}")
            sys.exit(1)
        if single_case < 1 or single_case > total:
            print(f"Error: case number must be between 1 and {total}")
            sys.exit(1)

    # ── Load existing results ────────────────────────────────────────
    results = load_results()
    loaded_count = len(results)

    if loaded_count > 0:
        print(f"Loaded {loaded_count} previous result(s) from {RESULTS_FILE}")
        print(f"  (delete that file and rerun to regenerate everything)")
        print()

    # ── Determine which cases to run ─────────────────────────────────
    if single_case is not None:
        cases_to_run = [single_case]
        if single_case in results:
            print(f"Case {single_case} exists in results — will recompute it.")
    else:
        all_cases = set(range(1, total + 1))
        done_cases = set(results.keys())
        cases_to_run = sorted(all_cases - done_cases)
        if not cases_to_run:
            print(f"All {total} cases already computed.")
        else:
            print(
                f"{len(cases_to_run)} case(s) to compute: "
                f"{', '.join(str(c) for c in cases_to_run)}"
            )
        print()

    # ── Run cases ────────────────────────────────────────────────────
    new_count = 0
    fail_count = 0
    failures = []
    t0 = time.time()

    for i, number in enumerate(cases_to_run):
        if i > 0 and DELAY_BETWEEN_MESHES > 0:
            time.sleep(DELAY_BETWEEN_MESHES)

        ix, iy = case_to_grid(number)
        bounds = make_bounds(ix, iy)
        label = (
            f"[{i + 1}/{len(cases_to_run)}]  case {number:3d}  "
            f"({bounds.xmin:.0f}, {bounds.ymin:.0f}) -> "
            f"({bounds.xmax:.0f}, {bounds.ymax:.0f})"
        )

        try:
            print(f"{label}  generating...", end="", flush=True)
            result = run_case(number)
            q = result["quality"]
            eq = q["element_quality"]
            print(
                f"  {q['num_cells']:>6} cells  "
                f"elemQ {eq['min']:.3f}/{eq['mean']:.3f}  "
                f"{result['time']:.1f}s"
            )
            results[number] = result
            new_count += 1

            # Save after each success so progress is not lost
            save_results(results)

        except Exception as e:
            print(f"  FAILED: {e}")
            traceback.print_exc()
            failures.append({"number": number, "error": str(e)})
            fail_count += 1

    elapsed = time.time() - t0

    # ── Summary ──────────────────────────────────────────────────────
    if cases_to_run:
        print()
        print(
            f"Finished in {elapsed:.0f}s  "
            f"({new_count} succeeded, {fail_count} failed, "
            f"{loaded_count} previously cached)"
        )

    if failures:
        print()
        print(f"Failed cases ({len(failures)})  — rerun the script to retry:")
        for f in failures:
            print(f"  case {f['number']:3d}: {f['error']}")

    # ── Present all results ──────────────────────────────────────────
    missing = total - len(results)
    print()
    if loaded_count > 0 and new_count > 0:
        print(
            f"Results: {loaded_count} from file + {new_count} newly computed"
            f" = {len(results)}/{total}"
        )
    elif loaded_count > 0:
        print(f"Results: {len(results)}/{total} (all from file)")
    else:
        print(f"Results: {len(results)}/{total}")

    if missing > 0:
        missing_cases = sorted(set(range(1, total + 1)) - set(results.keys()))
        print(
            f"Missing cases ({missing}): "
            f"{', '.join(str(c) for c in missing_cases)}"
        )

    if results:
        print_summary_table(results)


if __name__ == "__main__":
    main()
