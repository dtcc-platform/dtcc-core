"""
Survey of flat mesh quality across central Stockholm.

Generates a 10×10 grid of 500m bounding boxes spanning central Stockholm
(Kungsholmen – Östermalm, Södermalm – Vasastan), builds a flat mesh for
each, and reports mesh quality metrics.

Meshes are saved to sandbox/output/ numbered 001–100 so they can be
opened for visual inspection.

Coordinate system: SWEREF99 TM (EPSG:3006)
"""

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
#   x  673 000 → 678 000  (5 km, covers Kungsholmen – Östermalm)
#   y  6 578 500 → 6 583 500  (5 km, covers Södermalm – Vasastan)
X_MIN = 673_000
Y_MIN = 6_578_500
NX, NY = 10, 10
BOX_SIZE = 500  # metres

# Meshing parameters
MAX_MESH_SIZE = 10.0
MIN_MESH_ANGLE = 25.0

# Output
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

# ── Helpers ──────────────────────────────────────────────────────────────


def make_bounds(ix, iy):
    """Return a Bounds for grid cell (ix, iy)."""
    x0 = X_MIN + ix * BOX_SIZE
    y0 = Y_MIN + iy * BOX_SIZE
    return Bounds(x0, y0, x0 + BOX_SIZE, y0 + BOX_SIZE)


def print_summary_table(results):
    """Print a combined table of all mesh quality results."""

    sep = "-" * 115

    # Header
    print()
    print(f"{'#':>4}  {'Bounds (xmin, ymin)':>22}  {'Cells':>8}  "
          f"{'ElemQ min':>9}  {'ElemQ mean':>10}  "
          f"{'AR max':>8}  {'AR mean':>8}  "
          f"{'ER max':>8}  {'ER mean':>8}  "
          f"{'Skew max':>8}")
    print(sep)

    eq_mins, eq_means = [], []
    ar_maxs, ar_means = [], []
    er_maxs, er_means = [], []
    sk_maxs = []

    for r in results:
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
            f"{r['number']:>4}  ({b.xmin:10.0f}, {b.ymin:10.0f})  {q['num_cells']:>8}  "
            f"{eq['min']:>9.4f}  {eq['mean']:>10.4f}  "
            f"{ar['max']:>8.4f}  {ar['mean']:>8.4f}  "
            f"{er['max']:>8.4f}  {er['mean']:>8.4f}  "
            f"{sk['max']:>8.4f}"
        )

    print(sep)

    # Global statistics
    print()
    print("Global statistics across all meshes")
    print(sep)
    total_cells = sum(r["quality"]["num_cells"] for r in results)
    print(f"  Total meshes generated : {len(results)}")
    print(f"  Total cells            : {total_cells}")
    print()
    print(f"  Element quality  min   : {np.min(eq_mins):.4f}  (worst single element)")
    print(f"  Element quality  mean  : {np.mean(eq_means):.4f}  (avg of per-mesh means)")
    print(f"  Aspect ratio     max   : {np.max(ar_maxs):.4f}  (worst single element)")
    print(f"  Aspect ratio     mean  : {np.mean(ar_means):.4f}  (avg of per-mesh means)")
    print(f"  Edge ratio       max   : {np.max(er_maxs):.4f}  (worst single element)")
    print(f"  Edge ratio       mean  : {np.mean(er_means):.4f}  (avg of per-mesh means)")
    print(f"  Skewness         max   : {np.max(sk_maxs):.4f}  (worst single element)")
    print(sep)


# ── Main ─────────────────────────────────────────────────────────────────


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    results = []
    failures = []
    total = NX * NY
    number = 0
    t0 = time.time()

    for iy in range(NY):
        for ix in range(NX):
            number += 1
            bounds = make_bounds(ix, iy)
            label = (f"[{number:3d}/{total}]  "
                     f"({bounds.xmin:.0f}, {bounds.ymin:.0f}) – "
                     f"({bounds.xmax:.0f}, {bounds.ymax:.0f})")

            try:
                print(f"{label}  generating...", end="", flush=True)
                t1 = time.time()

                mesh = dtcc_core.datasets.city_flat_mesh(
                    bounds=bounds,
                    max_mesh_size=MAX_MESH_SIZE,
                    min_mesh_angle=MIN_MESH_ANGLE,
                )

                q = mesh.quality()
                dt = time.time() - t1

                # Save mesh
                filename = f"{number:03d}.vtu"
                mesh.save(os.path.join(OUTPUT_DIR, filename))

                results.append({
                    "number": number,
                    "bounds": bounds,
                    "quality": q,
                    "time": dt,
                    "file": filename,
                })

                eq = q["element_quality"]
                print(f"  {q['num_cells']:>6} cells  "
                      f"elemQ {eq['min']:.3f}/{eq['mean']:.3f}  "
                      f"{dt:.1f}s")

            except Exception as e:
                dt = time.time() - t1
                failures.append({
                    "number": number,
                    "bounds": bounds,
                    "error": str(e),
                })
                print(f"  FAILED ({dt:.1f}s): {e}")
                traceback.print_exc()

    elapsed = time.time() - t0

    # ── Summary ──────────────────────────────────────────────────────
    print()
    print(f"Completed in {elapsed:.0f}s  "
          f"({len(results)} succeeded, {len(failures)} failed)")

    if results:
        print_summary_table(results)

    if failures:
        print()
        print(f"Failed meshes ({len(failures)}):")
        for f in failures:
            b = f["bounds"]
            print(f"  #{f['number']:03d}  ({b.xmin:.0f}, {b.ymin:.0f}): {f['error']}")


if __name__ == "__main__":
    main()
