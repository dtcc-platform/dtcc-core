"""Measure what a cleaning floor costs the mesh, by meshing it; research only.

The nearest-feature cost model this replaces charged a whole canonical edge at
the closest distance found anywhere along it, so one short sampling edge billed
an entire wall as a narrow contact. Checked against real meshes it overstated a
genuine 5 cm wall contact by about 2.5 times and an isolated 5 cm edge on an
otherwise straight wall by about 8. Because the corpus's defect mass is mostly
sampling, the total was inflated by exactly what it should not have charged.
That model and the delta it recommended are withdrawn.

This measures the same quantity without a model. Clean the fixed corpus at a
ladder of floors, mesh the results with the pinned mesher, and count elements.
No extrapolation, no per-family coefficient and no contact-length assumption.

Comparability is the whole point, so the ground domain is derived from the
*raw* input and never from the cleaned output: the rectangle a group or a city
is meshed on is identical at every floor, and only the footprints inside it
differ. Group rows are compared only across floors where the group conforms at
every floor on the ladder; city rows only for cities complete at that floor.

  .venv/bin/python sandbox/cleaning_delta_cost_probe.py --output /tmp/delta-cost \
      --run benchmarks/runs/2026-09-17_100049_quick --delta 0.5
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from shapely import box, wkb
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning.contract import check_cleaning_contract  # noqa: E402
from dtcc_core.builder.meshing.dtcc_mesher_backend import (  # noqa: E402
    build_city_flat_mesh_with_dtcc_mesher,
)
from sandbox.cleaning_graph_prototype import polygon_parts  # noqa: E402
from sandbox.cleaning_staged_probe import construct  # noqa: E402

# The corrected city meshes used a 10 m maximum edge length. This build of the
# pinned mesher refuses a maximum edge length on a coverage whose ground has a
# building nested in it - "invalid mesh topology", the same defect as the two
# meshing integration tests that fail here and not on macOS. Quality-driven
# refinement alone is used instead. It works on every case here, it keeps the
# cost signal intact at every edge length tried, and it is arguably the cleaner
# measure: a global edge ceiling adds elements that have nothing to do with the
# cleaning floor. The absolute counts are therefore not comparable with the
# corrected city meshes; the comparison across floors, which is the point, is.
GROUND_PADDING = 10.0
MIN_MESH_ANGLE = 20.0
MAX_EDGE_LENGTH = None
GROUND_MARKER = -2
POOR_ELEMENT_QUALITY = 0.02
ATTEMPT_SECONDS = 90.0


def domain(raw):
    """The ground rectangle, from the raw input, so a floor cannot move it."""
    minx, miny, maxx, maxy = unary_union(list(raw)).bounds
    return box(
        minx - GROUND_PADDING,
        miny - GROUND_PADDING,
        maxx + GROUND_PADDING,
        maxy + GROUND_PADDING,
    )


def _mesh(regions, markers):
    from dtcc_core.model.mixins.mesh.quality import tri_element_quality

    started = time.perf_counter()
    try:
        mesh = build_city_flat_mesh_with_dtcc_mesher(
            region_polygons=regions,
            region_markers=markers,
            max_mesh_size=MAX_EDGE_LENGTH,
            min_mesh_angle=MIN_MESH_ANGLE,
        )
    except Exception as error:
        # An admissible coverage the mesher still refuses is a result about the
        # mesher, not a reason to drop the row. It is recorded and excluded
        # from the cost curve rather than counted as free.
        return {
            "outcome": "raised",
            "error_type": type(error).__name__,
            "error": str(error),
            "seconds": time.perf_counter() - started,
        }
    seconds = time.perf_counter() - started
    faces = np.asarray(mesh.faces, dtype=int)
    if not len(faces):
        return {"outcome": "empty", "seconds": seconds}
    quality = tri_element_quality(mesh.vertices, mesh.faces)
    return {
        "outcome": "meshed",
        "seconds": seconds,
        "faces": int(len(faces)),
        "vertices": int(len(mesh.vertices)),
        "element_quality_min": float(np.min(quality)),
        "degenerate_faces": int(np.count_nonzero(quality <= 0)),
        "element_quality_below_threshold": int(
            np.count_nonzero(quality < POOR_ELEMENT_QUALITY)
        ),
    }


def _worker(queue, regions, markers):
    queue.put(_mesh(regions, markers))


def mesh_on(ground, cleaned):
    """Mesh a cleaned coverage on a fixed ground rectangle, under a wall clock.

    Refinement around a very small feature need not terminate, and a mesher
    that does not terminate is a result rather than a hang to wait out.
    """
    import multiprocessing

    parts = [
        orient(polygon, sign=1.0)
        for polygon in (polygon_parts(unary_union(list(cleaned))) if cleaned else [])
    ]
    ground_parts = (
        polygon_parts(ground.difference(unary_union(parts))) if parts else [ground]
    )
    regions = [*ground_parts, *parts]
    markers = [GROUND_MARKER] * len(ground_parts) + list(range(len(parts)))
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    process = context.Process(target=_worker, args=(queue, regions, markers))
    started = time.perf_counter()
    process.start()
    try:
        return queue.get(timeout=ATTEMPT_SECONDS)
    except Exception:
        alive = process.is_alive()
        return {
            "outcome": "did_not_terminate" if alive else "crashed",
            "seconds": time.perf_counter() - started,
        }
    finally:
        process.terminate()
        process.join()


def survey(run, destination, delta, epsilon, budget=None):
    """Clean and mesh every group at one floor, appending rows as they finish."""
    from sandbox.cleaning_pipeline_probe import saved_city_groups

    started = time.perf_counter()
    done, cached = set(), {}
    if destination.exists():
        for line in destination.read_text().splitlines():
            row = json.loads(line)
            done.add((row["kind"], row["city"], row.get("group")))
            if row.get("geometry"):
                cached[row["city"], row["group"]] = wkb.loads(
                    bytes.fromhex(row["geometry"])
                )
    for record, _, _, groups in saved_city_groups(run):
        city = record["case"]["city"]
        assembled, complete = [], True
        for index, raw in enumerate(groups):
            if ("group", city, index) in done:
                # Each row carries its cleaned geometry, so a resumed pass can
                # still assemble the city. Without it, whether a city gets an
                # assembly check would depend on where the budget happened to
                # stop, which is not a property of the measurement.
                if (city, index) in cached:
                    assembled.extend(polygon_parts(cached[city, index]))
                else:
                    complete = False
                continue
            if budget is not None and time.perf_counter() - started > budget:
                print(f"budget reached before {city} group {index}", flush=True)
                return False
            output, report = construct(raw, delta=delta, epsilon=epsilon)
            row = {
                "kind": "group",
                "delta": delta,
                "epsilon": epsilon,
                "city": city,
                "group": index,
                "outcome": report["outcome"],
                "construction_seconds": round(report["seconds"], 6),
                "initial_vertices": report["initial"][2],
            }
            if output is None:
                complete = False
            else:
                assembled.extend(output)
                row["geometry"] = unary_union(output).wkb_hex
                row["mesh"] = mesh_on(domain(raw), output)
            with destination.open("a") as handle:
                handle.write(json.dumps(row) + "\n")
        if complete and ("city", city, None) not in done:
            whole = check_cleaning_contract(
                [p for group in groups for p in group],
                assembled,
                delta=delta,
                epsilon=epsilon,
            )
            row = {
                "kind": "city",
                "delta": delta,
                "epsilon": epsilon,
                "city": city,
                "group": None,
                "contract": whole["status"],
                "mesh": mesh_on(
                    domain([p for group in groups for p in group]), assembled
                ),
            }
            with destination.open("a") as handle:
                handle.write(json.dumps(row) + "\n")
            print(
                f"{city}: complete, city contract {whole['status']}, "
                f"{row['mesh'].get('faces')} faces",
                flush=True,
            )
        elif not complete:
            print(f"{city}: incomplete at delta {delta:g}", flush=True)
    return True


def provenance():
    import hashlib
    import platform

    import numpy
    import shapely
    from importlib.metadata import distribution

    try:
        mesher = json.loads(distribution("dtcc-mesher").read_text("direct_url.json"))
    except Exception:  # pragma: no cover - only when installed from a wheel
        mesher = None
    return {
        "scripts": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in (
                "sandbox/cleaning_delta_cost_probe.py",
                "sandbox/cleaning_staged_probe.py",
                "dtcc_core/builder/cleaning/contract.py",
            )
        },
        "versions": {
            "shapely": shapely.__version__,
            "geos": ".".join(str(n) for n in shapely.geos_version),
            "numpy": numpy.__version__,
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "mesher": mesher,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--delta", required=True, type=float)
    parser.add_argument("--epsilon", type=float, default=None)
    parser.add_argument("--budget", type=float, default=None)
    args = parser.parse_args()
    if not np.isfinite(args.delta) or args.delta <= 0:
        parser.error("delta must be finite and positive")
    args.output.mkdir(parents=True, exist_ok=True)
    epsilon = args.delta / 2 if args.epsilon is None else args.epsilon
    lines = args.output / f"rows-delta{args.delta:g}.jsonl"
    complete = survey(args.run, lines, args.delta, epsilon, budget=args.budget)
    rows = [json.loads(line) for line in lines.read_text().splitlines()]
    groups = [r for r in rows if r["kind"] == "group"]
    print(
        f"delta {args.delta:g}: "
        f"{sum(r['outcome'] in ('conforming', 'unchanged') for r in groups)}"
        f"/{len(groups)} conforming; complete={complete}",
        flush=True,
    )
    if not complete:
        return
    (args.output / f"cost-delta{args.delta:g}.json").write_text(
        json.dumps(
            {
                "scope": (
                    "Measured mesh cost of a cleaning floor. The ground domain "
                    "comes from the raw input so it is identical at every "
                    "floor. Research only."
                ),
                "delta": args.delta,
                "epsilon": epsilon,
                "ground_padding": GROUND_PADDING,
                "min_mesh_angle": MIN_MESH_ANGLE,
                "max_edge_length": MAX_EDGE_LENGTH,
                "run": str(args.run),
                **provenance(),
                "rows": rows,
            },
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
