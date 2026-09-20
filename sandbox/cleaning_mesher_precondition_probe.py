"""What the mesher actually requires; never called by production cleaning.

Milestone D of
`.agent/plans/2026-09-18-footprint-cleaning-feasibility-and-architecture-decision.md`.
The contract promises uniform delta separation everywhere while its own mesher
section says the mesher owns refinement. This feeds deliberately sub-delta
configurations to core's pinned mesher through the production flat-coverage
entry point and finds where it actually fails or degrades.

Each configuration is meshed at several requested spacings. A failure that a
smaller spacing removes is one refinement could absorb; a failure that survives
every spacing is a genuine precondition. A global spacing is a crude stand-in
for local grading, and it is clamped from below, so the finest features answer
"not established" rather than "absorbable".

  .venv/bin/python sandbox/cleaning_mesher_precondition_probe.py --output /tmp/preconditions
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning.contract import admissibility  # noqa: E402
from dtcc_core.builder.meshing.dtcc_mesher_backend import (  # noqa: E402
    build_city_flat_mesh_with_dtcc_mesher,
)
from sandbox.cleaning_graph_prototype import polygon_parts  # noqa: E402

DELTA = 0.5
MIN_MESH_ANGLE = 25.0  # the saved corpus parameter, not a research choice
GROUND_MARKER = -2  # the marker geometry_builders/meshes.py gives ground
GROUND_PADDING = 2.0
# Element quality below this is already the flagged band in benchmark metrics.
POOR_ELEMENT_QUALITY = 0.02
SEPARATIONS = (
    1.0,
    0.5,
    0.25,
    0.125,
    0.1,
    0.05,
    0.025,
    0.01,
    0.005,
    0.002,
    0.001,
    0.0005,
    0.0001,
)
MESH_SIZES = (10.0, 1.0, 0.25)
FINEST_REQUESTED_SPACING = 0.02
# A global spacing stands in for local grading, so the finest attempt can ask
# for an absurd number of elements. Refuse those rather than let one wedge
# decide how long the sweep takes; the row then reports that refinement was not
# established, which is the honest answer.
MAX_REQUESTED_ELEMENTS = 400_000
# Refinement around a very small input angle need not terminate, and a hung
# mesher is itself a result. Each attempt therefore runs in its own process
# under a wall clock. The limit is above every attempt that has finished here,
# so it decides nothing that finishing decides.
ATTEMPT_SECONDS = 90.0


def configurations(separation):
    """Four sub-delta feature families, each shrinking with one parameter.

    Every family is a legal, valid coverage at any positive separation: what
    changes is the size of one geometric feature, so a mesher failure can be
    attributed to that feature and not to malformed input.
    """
    s = separation
    return {
        # Two buildings whose facing walls are `s` apart.
        "near_walls": [box(0, 0, 10, 10), box(10 + s, 0, 20 + s, 10)],
        # A courtyard whose entrance passage is `s` wide.
        "thin_passage": [
            Polygon(
                [
                    (0, 0),
                    (10, 0),
                    (10, 5 - s / 2),
                    (3, 5 - s / 2),
                    (3, 5 + s / 2),
                    (10, 5 + s / 2),
                    (10, 10),
                    (0, 10),
                ]
            )
        ],
        # A square hole of side `s` inside one building.
        "small_hole": [
            Polygon(
                [(0, 0), (10, 0), (10, 10), (0, 10)],
                [[(5, 5), (5 + s, 5), (5 + s, 5 + s), (5, 5 + s)]],
            )
        ],
        # A wedge whose tip subtends atan(s / 20): a sharp corner and a short
        # edge of length `s` at the far end.
        "sharp_tip": [Polygon([(0, 0), (20, 0), (20, s)])],
    }


def coverage_regions(buildings):
    """Buildings plus surrounding ground, in the production region convention."""
    occupied = unary_union(buildings)
    parts = polygon_parts(occupied)
    minx, miny, maxx, maxy = occupied.bounds
    ground = box(
        minx - GROUND_PADDING,
        miny - GROUND_PADDING,
        maxx + GROUND_PADDING,
        maxy + GROUND_PADDING,
    ).difference(occupied)
    ground_parts = polygon_parts(ground)
    return (
        [*ground_parts, *parts],
        [GROUND_MARKER] * len(ground_parts) + list(range(len(parts))),
        parts,
    )


def triangle_angles(vertices, faces):
    points = np.asarray(vertices, dtype=float)[:, :2]
    corners = points[np.asarray(faces, dtype=int)]
    angles = []
    for i in range(3):
        a = corners[:, i]
        u = corners[:, (i + 1) % 3] - a
        v = corners[:, (i + 2) % 3] - a
        norms = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            cosine = np.clip((u * v).sum(axis=1) / norms, -1.0, 1.0)
        angles.append(np.degrees(np.arccos(cosine)))
    return np.concatenate(angles)


def _mesh_once(regions, markers, spacing):
    """The meshing call itself, to be run in a child process under a wall clock."""
    from dtcc_core.model.mixins.mesh.quality import tri_element_quality

    started = time.perf_counter()
    try:
        mesh = build_city_flat_mesh_with_dtcc_mesher(
            region_polygons=regions,
            region_markers=markers,
            max_mesh_size=spacing,
            min_mesh_angle=MIN_MESH_ANGLE,
        )
    except Exception as error:  # the mesher raises plain RuntimeError
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
    angles = triangle_angles(mesh.vertices, mesh.faces)
    return {
        "outcome": "meshed",
        "seconds": seconds,
        "faces": int(len(faces)),
        "vertices": int(len(mesh.vertices)),
        "element_quality_min": float(np.min(quality)),
        "element_quality_p01": float(np.percentile(quality, 1)),
        "element_quality_below_threshold": int(
            np.count_nonzero(quality < POOR_ELEMENT_QUALITY)
        ),
        "degenerate_faces": int(np.count_nonzero(quality <= 0)),
        "smallest_angle_degrees": float(np.min(angles)),
        "angles_below_requested_minimum": int(
            np.count_nonzero(angles < MIN_MESH_ANGLE)
        ),
    }


def _worker(queue, regions, markers, spacing):
    queue.put(_mesh_once(regions, markers, spacing))


def mesh_attempt(buildings, spacing):
    """Mesh one configuration through the production flat-coverage entry point."""
    import multiprocessing

    regions, markers, parts = coverage_regions(buildings)
    state = admissibility(parts, DELTA)
    record = {
        "requested_spacing": spacing,
        "contract_resolved": state["resolved"],
        "contract_subscale_pairs": state.get("subscale_pairs"),
        "contract_separation_capped_at_delta": state.get("separation_capped_at_delta"),
    }
    area = unary_union(regions).area
    requested_elements = 4 * area / spacing**2
    if requested_elements > MAX_REQUESTED_ELEMENTS:
        record.update(
            outcome="not_attempted",
            reason="requested spacing exceeds the element budget",
            estimated_elements=requested_elements,
        )
        return record
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    process = context.Process(target=_worker, args=(queue, regions, markers, spacing))
    started = time.perf_counter()
    process.start()
    try:
        record.update(queue.get(timeout=ATTEMPT_SECONDS))
    except Exception:
        # A child that is gone without a result crashed; one still running at
        # the limit did not terminate. Both are failures, and they are not the
        # same failure, so they are not given the same name.
        alive = process.is_alive()
        record.update(
            outcome="did_not_terminate" if alive else "crashed",
            attempt_limit_seconds=ATTEMPT_SECONDS,
            child_alive_at_limit=alive,
            child_exit_code=process.exitcode,
            seconds=time.perf_counter() - started,
        )
    finally:
        process.terminate()
        process.join()
    return record


def provenance():
    """Record what produced a number, in the form the earlier manifests use."""
    import hashlib
    import platform

    import shapely

    return {
        "scripts": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in (
                "sandbox/cleaning_mesher_precondition_probe.py",
                "dtcc_core/builder/cleaning/contract.py",
                "dtcc_core/builder/meshing/dtcc_mesher_backend.py",
            )
        },
        "versions": {
            "shapely": shapely.__version__,
            "geos": ".".join(str(n) for n in shapely.geos_version),
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "measurement_environment": platform.platform(),
    }


def mesher_revision():
    from importlib.metadata import distribution

    import dtcc_mesher

    try:
        direct = json.loads(distribution("dtcc-mesher").read_text("direct_url.json"))
    except Exception:  # pragma: no cover - only when installed from a wheel
        direct = None
    return {
        "module": dtcc_mesher.__file__,
        "version": getattr(dtcc_mesher, "__version__", None),
        "direct_url": direct,
    }


def requested_spacings(separation):
    """The fixed spacing ladder for one separation, coarsest first."""
    finest = max(FINEST_REQUESTED_SPACING, min(4 * separation, MESH_SIZES[0]))
    return sorted({*MESH_SIZES, finest}, reverse=True)


def survey(destination, separations, budget=None):
    """Append one record per attempt; never redo a finished attempt.

    The measurement harness cannot hold a process for the whole sweep and a
    single attempt can take a minute and a half, so the unit that lands on disk
    is the attempt, not the configuration.
    """
    started = time.perf_counter()
    done = {}
    if destination.exists():
        for line in destination.read_text().splitlines():
            record = json.loads(line)
            done[
                (record["family"], record["separation"], record["requested_spacing"])
            ] = record
    complete = True
    for separation in separations:
        for name, buildings in configurations(separation).items():
            for spacing in requested_spacings(separation):
                key = (name, separation, spacing)
                if key in done:
                    continue
                if budget is not None and time.perf_counter() - started > budget:
                    print(f"budget reached before {name} s={separation} h={spacing}")
                    return done, False
                record = mesh_attempt(buildings, spacing)
                done[key] = {"family": name, "separation": separation, **record}
                with destination.open("a") as handle:
                    handle.write(json.dumps(done[key]) + "\n")
                print(
                    f"{name:14s} s={separation:<8g} h={spacing:<6g} "
                    f"{record['outcome']:18s} {record.get('seconds', 0):7.2f}s "
                    f"q={record.get('element_quality_min', float('nan')):.4f}",
                    flush=True,
                )
    return done, complete


def assemble(done, separations):
    """Group finished attempts into one verdict per family and separation."""
    rows = []
    for separation in separations:
        for name in configurations(separation):
            attempts = [
                done[(name, separation, spacing)]
                for spacing in requested_spacings(separation)
                if (name, separation, spacing) in done
            ]
            if len(attempts) != len(requested_spacings(separation)):
                continue
            tried = [a for a in attempts if a["outcome"] != "not_attempted"]
            failed = [a for a in tried if a["outcome"] != "meshed"]
            meshed = [a for a in tried if a["outcome"] == "meshed"]
            finest = min(a["requested_spacing"] for a in tried)
            rows.append(
                {
                    "family": name,
                    "separation": separation,
                    "verdict": (
                        "meshed at every requested spacing"
                        if not failed
                        else (
                            "refinement removed the failure"
                            if meshed
                            else "failed at every requested spacing"
                        )
                    ),
                    "coarsest_failure_spacing": (
                        max(a["requested_spacing"] for a in failed) if failed else None
                    ),
                    "worst_element_quality": min(
                        (a["element_quality_min"] for a in meshed), default=None
                    ),
                    "smallest_angle_degrees": min(
                        (a["smallest_angle_degrees"] for a in meshed), default=None
                    ),
                    "finest_requested_spacing": finest,
                    "grading_conclusion_established": finest <= separation
                    or not failed,
                    "attempts": attempts,
                }
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--budget", type=float, default=None)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("output directory is not empty; use a fresh directory")
    args.output.mkdir(parents=True, exist_ok=True)
    done, complete = survey(
        args.output / "attempts.jsonl", SEPARATIONS, budget=args.budget
    )
    rows = assemble(done, SEPARATIONS)
    print(f"{len(done)} attempts, {len(rows)} rows; complete={complete}", flush=True)
    if not complete:
        return
    (args.output / "preconditions.json").write_text(
        json.dumps(
            {
                "scope": (
                    "Milestone D: where core's pinned mesher actually fails on "
                    "sub-delta features, measured through the production "
                    "flat-coverage entry point. Not a cleaning result."
                ),
                "delta": DELTA,
                "min_mesh_angle": MIN_MESH_ANGLE,
                "ground_marker": GROUND_MARKER,
                "ground_padding": GROUND_PADDING,
                "poor_element_quality": POOR_ELEMENT_QUALITY,
                "max_requested_elements": MAX_REQUESTED_ELEMENTS,
                "attempt_limit_seconds": ATTEMPT_SECONDS,
                "separations": list(SEPARATIONS),
                "requested_spacings": list(MESH_SIZES),
                "finest_requested_spacing": FINEST_REQUESTED_SPACING,
                "mesher": mesher_revision(),
                **provenance(),
                "rows": rows,
            },
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
