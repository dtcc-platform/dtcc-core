"""Replay frozen survey inputs through the shared public cleaning/flat workflow.

This is an evidence client, not an alternate adapter: construction, selection,
source attribution, flat-city preparation and mesh measurements all come from
the production builder and benchmark modules.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.io.model import load_model
from dtcc_core.model import Bounds, City, GeometryType

from benchmarks.benchmark_catalog import CLEANING_PARAMETERS
from benchmarks.benchmark_datasets import _mesh_metrics, json_ready
from benchmarks.benchmark_phases import (
    PREPARATION_PARAMETERS,
    _flat_city,
    footprint_metrics,
    prepare_city_from_bounds,
)


def _raw_polygons(city: City):
    return [
        polygon
        for building in city.buildings
        if (geometry := building.flatten_geometry(GeometryType.LOD0)) is not None
        if (polygon := geometry.to_polygon(simplify=0.0)) is not None
        and not polygon.is_empty
    ]


def _raw_hash(polygons) -> str:
    return hashlib.sha256(b"".join(polygon.wkb for polygon in polygons)).hexdigest()


def _write_rows(path: Path, rows: list[dict]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        "".join(json.dumps(json_ready(row)) + "\n" for row in rows)
    )
    temporary.replace(path)


def replay(
    run: Path,
    staged: Path,
    output: Path,
    *,
    cleaning_only: bool,
    surface: bool,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = (
        [json.loads(line) for line in output.read_text().splitlines()]
        if output.exists()
        else []
    )
    completed = {row["task_id"] for row in rows}
    expected = {
        row["task_id"]: row
        for line in staged.read_text().splitlines()
        if (row := json.loads(line))
    }
    records = [
        row
        for row in json.loads((run / "results.json").read_text())["results"]
        if row["dataset"] == "city_flat_mesh"
    ]
    for record in sorted(records, key=lambda row: row["task_id"]):
        task_id = record["task_id"]
        if task_id in completed:
            continue
        row = {"task_id": task_id, "case": record["case"]}
        phase = "load"
        started = time.perf_counter()
        try:
            artifact = (run / record["artifacts"]["cleaning_input"]["path"]).resolve()
            if not artifact.is_relative_to(run.resolve()):
                raise ValueError("saved cleaning artifact escaped its run")
            city = load_model(artifact.parent / "raw.dtcc", expected_type=City)
            raw = _raw_polygons(city)
            raw_hash = _raw_hash(raw)
            expected_hash = expected[task_id]["raw_geometry_hash"]
            if raw_hash != expected_hash:
                raise ValueError(
                    f"raw geometry hash mismatch: {raw_hash} != {expected_hash}"
                )
            row.update(
                raw_geometry_hash=raw_hash,
                raw_identity_verified=True,
                original_input_count=len(raw),
            )

            phase = "cleaning"
            phase_started = time.perf_counter()
            parameters = record["parameters"]
            conditioned = meshes.build_conditioned_footprints(
                city.buildings,
                lod=GeometryType.LOD0,
                **{key: parameters[key] for key in CLEANING_PARAMETERS},
                max_mesh_size=parameters["max_mesh_size"],
                cleaning_diagnostics=False,
            )
            cleaned = [
                surface.to_polygon(simplify=0.0)
                for surface in conditioned.surfaces
            ]
            row.update(
                outcome="accepted",
                cleaning_seconds=time.perf_counter() - phase_started,
                cleaning=footprint_metrics(
                    raw,
                    cleaned,
                    conditioned.source_map,
                    conditioned.declared_scale,
                    epsilon=conditioned.declared_scale / 2,
                ),
                before_selection_contract=conditioned.diagnostics.get(
                    "before_selection_contract"
                ),
                selection=conditioned.diagnostics.get("selection"),
                final_handoff_contract=conditioned.diagnostics.get(
                    "final_handoff_contract"
                ),
                policy_exclusions=conditioned.diagnostics.get(
                    "policy_exclusions", []
                ),
                source_map=conditioned.source_map,
            )

            if not cleaning_only:
                phase = "surface_preparation" if surface else "flat_meshing"
                phase_started = time.perf_counter()
                bounds = record["case"]["bounds"]
                if surface:
                    city = prepare_city_from_bounds(
                        Bounds(*bounds),
                        buildings=city.buildings,
                        **{
                            key: parameters[key]
                            for key in PREPARATION_PARAMETERS
                        },
                    )
                    row["surface_preparation_seconds"] = (
                        time.perf_counter() - phase_started
                    )
                    phase = "surface_meshing"
                    phase_started = time.perf_counter()
                    build = meshes.build_city_surface_mesh
                else:
                    city = _flat_city(city, Bounds(*bounds))
                    build = meshes.build_city_flat_mesh
                mesh = build(
                    city,
                    lod=GeometryType.LOD0,
                    conditioned_footprints=conditioned,
                    report_mesh_quality=False,
                    cleaning_diagnostics=False,
                    stage_audit={},
                    **{
                        key: parameters[key]
                        for key in (
                            "max_mesh_size",
                            "min_mesh_angle",
                            *CLEANING_PARAMETERS,
                        )
                    },
                )
                row["mesh"] = {
                    "outcome": "meshed",
                    "dataset": "surface" if surface else "flat",
                    "seconds": time.perf_counter() - phase_started,
                    **_mesh_metrics(mesh),
                }
        except Exception as error:
            row.update(
                outcome="unresolved" if phase == "cleaning" else "failed",
                failure_phase=phase,
                error_type=type(error).__name__,
                error=str(error),
            )
            diagnostics = getattr(error, "diagnostics", None)
            if diagnostics is not None:
                row["diagnostics"] = diagnostics
        row["wall_seconds"] = time.perf_counter() - started
        rows.append(row)
        _write_rows(output, rows)
        print(
            f"{task_id}: {row['outcome']} {row['wall_seconds']:.3f}s",
            flush=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--staged", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cleaning-only", action="store_true")
    parser.add_argument("--surface", action="store_true")
    args = parser.parse_args()
    replay(
        args.run,
        args.staged,
        args.output,
        cleaning_only=args.cleaning_only,
        surface=args.surface,
    )


if __name__ == "__main__":
    main()
