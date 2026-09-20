"""Observational probe for docs/design/footprint-cleaning-contract.md.

No acceptance gates and no new cleaning implementation. Run from dtcc-core:
  .venv/bin/python sandbox/cleaning_contract_probe.py --output /tmp/contract-probe
Optionally pass --run <saved-benchmark-run> to inspect its flat/cleaning tasks.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from shapely.geometry import Polygon, box

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from benchmarks.benchmark_catalog import CLEANING_PARAMETERS, DEFAULT_PARAMETERS
from benchmarks.benchmark_phases import footprint_metrics, load_cleaning, write_json
from dtcc_core.builder.geometry_builders.meshes import build_conditioned_footprints
from dtcc_core.builder.cleaning.contract import admissibility
from dtcc_core.model import Building, GeometryType, Surface


def examples():
    return {
        "redundant_wall_vertex": [
            Polygon([(0, 0), (0.1, 0), (10, 0), (10, 6), (0, 6)])
        ],
        "near_buildings": [box(0, 0, 10, 6), box(10.2, 0, 20.2, 6)],
        "courtyard_passage": [
            box(0, 0, 12, 10).difference(box(4, 3, 8, 7).union(box(5.9, 7, 6.1, 10)))
        ],
        "tiny_hole": [
            Polygon(
                box(0, 0, 10, 10).exterior.coords,
                [box(4.9, 4.9, 5.1, 5.1).exterior.coords],
            )
        ],
        "sharp_corner": [Polygon([(0, 0), (30, -2), (30, 2)])],
        "small_resolved_building": [box(0, 0, 3, 3)],
    }


def measure(name, raw, cleaned, source_map, delta, epsilon, existing_contract):
    metrics = footprint_metrics(raw, cleaned, source_map, delta, epsilon=epsilon)
    contract = metrics["geometric_contract"]
    return {
        "name": name,
        "delta": delta,
        "epsilon": epsilon,
        "input": admissibility(raw, delta),
        "output": contract["admissibility"],
        "fidelity": contract["fidelity"],
        "existing_contract": existing_contract,
        "metrics": metrics,
    }


def plot_examples(cases, rows, destination):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from dtcc_core.plotting.style import plot_polygon_geometries, set_axes_extent

    fig, axes = plt.subplots(
        len(cases), 2, figsize=(12, 3 * len(cases)), layout="constrained"
    )
    for (name, (raw, cleaned)), row, pair in zip(cases.items(), rows, axes):
        for ax, polygons in zip(pair, (raw, cleaned)):
            plot_polygon_geometries(
                ax, polygons, palette=["#66a8bd"], edgecolor="#234653"
            )
            for polygon in polygons:
                for ring in [polygon.exterior, *polygon.interiors]:
                    xy = np.asarray(ring.coords)
                    ax.scatter(xy[:, 0], xy[:, 1], s=8, color="#234653", zorder=3)
        set_axes_extent(pair, raw + cleaned)
        pair[0].set_title(name.replace("_", " ") + " — input")
        pair[1].set_title(
            f"Current cleaner — resolved: {row['output']['resolved']}; "
            f"fidelity: {row['fidelity']['status']}"
        )
        if row["input"]["resolved"]:
            pair[0].text(0.01, 0.03, "Already admissible", transform=pair[0].transAxes)
    fig.suptitle(
        "Cleaning contract examples · δ = 0.5 m, ε = 0.25 m\nArea selection disabled"
    )
    fig.savefig(destination, dpi=150)
    plt.close(fig)


def saved_cases(run):
    payload = json.loads((run / "results.json").read_text())
    seen = set()
    for result in payload["results"]:
        if result["dataset"] not in {"city_footprints", "city_flat_mesh"}:
            continue
        key = result["case"]["id"], result["scenario"]["id"]
        artifact = result.get("artifacts", {}).get("cleaning_input")
        if key in seen or not artifact:
            continue
        seen.add(key)
        path = (run / artifact["path"]).resolve()
        if not path.is_relative_to(run.resolve()):
            raise ValueError("Saved cleaning artifact must be inside its run")
        city, conditioned = load_cleaning(path, result)
        raw = [
            b.flatten_geometry(GeometryType.LOD0).to_polygon(simplify=0.0)
            for b in city.buildings
        ]
        cleaned = [s.to_polygon(simplify=0.0) for s in conditioned.surfaces]
        delta = result["parameters"]["min_building_detail"]
        row = measure(
            "/".join(key),
            raw,
            cleaned,
            conditioned.source_map,
            delta,
            delta / 2,
            result["metrics"]["cleaning"]
            .get("contract", {})
            .get("status", "unavailable"),
        )
        row["existing_clearance_tolerance"] = (
            result["metrics"]["cleaning"]
            .get("contract", {})
            .get("metrics", {})
            .get("contract_tolerance")
        )
        yield row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows, plotted = [], {}
    delta, epsilon = 0.5, 0.25
    parameters = {k: DEFAULT_PARAMETERS[k] for k in CLEANING_PARAMETERS}
    parameters["min_building_area"] = 0.0
    for name, raw in examples().items():
        buildings = []
        for polygon in raw:
            surface = Surface()
            surface.from_polygon(polygon)
            building = Building()
            building.add_geometry(surface, GeometryType.LOD0)
            buildings.append(building)
        conditioned = build_conditioned_footprints(
            buildings,
            lod=GeometryType.LOD0,
            **parameters,
            max_mesh_size=10.0,
            cleaning_diagnostics=False,
            raise_on_contract_error=False,
        )
        cleaned = [s.to_polygon(simplify=0.0) for s in conditioned.surfaces]
        row = measure(
            name,
            raw,
            cleaned,
            conditioned.source_map,
            delta,
            epsilon,
            conditioned.contract["status"],
        )
        rows.append(row)
        plotted[name] = raw, cleaned
        print(
            name,
            row["input"]["resolved"],
            row["output"]["resolved"],
            row["fidelity"],
            flush=True,
        )
    plot_examples(plotted, rows, args.output / "examples.png")
    cities = list(saved_cases(args.run)) if args.run else []
    write_json(
        args.output / "observations.json",
        {
            "core_revision": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
            "saved_run": str(args.run.resolve()) if args.run else None,
            "examples": rows,
            "cities": cities,
        },
    )
    for row in cities:
        print(row["name"], row["output"], row["fidelity"], flush=True)


if __name__ == "__main__":
    main()
