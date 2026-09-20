"""Run the predeclared east-adjacent validation tiles via the benchmark path."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.benchmark_catalog import CITIES, DEFAULT_PARAMETERS
from benchmarks.benchmark_datasets import json_ready, run_dataset


VALIDATION_CITIES = ("lund", "gothenburg", "stockholm")


def validation_case(city_name: str) -> dict:
    city = CITIES[city_name]
    xmin = city.grid_xmin + city.grid_nx * city.grid_box_size
    ymin = city.grid_ymin + city.grid_center_iy * city.grid_box_size
    return {
        "id": f"fresh_east:{city_name}:central_row",
        "kind": "fresh_east_adjacent",
        "city": city_name,
        "label": f"{city.label} east-adjacent central-row validation tile",
        "bounds": [xmin, ymin, xmin + city.grid_box_size, ymin + city.grid_box_size],
        "tags": ["fresh_validation", city_name],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; use a fresh directory")
    args.output.mkdir(parents=True)
    results = []
    for city_name in VALIDATION_CITIES:
        case = validation_case(city_name)
        task = {
            "id": f"fresh:city_flat_mesh:{case['id']}:baseline",
            "suite": "fresh_validation",
            "dataset": "city_flat_mesh",
            "case": case,
            "scenario": {"id": "baseline", "parameters": {}, "tags": []},
            "parameters": dict(DEFAULT_PARAMETERS),
            "timeout_seconds": 420,
            "phase": "both",
            "task_dir": str(args.output / city_name),
        }
        result = run_dataset(task)
        results.append(result)
        print(city_name, result["status"], flush=True)
    (args.output / "results.json").write_text(
        json.dumps(json_ready(results), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
