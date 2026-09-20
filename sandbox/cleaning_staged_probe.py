"""Research CLI for the shared production footprint constructor.

Construction and acceptance live in ``dtcc_core.builder.cleaning.construction``;
this module remains as the documented survey entry point and compatibility
import for research scripts.
"""

from dtcc_core.builder.cleaning.construction import *  # noqa: F401,F403


def main():
    """Retain the historical group-survey command as thin orchestration."""
    import argparse
    import json
    import time
    from pathlib import Path

    import numpy as np

    from sandbox.cleaning_pipeline_probe import saved_city_groups, topology_examples

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--delta", type=float, default=0.5)
    parser.add_argument("--epsilon", type=float, default=None)
    parser.add_argument("--cities", type=str, default="")
    parser.add_argument("--budget", type=float, default=None)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error(
            "output directory is not empty; use a fresh directory (resume is unsupported)"
        )
    if not np.isfinite(args.delta) or args.delta <= 0:
        parser.error("delta must be finite and positive")
    if args.epsilon is not None and (
        not np.isfinite(args.epsilon) or args.epsilon < 0
    ):
        parser.error("epsilon must be finite and nonnegative")
    args.output.mkdir(parents=True, exist_ok=True)
    epsilon = args.delta / 2 if args.epsilon is None else args.epsilon

    for name, (raw, _) in topology_examples().items():
        _, report = construct(raw, delta=args.delta, epsilon=epsilon)
        print(
            f"{name:20s} {report['outcome']:12s} "
            f"init={report['initial']} final={report.get('final')}",
            flush=True,
        )
    if args.run is None:
        return

    selected_cities = {city for city in args.cities.split(",") if city}
    destination = args.output / f"groups-delta{args.delta:g}.jsonl"
    started = time.perf_counter()
    for record, _, _, groups in saved_city_groups(args.run):
        city = record["case"]["city"]
        if selected_cities and city not in selected_cities:
            continue
        for group_index, raw in enumerate(groups):
            if args.budget is not None and time.perf_counter() - started > args.budget:
                return
            _, report = construct(raw, delta=args.delta, epsilon=epsilon)
            row = {
                "city": city,
                "group": group_index,
                "polygon_count": len(raw),
                "delta": args.delta,
                "epsilon": epsilon,
                "outcome": report["outcome"],
                "reason": report["reason"],
                "initial": report["initial"],
                "final": report.get("final"),
                "edits": report["edits"],
                "evaluations": report["evaluations"],
                "global_evaluations": report["global_evaluations"],
                "seconds": report["seconds"],
            }
            with destination.open("a") as handle:
                handle.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
