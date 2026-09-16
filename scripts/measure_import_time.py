#!/usr/bin/env python3
"""
Measure the wall-clock cost of importing dtcc_core and report the most
expensive modules pulled in along the way.

Usage:
  python scripts/measure_import_time.py
  python scripts/measure_import_time.py --module dtcc_core.io --runs 7
  python scripts/measure_import_time.py --top 25

The reported time is the minimum across runs, which is the most stable
estimate available: it discards scheduler noise and cold-cache outliers.
Run it twice in a row if the .pyc cache may be cold.

Exit codes:
  0: Measurement completed
  2: Invalid input or unexpected error
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from collections import defaultdict

IMPORTTIME_LINE = re.compile(r"import time:\s+(\d+) \|\s+(\d+) \|(\s*)(\S+)")


def measure_wall_clock(statement: str, runs: int) -> tuple[float, float]:
    """Return (minimum, median) wall-clock seconds for running `statement`."""
    timings = []
    for _ in range(runs):
        start = time.perf_counter()
        result = subprocess.run(
            [sys.executable, "-c", statement],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        timings.append(time.perf_counter() - start)
        if result.returncode != 0:
            raise SystemExit(f"error: `{statement}` failed with {result.returncode}")
    timings.sort()
    return timings[0], timings[len(timings) // 2]


def collect_importtime(statement: str) -> list[tuple[int, int, int, str]]:
    """Run `statement` under -X importtime and parse the trace."""
    result = subprocess.run(
        [sys.executable, "-X", "importtime", "-c", statement],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise SystemExit(f"error: `{statement}` failed with {result.returncode}")

    rows = []
    for line in result.stderr.splitlines():
        match = IMPORTTIME_LINE.match(line)
        if match:
            self_us, cumulative_us, indent, name = match.groups()
            rows.append((int(self_us), int(cumulative_us), len(indent) // 2, name))
    return rows


def report(rows, top: int) -> None:
    """Print the cost breakdown for a parsed importtime trace."""
    total_us = sum(row[0] for row in rows)
    print(f"modules imported: {len(rows)}")
    print(f"total self time:  {total_us / 1000:.1f} ms")

    by_root = defaultdict(int)
    for self_us, _, _, name in rows:
        by_root[name.split(".")[0]] += self_us

    print(f"\ncumulative cost by top-level package (top {top}):")
    ranked = sorted(by_root.items(), key=lambda item: -item[1])[:top]
    for name, self_us in ranked:
        share = 100 * self_us / total_us if total_us else 0
        print(f"  {self_us / 1000:8.1f} ms  {share:4.1f}%  {name}")

    print(f"\nmost expensive single imports, cumulative (top {top}):")
    for self_us, cumulative_us, _, name in sorted(rows, key=lambda r: -r[1])[:top]:
        print(f"  {cumulative_us / 1000:8.1f} ms cum  {self_us / 1000:7.1f} ms self  {name}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--module",
        default="dtcc_core",
        help="module to import (default: dtcc_core)",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=5,
        help="wall-clock runs to time (default: 5)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=15,
        help="entries to show per table (default: 15)",
    )
    args = parser.parse_args()

    if args.runs < 1:
        parser.error("--runs must be at least 1")
    if args.top < 1:
        parser.error("--top must be at least 1")

    statement = f"import {args.module}"

    baseline_min, _ = measure_wall_clock("pass", args.runs)
    module_min, module_median = measure_wall_clock(statement, args.runs)

    print(f"$ python -c '{statement}'")
    print(f"interpreter baseline: {baseline_min:.3f} s")
    print(f"minimum of {args.runs}:      {module_min:.3f} s")
    print(f"median of {args.runs}:       {module_median:.3f} s")
    print(f"attributable to import: {module_min - baseline_min:.3f} s\n")

    report(collect_importtime(statement), args.top)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as exc:  # pragma: no cover - defensive
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
