"""Replay the frozen full-survey focus fixtures through a cleaner.

This is a thin research client of the staged constructor.  It does not contain
repair or acceptance logic; its purpose is to preserve comparable case-level
evidence while that implementation moves into the cleaning package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from shapely import from_wkb

from dtcc_core.builder.cleaning import construction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURES = ROOT / "tests/data/cleaning/full-survey-focused-groups.json"


def _source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(fixtures: Path) -> dict:
    payload = json.loads(fixtures.read_text())
    rows = []
    started = time.perf_counter()
    for case in payload["cases"]:
        for group in case["groups"]:
            raw = [from_wkb(bytes.fromhex(value)) for value in group["polygon_wkb_hex"]]
            group_started = time.perf_counter()
            output, report = construction.construct(
                raw, delta=payload["delta"], epsilon=payload["epsilon"]
            )
            elapsed = time.perf_counter() - group_started
            rows.append(
                {
                    "case_id": case["case_id"],
                    "group": group["group"],
                    "polygon_count": len(raw),
                    "outcome": report["outcome"],
                    "reason": report["reason"],
                    "initial": report["initial"],
                    "final": report.get("final"),
                    "sampling_attempts": report.get("sampling_attempts"),
                    "edits": report["edits"],
                    "sites": report["sites"],
                    "evaluations": report["evaluations"],
                    "admissions": report["admissions"],
                    "global_evaluations": report["global_evaluations"],
                    "total_work": report.get("total_work", {}),
                    "contract": report["contract"],
                    "mesher_profile": report.get("mesher_profile"),
                    "fallback_attempts": report.get("fallback_attempts", []),
                    "seconds": elapsed,
                }
            )
            print(
                f"{case['case_id']} group {group['group']}: "
                f"{report['outcome']} {report.get('final')} {elapsed:.3f}s",
                flush=True,
            )
    return {
        "scope": "Frozen focus-fixture replay through production construction",
        "fixtures": str(fixtures.relative_to(ROOT)),
        "fixture_sha256": _source_hash(fixtures),
        "constructor": "dtcc_core/builder/cleaning/construction.py",
        "constructor_sha256": _source_hash(
            ROOT / "dtcc_core/builder/cleaning/construction.py"
        ),
        "contract_sha256": _source_hash(
            ROOT / "dtcc_core/builder/cleaning/contract.py"
        ),
        "delta": payload["delta"],
        "epsilon": payload["epsilon"],
        "seconds": time.perf_counter() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; use a fresh path")
    result = replay(args.fixtures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
