"""Aggregate the frozen legacy/staged/public survey rows into promotion evidence."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np


LEGACY_RUN = {
    "gothenburg": "legacy2",
    "helsingborg": "legacy2",
    "linkoping": "legacy2",
    "stockholm": "legacy2",
    "uppsala": "legacy2",
    "lund": "legacy3",
    "malmo": "legacy3",
    "norrkoping": "legacy3",
    "orebro": "legacy3",
    "vasteras": "legacy3",
}


def distribution(values) -> dict:
    values = np.asarray(list(values), dtype=float)
    if not len(values):
        return {"count": 0}
    return {
        "count": int(len(values)),
        "sum": float(values.sum()),
        "min": float(values.min()),
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "max": float(values.max()),
    }


def passed(contract) -> bool:
    return isinstance(contract, dict) and contract.get("status") == "pass"


def turnover(ids, old_status, new_status) -> dict:
    result = {name: [] for name in ("both_pass", "gain", "loss", "both_fail")}
    for task_id in ids:
        old = old_status(task_id)
        new = new_status(task_id)
        key = (
            "both_pass" if old and new else "gain" if new else "loss" if old else "both_fail"
        )
        result[key].append(task_id)
    return {key: {"count": len(value), "task_ids": value} for key, value in result.items()}


def surface_summary(path: Path, flat: dict, ids: list[str]) -> dict:
    rows = {}
    for city in LEGACY_RUN:
        for line in (path / f"{city}.jsonl").read_text().splitlines():
            row = json.loads(line)
            rows[row["task_id"]] = row
    if len(rows) != 1000 or set(rows) != set(ids):
        raise ValueError(f"expected 1000 aligned surface rows, got {len(rows)}")
    meshed = [
        task_id
        for task_id in ids
        if rows[task_id].get("mesh", {}).get("outcome") == "meshed"
    ]
    failures = [task_id for task_id in ids if rows[task_id]["outcome"] == "failed"]
    return {
        "cases": len(rows),
        "raw_identity_verified": sum(
            row.get("raw_identity_verified") is True for row in rows.values()
        ),
        "cleaning_outcome_matches_flat": sum(
            (
                "unresolved"
                if rows[task_id].get("failure_phase") == "cleaning"
                else "accepted"
            )
            == flat[task_id]["outcome"]
            for task_id in ids
        ),
        "outcomes": dict(Counter(row["outcome"] for row in rows.values())),
        "meshed": len(meshed),
        "unattempted_after_cleaning_rejection": sum(
            row["outcome"] == "unresolved" for row in rows.values()
        ),
        "failures_by_phase": dict(
            Counter(rows[task_id].get("failure_phase", "unknown") for task_id in failures)
        ),
        "failure_task_ids": failures,
        "faces": sum(rows[task_id]["mesh"]["num_faces"] for task_id in meshed),
        "quality_min": distribution(
            rows[task_id]["mesh"]["quality"]["element_quality"]["min"]
            for task_id in meshed
        ),
        "quality_p01": distribution(
            rows[task_id]["mesh"]["element_quality_p01"] for task_id in meshed
        ),
        "q_below_0_02_cells": sum(
            rows[task_id]["mesh"]["element_quality_below_0_02_count"]
            for task_id in meshed
        ),
        "degenerate_cells": sum(
            rows[task_id]["mesh"]["degenerate_cell_count"] for task_id in meshed
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--production", required=True, type=Path)
    parser.add_argument("--surface", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    new = {}
    staged = {}
    legacy = {}
    for city, legacy_name in LEGACY_RUN.items():
        for line in (args.production / f"{city}.jsonl").read_text().splitlines():
            row = json.loads(line)
            new[row["task_id"]] = row
        for line in (args.audit / "staged" / f"{city}.jsonl").read_text().splitlines():
            row = json.loads(line)
            staged[row["task_id"]] = row
        rows = json.loads(
            (args.audit / legacy_name / city / "results.json").read_text()
        )["results"]
        legacy.update({row["task_id"]: row for row in rows})
    ids = sorted(new)
    if not (len(ids) == len(staged) == len(legacy) == 1000):
        raise ValueError(
            f"expected 1000 aligned rows, got {len(new)}/{len(staged)}/{len(legacy)}"
        )
    if set(ids) != set(staged) or set(ids) != set(legacy):
        raise ValueError("task ID sets differ")

    def legacy_pre(task_id):
        cleaning = legacy[task_id]["metrics"]["cleaning"]
        return not cleaning.get("input_count", 0) or passed(
            cleaning.get("before_selection_contract")
        )

    def staged_pre(task_id):
        return passed(staged[task_id].get("contract"))

    def new_pre(task_id):
        return passed(new[task_id].get("before_selection_contract"))

    accepted = [task_id for task_id in ids if new[task_id]["outcome"] == "accepted"]
    meshed = [
        task_id
        for task_id in accepted
        if new[task_id].get("mesh", {}).get("outcome") == "meshed"
    ]
    nonempty = [task_id for task_id in ids if new[task_id]["original_input_count"]]
    exclusions = [
        exclusion
        for task_id in accepted
        for exclusion in new[task_id].get("selection", {}).get("excluded_regions", [])
    ]
    policy_exclusions = [
        exclusion
        for task_id in accepted
        for exclusion in new[task_id].get("policy_exclusions", [])
    ]
    failure_reasons = Counter()
    for task_id in ids:
        row = new[task_id]
        if row["outcome"] == "accepted":
            continue
        diagnostics = row.get("diagnostics", {})
        groups = diagnostics.get("unresolved_groups", [])
        if groups:
            failure_reasons.update(group.get("reason", "unknown") for group in groups)
        else:
            failure_reasons[row.get("failure_phase", "unknown")] += 1

    quality_min = [
        new[task_id]["mesh"]["quality"]["element_quality"]["min"]
        for task_id in meshed
    ]
    quality_p01 = [new[task_id]["mesh"]["element_quality_p01"] for task_id in meshed]
    report = {
        "scope": "Exact saved-input production cleaning and flat-mesh replay",
        "cases": len(ids),
        "nonempty": len(nonempty),
        "empty": len(ids) - len(nonempty),
        "raw_identity_verified": sum(
            row.get("raw_identity_verified") is True for row in new.values()
        ),
        "outcomes": dict(Counter(row["outcome"] for row in new.values())),
        "accepted": len(accepted),
        "flat_meshed": len(meshed),
        "failure_reasons": dict(failure_reasons),
        "pre_selection_contract_pass": {
            "legacy": sum(legacy_pre(task_id) for task_id in ids),
            "staged": sum(staged_pre(task_id) for task_id in ids),
            "production": sum(new_pre(task_id) for task_id in ids),
        },
        "production_handoff_profile_pass": {
            "before_selection": sum(
                passed(
                    new[task_id]
                    .get("before_selection_contract", {})
                    .get("mesher_profile")
                )
                for task_id in accepted
            ),
            "final_subdivision": sum(
                passed(
                    new[task_id]
                    .get("final_handoff_contract", {})
                    .get("mesher_profile")
                )
                for task_id in accepted
            ),
        },
        "legacy_to_production_turnover": turnover(ids, legacy_pre, new_pre),
        "staged_to_production_turnover": turnover(ids, staged_pre, new_pre),
        "selection": {
            "excluded_regions": len(exclusions),
            "excluded_area": sum(row["area"] for row in exclusions),
            "policy_exclusions": len(policy_exclusions),
            "policy_excluded_area": sum(row["area"] for row in policy_exclusions),
            "raw_reference_final_fidelity_pass": sum(
                passed(new[task_id].get("final_handoff_contract"))
                for task_id in accepted
            ),
        },
        "timing": {
            "cleaning_seconds": distribution(
                new[task_id]["cleaning_seconds"] for task_id in accepted
            ),
            "flat_mesh_seconds": distribution(
                new[task_id]["mesh"]["seconds"] for task_id in meshed
            ),
        },
        "drift": distribution(
            new[task_id]["cleaning"]["removed_area"]
            + new[task_id]["cleaning"]["added_area"]
            for task_id in accepted
        ),
        "flat_mesh": {
            "faces": sum(new[task_id]["mesh"]["num_faces"] for task_id in meshed),
            "quality_min": distribution(quality_min),
            "quality_p01": distribution(quality_p01),
            "q_below_0_02_cells": sum(
                new[task_id]["mesh"]["element_quality_below_0_02_count"]
                for task_id in meshed
            ),
            "degenerate_cells": sum(
                new[task_id]["mesh"]["degenerate_cell_count"]
                for task_id in meshed
            ),
            "qmin_below_0_02_task_ids": [
                task_id
                for task_id in meshed
                if new[task_id]["mesh"]["quality"]["element_quality"]["min"] < 0.02
            ],
        },
    }
    if args.surface is not None:
        report["surface_mesh"] = surface_summary(args.surface, new, ids)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
