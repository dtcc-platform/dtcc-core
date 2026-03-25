"""Run and summarize a 100-tile Stockholm flat-mesh comparison.

This script compares the current Shapely-based cleaner against the legacy
pre-refactor implementation across the full 10x10 Stockholm grid. It:

1. ensures a legacy git worktree exists at the requested commit
2. runs the branch-compatible comparison harness in summary-only mode for all
   missing cases in both legacy and new modes
3. aggregates all per-case JSON summaries into one combined JSON file
4. writes summary tables
5. writes comparison plots (distribution panels, delta histograms, grid heatmaps)

The batch is resumable: if a per-case ``summary.json`` already exists for a
case/mode, that case is skipped unless ``--force-rerun`` is used.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import sysconfig
from datetime import datetime, UTC
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


REPO_ROOT = Path(__file__).resolve().parent.parent
COMPARE_SCRIPT = REPO_ROOT / "sandbox" / "compare_stockholm_flat_mesh.py"
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "sandbox" / "output" / "stockholm_flat_mesh_compare" / "stockholm_100_tile_batch"
)
DEFAULT_LEGACY_COMMIT = "0493c62"
DEFAULT_LEGACY_ROOT = REPO_ROOT.parent / f"dtcc-core-legacy-{DEFAULT_LEGACY_COMMIT}"
CURRENT_DTCC_BUILDER = REPO_ROOT / "dtcc_core" / "builder" / "_dtcc_builder.cpython-312-darwin.so"
LEGACY_DTCC_BUILDER_REL = Path("dtcc_core") / "builder" / "_dtcc_builder.cpython-312-darwin.so"

NX = 10
NY = 10
CASES = list(range(1, NX * NY + 1))

LEGACY_LABEL = "legacy-100"
NEW_LABEL = "new-100"


METRIC_SPECS = [
    {
        "key": "eq_worst",
        "label": "Element Quality Worst",
        "path": ("flat_mesh_quality", "element_quality", "min"),
        "better": "higher",
        "hist_log": False,
        "heatmap": True,
    },
    {
        "key": "eq_mean",
        "label": "Element Quality Mean",
        "path": ("flat_mesh_quality", "element_quality", "mean"),
        "better": "higher",
        "hist_log": False,
        "heatmap": False,
    },
    {
        "key": "ar_worst",
        "label": "Aspect Ratio Worst",
        "path": ("flat_mesh_quality", "aspect_ratio", "max"),
        "better": "lower",
        "hist_log": False,
        "heatmap": True,
    },
    {
        "key": "er_worst",
        "label": "Edge Ratio Worst",
        "path": ("flat_mesh_quality", "edge_ratio", "max"),
        "better": "lower",
        "hist_log": False,
        "heatmap": True,
    },
    {
        "key": "skew_worst",
        "label": "Skewness Worst",
        "path": ("flat_mesh_quality", "skewness", "max"),
        "better": "lower",
        "hist_log": False,
        "heatmap": False,
    },
    {
        "key": "symdiff_area",
        "label": "Polygon SymDiff Area (m^2)",
        "path": ("raw_to_conditioned_difference_metrics", "symmetric_difference_area"),
        "better": "lower",
        "hist_log": True,
        "heatmap": True,
    },
    {
        "key": "missing_area",
        "label": "Polygon Missing Area (m^2)",
        "path": ("raw_to_conditioned_difference_metrics", "reference_minus_candidate_area"),
        "better": "lower",
        "hist_log": True,
        "heatmap": True,
    },
    {
        "key": "extra_area",
        "label": "Polygon Extra Area (m^2)",
        "path": ("raw_to_conditioned_difference_metrics", "candidate_minus_reference_area"),
        "better": "lower",
        "hist_log": True,
        "heatmap": False,
    },
    {
        "key": "conditioned_polygon_count",
        "label": "Conditioned Polygon Count",
        "path": ("conditioned_polygon_count",),
        "better": "lower",
        "hist_log": False,
        "heatmap": True,
    },
    {
        "key": "conditioning_seconds",
        "label": "Conditioning Time (s)",
        "path": ("timing_summary", "conditioning_seconds"),
        "better": "lower",
        "hist_log": False,
        "heatmap": False,
    },
    {
        "key": "core_seconds",
        "label": "Core Time (s)",
        "path": ("timing_summary", "core_seconds"),
        "better": "lower",
        "hist_log": False,
        "heatmap": False,
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
    )
    parser.add_argument(
        "--legacy-commit",
        default=DEFAULT_LEGACY_COMMIT,
    )
    parser.add_argument(
        "--legacy-root",
        type=Path,
        default=DEFAULT_LEGACY_ROOT,
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.0,
        help="Delay forwarded to the per-case harness.",
    )
    parser.add_argument(
        "--force-rerun",
        action="store_true",
        help="Recompute all cases even if per-case summaries already exist.",
    )
    parser.add_argument(
        "--skip-runs",
        action="store_true",
        help="Skip the legacy/new runs and only rebuild aggregated outputs from existing summaries.",
    )
    parser.add_argument(
        "--disable-cleaning-diagnostics",
        action="store_true",
        help=(
            "Forwarded to the new cleaner runs to disable cleaner stage metrics "
            "and cleaner log output for runtime comparisons."
        ),
    )
    parser.add_argument(
        "--legacy-pythonpath-prefix",
        type=Path,
        default=None,
        help=(
            "Optional path prepended to PYTHONPATH for legacy runs, for example "
            "a vendored Polyforge checkout."
        ),
    )
    return parser.parse_args()


def _git(args: list[str], *, cwd: Path) -> str:
    result = subprocess.run(
        args,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def ensure_legacy_worktree(legacy_root: Path, legacy_commit: str) -> None:
    if not legacy_root.exists():
        subprocess.run(
            ["git", "worktree", "add", str(legacy_root), legacy_commit],
            cwd=REPO_ROOT,
            check=True,
        )

    legacy_builder = legacy_root / LEGACY_DTCC_BUILDER_REL
    if not legacy_builder.exists():
        legacy_builder.symlink_to(CURRENT_DTCC_BUILDER)


def _dyld_fallback_library_path() -> str:
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if not conda_prefix:
        conda_prefix = str(Path(sys.executable).resolve().parent.parent)
    pieces = []
    if conda_prefix:
        pieces.append(str(Path(conda_prefix) / "lib"))
        pieces.append(
            str(Path(conda_prefix) / "lib" / "python3.12" / "site-packages" / "pyspade_native" / "lib")
        )
    existing = os.environ.get("DYLD_FALLBACK_LIBRARY_PATH")
    if existing:
        pieces.append(existing)
    return ":".join(piece for piece in pieces if piece)


def _python_site_packages() -> str:
    return sysconfig.get_paths()["purelib"]


def _isolated_python_command(
    *,
    script_path: Path,
    argv: list[str],
    pythonpath_entries: list[str],
) -> list[str]:
    launcher = (
        "import runpy, sys\n"
        f"sys.path[:0] = {pythonpath_entries!r}\n"
        f"sys.argv = {argv!r}\n"
        f"runpy.run_path({str(script_path)!r}, run_name='__main__')\n"
    )
    return [sys.executable, "-S", "-c", launcher]


def _mode_root(output_root: Path, label: str, mode: str) -> Path:
    return output_root / label / mode


def _case_summary_path(output_root: Path, label: str, mode: str, case: int) -> Path:
    return _mode_root(output_root, label, mode) / f"case_{case:03d}" / "summary.json"


def missing_cases(output_root: Path, label: str, mode: str, *, force_rerun: bool) -> list[int]:
    if force_rerun:
        return CASES
    return [
        case
        for case in CASES
        if not _case_summary_path(output_root, label, mode, case).exists()
    ]


def run_harness(
    *,
    repo_root: Path,
    git_root: Path,
    output_root: Path,
    label: str,
    mode: str,
    cases: list[int],
    delay: float,
    disable_cleaning_diagnostics: bool,
    pythonpath_prefix: Path | None = None,
) -> None:
    if not cases:
        return

    (repo_root / "sandbox" / "output").mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    pythonpath_entries: list[str] = []
    if pythonpath_prefix is not None:
        pythonpath_entries.append(str(pythonpath_prefix))
    pythonpath_entries.append(str(repo_root))
    pythonpath_entries.append(_python_site_packages())
    env["DYLD_FALLBACK_LIBRARY_PATH"] = _dyld_fallback_library_path()

    argv = [
        str(COMPARE_SCRIPT),
        "--output-root",
        str(output_root),
        "--git-root",
        str(git_root),
        "--label",
        label,
        "--mode",
        mode,
        "--delay",
        str(delay),
        "--summary-only",
        "--cases",
        *[str(case) for case in cases],
    ]
    if disable_cleaning_diagnostics and mode == "new":
        argv.append("--disable-cleaning-diagnostics")

    cmd = _isolated_python_command(
        script_path=COMPARE_SCRIPT,
        argv=argv,
        pythonpath_entries=pythonpath_entries,
    )
    subprocess.run(cmd, cwd=repo_root, env=env, check=True)


def _get_nested(mapping: dict[str, Any], path: tuple[str, ...]) -> float | int | None:
    value: Any = mapping
    for part in path:
        if not isinstance(value, dict):
            return None
        value = value.get(part)
    return value


def case_to_grid(case: int) -> tuple[int, int]:
    iy, ix = divmod(case - 1, NX)
    return ix, iy


def extract_metrics(summary: dict[str, Any]) -> dict[str, float | int | None]:
    metrics = {}
    for spec in METRIC_SPECS:
        metrics[spec["key"]] = _get_nested(summary, spec["path"])
    return metrics


def _delta_value(new_value: float | int | None, legacy_value: float | int | None) -> float | None:
    if new_value is None or legacy_value is None:
        return None
    return float(new_value) - float(legacy_value)


def _aggregate(values: list[float]) -> dict[str, float]:
    if not values:
        return {
            "count": 0,
            "mean": None,
            "median": None,
            "min": None,
            "max": None,
            "p10": None,
            "p90": None,
        }
    array = np.asarray(values, dtype=float)
    return {
        "count": int(array.size),
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "p10": float(np.percentile(array, 10)),
        "p90": float(np.percentile(array, 90)),
    }


def build_case_records(output_root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for case in CASES:
        legacy_summary_path = _case_summary_path(output_root, LEGACY_LABEL, "legacy", case)
        new_summary_path = _case_summary_path(output_root, NEW_LABEL, "new", case)
        if not legacy_summary_path.exists():
            raise FileNotFoundError(f"Missing legacy summary for case {case}: {legacy_summary_path}")
        if not new_summary_path.exists():
            raise FileNotFoundError(f"Missing new summary for case {case}: {new_summary_path}")

        legacy_summary = json.loads(legacy_summary_path.read_text())
        new_summary = json.loads(new_summary_path.read_text())
        legacy_metrics = extract_metrics(legacy_summary)
        new_metrics = extract_metrics(new_summary)
        deltas = {
            key: _delta_value(new_metrics[key], legacy_metrics[key])
            for key in legacy_metrics
        }
        ix, iy = case_to_grid(case)
        records.append(
            {
                "case": case,
                "ix": ix,
                "iy": iy,
                "bounds": legacy_summary["bounds"],
                "legacy_status": legacy_summary.get("status", "completed"),
                "new_status": new_summary.get("status", "completed"),
                "legacy_error": legacy_summary.get("error"),
                "new_error": new_summary.get("error"),
                "legacy_summary_path": str(legacy_summary_path),
                "new_summary_path": str(new_summary_path),
                "legacy_metrics": legacy_metrics,
                "new_metrics": new_metrics,
                "delta_new_minus_legacy": deltas,
                "legacy_summary": legacy_summary,
                "new_summary": new_summary,
            }
        )
    return records


def build_aggregate_summary(case_records: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate: dict[str, Any] = {"metrics": {}}
    for spec in METRIC_SPECS:
        key = spec["key"]
        legacy_values = [
            float(record["legacy_metrics"][key])
            for record in case_records
            if record["legacy_metrics"][key] is not None
        ]
        new_values = [
            float(record["new_metrics"][key])
            for record in case_records
            if record["new_metrics"][key] is not None
        ]
        delta_values = [
            float(record["delta_new_minus_legacy"][key])
            for record in case_records
            if record["delta_new_minus_legacy"][key] is not None
        ]

        better = spec["better"]
        epsilon = 1e-9
        if better == "higher":
            new_wins = sum(delta > epsilon for delta in delta_values)
            legacy_wins = sum(delta < -epsilon for delta in delta_values)
        else:
            new_wins = sum(delta < -epsilon for delta in delta_values)
            legacy_wins = sum(delta > epsilon for delta in delta_values)
        ties = len(delta_values) - new_wins - legacy_wins

        aggregate["metrics"][key] = {
            "label": spec["label"],
            "better": better,
            "legacy": _aggregate(legacy_values),
            "new": _aggregate(new_values),
            "delta_new_minus_legacy": _aggregate(delta_values),
            "new_wins": new_wins,
            "legacy_wins": legacy_wins,
            "ties": ties,
        }
    aggregate["status_summary"] = {
        "legacy": {
            "completed": sum(record["legacy_status"] == "completed" for record in case_records),
            "failed": sum(record["legacy_status"] != "completed" for record in case_records),
        },
        "new": {
            "completed": sum(record["new_status"] == "completed" for record in case_records),
            "failed": sum(record["new_status"] != "completed" for record in case_records),
        },
    }
    aggregate["failed_cases"] = [
        {
            "case": record["case"],
            "legacy_status": record["legacy_status"],
            "new_status": record["new_status"],
            "legacy_error": record["legacy_error"],
            "new_error": record["new_error"],
        }
        for record in case_records
        if record["legacy_status"] != "completed" or record["new_status"] != "completed"
    ]
    return aggregate


def write_combined_json(
    output_root: Path,
    *,
    legacy_commit: str,
    legacy_root: Path,
    legacy_pythonpath_prefix: Path | None,
    case_records: list[dict[str, Any]],
    aggregate_summary: dict[str, Any],
) -> Path:
    payload = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "repo_root": str(REPO_ROOT),
        "legacy_root": str(legacy_root),
        "legacy_commit": legacy_commit,
        "legacy_pythonpath_prefix": (
            str(legacy_pythonpath_prefix) if legacy_pythonpath_prefix is not None else None
        ),
        "new_commit": _git(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT),
        "new_branch": _git(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=REPO_ROOT),
        "case_count": len(case_records),
        "cases": case_records,
        "aggregate_summary": aggregate_summary,
    }
    path = output_root / "stockholm_100_tile_combined_results.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def write_case_table_csv(output_root: Path, case_records: list[dict[str, Any]]) -> Path:
    path = output_root / "stockholm_100_tile_case_table.csv"
    fieldnames = [
        "case",
        "ix",
        "iy",
        "legacy_status",
        "new_status",
        "legacy_error_stage",
        "new_error_stage",
        "eq_worst_legacy",
        "eq_worst_new",
        "eq_worst_delta",
        "eq_mean_legacy",
        "eq_mean_new",
        "eq_mean_delta",
        "ar_worst_legacy",
        "ar_worst_new",
        "ar_worst_delta",
        "er_worst_legacy",
        "er_worst_new",
        "er_worst_delta",
        "skew_worst_legacy",
        "skew_worst_new",
        "skew_worst_delta",
        "symdiff_legacy",
        "symdiff_new",
        "symdiff_delta",
        "missing_legacy",
        "missing_new",
        "missing_delta",
        "extra_legacy",
        "extra_new",
        "extra_delta",
        "conditioned_polygon_count_legacy",
        "conditioned_polygon_count_new",
        "conditioned_polygon_count_delta",
        "conditioning_seconds_legacy",
        "conditioning_seconds_new",
        "conditioning_seconds_delta",
        "core_seconds_legacy",
        "core_seconds_new",
        "core_seconds_delta",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in case_records:
            legacy = record["legacy_metrics"]
            new = record["new_metrics"]
            delta = record["delta_new_minus_legacy"]
            writer.writerow(
                {
                    "case": record["case"],
                    "ix": record["ix"],
                    "iy": record["iy"],
                    "legacy_status": record["legacy_status"],
                    "new_status": record["new_status"],
                    "legacy_error_stage": (record["legacy_error"] or {}).get("stage"),
                    "new_error_stage": (record["new_error"] or {}).get("stage"),
                    "eq_worst_legacy": legacy["eq_worst"],
                    "eq_worst_new": new["eq_worst"],
                    "eq_worst_delta": delta["eq_worst"],
                    "eq_mean_legacy": legacy["eq_mean"],
                    "eq_mean_new": new["eq_mean"],
                    "eq_mean_delta": delta["eq_mean"],
                    "ar_worst_legacy": legacy["ar_worst"],
                    "ar_worst_new": new["ar_worst"],
                    "ar_worst_delta": delta["ar_worst"],
                    "er_worst_legacy": legacy["er_worst"],
                    "er_worst_new": new["er_worst"],
                    "er_worst_delta": delta["er_worst"],
                    "skew_worst_legacy": legacy["skew_worst"],
                    "skew_worst_new": new["skew_worst"],
                    "skew_worst_delta": delta["skew_worst"],
                    "symdiff_legacy": legacy["symdiff_area"],
                    "symdiff_new": new["symdiff_area"],
                    "symdiff_delta": delta["symdiff_area"],
                    "missing_legacy": legacy["missing_area"],
                    "missing_new": new["missing_area"],
                    "missing_delta": delta["missing_area"],
                    "extra_legacy": legacy["extra_area"],
                    "extra_new": new["extra_area"],
                    "extra_delta": delta["extra_area"],
                    "conditioned_polygon_count_legacy": legacy["conditioned_polygon_count"],
                    "conditioned_polygon_count_new": new["conditioned_polygon_count"],
                    "conditioned_polygon_count_delta": delta["conditioned_polygon_count"],
                    "conditioning_seconds_legacy": legacy["conditioning_seconds"],
                    "conditioning_seconds_new": new["conditioning_seconds"],
                    "conditioning_seconds_delta": delta["conditioning_seconds"],
                    "core_seconds_legacy": legacy["core_seconds"],
                    "core_seconds_new": new["core_seconds"],
                    "core_seconds_delta": delta["core_seconds"],
                }
            )
    return path


def _fmt_md_value(value: float | int | None, digits: int = 3) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def write_case_table_md(output_root: Path, case_records: list[dict[str, Any]]) -> Path:
    path = output_root / "stockholm_100_tile_case_table.md"
    lines = [
        "| Case | Legacy Status | New Status | EQ Worst Old | EQ Worst New | AR Worst Old | AR Worst New | ER Worst Old | ER Worst New | SymDiff Old | SymDiff New | Missing Old | Missing New |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for record in case_records:
        legacy = record["legacy_metrics"]
        new = record["new_metrics"]
        lines.append(
            "| {case} | {legacy_status} | {new_status} | {leq} | {neq} | {lar} | {nar} | {ler} | {ner} | {ls} | {ns} | {lm} | {nm} |".format(
                case=record["case"],
                legacy_status=record["legacy_status"],
                new_status=record["new_status"],
                leq=_fmt_md_value(legacy["eq_worst"]),
                neq=_fmt_md_value(new["eq_worst"]),
                lar=_fmt_md_value(legacy["ar_worst"]),
                nar=_fmt_md_value(new["ar_worst"]),
                ler=_fmt_md_value(legacy["er_worst"]),
                ner=_fmt_md_value(new["er_worst"]),
                ls=_fmt_md_value(legacy["symdiff_area"], digits=1),
                ns=_fmt_md_value(new["symdiff_area"], digits=1),
                lm=_fmt_md_value(legacy["missing_area"], digits=1),
                nm=_fmt_md_value(new["missing_area"], digits=1),
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_overview_md(output_root: Path, aggregate_summary: dict[str, Any]) -> Path:
    path = output_root / "stockholm_100_tile_overview.md"
    legacy_status = aggregate_summary["status_summary"]["legacy"]
    new_status = aggregate_summary["status_summary"]["new"]
    lines = [
        "# Stockholm 100-tile comparison",
        "",
        f"Legacy completed: {legacy_status['completed']} | Legacy failed: {legacy_status['failed']}",
        f"New completed: {new_status['completed']} | New failed: {new_status['failed']}",
        "",
        "| Metric | Better | Legacy N | Legacy Mean | New N | New Mean | Mean Delta (new-old) | Median Delta | New Wins | Legacy Wins | Ties |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    failed_cases = aggregate_summary.get("failed_cases", [])
    if failed_cases:
        lines.append("")
        lines.append("## Failed cases")
        lines.append("")
        for failed_case in failed_cases:
            lines.append(
                f"- case {failed_case['case']}: legacy={failed_case['legacy_status']} new={failed_case['new_status']}"
            )
        lines.append("")
    for spec in METRIC_SPECS:
        key = spec["key"]
        aggregate = aggregate_summary["metrics"][key]
        lines.append(
            "| {label} | {better} | {ln} | {lmean} | {nn} | {nmean} | {dmean} | {dmedian} | {new_wins} | {legacy_wins} | {ties} |".format(
                label=spec["label"],
                better=aggregate["better"],
                ln=aggregate["legacy"]["count"],
                lmean=_fmt_md_value(aggregate["legacy"]["mean"], digits=4),
                nn=aggregate["new"]["count"],
                nmean=_fmt_md_value(aggregate["new"]["mean"], digits=4),
                dmean=_fmt_md_value(aggregate["delta_new_minus_legacy"]["mean"], digits=4),
                dmedian=_fmt_md_value(aggregate["delta_new_minus_legacy"]["median"], digits=4),
                new_wins=aggregate["new_wins"],
                legacy_wins=aggregate["legacy_wins"],
                ties=aggregate["ties"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _values(case_records: list[dict[str, Any]], key: str, side: str) -> np.ndarray:
    return np.asarray(
        [
            float(value)
            for record in case_records
            for value in [record[f"{side}_metrics"][key]]
            if value is not None
        ],
        dtype=float,
    )


def plot_distribution_panels(output_root: Path, case_records: list[dict[str, Any]]) -> Path:
    specs = [spec for spec in METRIC_SPECS if spec["key"] in {"eq_worst", "eq_mean", "ar_worst", "er_worst", "skew_worst", "symdiff_area", "missing_area", "conditioned_polygon_count"}]
    fig, axes = plt.subplots(2, 4, figsize=(18, 10), constrained_layout=True)
    axes = axes.ravel()
    colors = {"legacy": "#4C78A8", "new": "#F58518"}

    for ax, spec in zip(axes, specs):
        legacy_values = _values(case_records, spec["key"], "legacy")
        new_values = _values(case_records, spec["key"], "new")
        if legacy_values.size == 0 and new_values.size == 0:
            ax.set_title(spec["label"])
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue
        combined = np.concatenate([legacy_values, new_values])

        if spec["hist_log"]:
            minimum = max(np.min(combined[combined > 0]), 1e-6)
            maximum = np.max(combined)
            bins = np.logspace(np.log10(minimum), np.log10(maximum), 20)
            ax.set_xscale("log")
        else:
            bins = 20

        ax.hist(legacy_values, bins=bins, alpha=0.55, color=colors["legacy"], label="Legacy")
        ax.hist(new_values, bins=bins, alpha=0.55, color=colors["new"], label="New")
        if legacy_values.size:
            ax.axvline(np.mean(legacy_values), color=colors["legacy"], linestyle="--", linewidth=1.2)
        if new_values.size:
            ax.axvline(np.mean(new_values), color=colors["new"], linestyle="--", linewidth=1.2)
        ax.set_title(f"{spec['label']}\nlegacy n={legacy_values.size}, new n={new_values.size}")
        ax.grid(alpha=0.25)

    axes[0].legend(loc="best")
    fig.suptitle("Stockholm 100-tile metric distributions")
    path = output_root / "stockholm_100_tile_distribution_panels.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_delta_histograms(output_root: Path, case_records: list[dict[str, Any]]) -> Path:
    specs = [spec for spec in METRIC_SPECS if spec["key"] in {"eq_worst", "eq_mean", "ar_worst", "er_worst", "symdiff_area", "missing_area", "conditioning_seconds", "core_seconds"}]
    fig, axes = plt.subplots(2, 4, figsize=(18, 10), constrained_layout=True)
    axes = axes.ravel()

    for ax, spec in zip(axes, specs):
        deltas = np.asarray(
            [
                float(delta)
                for record in case_records
                for delta in [record["delta_new_minus_legacy"][spec["key"]]]
                if delta is not None
            ],
            dtype=float,
        )
        if deltas.size == 0:
            ax.set_title(spec["label"])
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue
        ax.hist(deltas, bins=20, color="#72B7B2", alpha=0.85)
        ax.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
        ax.axvline(np.mean(deltas), color="#E45756", linestyle="-", linewidth=1.2)
        direction = "higher is better" if spec["better"] == "higher" else "lower is better"
        ax.set_title(f"{spec['label']}\nnew - legacy ({direction}), n={deltas.size}")
        ax.grid(alpha=0.25)

    fig.suptitle("Stockholm 100-tile delta histograms")
    path = output_root / "stockholm_100_tile_delta_histograms.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_delta_heatmaps(output_root: Path, case_records: list[dict[str, Any]]) -> Path:
    specs = [spec for spec in METRIC_SPECS if spec["heatmap"]][:6]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes = axes.ravel()
    cmap = plt.get_cmap("RdYlGn").copy()
    cmap.set_bad(color="#e6e6e6")

    for ax, spec in zip(axes, specs):
        grid = np.full((NY, NX), np.nan, dtype=float)
        values = []
        for record in case_records:
            raw_delta = record["delta_new_minus_legacy"][spec["key"]]
            if raw_delta is None:
                continue
            signed_advantage = float(raw_delta)
            if spec["better"] == "lower":
                signed_advantage = -signed_advantage
            grid[record["iy"], record["ix"]] = signed_advantage
            values.append(signed_advantage)

        vmax = max(abs(np.min(values)), abs(np.max(values))) if values else 1.0
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        image = ax.imshow(grid, origin="lower", cmap=cmap, norm=norm)
        ax.set_title(f"{spec['label']}\nSigned advantage: + new better, - legacy better")
        ax.set_xticks(range(NX))
        ax.set_yticks(range(NY))
        ax.set_xlabel("ix")
        ax.set_ylabel("iy")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Green = new better, red = legacy better")

    fig.suptitle(
        "Stockholm 100-tile heatmaps\n"
        "Sign normalized so positive/green always means the new cleaner is better"
    )
    path = output_root / "stockholm_100_tile_delta_heatmaps.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    ensure_legacy_worktree(args.legacy_root, args.legacy_commit)

    if not args.skip_runs:
        legacy_missing = missing_cases(
            args.output_root,
            LEGACY_LABEL,
            "legacy",
            force_rerun=args.force_rerun,
        )
        new_missing = missing_cases(
            args.output_root,
            NEW_LABEL,
            "new",
            force_rerun=args.force_rerun,
        )

        run_harness(
            repo_root=args.legacy_root,
            git_root=args.legacy_root,
            output_root=args.output_root,
            label=LEGACY_LABEL,
            mode="legacy",
            cases=legacy_missing,
            delay=args.delay,
            disable_cleaning_diagnostics=args.disable_cleaning_diagnostics,
            pythonpath_prefix=args.legacy_pythonpath_prefix,
        )
        run_harness(
            repo_root=REPO_ROOT,
            git_root=REPO_ROOT,
            output_root=args.output_root,
            label=NEW_LABEL,
            mode="new",
            cases=new_missing,
            delay=args.delay,
            disable_cleaning_diagnostics=args.disable_cleaning_diagnostics,
            pythonpath_prefix=None,
        )

    case_records = build_case_records(args.output_root)
    case_records.sort(key=lambda record: record["case"])
    aggregate_summary = build_aggregate_summary(case_records)

    combined_json = write_combined_json(
        args.output_root,
        legacy_commit=args.legacy_commit,
        legacy_root=args.legacy_root,
        legacy_pythonpath_prefix=args.legacy_pythonpath_prefix,
        case_records=case_records,
        aggregate_summary=aggregate_summary,
    )
    case_table_csv = write_case_table_csv(args.output_root, case_records)
    case_table_md = write_case_table_md(args.output_root, case_records)
    overview_md = write_overview_md(args.output_root, aggregate_summary)
    dist_plot = plot_distribution_panels(args.output_root, case_records)
    delta_hist_plot = plot_delta_histograms(args.output_root, case_records)
    heatmap_plot = plot_delta_heatmaps(args.output_root, case_records)

    print(combined_json)
    print(case_table_csv)
    print(case_table_md)
    print(overview_md)
    print(dist_plot)
    print(delta_hist_plot)
    print(heatmap_plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
