"""
Run a quick single-tile benchmark comparison for footprints, mesh_2d, and mesh_3d.

Edit ``TILE`` below to choose the Stockholm tile. This wrapper intentionally
takes no command-line arguments; it runs the three benchmark scripts for both
``legacy`` and ``new`` conditioning, stores all outputs under
``benchmarks/bench_test/``, and prints one combined summary at the end.

Typical usage:
    /Users/logg/scratch/dtcc/venv/bin/python benchmarks/bench_test.py
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import bench_footprints
    import bench_mesh_2d
    import bench_mesh_3d
    from _stockholm_common import format_console_table, json_ready, load_results
except ImportError:
    from benchmarks import bench_footprints, bench_mesh_2d, bench_mesh_3d
    from benchmarks._stockholm_common import format_console_table, json_ready, load_results


TILE = 54
MODES = ("legacy", "new")
BENCH_2D_MESHERS = ("dtcc_mesher", "spade")

BENCHMARKS_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = BENCHMARKS_DIR / "bench_test"
RUN_OUTPUT_DIR = OUTPUT_ROOT / f"{TILE:03d}"
COMBINED_RESULTS_PATH = RUN_OUTPUT_DIR / "results.json"
COMBINED_SUMMARY_PATH = RUN_OUTPUT_DIR / "summary.txt"

MESH_2D_MAX_MESH_SIZE = bench_mesh_2d.DEFAULT_MAX_MESH_SIZE
MESH_2D_MIN_MESH_ANGLE = bench_mesh_2d.MIN_MESH_ANGLE
MESH_2D_PIPELINE_MODE = "strict"
MESH_2D_STAGE_AUDIT = True

MESH_3D_MAX_MESH_SIZE = bench_mesh_3d.DEFAULT_MAX_MESH_SIZE
MESH_3D_MIN_MESH_ANGLE = bench_mesh_3d.MIN_MESH_ANGLE
MESH_3D_DOMAIN_HEIGHT = bench_mesh_3d.DEFAULT_DOMAIN_HEIGHT
MESH_3D_QUALITY_RATIO = bench_mesh_3d.DEFAULT_QUALITY_RATIO
MESH_3D_QUALITY_ENABLED = True
MESH_3D_PRESERVE_SURFACE = False
MESH_3D_MAX_ADDED_POINTS = bench_mesh_3d.DEFAULT_MAX_ADDED_POINTS
MESH_3D_LOD = bench_mesh_3d.DEFAULT_LOD
MESH_3D_MERGE_BUILDINGS = bench_mesh_3d.MERGE_BUILDINGS
MESH_3D_MESHER = bench_mesh_3d.DEFAULT_MESHER
MESH_3D_PIPELINE_MODE = "strict"
MESH_3D_SAVE_TETGEN_INPUT = True
MESH_3D_STAGE_AUDIT = True


@dataclass
class RunResult:
    benchmark: str
    mode: str
    command: list[str]
    output_dir: Path
    case_dir: Path | None
    results_path: Path | None
    summary_path: Path | None
    returncode: int
    elapsed_seconds: float
    status: str
    note: str | None
    case_record: dict[str, Any] | None

    @property
    def name(self) -> str:
        return f"{self.mode}:{self.benchmark}"


def mode_output_dir(mode: str, benchmark: str) -> Path:
    return RUN_OUTPUT_DIR / mode / benchmark


def case_output_dir(output_dir: Path) -> Path:
    return output_dir / f"{TILE:03d}"


def display_path(path: Path | None) -> str:
    if path is None:
        return "-"
    try:
        return str(path.relative_to(RUN_OUTPUT_DIR))
    except ValueError:
        return str(path)


def fmt(value: Any, digits: int = 2) -> str:
    if value is None:
        return "-"
    if isinstance(value, (int, float)):
        return f"{float(value):.{digits}f}"
    return str(value)


def derive_status(returncode: int, case_record: dict[str, Any] | None) -> str:
    if returncode != 0:
        return f"failed({returncode})"
    if case_record is None:
        return "no-results"
    benchmark_status = str(case_record.get("status", "")).lower()
    if benchmark_status == "completed":
        return "ok"
    result = case_record.get("result", {})
    if isinstance(result, dict) and result.get("status") == "success":
        return "ok"
    if benchmark_status:
        return benchmark_status
    if isinstance(result, dict) and result.get("status"):
        return str(result.get("status"))
    return "partial"


def summarize_error(error: Any) -> str | None:
    if error is None:
        return None
    if isinstance(error, dict):
        message = error.get("message")
        if isinstance(message, str) and message.strip():
            return message.strip()
        error_type = error.get("type")
        if isinstance(error_type, str) and error_type.strip():
            return error_type.strip()
        return None
    text = str(error).strip()
    return text or None


def shorten_note(note: str | None, limit: int = 72) -> str:
    if note is None:
        return "-"
    compact = " ".join(note.split())
    if len(compact) <= limit:
        return compact
    return compact[: limit - 3].rstrip() + "..."


def status_badge(status: str) -> str:
    normalized = status.strip().lower()
    if normalized == "ok":
        return "✅ ok"
    if normalized in {"fail", "failed"} or normalized.startswith("failed("):
        return "❌ fail"
    if normalized == "partial":
        return "⚠️ partial"
    if normalized == "no-results":
        return "⚪ no-results"
    return status


def section_heading(title: str) -> str:
    return f"{title}\n" + ("=" * len(title))


def derive_mesh_2d_status_and_note(
    returncode: int,
    case_record: dict[str, Any] | None,
    meshers: list[str],
) -> tuple[str, str | None]:
    if returncode != 0:
        return f"failed({returncode})", None
    if case_record is None:
        return "no-results", None

    mesher_results = case_record.get("meshers", {})
    if not isinstance(mesher_results, dict):
        return derive_status(returncode, case_record), None

    statuses = [mesher_results.get(mesher, {}).get("status") for mesher in meshers]
    if all(status == "success" for status in statuses) and bench_mesh_2d.case_complete(
        case_record, meshers
    ):
        return "ok", None

    failed_messages = [
        summarize_error(mesher_results.get(mesher, {}).get("error"))
        for mesher in meshers
        if mesher_results.get(mesher, {}).get("status") == "failed"
    ]
    note = next((message for message in failed_messages if message), None)

    if all(status == "failed" for status in statuses):
        return "fail", note
    if any(status == "failed" for status in statuses):
        return "partial", note
    return "partial", note


def derive_mesh_3d_status_and_note(
    returncode: int,
    case_record: dict[str, Any] | None,
) -> tuple[str, str | None]:
    if returncode != 0:
        return f"failed({returncode})", None
    if case_record is None:
        return "no-results", None

    result = case_record.get("result", {})
    status = result.get("status")
    if status == "success" and bench_mesh_3d.case_complete(case_record):
        return "ok", None
    if status == "failed":
        return "fail", summarize_error(result.get("error"))
    return derive_status(returncode, case_record), summarize_error(result.get("error"))


def run_subprocess(name: str, command: list[str]) -> tuple[int, float]:
    print("=" * 80)
    print(f"{name}: {shlex.join(command)}")
    print("=" * 80)
    print()

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    started = time.perf_counter()
    completed = subprocess.run(command, cwd=BENCHMARKS_DIR, env=env)
    elapsed = time.perf_counter() - started

    print()
    return completed.returncode, elapsed


def build_footprints_command(mode: str) -> list[str]:
    output_dir = mode_output_dir(mode, "footprints")
    return [
        sys.executable,
        "-u",
        str(BENCHMARKS_DIR / "bench_footprints.py"),
        "--cases",
        str(TILE),
        "--mode",
        mode,
        "--label",
        mode,
        "--delay",
        "0",
        "--output-dir",
        str(output_dir),
        "--no-plots",
    ]


def build_mesh_2d_command(mode: str) -> list[str]:
    output_dir = mode_output_dir(mode, "mesh_2d")
    command = [
        sys.executable,
        "-u",
        str(BENCHMARKS_DIR / "bench_mesh_2d.py"),
        "--cases",
        str(TILE),
        "--meshers",
        *BENCH_2D_MESHERS,
        "--conditioning-mode",
        mode,
        "--delay",
        "0",
        "--max-mesh-size",
        "none" if MESH_2D_MAX_MESH_SIZE is None else str(MESH_2D_MAX_MESH_SIZE),
        "--output-dir",
        str(output_dir),
        "--no-plots",
    ]
    if MESH_2D_STAGE_AUDIT:
        command.append("--stage-audit")
    return command


def build_mesh_3d_command(mode: str) -> list[str]:
    output_dir = mode_output_dir(mode, "mesh_3d")
    command = [
        sys.executable,
        "-u",
        str(BENCHMARKS_DIR / "bench_mesh_3d.py"),
        "--cases",
        str(TILE),
        "--conditioning-mode",
        mode,
        "--delay",
        "0",
        "--max-mesh-size",
        str(MESH_3D_MAX_MESH_SIZE),
        "--min-mesh-angle",
        str(MESH_3D_MIN_MESH_ANGLE),
        "--domain-height",
        str(MESH_3D_DOMAIN_HEIGHT),
        "--quality-ratio",
        str(MESH_3D_QUALITY_RATIO),
        "--lod",
        MESH_3D_LOD.name,
        "--mesher",
        MESH_3D_MESHER,
        "--output-dir",
        str(output_dir),
        "--no-plots",
    ]
    if MESH_3D_MERGE_BUILDINGS:
        command.append("--merge-buildings")
    else:
        command.append("--no-merge-buildings")
    if not MESH_3D_QUALITY_ENABLED:
        command.append("--size-only")
    if MESH_3D_PRESERVE_SURFACE:
        command.append("--preserve-surface")
    if MESH_3D_SAVE_TETGEN_INPUT:
        command.append("--save-tetgen-input")
    else:
        command.append("--no-save-tetgen-input")
    if MESH_3D_STAGE_AUDIT:
        command.append("--stage-audit")
    if MESH_3D_MAX_ADDED_POINTS is not None:
        command.extend(["--max-added-points", str(MESH_3D_MAX_ADDED_POINTS)])
    return command


def collect_footprints_result(
    mode: str,
    returncode: int,
    elapsed_seconds: float,
    command: list[str],
) -> RunResult:
    output_dir = mode_output_dir(mode, "footprints")
    results_path = output_dir / "results.json"
    summary_path = output_dir / "summary.txt"
    case_dir = case_output_dir(output_dir)
    results = load_results(results_path) if results_path.exists() else {}
    case_record = results.get(TILE)
    return RunResult(
        benchmark="footprints",
        mode=mode,
        command=command,
        output_dir=output_dir,
        case_dir=case_dir if case_dir.exists() else None,
        results_path=results_path if results_path.exists() else None,
        summary_path=summary_path if summary_path.exists() else None,
        returncode=returncode,
        elapsed_seconds=elapsed_seconds,
        status=derive_status(returncode, case_record),
        note=None,
        case_record=case_record,
    )


def collect_mesh_2d_result(
    mode: str,
    returncode: int,
    elapsed_seconds: float,
    command: list[str],
) -> RunResult:
    output_dir = mode_output_dir(mode, "mesh_2d")
    meshers = bench_mesh_2d.resolve_requested_meshers(list(BENCH_2D_MESHERS))
    results_path = bench_mesh_2d.results_file_path(
        output_dir,
        meshers,
        max_mesh_size=MESH_2D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
        pipeline_mode=MESH_2D_PIPELINE_MODE,
        stage_audit_enabled=MESH_2D_STAGE_AUDIT,
    )
    summary_path = bench_mesh_2d.summary_text_path(
        output_dir,
        meshers,
        max_mesh_size=MESH_2D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
        pipeline_mode=MESH_2D_PIPELINE_MODE,
        stage_audit_enabled=MESH_2D_STAGE_AUDIT,
    )
    case_dir = bench_mesh_2d.case_output_dir(
        output_dir,
        TILE,
        meshers,
        max_mesh_size=MESH_2D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
        pipeline_mode=MESH_2D_PIPELINE_MODE,
        stage_audit_enabled=MESH_2D_STAGE_AUDIT,
    )
    results = load_results(results_path) if results_path.exists() else {}
    case_record = results.get(TILE)
    status, note = derive_mesh_2d_status_and_note(returncode, case_record, meshers)
    return RunResult(
        benchmark="mesh_2d",
        mode=mode,
        command=command,
        output_dir=output_dir,
        case_dir=case_dir if case_dir.exists() else None,
        results_path=results_path if results_path.exists() else None,
        summary_path=summary_path if summary_path.exists() else None,
        returncode=returncode,
        elapsed_seconds=elapsed_seconds,
        status=status,
        note=note,
        case_record=case_record,
    )


def collect_mesh_3d_result(
    mode: str,
    returncode: int,
    elapsed_seconds: float,
    command: list[str],
) -> RunResult:
    output_dir = mode_output_dir(mode, "mesh_3d")
    results_path = bench_mesh_3d.results_file_path(
        output_dir,
        max_mesh_size=MESH_3D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_3D_MIN_MESH_ANGLE,
        domain_height=MESH_3D_DOMAIN_HEIGHT,
        quality_ratio=MESH_3D_QUALITY_RATIO,
        quality_enabled=MESH_3D_QUALITY_ENABLED,
        preserve_surface=MESH_3D_PRESERVE_SURFACE,
        max_added_points=MESH_3D_MAX_ADDED_POINTS,
        lod=MESH_3D_LOD,
        merge_buildings=MESH_3D_MERGE_BUILDINGS,
        mesher=MESH_3D_MESHER,
        pipeline_mode=MESH_3D_PIPELINE_MODE,
        stage_audit_enabled=MESH_3D_STAGE_AUDIT,
    )
    summary_path = bench_mesh_3d.summary_text_path(
        output_dir,
        max_mesh_size=MESH_3D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_3D_MIN_MESH_ANGLE,
        domain_height=MESH_3D_DOMAIN_HEIGHT,
        quality_ratio=MESH_3D_QUALITY_RATIO,
        quality_enabled=MESH_3D_QUALITY_ENABLED,
        preserve_surface=MESH_3D_PRESERVE_SURFACE,
        max_added_points=MESH_3D_MAX_ADDED_POINTS,
        lod=MESH_3D_LOD,
        merge_buildings=MESH_3D_MERGE_BUILDINGS,
        mesher=MESH_3D_MESHER,
        pipeline_mode=MESH_3D_PIPELINE_MODE,
        stage_audit_enabled=MESH_3D_STAGE_AUDIT,
    )
    case_dir = bench_mesh_3d.case_output_dir(
        output_dir,
        TILE,
        max_mesh_size=MESH_3D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_3D_MIN_MESH_ANGLE,
        domain_height=MESH_3D_DOMAIN_HEIGHT,
        quality_ratio=MESH_3D_QUALITY_RATIO,
        quality_enabled=MESH_3D_QUALITY_ENABLED,
        preserve_surface=MESH_3D_PRESERVE_SURFACE,
        max_added_points=MESH_3D_MAX_ADDED_POINTS,
        lod=MESH_3D_LOD,
        merge_buildings=MESH_3D_MERGE_BUILDINGS,
        mesher=MESH_3D_MESHER,
        pipeline_mode=MESH_3D_PIPELINE_MODE,
        stage_audit_enabled=MESH_3D_STAGE_AUDIT,
    )
    results = load_results(results_path) if results_path.exists() else {}
    case_record = results.get(TILE)
    status, note = derive_mesh_3d_status_and_note(returncode, case_record)
    return RunResult(
        benchmark="mesh_3d",
        mode=mode,
        command=command,
        output_dir=output_dir,
        case_dir=case_dir if case_dir.exists() else None,
        results_path=results_path if results_path.exists() else None,
        summary_path=summary_path if summary_path.exists() else None,
        returncode=returncode,
        elapsed_seconds=elapsed_seconds,
        status=status,
        note=note,
        case_record=case_record,
    )


def footprints_summary_rows(runs: list[RunResult]) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for mode in MODES:
        run = next(run for run in runs if run.benchmark == "footprints" and run.mode == mode)
        record = run.case_record or {}
        timings = record.get("timing_summary", {})
        mesh = record.get("flat_mesh_quality_summary", {})
        boundary = record.get("conditioned_polygon_boundary_metrics", {})
        delta = record.get("raw_to_conditioned_difference_metrics", {})
        rows.append(
            [
                mode,
                status_badge(run.status),
                record.get("raw_polygon_count", "-"),
                record.get("conditioned_polygon_count", "-"),
                fmt(boundary.get("min_clearance"), 2),
                fmt(delta.get("symmetric_difference_area"), 2),
                fmt(mesh.get("element_quality_worst"), 3),
                fmt(mesh.get("aspect_ratio_worst"), 2),
                f"{fmt(timings.get('conditioning_seconds'), 2)}s",
                f"{fmt(timings.get('meshing_seconds'), 2)}s",
                f"{fmt(timings.get('core_seconds'), 2)}s",
            ]
        )
    return rows


def mesh_2d_summary_rows(runs: list[RunResult]) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for mode in MODES:
        run = next(run for run in runs if run.benchmark == "mesh_2d" and run.mode == mode)
        case_record = run.case_record or {}
        mesher_results = case_record.get("meshers", {})
        for mesher in bench_mesh_2d.resolve_requested_meshers(list(BENCH_2D_MESHERS)):
            result = mesher_results.get(mesher, {})
            metrics = result.get("metrics", {})
            rows.append(
                [
                    mode,
                    mesher,
                    "✅ ok" if result.get("status") == "success" else "❌ fail",
                    metrics.get("num_cells", "-"),
                    fmt(metrics.get("element_quality_min"), 4),
                    fmt(metrics.get("aspect_ratio_max"), 2),
                    fmt(metrics.get("edge_length_p01"), 2),
                    metrics.get("short_edges_lt_0_5_count", "-"),
                    f"{fmt(result.get('time'), 2)}s",
                ]
            )
    return rows


def mesh_3d_summary_rows(runs: list[RunResult]) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for mode in MODES:
        run = next(run for run in runs if run.benchmark == "mesh_3d" and run.mode == mode)
        result = (run.case_record or {}).get("result", {})
        metrics = result.get("metrics", {})
        rows.append(
            [
                mode,
                "✅ ok" if result.get("status") == "success" else "❌ fail",
                metrics.get("num_cells", "-"),
                fmt(metrics.get("element_quality_min"), 4),
                fmt(metrics.get("aspect_ratio_max"), 2),
                fmt(metrics.get("edge_length_p01"), 2),
                fmt(metrics.get("volume_p01"), 2),
                metrics.get("low_quality_lt_0_05_count", "-"),
                f"{fmt(result.get('time'), 2)}s",
            ]
        )
    return rows


def build_combined_report(runs: list[RunResult], total_elapsed: float) -> str:
    run_rows = [
        [
            run.mode,
            run.benchmark,
            status_badge(run.status),
            f"{run.elapsed_seconds:.1f}s",
            display_path(run.output_dir),
            display_path(run.case_dir),
            shorten_note(run.note),
        ]
        for run in runs
    ]

    new_runs_ok = all(run.status == "ok" for run in runs if run.mode == "new")
    overall_status = "✅ ok" if new_runs_ok else "❌ fail"
    commands = [f"{run.name}: {shlex.join(run.command)}" for run in runs]

    lines = [
        "bench_test",
        f"Tile: {TILE:03d}",
        f"Output directory: {RUN_OUTPUT_DIR}",
        f"Status: {overall_status} (based on new pipeline)",
        f"Elapsed time: {total_elapsed:.1f}s",
        "",
        section_heading("▶ COMMANDS"),
    ]
    lines.extend(commands)
    lines.extend(
        [
            "",
        format_console_table(
            ["Mode", "Benchmark", "Status", "Elapsed", "Output dir", "Case artifacts", "Note"],
            run_rows,
            title="📋 RUNS",
        ),
        "",
        format_console_table(
            [
                "Mode",
                "Status",
                "Raw",
                "Cond",
                "Clear min",
                "SymDiff",
                "EQ worst",
                "AR worst",
                "Cond s",
                "Mesh s",
                "Core s",
            ],
            footprints_summary_rows(runs),
            title="🧹 FOOTPRINTS",
        ),
        "",
        format_console_table(
            [
                "Mode",
                "Mesher",
                "Status",
                "Cells",
                "EQ min",
                "AR max",
                "Edge p01",
                "E<0.5m",
                "Time",
            ],
            mesh_2d_summary_rows(runs),
            title="🔺 MESH 2D",
        ),
        "",
        format_console_table(
            [
                "Mode",
                "Status",
                "Cells",
                "EQ min",
                "AR max",
                "Edge p01",
                "Vol p01",
                "Q<0.05",
                "Time",
            ],
            mesh_3d_summary_rows(runs),
            title="🧱 MESH 3D",
        ),
        ]
    )
    return "\n".join(lines)


def save_combined_results(runs: list[RunResult], total_elapsed: float) -> None:
    payload = {
        "tile": TILE,
        "output_dir": str(RUN_OUTPUT_DIR),
        "elapsed_seconds": total_elapsed,
        "runs": [
            {
                "mode": run.mode,
                "benchmark": run.benchmark,
                "command": run.command,
                "returncode": run.returncode,
                "elapsed_seconds": run.elapsed_seconds,
                "status": run.status,
                "note": run.note,
                "output_dir": str(run.output_dir),
                "case_dir": str(run.case_dir) if run.case_dir is not None else None,
                "results_path": str(run.results_path) if run.results_path is not None else None,
                "summary_path": str(run.summary_path) if run.summary_path is not None else None,
                "case": run.case_record,
            }
            for run in runs
        ],
    }
    COMBINED_RESULTS_PATH.write_text(
        json.dumps(json_ready(payload), indent=2) + "\n",
        encoding="utf-8",
    )


def run_mode(mode: str) -> list[RunResult]:
    footprints_command = build_footprints_command(mode)
    footprints_returncode, footprints_elapsed = run_subprocess(
        f"{mode}:bench_footprints",
        footprints_command,
    )
    footprints_run = collect_footprints_result(
        mode,
        footprints_returncode,
        footprints_elapsed,
        footprints_command,
    )

    mesh_2d_command = build_mesh_2d_command(mode)
    mesh_2d_returncode, mesh_2d_elapsed = run_subprocess(
        f"{mode}:bench_mesh_2d",
        mesh_2d_command,
    )
    mesh_2d_run = collect_mesh_2d_result(
        mode,
        mesh_2d_returncode,
        mesh_2d_elapsed,
        mesh_2d_command,
    )

    mesh_3d_command = build_mesh_3d_command(mode)
    mesh_3d_returncode, mesh_3d_elapsed = run_subprocess(
        f"{mode}:bench_mesh_3d",
        mesh_3d_command,
    )
    mesh_3d_run = collect_mesh_3d_result(
        mode,
        mesh_3d_returncode,
        mesh_3d_elapsed,
        mesh_3d_command,
    )

    return [footprints_run, mesh_2d_run, mesh_3d_run]


def main() -> int:
    RUN_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    runs: list[RunResult] = []
    for mode in MODES:
        runs.extend(run_mode(mode))

    total_elapsed = time.perf_counter() - started
    report = build_combined_report(runs, total_elapsed)
    save_combined_results(runs, total_elapsed)
    COMBINED_SUMMARY_PATH.write_text(report + "\n", encoding="utf-8")

    print(report)
    print()
    print(f"Combined results: {COMBINED_RESULTS_PATH}")
    print(f"Combined summary: {COMBINED_SUMMARY_PATH}")

    return 0 if all(run.status == "ok" for run in runs if run.mode == "new") else 1


if __name__ == "__main__":
    raise SystemExit(main())
