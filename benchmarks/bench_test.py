"""
Run a quick single-tile 2D + 3D benchmark smoke test.

Edit ``TILE`` below to choose the Stockholm tile. This wrapper intentionally
takes no command-line arguments; it runs ``bench_mesh_2d.py`` and
``bench_mesh_3d.py`` with a small fixed configuration, stores all outputs under
``benchmarks/bench_test/``, and prints a combined summary.

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
    import bench_mesh_2d
    import bench_mesh_3d
    from _stockholm_common import format_console_table, json_ready
except ImportError:
    from benchmarks import bench_mesh_2d, bench_mesh_3d
    from benchmarks._stockholm_common import format_console_table, json_ready


TILE = 54

BENCHMARKS_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = BENCHMARKS_DIR / "bench_test"
RUN_OUTPUT_DIR = OUTPUT_ROOT / f"tile_{TILE:03d}"
MESH_2D_OUTPUT_DIR = RUN_OUTPUT_DIR / "mesh_2d"
MESH_3D_OUTPUT_DIR = RUN_OUTPUT_DIR / "mesh_3d"
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
    name: str
    command: list[str]
    output_dir: Path
    case_dir: Path | None
    results_path: Path | None
    summary_path: Path | None
    returncode: int
    elapsed_seconds: float
    status: str
    report: str
    selected_results: dict[int, dict[str, Any]]


def display_path(path: Path | None) -> str:
    if path is None:
        return "-"
    try:
        return str(path.relative_to(RUN_OUTPUT_DIR))
    except ValueError:
        return str(path)


def derive_status(*, returncode: int, has_results: bool, complete: bool) -> str:
    if returncode != 0:
        return f"failed({returncode})"
    if complete:
        return "ok"
    if has_results:
        return "partial"
    return "no-results"


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


def build_mesh_2d_command() -> list[str]:
    command = [
        sys.executable,
        "-u",
        str(BENCHMARKS_DIR / "bench_mesh_2d.py"),
        "--cases",
        str(TILE),
        "--delay",
        "0",
        "--max-mesh-size",
        "none" if MESH_2D_MAX_MESH_SIZE is None else str(MESH_2D_MAX_MESH_SIZE),
        "--output-dir",
        str(MESH_2D_OUTPUT_DIR),
    ]
    if MESH_2D_STAGE_AUDIT:
        command.append("--stage-audit")
    return command


def build_mesh_3d_command() -> list[str]:
    command = [
        sys.executable,
        "-u",
        str(BENCHMARKS_DIR / "bench_mesh_3d.py"),
        "--cases",
        str(TILE),
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
        str(MESH_3D_OUTPUT_DIR),
    ]
    if MESH_3D_MERGE_BUILDINGS:
        command.append("--merge-buildings")
    else:
        command.append("--no-merge-buildings")
    if MESH_3D_QUALITY_ENABLED:
        pass
    else:
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


def collect_mesh_2d_result(returncode: int, elapsed_seconds: float, command: list[str]) -> RunResult:
    meshers = bench_mesh_2d.resolve_requested_meshers(None)
    results_path = bench_mesh_2d.results_file_path(
        MESH_2D_OUTPUT_DIR,
        meshers,
        max_mesh_size=MESH_2D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
        pipeline_mode=MESH_2D_PIPELINE_MODE,
        stage_audit_enabled=MESH_2D_STAGE_AUDIT,
    )
    summary_path = bench_mesh_2d.summary_text_path(
        MESH_2D_OUTPUT_DIR,
        meshers,
        max_mesh_size=MESH_2D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
        pipeline_mode=MESH_2D_PIPELINE_MODE,
        stage_audit_enabled=MESH_2D_STAGE_AUDIT,
    )
    case_dir = bench_mesh_2d.case_output_dir(
        MESH_2D_OUTPUT_DIR,
        TILE,
        meshers,
        max_mesh_size=MESH_2D_MAX_MESH_SIZE,
        min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
        pipeline_mode=MESH_2D_PIPELINE_MODE,
        stage_audit_enabled=MESH_2D_STAGE_AUDIT,
    )

    results = bench_mesh_2d.load_results(results_path) if results_path.exists() else {}
    selected_results = {TILE: results[TILE]} if TILE in results else {}
    complete = False
    if selected_results:
        complete = bench_mesh_2d.case_complete(selected_results[TILE], meshers)
        report = bench_mesh_2d.build_summary_report(
            selected_results,
            meshers,
            max_mesh_size=MESH_2D_MAX_MESH_SIZE,
            min_mesh_angle=MESH_2D_MIN_MESH_ANGLE,
            pipeline_mode=MESH_2D_PIPELINE_MODE,
            stage_audit_enabled=MESH_2D_STAGE_AUDIT,
            elapsed_seconds=elapsed_seconds,
        )
    else:
        report = "bench_mesh_2d\nNo result record was written for the requested tile."

    status = derive_status(
        returncode=returncode,
        has_results=bool(selected_results),
        complete=complete,
    )
    return RunResult(
        name="bench_mesh_2d",
        command=command,
        output_dir=MESH_2D_OUTPUT_DIR,
        case_dir=case_dir,
        results_path=results_path,
        summary_path=summary_path,
        returncode=returncode,
        elapsed_seconds=elapsed_seconds,
        status=status,
        report=report,
        selected_results=selected_results,
    )


def collect_mesh_3d_result(returncode: int, elapsed_seconds: float, command: list[str]) -> RunResult:
    results_path = bench_mesh_3d.results_file_path(
        MESH_3D_OUTPUT_DIR,
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
        MESH_3D_OUTPUT_DIR,
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
        MESH_3D_OUTPUT_DIR,
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

    results = bench_mesh_3d.load_results(results_path) if results_path.exists() else {}
    selected_results = {TILE: results[TILE]} if TILE in results else {}
    complete = False
    if selected_results:
        complete = bench_mesh_3d.case_complete(selected_results[TILE])
        report = bench_mesh_3d.build_summary_report(
            selected_results,
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
            elapsed_seconds=elapsed_seconds,
        )
    else:
        report = "bench_mesh_3d\nNo result record was written for the requested tile."

    status = derive_status(
        returncode=returncode,
        has_results=bool(selected_results),
        complete=complete,
    )
    return RunResult(
        name="bench_mesh_3d",
        command=command,
        output_dir=MESH_3D_OUTPUT_DIR,
        case_dir=case_dir,
        results_path=results_path,
        summary_path=summary_path,
        returncode=returncode,
        elapsed_seconds=elapsed_seconds,
        status=status,
        report=report,
        selected_results=selected_results,
    )


def build_combined_report(runs: list[RunResult], total_elapsed: float) -> str:
    overview_rows = [
        [
            run.name,
            run.status,
            f"{run.elapsed_seconds:.1f}s",
            display_path(run.output_dir),
            display_path(run.case_dir),
        ]
        for run in runs
    ]

    lines = [
        "bench_test",
        f"Tile: {TILE:03d}",
        f"Output directory: {RUN_OUTPUT_DIR}",
        f"Elapsed time: {total_elapsed:.1f}s",
        "",
        format_console_table(
            ["Benchmark", "Status", "Elapsed", "Output dir", "Case artifacts"],
            overview_rows,
            title="Runs",
        ),
        "",
        "Commands",
    ]
    lines.extend(f"{run.name}: {shlex.join(run.command)}" for run in runs)

    for run in runs:
        lines.extend(
            [
                "",
                "=" * 80,
                run.name,
                "=" * 80,
                run.report.strip(),
                "",
                f"Output dir: {run.output_dir}",
                f"Case artifacts: {run.case_dir if run.case_dir is not None else '-'}",
                f"Results file: {run.results_path if run.results_path is not None else '-'}",
                f"Summary text: {run.summary_path if run.summary_path is not None else '-'}",
            ]
        )

    return "\n".join(lines)


def save_combined_results(runs: list[RunResult], total_elapsed: float) -> None:
    payload = {
        "tile": TILE,
        "output_dir": str(RUN_OUTPUT_DIR),
        "elapsed_seconds": total_elapsed,
        "benchmarks": {
            run.name: {
                "command": run.command,
                "returncode": run.returncode,
                "elapsed_seconds": run.elapsed_seconds,
                "status": run.status,
                "output_dir": str(run.output_dir),
                "case_dir": str(run.case_dir) if run.case_dir is not None else None,
                "results_path": str(run.results_path) if run.results_path is not None else None,
                "summary_path": str(run.summary_path) if run.summary_path is not None else None,
                "cases": run.selected_results,
            }
            for run in runs
        },
    }
    COMBINED_RESULTS_PATH.write_text(
        json.dumps(json_ready(payload), indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    RUN_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()

    mesh_2d_command = build_mesh_2d_command()
    mesh_2d_returncode, mesh_2d_elapsed = run_subprocess("bench_mesh_2d", mesh_2d_command)
    mesh_2d_run = collect_mesh_2d_result(
        mesh_2d_returncode,
        mesh_2d_elapsed,
        mesh_2d_command,
    )

    mesh_3d_command = build_mesh_3d_command()
    mesh_3d_returncode, mesh_3d_elapsed = run_subprocess("bench_mesh_3d", mesh_3d_command)
    mesh_3d_run = collect_mesh_3d_result(
        mesh_3d_returncode,
        mesh_3d_elapsed,
        mesh_3d_command,
    )

    runs = [mesh_2d_run, mesh_3d_run]
    total_elapsed = time.perf_counter() - started

    report = build_combined_report(runs, total_elapsed)
    save_combined_results(runs, total_elapsed)
    COMBINED_SUMMARY_PATH.write_text(report + "\n", encoding="utf-8")

    print(report)
    print()
    print(f"Combined results: {COMBINED_RESULTS_PATH}")
    print(f"Combined summary: {COMBINED_SUMMARY_PATH}")

    return 0 if all(run.status == "ok" for run in runs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
