"""Opt-in geographic acceptance test: execute the same CLI as manual benchmarks."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from benchmarks.benchmark_catalog import build_tasks


@pytest.mark.skipif(
    os.environ.get("DTCC_RUN_QUICK_BENCHMARKS") != "1",
    reason="Set DTCC_RUN_QUICK_BENCHMARKS=1 for the real ten-tile suite",
)
def test_quick_lod2_central_tiles(tmp_path):
    root = Path(__file__).resolve().parents[2]
    output_root = Path(os.environ.get("DTCC_QUICK_BENCHMARK_OUTPUT", tmp_path))
    output_root.mkdir(parents=True, exist_ok=True)
    saved_input = os.environ.get("DTCC_QUICK_BENCHMARK_INPUT")
    payloads = {}
    errors = []
    for method in ("catalogue", "ransac"):
        output = output_root / method
        tasks = build_tasks(
            "quick", datasets=["city_surface_mesh"], lod="lod2", lod2_method=method
        )
        # Allow the runner's SIGINT/TERM/KILL grace periods as well as startup.
        timeout = sum(task.timeout_seconds + 30 for task in tasks) + 120
        command = [
            sys.executable,
            str(root / "benchmarks" / "bench"),
            "quick",
            "--dataset",
            "city_surface_mesh",
            "--lod",
            "lod2",
            "--lod2-method",
            method,
            "--save-artifacts",
            "--output",
            str(output),
        ]
        if method == "ransac":
            command += ["--phase", "meshing", "--input", str(output_root / "catalogue")]
        elif saved_input:
            command += [
                "--phase",
                "replay",
                "--input",
                str(Path(saved_input).resolve()),
            ]
        with (output_root / f"{method}.log").open("w") as log:
            result = subprocess.run(
                command, cwd=root, stdout=log, stderr=subprocess.STDOUT, timeout=timeout
            )
        if result.returncode:
            errors.append(f"{method}: CLI exited {result.returncode}; see {output}")
        path = output / "results.json"
        if not path.exists():
            errors.append(f"{method}: no results.json")
            continue
        payload = json.loads(path.read_text())
        rows = payload["results"]
        assert len(rows) == len(tasks) == 10
        payloads[method] = {r["case"]["id"]: r for r in rows}
        for row in rows:
            mesh = row.get("metrics", {}).get("meshing", {})
            if (
                mesh.get("status") not in ("success", "warning")
                or mesh.get("degenerate_cell_count", 1) != 0
                or not mesh.get("lod2", {}).get("retained_nonflat_fits", 0)
            ):
                errors.append(
                    f"{row['task_id']}: no valid nonflat LOD2 mesh; {row.get('error')}"
                )
    if set(payloads) == {"catalogue", "ransac"}:
        assert payloads["catalogue"].keys() == payloads["ransac"].keys()
        for key, left in payloads["catalogue"].items():
            right = payloads["ransac"][key]
            before = left.get("metrics", {}).get("mesh_input", {})
            after = right.get("metrics", {}).get("mesh_input", {})
            if (
                not before.get("roof_evidence")
                or after.get("source") != "saved"
                or before["roof_evidence"] != after.get("roof_evidence")
            ):
                errors.append(f"{key}: roof evidence was not replayed exactly")
    assert not errors, "\n".join(errors)
