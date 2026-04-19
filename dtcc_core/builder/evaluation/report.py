from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional

from dtcc_core.builder.evaluation.runner import BuildingResult


_BASE_HEADERS = [
    "building_id", "stage_outcome", "roof_type_predicted", "roof_type_truth",
    "roof_type_correct", "confidence", "classifier_source", "fallback_reason_attr",
    "plane_count", "plane_count_delta", "coverage_ratio", "top2_area_share",
    "symmetry_error_deg", "watertight", "ridge_height_error", "eave_height_error",
    "total_ms",
]


def write_per_building_csv(results: Iterable[BuildingResult], path: Path | str) -> None:
    """Write one CSV row per result with flattened timing + semantic-IoU columns."""
    results = list(results)
    path_obj = Path(path)
    path_obj.parent.mkdir(parents=True, exist_ok=True)

    if not results:
        path_obj.write_text("")
        return

    iou_keys = sorted({k for r in results for k in r.semantic_iou})
    timing_keys = sorted({k for r in results for k in r.timings})

    headers = _BASE_HEADERS + [f"iou_{k}" for k in iou_keys] + [f"timing_{k}_ms" for k in timing_keys]

    with path_obj.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        for r in results:
            row = {h: getattr(r, h, None) for h in _BASE_HEADERS}
            for k in iou_keys:
                row[f"iou_{k}"] = r.semantic_iou.get(k)
            for k in timing_keys:
                row[f"timing_{k}_ms"] = r.timings.get(k)
            w.writerow(row)


def summarize(results: List[BuildingResult]) -> Dict[str, Any]:
    """Aggregate dataset-level statistics from per-building results."""
    n = len(results)
    if n == 0:
        return {"count": 0}

    stage_counts = Counter(r.stage_outcome for r in results)
    success = [r for r in results if r.stage_outcome == "success"]
    rt_acc = (
        sum(r.roof_type_correct for r in success) / len(success)
        if success else None
    )
    wt_rate = sum(1 for r in results if r.watertight) / n

    def _pct(vals: List[float], q: float) -> float:
        s = sorted(vals)
        k = max(0, min(len(s) - 1, int(round(q * (len(s) - 1)))))
        return s[k]

    total_times = [r.total_ms for r in results]
    timing = {
        "min": min(total_times),
        "median": median(total_times),
        "p95": _pct(total_times, 0.95),
        "max": max(total_times),
    }

    pc_deltas = [r.plane_count_delta for r in results if r.plane_count_delta is not None]
    coverage = [r.coverage_ratio for r in results if r.coverage_ratio is not None]
    sym = [r.symmetry_error_deg for r in results if r.symmetry_error_deg is not None]

    return {
        "count": n,
        "stage_outcome_counts": dict(stage_counts),
        "roof_type_accuracy_on_success": rt_acc,
        "watertight_rate": wt_rate,
        "total_ms": timing,
        "plane_count_delta": {
            "median": median(pc_deltas) if pc_deltas else None,
            "max": max(pc_deltas) if pc_deltas else None,
        },
        "coverage_ratio": {
            "median": median(coverage) if coverage else None,
        },
        "slope_symmetry_deg": {
            "median": median(sym) if sym else None,
        },
    }


def write_summary_json(results: List[BuildingResult], path: Path | str) -> None:
    """Write aggregate stats as JSON to the given path."""
    summary = summarize(results)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w") as f:
        json.dump(summary, f, indent=2, default=str)


@dataclass
class FailureThresholds:
    max_plane_count_delta: Optional[int] = 2
    min_coverage: Optional[float] = 0.7
    min_semantic_iou: Optional[float] = 0.7
    require_watertight: bool = True


def _is_failure(r: BuildingResult, t: FailureThresholds) -> bool:
    if r.stage_outcome != "success":
        return True
    if t.max_plane_count_delta is not None and r.plane_count_delta is not None:
        if r.plane_count_delta > t.max_plane_count_delta:
            return True
    if t.min_coverage is not None and r.coverage_ratio is not None:
        if r.coverage_ratio < t.min_coverage:
            return True
    if t.min_semantic_iou is not None and r.semantic_iou:
        if any(v < t.min_semantic_iou for v in r.semantic_iou.values()):
            return True
    if t.require_watertight and not r.watertight:
        return True
    return False


def write_failure_gallery(
    dataset_root,
    results: List[BuildingResult],
    thresholds: FailureThresholds,
    out_dir: Path | str,
) -> List[Path]:
    """For each failing building, re-run the pipeline and write a VTK file."""
    from dtcc_core.builder.evaluation.dataset import EvalDataset
    from dtcc_core.builder.evaluation.runner import Runner
    from dtcc_core.builder.evaluation.vtk_export import write_building_vtk
    from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings

    failing_ids = {r.building_id for r in results if _is_failure(r, thresholds)}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths: List[Path] = []
    if not failing_ids:
        return paths

    ds = EvalDataset(dataset_root)
    runner = Runner()
    for sample in ds:
        if sample.building_id not in failing_ids:
            continue
        building = runner._sample_to_building(sample)
        build_lod2_buildings([building])
        if building.lod2 is None:
            continue
        path = out_dir / f"{sample.building_id}.vtk"
        write_building_vtk(building, path)
        paths.append(path)
    return paths
