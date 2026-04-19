#!/usr/bin/env python3
"""Evaluate the rule-based LoD2 pipeline against a labeled dataset.

Usage:
    python sandbox/evaluate_lod2.py --dataset path/to/dataset --out path/to/output
"""
import argparse
import sys
from pathlib import Path

from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.runner import Runner
from dtcc_core.builder.evaluation.report import (
    FailureThresholds,
    summarize,
    write_failure_gallery,
    write_per_building_csv,
    write_summary_json,
)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--max-plane-count-delta", type=int, default=2)
    p.add_argument("--min-coverage", type=float, default=0.7)
    p.add_argument("--min-semantic-iou", type=float, default=0.7)
    p.add_argument("--no-watertight-required", action="store_true")
    args = p.parse_args(argv)

    ds = EvalDataset(args.dataset)
    runner = Runner()
    results = runner.run(ds)

    args.out.mkdir(parents=True, exist_ok=True)
    csv_path = args.out / "per_building.csv"
    summary_path = args.out / "summary.json"
    gallery_dir = args.out / "failures"

    write_per_building_csv(results, csv_path)
    write_summary_json(results, summary_path)
    thresholds = FailureThresholds(
        max_plane_count_delta=args.max_plane_count_delta,
        min_coverage=args.min_coverage,
        min_semantic_iou=args.min_semantic_iou,
        require_watertight=not args.no_watertight_required,
    )
    gallery_paths = write_failure_gallery(
        dataset_root=args.dataset,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )

    summary = summarize(results)
    print(f"Evaluated {summary['count']} buildings")
    for reason, count in sorted(summary["stage_outcome_counts"].items()):
        print(f"  {reason}: {count}")
    if summary.get("roof_type_accuracy_on_success") is not None:
        print(
            f"Roof-type accuracy on success subset: "
            f"{summary['roof_type_accuracy_on_success']:.3f}"
        )
    print(f"Watertight rate: {summary['watertight_rate']:.3f}")
    print(f"Wrote: {csv_path}, {summary_path}")
    print(f"Failure gallery: {len(gallery_paths)} VTK files under {gallery_dir}")


if __name__ == "__main__":
    main(sys.argv[1:])
