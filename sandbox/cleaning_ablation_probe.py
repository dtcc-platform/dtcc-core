"""Run a small performance/quality ablation study for the new cleaner.

This probe focuses on a few representative slow/important Stockholm tiles and
compares the current cleaner against several principled stage ablations. It is
meant to answer two questions:

1. Which stages dominate runtime?
2. Which stages are actually buying geometry / mesh quality?
"""

from __future__ import annotations

import argparse
import json
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.bench_footprints import (
    build_mesh_from_conditioned_footprints,
    case_to_grid,
    coverage_difference_metrics,
    extract_raw_footprints,
    make_bounds,
    mesh_quality_summary,
    prepare_city,
)

import dtcc_core.builder.cleaning as cleaning
from dtcc_core.model import GeometryType


OUTPUT_ROOT = (
    REPO_ROOT / "benchmarks" / "output_footprints" / "cleaning_ablation_probe"
)


@dataclass(frozen=True)
class Variant:
    name: str
    description: str
    patches: tuple[tuple[str, Callable[..., Any]], ...]


def _identity_stage(
    polygons,
    source_map,
    **_: Any,
):
    return list(polygons), [list(indices) for indices in source_map]


def _no_local_coverage(*args: Any, **kwargs: Any):
    return None


VARIANTS: tuple[Variant, ...] = (
    Variant(
        name="baseline",
        description="Current implementation.",
        patches=(),
    ),
    Variant(
        name="no_local_coverage",
        description="Disable local shared-boundary coverage simplification branch.",
        patches=(("_simplify_coverage_locally", _no_local_coverage),),
    ),
    Variant(
        name="no_coverage_meshing_regularization",
        description="Skip final coverage meshing regularization.",
        patches=(("_regularize_coverage_for_meshing", _identity_stage),),
    ),
    Variant(
        name="no_source_recovery",
        description="Skip source-coordinate recovery.",
        patches=(("_recover_source_supported_coordinates", _identity_stage),),
    ),
    Variant(
        name="post_refinement_off",
        description="Skip all post-coverage refinement stages after branch selection.",
        patches=(
            ("_reclaim_source_supported_area", _identity_stage),
            ("_absorb_small_supported_components", _identity_stage),
            ("_simplify_polygons_for_meshing", _identity_stage),
            ("_regularize_low_clearance_polygons", _identity_stage),
            ("_recover_source_supported_coordinates", _identity_stage),
            ("_regularize_coverage_for_meshing", _identity_stage),
        ),
    ),
)


@contextmanager
def patched_cleaning(variant: Variant):
    originals: list[tuple[str, Any]] = []
    try:
        for name, replacement in variant.patches:
            originals.append((name, getattr(cleaning.footprints, name)))
            setattr(cleaning.footprints, name, replacement)
        yield
    finally:
        for name, original in reversed(originals):
            setattr(cleaning.footprints, name, original)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        nargs="+",
        type=int,
        default=[54, 55, 63],
        help="Representative Stockholm tiles to probe.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
    )
    parser.add_argument(
        "--max-mesh-size",
        type=float,
        default=10.0,
    )
    parser.add_argument(
        "--min-mesh-angle",
        type=float,
        default=25.0,
    )
    parser.add_argument(
        "--merge-tolerance",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--min-building-detail",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--min-building-area",
        type=float,
        default=15.0,
    )
    parser.add_argument(
        "--raster-cell-size",
        type=float,
        default=2.0,
    )
    parser.add_argument(
        "--raster-radius",
        type=float,
        default=3.0,
    )
    return parser.parse_args()


def _mean(values: list[float | None]) -> float | None:
    filtered = [value for value in values if value is not None]
    if not filtered:
        return None
    return sum(filtered) / len(filtered)


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    case_inputs: dict[int, dict[str, Any]] = {}
    for case in args.cases:
        ix, iy = case_to_grid(case)
        bounds = make_bounds(ix, iy)
        terrain_raster, buildings, _ = prepare_city(
            bounds,
            raster_cell_size=args.raster_cell_size,
            raster_radius=args.raster_radius,
        )
        raw_polygons, raw_source_map = extract_raw_footprints(buildings)
        case_inputs[case] = {
            "bounds": bounds,
            "terrain_raster": terrain_raster,
            "buildings": buildings,
            "raw_polygons": raw_polygons,
            "raw_source_map": raw_source_map,
        }

    results: list[dict[str, Any]] = []
    for case in args.cases:
        case_input = case_inputs[case]
        for variant in VARIANTS:
            with patched_cleaning(variant):
                t0 = time.perf_counter()
                conditioning_result = cleaning.condition_building_footprints(
                    case_input["buildings"],
                    lod=GeometryType.LOD0,
                    options=cleaning.ConditioningOptions(
                        precision_grid=None,
                        min_feature_size=args.min_building_detail,
                        merge_distance=args.merge_tolerance,
                        min_area=args.min_building_area,
                        min_hole_area=args.min_building_detail**2,
                        collect_stage_metrics=False,
                        enable_logging=False,
                    ),
                )
                conditioning_seconds = time.perf_counter() - t0

                t0 = time.perf_counter()
                mesh = build_mesh_from_conditioned_footprints(
                    case_input["terrain_raster"],
                    conditioning_result.polygons,
                    conditioning_result.source_map,
                    case_input["buildings"],
                    max_mesh_size=args.max_mesh_size,
                    min_mesh_angle=args.min_mesh_angle,
                    disable_cleaning_diagnostics=True,
                )
                meshing_seconds = time.perf_counter() - t0

            difference = coverage_difference_metrics(
                case_input["raw_polygons"],
                conditioning_result.polygons,
            )
            mesh_metrics = mesh_quality_summary(mesh.quality())

            results.append(
                {
                    "case": case,
                    "variant": variant.name,
                    "description": variant.description,
                    "conditioning_seconds": conditioning_seconds,
                    "meshing_seconds": meshing_seconds,
                    "core_seconds": conditioning_seconds + meshing_seconds,
                    "conditioned_polygon_count": len(conditioning_result.polygons),
                    "symdiff_area": difference["symmetric_difference_area"],
                    "missing_area": difference["reference_minus_candidate_area"],
                    "extra_area": difference["candidate_minus_reference_area"],
                    "eq_worst": mesh_metrics["element_quality_worst"],
                    "eq_mean": mesh_metrics["element_quality_mean"],
                    "ar_worst": mesh_metrics["aspect_ratio_worst"],
                    "er_worst": mesh_metrics["edge_ratio_worst"],
                    "skew_worst": mesh_metrics["skewness_worst"],
                }
            )

    by_variant: dict[str, list[dict[str, Any]]] = {}
    for result in results:
        by_variant.setdefault(result["variant"], []).append(result)

    summary_rows: list[dict[str, Any]] = []
    baseline_rows = {row["case"]: row for row in by_variant["baseline"]}
    for variant in VARIANTS:
        rows = by_variant[variant.name]
        conditioning_mean = _mean([row["conditioning_seconds"] for row in rows])
        core_mean = _mean([row["core_seconds"] for row in rows])
        summary_rows.append(
            {
                "variant": variant.name,
                "description": variant.description,
                "conditioning_seconds_mean": conditioning_mean,
                "core_seconds_mean": core_mean,
                "symdiff_area_mean": _mean([row["symdiff_area"] for row in rows]),
                "missing_area_mean": _mean([row["missing_area"] for row in rows]),
                "extra_area_mean": _mean([row["extra_area"] for row in rows]),
                "eq_worst_mean": _mean([row["eq_worst"] for row in rows]),
                "ar_worst_mean": _mean([row["ar_worst"] for row in rows]),
                "er_worst_mean": _mean([row["er_worst"] for row in rows]),
                "conditioning_speedup_vs_baseline": (
                    baseline_rows and conditioning_mean is not None
                    and _mean(
                        [baseline_rows[row["case"]]["conditioning_seconds"] for row in rows]
                    )
                    / conditioning_mean
                    if variant.name != "baseline"
                    else 1.0
                ),
            }
        )

    payload = {
        "cases": args.cases,
        "results": results,
        "summary": summary_rows,
    }
    json_path = args.output_root / "cleaning_ablation_probe.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    md_lines = [
        "# Cleaning Ablation Probe",
        "",
        f"Cases: {', '.join(str(case) for case in args.cases)}",
        "",
        "| Variant | Cond Mean (s) | Core Mean (s) | Speedup vs Baseline | SymDiff Mean | Missing Mean | Extra Mean | EQ Worst Mean | AR Worst Mean | ER Worst Mean |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        speedup = row["conditioning_speedup_vs_baseline"]
        md_lines.append(
            "| {variant} | {conditioning:.3f} | {core:.3f} | {speedup:.2f} | {symdiff:.1f} | {missing:.1f} | {extra:.1f} | {eq:.3f} | {ar:.3f} | {er:.3f} |".format(
                variant=row["variant"],
                conditioning=row["conditioning_seconds_mean"] or 0.0,
                core=row["core_seconds_mean"] or 0.0,
                speedup=speedup or 0.0,
                symdiff=row["symdiff_area_mean"] or 0.0,
                missing=row["missing_area_mean"] or 0.0,
                extra=row["extra_area_mean"] or 0.0,
                eq=row["eq_worst_mean"] or 0.0,
                ar=row["ar_worst_mean"] or 0.0,
                er=row["er_worst_mean"] or 0.0,
            )
        )
    md_lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Lower is better for runtime, `symdiff`, `missing`, `extra`, `AR worst`, and `ER worst`.",
            "- Higher is better for `EQ worst`.",
            "- This probe isolates a few representative slow/important tiles instead of the full grid.",
        ]
    )
    md_path = args.output_root / "cleaning_ablation_probe.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
