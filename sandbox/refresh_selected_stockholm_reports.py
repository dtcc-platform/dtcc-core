from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
from PIL import Image, ImageDraw, ImageFont


OUTPUT_ROOT = (
    Path(__file__).resolve().parent / "output" / "stockholm_flat_mesh_compare"
)
LEGACY_ROOT = OUTPUT_ROOT / "selected-old-legacy"
NEW_ROOT = OUTPUT_ROOT / "selected-new-final"
TABLE_PATH = OUTPUT_ROOT / "worst_case_comparison_table.md"

SHORT_EDGE_THRESHOLD = 0.5
OLD_PAD_BOTTOM = 180
NEW_PAD_BOTTOM = 230
PANEL_MARGIN = 18
FONT = ImageFont.load_default()


def polygon_boundary_metrics(
    polygons: list,
    *,
    short_edge_threshold: float,
) -> dict[str, float | int | None]:
    if not polygons:
        return {
            "polygon_count": 0,
            "vertex_count": 0,
            "mean_edge_length": None,
            "min_edge_length": None,
            "short_edge_count": 0,
            "min_clearance": None,
        }

    vertex_count = 0
    edge_lengths: list[float] = []
    short_edge_count = 0
    clearances: list[float] = []

    for polygon in polygons:
        rings = [polygon.exterior, *polygon.interiors]
        for ring in rings:
            coords = list(ring.coords)
            vertex_count += len(coords) - 1
            for start, end in zip(coords, coords[1:]):
                length = float(np.hypot(end[0] - start[0], end[1] - start[1]))
                edge_lengths.append(length)
                if short_edge_threshold > 0 and length + 1e-12 < short_edge_threshold:
                    short_edge_count += 1
        try:
            clearance = polygon.minimum_clearance
        except Exception:
            continue
        if np.isfinite(clearance):
            clearances.append(float(clearance))

    return {
        "polygon_count": len(polygons),
        "vertex_count": vertex_count,
        "mean_edge_length": float(np.mean(edge_lengths)) if edge_lengths else None,
        "min_edge_length": float(np.min(edge_lengths)) if edge_lengths else None,
        "short_edge_count": short_edge_count,
        "min_clearance": float(np.min(clearances)) if clearances else None,
    }


def timing_summary(timings_seconds: dict[str, float]) -> dict[str, float]:
    conditioning = float(timings_seconds.get("conditioning", 0.0))
    meshing = float(timings_seconds.get("meshing", 0.0))
    return {
        "conditioning_seconds": conditioning,
        "meshing_seconds": meshing,
        "core_seconds": conditioning + meshing,
    }


def fmt(value: float | int | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def polygon_text(
    metrics: dict[str, float | int | None],
    *,
    delta_metrics: dict[str, float] | None = None,
) -> str:
    lines = [
        f"polys {metrics['polygon_count']}  verts {metrics['vertex_count']}",
        f"edge mean {fmt(metrics['mean_edge_length'])}  min {fmt(metrics['min_edge_length'])}",
        f"short < {SHORT_EDGE_THRESHOLD:.2f}m: {metrics['short_edge_count']}  clear {fmt(metrics['min_clearance'])}",
    ]
    if delta_metrics is not None:
        lines.extend(
            [
                f"symdiff {fmt(delta_metrics['symmetric_difference_area'])} m2",
                f"missing {fmt(delta_metrics['reference_minus_candidate_area'])}  extra {fmt(delta_metrics['candidate_minus_reference_area'])}",
            ]
        )
    return "\n".join(lines)


def mesh_text(summary: dict[str, object]) -> str:
    mesh_quality = summary["flat_mesh_quality"]
    timings = summary.get("timing_summary") or timing_summary(summary["timings_seconds"])
    return "\n".join(
        [
            f"EQ mean {mesh_quality['element_quality']['mean']:.3f}  worst {mesh_quality['element_quality']['min']:.3f}",
            f"AR mean {mesh_quality['aspect_ratio']['mean']:.3f}  worst {mesh_quality['aspect_ratio']['max']:.3f}",
            f"ER mean {mesh_quality['edge_ratio']['mean']:.3f}  worst {mesh_quality['edge_ratio']['max']:.3f}",
            f"Skew mean {mesh_quality['skewness']['mean']:.3f}  worst {mesh_quality['skewness']['max']:.3f}",
            f"Cond {timings['conditioning_seconds']:.3f}s  Mesh {timings['meshing_seconds']:.3f}s",
            f"Core {timings['core_seconds']:.3f}s",
        ]
    )


def update_case_image(case_dir: Path) -> None:
    summary_path = case_dir / "summary.json"
    summary = json.loads(summary_path.read_text())
    if "timing_summary" not in summary:
        summary["timing_summary"] = timing_summary(summary["timings_seconds"])
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    raw = list(gpd.read_file(case_dir / "footprints.gpkg", layer="raw").geometry)
    conditioned = list(
        gpd.read_file(case_dir / "footprints.gpkg", layer="conditioned").geometry
    )

    raw_metrics = polygon_boundary_metrics(
        raw,
        short_edge_threshold=SHORT_EDGE_THRESHOLD,
    )
    conditioned_metrics = polygon_boundary_metrics(
        conditioned,
        short_edge_threshold=SHORT_EDGE_THRESHOLD,
    )
    delta_metrics = summary["raw_to_conditioned_difference_metrics"]

    image_path = case_dir / "comparison.png"
    image = Image.open(image_path).convert("RGB")
    width, height = image.size

    # Current selected reports already have a metrics band appended.
    if height > 1250:
        image = image.crop((0, 0, width, height - OLD_PAD_BOTTOM))
        width, height = image.size

    canvas = Image.new("RGB", (width, height + NEW_PAD_BOTTOM), "white")
    canvas.paste(image, (0, 0))
    draw = ImageDraw.Draw(canvas)

    panel_width = width // 3
    boxes = [
        polygon_text(raw_metrics),
        polygon_text(conditioned_metrics, delta_metrics=delta_metrics),
        mesh_text(summary),
    ]

    for index, text in enumerate(boxes):
        x0 = index * panel_width + PANEL_MARGIN
        y0 = height + 14
        x1 = (index + 1) * panel_width - PANEL_MARGIN
        y1 = height + NEW_PAD_BOTTOM - 14
        draw.rounded_rectangle(
            (x0, y0, x1, y1),
            radius=10,
            fill="white",
            outline="black",
            width=1,
        )
        draw.multiline_text(
            (x0 + 12, y0 + 10),
            text,
            fill="black",
            font=FONT,
            spacing=4,
        )

    canvas.save(image_path)


def table_row(case: str, legacy_summary: dict, new_summary: dict) -> str:
    legacy_delta = legacy_summary["raw_to_conditioned_difference_metrics"]
    new_delta = new_summary["raw_to_conditioned_difference_metrics"]
    legacy_mesh = legacy_summary["flat_mesh_quality"]
    new_mesh = new_summary["flat_mesh_quality"]
    legacy_timing = legacy_summary.get("timing_summary") or timing_summary(
        legacy_summary["timings_seconds"]
    )
    new_timing = new_summary.get("timing_summary") or timing_summary(
        new_summary["timings_seconds"]
    )

    return (
        "| {case} | {ls:.1f} | {ns:.1f} | {lm:.1f} | {nm:.1f} | "
        "{leqmin:.3f} | {neqmin:.3f} | {leqmean:.3f} | {neqmean:.3f} | "
        "{lar:.3f} | {nar:.3f} | {ler:.3f} | {ner:.3f} | "
        "{lcond:.3f} | {ncond:.3f} | {lmesh:.3f} | {nmesh:.3f} | {lcore:.3f} | {ncore:.3f} |"
    ).format(
        case=case,
        ls=legacy_delta["symmetric_difference_area"],
        ns=new_delta["symmetric_difference_area"],
        lm=legacy_delta["reference_minus_candidate_area"],
        nm=new_delta["reference_minus_candidate_area"],
        leqmin=legacy_mesh["element_quality"]["min"],
        neqmin=new_mesh["element_quality"]["min"],
        leqmean=legacy_mesh["element_quality"]["mean"],
        neqmean=new_mesh["element_quality"]["mean"],
        lar=legacy_mesh["aspect_ratio"]["max"],
        nar=new_mesh["aspect_ratio"]["max"],
        ler=legacy_mesh["edge_ratio"]["max"],
        ner=new_mesh["edge_ratio"]["max"],
        lcond=legacy_timing["conditioning_seconds"],
        ncond=new_timing["conditioning_seconds"],
        lmesh=legacy_timing["meshing_seconds"],
        nmesh=new_timing["meshing_seconds"],
        lcore=legacy_timing["core_seconds"],
        ncore=new_timing["core_seconds"],
    )


def build_table() -> str:
    lines = [
        "| Case | Poly SymDiff Old | Poly SymDiff New | Missing Area Old | Missing Area New | EQ Worst Old | EQ Worst New | EQ Mean Old | EQ Mean New | AR Worst Old | AR Worst New | ER Worst Old | ER Worst New | Cond Old (s) | Cond New (s) | Mesh Old (s) | Mesh New (s) | Core Old (s) | Core New (s) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    cases = sorted(case_dir.name for case_dir in LEGACY_ROOT.glob("case_*"))
    for case_name in cases:
        legacy_summary = json.loads((LEGACY_ROOT / case_name / "summary.json").read_text())
        new_summary = json.loads((NEW_ROOT / case_name / "summary.json").read_text())
        lines.append(table_row(case_name.split("_")[1], legacy_summary, new_summary))

    return "\n".join(lines) + "\n"


def main() -> int:
    for root in (LEGACY_ROOT, NEW_ROOT):
        for case_dir in sorted(root.glob("case_*")):
            update_case_image(case_dir)

    TABLE_PATH.write_text(build_table(), encoding="utf-8")
    print(TABLE_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
