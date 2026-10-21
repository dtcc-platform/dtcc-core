"""Compare LOD2 geometry, quality and cost against independent saved baselines.

Run with ``python -m benchmarks.lod2_comparison CURRENT --production V7
--original-catalogue ORIGINAL_CAT --original-ransac ORIGINAL_RANSAC``.
No benchmark execution or source geometry is modified by this reporting tool.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import meshio
import numpy as np
import shapely

from dtcc_core.builder.meshing.compact_surface import quality_by_kind
from dtcc_core.model import Field, Mesh


OREBRO_TARGETS = {
    "catalogue": dict(q_p01=0.0928, below_q_0_02=1),
    "ransac": dict(q_p01=0.0595, below_q_0_02=4),
}


def with_fractions(metrics):
    """Retain absolute defect counts as well as fractions of each face kind."""
    result = dict(metrics)
    count = metrics["triangles"]
    for key, value in metrics.items():
        if key.startswith("below_"):
            result["fraction_" + key] = value / count if count else 0.0
    return result


def measured_quality(mesh, kinds, mask=None):
    """Measure actual exported 3D triangles with the production quality metric."""
    selected = np.ones(len(mesh.faces), dtype=bool) if mask is None else mask
    part = Mesh(
        vertices=mesh.vertices,
        faces=mesh.faces[selected],
        markers=np.asarray(mesh.markers)[selected],
    )
    part.fields = [
        Field(
            name="catalogue_surface_face",
            association="face",
            values=np.arange(len(part.faces)),
        )
    ]
    return {
        kind: with_fractions(row)
        for kind, row in quality_by_kind(part, np.asarray(kinds)[selected]).items()
    }


def read_surface(path, catalogue):
    """Load VTU markers explicitly; a geometry-only loader drops this provenance."""
    raw = meshio.read(path)
    mesh = Mesh(
        vertices=raw.points,
        faces=raw.cells_dict["triangle"],
        markers=raw.cell_data_dict["markers"]["triangle"],
    )
    markers = np.asarray(mesh.markers, dtype=int)
    patches = catalogue["patch_markers"]
    n = len(patches)
    kinds = np.where(markers < 0, "terrain", np.where(markers >= n, "roof", "wall"))
    owners = np.array(
        [
            "__terrain__" if marker < 0 else patches[marker % n]["source_key"]
            for marker in markers
        ]
    )
    return mesh, kinds, owners


def load_campaign(path, method=None):
    """Read completed records through the existing benchmark output schema."""
    path = Path(path).resolve()
    records = {}
    for roof_method in [method] if method else ["catalogue", "ransac"]:
        folder = path if method else path / roof_method
        payload = json.loads((folder / "results.json").read_text())
        for row in payload["results"]:
            metrics = row.get("metrics", {})
            catalogue = metrics.get(
                "stage_audit", metrics.get("meshing", {}).get("stage_audit", {})
            ).get("roof_catalogue", {})
            artifact = row.get("artifacts", {}).get("mesh", {}).get("path")
            target = Path(artifact) if artifact else None
            if target is not None and not target.is_absolute():
                target = folder / target
            records[(row["case"]["city"], roof_method)] = dict(
                catalogue=catalogue,
                meshing=metrics.get("meshing", {}),
                source_evidence=metrics.get("mesh_input", {}).get("roof_evidence"),
                elapsed_seconds=row["elapsed_seconds"],
                status=row["status"],
                cleaning_seconds=metrics.get("cleaning", {}).get("seconds"),
                artifact=target if target is not None and target.is_file() else None,
            )
    return records


def roof_signature(catalogue, source, xy):
    rows = []
    for p in catalogue["patch_markers"]:
        if p["source_key"] != source:
            continue
        coefficients = np.asarray(p["coefficients"])
        z = float((xy - p["origin"]) @ coefficients[:2] + coefficients[2])
        rows.append(
            (tuple(np.round(coefficients[:2], 10)), round(z, 7), p["fit_retained"])
        )
    return sorted(set(rows))


def source_change(before, after, source):
    """Separate changed roof retention/planes from changes of source support."""
    bm, bk, bo, bc = before
    am, ak, ao, ac = after
    selections = [(bm, bk, bo), (am, ak, ao)]
    domains = []
    for mesh, kinds, owners in selections:
        xy = mesh.vertices[mesh.faces[(kinds == "roof") & (owners == source)], :2]
        domains.append(
            shapely.union_all(shapely.polygons(xy)) if len(xy) else shapely.Polygon()
        )
    xy = (
        np.floor(np.asarray(domains[1].centroid.coords[0]))
        if not domains[1].is_empty
        else np.zeros(2)
    )
    planes_same = roof_signature(bc, source, xy) == roof_signature(ac, source, xy)
    area = float(domains[0].symmetric_difference(domains[1]).area)
    boundary = (
        float(domains[0].boundary.hausdorff_distance(domains[1].boundary))
        if all(not p.is_empty for p in domains)
        else None
    )
    budget = 1e-7 * sum(p.length for p in domains)
    return dict(
        same_roof_planes_and_retention=planes_same,
        support_area_change_m2=area,
        support_boundary_distance_m=boundary,
        classification="roof_model_or_retention_changed"
        if not planes_same
        else "source_support_changed"
        if area > budget
        else "unchanged_roof_geometry",
    )


def junction_mask(mesh, xy, radius=2.0):
    """Select every triangle intersecting a fixed XY neighborhood, including walls."""
    points = mesh.vertices[mesh.faces, :2]
    low = np.asarray(xy) - radius
    high = np.asarray(xy) + radius
    near = np.all(points.max(axis=1) >= low, axis=1) & np.all(
        points.min(axis=1) <= high, axis=1
    )
    ids = np.flatnonzero(near)
    if len(ids):
        geometry = shapely.convex_hull(shapely.multipoints(points[ids]))
        near[ids] = shapely.intersects(geometry, shapely.box(*low, *high))
    return near


def quality_flags(before, after):
    """Historical regressions are independent of the current minimal mesh guards."""
    if not before["triangles"] or not after["triangles"]:
        return []
    flags = []
    for key in ("minimum_q", "q_p01", "q_p05"):
        if after[key] < 0.9 * before[key]:
            flags.append(key + " decreased more than 10%")
    for threshold in (0.02, 0.1, 0.25):
        key = f"below_q_{threshold:g}"
        if after[key] > before[key]:
            flags.append(key + " count increased")
    return flags


def orebro_acceptance(quality, method):
    target = OREBRO_TARGETS[method]
    return dict(
        target=target,
        p01_met=quality["q_p01"] >= target["q_p01"],
        severe_count_met=quality["below_q_0.02"] <= target["below_q_0_02"],
    )


def compare_campaigns(current, production, originals, minimal=None):
    """Build three independent comparisons plus matched source/junction records."""
    cases = []
    quality = []
    sources = []
    junctions = []
    for (city, method), record in sorted(current.items()):
        cat = record["catalogue"]
        compact = cat.get("compact_surface", {})
        old_records = {
            "current": record,
            "v7": production[(city, method)],
            "original": originals[(city, method)],
        }
        if minimal and (city, method) in minimal:
            old_records["minimal"] = minimal[(city, method)]
        loaded = {}
        by_stage = {}
        for stage, entry in old_records.items():
            if entry["artifact"] is None:
                continue
            mesh, kinds, owners = read_surface(entry["artifact"], entry["catalogue"])
            loaded[stage] = (mesh, kinds, owners, entry["catalogue"])
            by_stage[stage] = measured_quality(mesh, kinds)
        if "current" not in loaded:
            cases.append(
                dict(city=city, method=method, valid=False, status=record["status"])
            )
            continue
        if "minimal" not in by_stage:
            for kind, metrics in compact["quality_before"].items():
                by_stage.setdefault("minimal", {})[kind] = with_fractions(metrics)
        for stage, measured in by_stage.items():
            for kind, metrics in measured.items():
                interior = (
                    sum(
                        p["kind"] == "interior"
                        and compact["faces"][p["proposer"]]["kind"] == kind
                        for p in compact.get("optional_vertices", [])
                    )
                    if stage == "current"
                    else None
                )
                shared = (
                    sum(
                        p["kind"] == "segment"
                        and any(
                            compact["faces"][i]["kind"] == kind
                            for i in p["affected_faces"]
                        )
                        for p in compact.get("optional_vertices", [])
                    )
                    if stage == "current"
                    else None
                )
                quality.append(
                    dict(
                        city=city,
                        method=method,
                        baseline=stage,
                        kind=kind,
                        optional_interior_points=interior,
                        incident_shared_points=shared,
                        **metrics,
                    )
                )
        groups = compact.get("optimization", {}).get("groups", [])
        stops = Counter(g["stopping"]["primary_reason"] for g in groups)
        bindings = Counter(
            name for g in groups for name in g["stopping"]["binding_limits"]
        )
        flags = {
            stage: {
                kind: quality_flags(metrics, by_stage["current"][kind])
                for kind, metrics in measured.items()
            }
            for stage, measured in by_stage.items()
            if stage != "current"
        }
        building = compact["building_triangles"]
        t0 = compact["T0"]
        bm = record["meshing"]
        valid = (
            bool(compact.get("world_topology", {}).get("valid"))
            and bm.get("status") in ("success", "warning")
            and bm.get("degenerate_cell_count", 1) == 0
        )
        row = dict(
            city=city,
            method=method,
            valid=valid,
            status=record["status"],
            triangles={stage: len(data[0].faces) for stage, data in loaded.items()},
            roof_quality=by_stage["current"]["roof"],
            quality_flags=flags,
            building_triangles=building,
            T0=t0,
            building_growth=building / t0,
            all_group_ceilings=all(
                g["end_triangles"] <= int(1.5 * g["T0"]) for g in groups
            ),
            all_point_caps=all(
                v <= 8 for v in compact["source_optional_vertices"].values()
            )
            and all(g["interior_points"] + g["shared_points"] <= 32 for g in groups)
            and all(
                sum(p["kind"] == "interior" for p in f["optional_vertices"]) <= 2
                for f in compact["faces"]
            ),
            retained_fits=cat["retained_fits"],
            retained_nonflat=cat["retained_nonflat_fits"],
            source_count=sum("source_key" in s for s in cat["sources"]),
            block_retries=len(cat["block_fallbacks"]),
            junction_changes=[
                dict(block_id=b["block_id"], **j)
                for b in cat["blocks"]
                for j in b["subdivision"]["changes"]
            ],
            stopping_reasons=dict(stops),
            binding_limits=dict(bindings),
            elapsed_seconds=record["elapsed_seconds"],
            meshing_seconds=bm["seconds"],
            stage_seconds=cat["assemblies"][-1].get("stage_seconds", {}),
            compact_stage_seconds=compact.get("stage_seconds", {}),
            optional_points=compact["optional_counts"],
            artifact=str(record["artifact"]),
            max_boundary_deviation_m=max(
                f["geometry"].get("boundary_deviation_m", 0.0) for f in compact["faces"]
            ),
        )
        mesh = loaded["current"][0]
        residuals = []
        for p in cat["patch_markers"]:
            xyz = mesh.vertices[
                mesh.faces[np.asarray(mesh.markers) == p["roof_marker"]]
            ].reshape(-1, 3)
            if len(xyz):
                co = np.asarray(p["coefficients"])
                expected = (xyz[:, :2] - p["origin"]) @ co[:2] + co[2]
                residuals.append(
                    float(
                        np.max(abs(xyz[:, 2] - expected)) / np.sqrt(1 + co[:2] @ co[:2])
                    )
                )
        row["world_roof_normal_deviation_m"] = max(residuals, default=0.0)
        row["world_roof_planes_valid"] = bool(
            row["world_roof_normal_deviation_m"]
            <= 1e-7 + 2 * np.spacing(np.max(abs(mesh.vertices)))
        )
        violations = []
        for point in compact.get("optional_vertices", []):
            if (
                point["before"]["penalty"] - point["after"]["penalty"]
                <= 0.015 * point["triangle_cost"]
            ):
                violations.append("insufficient accepted gain")
            pairs = (
                [(point["before"], point["after"])]
                if point["kind"] == "interior"
                else [
                    (point["per_face_before"][k], point["per_face_after"][k])
                    for k in point["per_face_before"]
                ]
            )
            for before, after in pairs:
                if after["minimum_q"] < before["minimum_q"] - 1e-10:
                    violations.append("incident worst quality")
                for angle in (5, 10):
                    if after[f"poor_{angle}_count"] > before[
                        f"poor_{angle}_count"
                    ] or after[f"poor_{angle}_area"] > before[
                        f"poor_{angle}_area"
                    ] + 1e-10 * max(1.0, before[f"poor_{angle}_area"]):
                        violations.append(f"incident {angle}-degree defects")
        row["transaction_guard_violations"] = violations
        previous_sources = {
            s["source_key"]
            for s in production[(city, method)]["catalogue"]["sources"]
            if "source_key" in s
        }
        selected_sources = {
            s["source_key"] for s in cat["sources"] if "source_key" in s
        }
        row["selected_sources_added"] = sorted(selected_sources - previous_sources)
        row["selected_sources_removed"] = sorted(previous_sources - selected_sources)
        if city == "orebro":
            row["orebro_targets"] = orebro_acceptance(row["roof_quality"], method)
        row["evidence_matches_v7"] = (
            record["source_evidence"] == production[(city, method)]["source_evidence"]
        )
        other = current.get((city, "ransac" if method == "catalogue" else "catalogue"))
        row["evidence_matches_other_method"] = (
            other is not None and record["source_evidence"] == other["source_evidence"]
        )
        reference_points = {}
        for stage in ("v7", "current"):
            for block in old_records[stage]["catalogue"].get("blocks", []):
                subdivision = block["subdivision"]
                for point in subdivision.get("changes", []) + subdivision.get(
                    "checked_junctions", []
                ):
                    reference_points[tuple(np.round(point["xy"], 6))] = point["xy"]
        owners = sorted(set(loaded["current"][2]) - {"__terrain__"})
        for source in owners:
            after = loaded["current"]
            for stage in ("v7", "original", "minimal"):
                if stage not in loaded:
                    continue
                before = loaded[stage]
                changes = source_change(before, after, source)
                previous = measured_quality(before[0], before[1], before[2] == source)
                new = measured_quality(after[0], after[1], after[2] == source)
                sources.append(
                    dict(
                        city=city,
                        method=method,
                        source_key=source,
                        baseline=stage,
                        geometry=changes,
                        before=previous,
                        after=new,
                        quality_flags={
                            kind: quality_flags(previous[kind], new[kind])
                            for kind in ("roof", "wall")
                        },
                    )
                )
        for index, xy in enumerate(reference_points.values()):
            measured = {}
            for stage, (mesh, kinds, owned, _) in loaded.items():
                mask = junction_mask(mesh, xy)
                measured[stage] = dict(
                    quality=measured_quality(mesh, kinds, mask),
                    sources={
                        source: measured_quality(mesh, kinds, mask & (owned == source))
                        for source in sorted(set(owned[mask]) - {"__terrain__"})
                    },
                )
            junctions.append(
                dict(
                    city=city,
                    method=method,
                    junction=index,
                    xy=xy,
                    radius_m=2.0,
                    measurements=measured,
                )
            )
        cases.append(row)
    return dict(cases=cases, quality=quality, sources=sources, junctions=junctions)


def write_comparison(output, data):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for name, key in [
        ("comparison", "cases"),
        ("sources", "sources"),
        ("junctions", "junctions"),
    ]:
        (output / f"{name}.json").write_text(json.dumps(data[key], indent=2) + "\n")
    with (output / "quality.csv").open("w") as f:
        rows = data["quality"]
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    cases = data["cases"]
    summary = dict(
        valid_surfaces=sum(r["valid"] for r in cases),
        cases=len(cases),
        all_group_ceilings=all(r.get("all_group_ceilings", False) for r in cases),
        all_point_caps=all(r.get("all_point_caps", False) for r in cases),
        evidence_preserved=all(
            r.get("evidence_matches_v7", False)
            and r.get("evidence_matches_other_method", False)
            for r in cases
        ),
        orebro_targets={
            r["method"]: r["orebro_targets"] for r in cases if r["city"] == "orebro"
        },
        max_city_building_growth=max(r.get("building_growth", 0.0) for r in cases),
    )
    summary["transaction_guard_violations"] = [
        f"{r['city']}/{r['method']}"
        for r in cases
        if r.get("transaction_guard_violations")
    ]
    summary["roof_plane_violations"] = [
        f"{r['city']}/{r['method']}"
        for r in cases
        if not r.get("world_roof_planes_valid")
    ]
    summary["source_selection_preserved"] = all(
        not r.get("selected_sources_added") and not r.get("selected_sources_removed")
        for r in cases
    )
    (output / "acceptance.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# LOD2 junction and stopping comparison",
        "",
        f"Valid exports: **{summary['valid_surfaces']}/{summary['cases']}**. Source evidence preserved: **{summary['evidence_preserved']}**.",
        "",
        "Minimal-to-current refinement, v7-to-current production changes, and original-to-current changes are independent comparisons. See [quality.csv](quality.csv), [per-source attribution](sources.json), [matched junction neighborhoods](junctions.json), and [acceptance targets](acceptance.json). Historical failure cases have no original exported mesh.",
        "",
        "| City / method | Original / v7 / current triangles | Roof p01 original / v7 / current | Roof q < .02 original / v7 / current | Building / T0 | Seconds | Stopping reasons |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in cases:
        if not row["valid"]:
            lines.append(f"| {row['city']} / {row['method']} | Failed | | | | | |")
            continue
        subset = [
            r
            for r in data["quality"]
            if r["city"] == row["city"]
            and r["method"] == row["method"]
            and r["kind"] == "roof"
        ]
        lookup = {r["baseline"]: r for r in subset}
        order = ["original", "v7", "current"]
        p01 = " / ".join(
            f"{lookup[k]['q_p01']:.4f}" if k in lookup else "-" for k in order
        )
        poor = " / ".join(
            str(lookup[k]["below_q_0.02"]) if k in lookup else "-" for k in order
        )
        triangles = " / ".join(str(row["triangles"].get(k, "-")) for k in order)
        reasons = ", ".join(f"{k}: {v}" for k, v in row["stopping_reasons"].items())
        lines.append(
            f"| [{row['city']} / {row['method']}]({row['artifact']}) | {triangles} | {p01} | {poor} | {row['building_growth']:.3f} | {row['elapsed_seconds']:.1f} | {reasons} |"
        )
    lines += ["", "## Örebro targets", ""]
    for method, acceptance in summary["orebro_targets"].items():
        lines.append(
            f"- {method}: p01 target {acceptance['target']['q_p01']}: **{'met' if acceptance['p01_met'] else 'UNMET'}**; severe-count target {acceptance['target']['below_q_0_02']}: **{'met' if acceptance['severe_count_met'] else 'UNMET'}**."
        )
    lines += [
        "",
        "Historical quality flags remain visible per source and face kind in the JSON; improvement over a changed minimal baseline does not clear them. Counts and fractions, severe-angle counts and areas, optional point costs, stage times, source support changes and roof retention are retained separately.",
        "",
    ]
    (output / "report.md").write_text("\n".join(lines))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current", type=Path)
    parser.add_argument("--production", type=Path, required=True)
    parser.add_argument("--original-catalogue", type=Path, required=True)
    parser.add_argument("--original-ransac", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--minimal",
        type=Path,
        help="Optional exported minimal surfaces, including partial city selections",
    )
    args = parser.parse_args()
    originals = load_campaign(args.original_catalogue, "catalogue") | load_campaign(
        args.original_ransac, "ransac"
    )
    minimal = load_campaign(args.minimal) if args.minimal else None
    data = compare_campaigns(
        load_campaign(args.current), load_campaign(args.production), originals, minimal
    )
    summary = write_comparison(args.output or args.current / "comparison", data)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
