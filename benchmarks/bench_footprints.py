"""
Benchmark Stockholm footprint conditioning and flat-mesh preparation.

This script downloads Stockholm tiles, extracts raw building footprints,
conditions them with either the legacy or new pipeline, builds a flat mesh, and
stores both visual and JSON artifacts for side-by-side inspection.

Typical usage:
    python benchmarks/bench_footprints.py --cases 45 46 55 56 --label current --mode legacy
    python benchmarks/bench_footprints.py --cases 45 46 55 56 --label new --mode new
    python benchmarks/bench_footprints.py --cases 55 --mode new --show-plot
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from polyforge import MergeStrategy, fix_clearance, merge_close_polygons, simplify_vwp
from polyforge.ops.clearance.protrusions import remove_narrow_wedges
from shapely import get_parts
from shapely.geometry import GeometryCollection, Polygon, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union
from shapely.validation import make_valid

import dtcc_core
from dtcc_core.builder import (
    build_terrain_raster,
    compute_building_heights,
    extract_roof_points,
)
from dtcc_core.builder.geometry_builders.meshes import (
    _build_ground_mesh_from_coverage,
    _condition_flat_mesh_coverage_regions,
)
from dtcc_core.model import Bounds, Building, City, GeometryType, Surface
try:
    from _stockholm_common import (
        BOX_SIZE,
        DEFAULT_DELAY_BETWEEN_CASES,
        MIN_BUILDING_AREA,
        MIN_BUILDING_DETAIL,
        MIN_MESH_ANGLE,
        NX,
        NY,
        add_cases_argument,
        annotate_heatmap,
        bounds_to_dict,
        case_to_grid,
        format_console_table,
        json_ready,
        load_plot_modules,
        make_bounds,
        repo_root,
        resolve_case_numbers,
        save_results,
        stockholm_output_dir,
    )
except ImportError:
    from benchmarks._stockholm_common import (
        BOX_SIZE,
        DEFAULT_DELAY_BETWEEN_CASES,
        MIN_BUILDING_AREA,
        MIN_BUILDING_DETAIL,
        MIN_MESH_ANGLE,
        NX,
        NY,
        add_cases_argument,
        annotate_heatmap,
        bounds_to_dict,
        case_to_grid,
        format_console_table,
        json_ready,
        load_plot_modules,
        make_bounds,
        repo_root,
        resolve_case_numbers,
        save_results,
        stockholm_output_dir,
    )

DEFAULT_MAX_MESH_SIZE = 10.0
DEFAULT_MIN_MESH_ANGLE = MIN_MESH_ANGLE
DEFAULT_MIN_BUILDING_DETAIL = MIN_BUILDING_DETAIL
DEFAULT_MIN_BUILDING_AREA = MIN_BUILDING_AREA
DEFAULT_MERGE_TOLERANCE = 0.5
DEFAULT_RASTER_CELL_SIZE = 2.0
DEFAULT_RASTER_RADIUS = 3.0

OUTPUT_DIR = stockholm_output_dir("output_footprints")
EPSG = "EPSG:3006"
CACHE_ROOT = Path.home() / "Library" / "Caches" / "dtcc-data"
CACHED_FOOTPRINTS_DIR = CACHE_ROOT / "downloaded-gpkg"

OVERVIEW_METRICS = (
    (("conditioned_polygon_count",), "Cond polys", "viridis", ".0f"),
    (("conditioned_polygon_boundary_metrics", "min_clearance"), "Min clearance", "RdYlGn", ".2f"),
    (
        ("raw_to_conditioned_difference_metrics", "symmetric_difference_area"),
        "SymDiff area",
        "RdYlGn_r",
        ".1f",
    ),
    (
        ("raw_to_conditioned_difference_metrics", "candidate_minus_reference_area"),
        "Extra area",
        "RdYlGn_r",
        ".1f",
    ),
    (("flat_mesh_quality_summary", "element_quality_worst"), "Mesh EQ worst", "RdYlGn", ".3f"),
    (("flat_mesh_quality_summary", "aspect_ratio_worst"), "Mesh AR worst", "RdYlGn_r", ".1f"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_cases_argument(parser)
    parser.add_argument(
        "--mode",
        choices=["auto", "legacy", "new"],
        default="new",
        help="Conditioning pipeline to use. 'new' is the default strict cleaner; 'legacy' runs the Polyforge pipeline for comparison.",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Output label. Defaults to the current git branch when available.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=DEFAULT_DELAY_BETWEEN_CASES,
        help="Delay in seconds between cases to avoid dataset rate limiting.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON, VTU, geopackage, and PNG outputs.",
    )
    parser.add_argument(
        "--git-root",
        type=Path,
        default=None,
        help="Repo root used for git metadata. Defaults to this script's repo root.",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Write per-case summary JSON only; skip mesh/plot/geopackage artifacts.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip PNG plot generation while keeping JSON/mesh artifacts.",
    )
    parser.add_argument(
        "--show-plot",
        action="store_true",
        help=(
            "Show the matplotlib comparison figure interactively after saving it. "
            "Useful for zooming into a case."
        ),
    )
    parser.add_argument(
        "--max-mesh-size",
        type=float,
        default=DEFAULT_MAX_MESH_SIZE,
    )
    parser.add_argument(
        "--min-mesh-angle",
        type=float,
        default=DEFAULT_MIN_MESH_ANGLE,
    )
    parser.add_argument(
        "--min-building-detail",
        type=float,
        default=DEFAULT_MIN_BUILDING_DETAIL,
    )
    parser.add_argument(
        "--min-building-area",
        type=float,
        default=DEFAULT_MIN_BUILDING_AREA,
    )
    parser.add_argument(
        "--merge-tolerance",
        type=float,
        default=DEFAULT_MERGE_TOLERANCE,
    )
    parser.add_argument(
        "--raster-cell-size",
        type=float,
        default=DEFAULT_RASTER_CELL_SIZE,
    )
    parser.add_argument(
        "--raster-radius",
        type=float,
        default=DEFAULT_RASTER_RADIUS,
    )
    parser.add_argument(
        "--no-merge-buildings",
        action="store_true",
        help="Disable gap-closing/merge behavior in the conditioning stage.",
    )
    parser.add_argument(
        "--disable-cleaning-diagnostics",
        action="store_true",
        help=(
            "Disable cleaner stage metrics and cleaner log output. "
            "Useful for runtime comparisons."
        ),
    )
    args = parser.parse_args()
    args.cases_explicit = args.cases is not None
    args.cases = resolve_case_numbers(args.cases)
    return args


def run_git_command(args: list[str], cwd: Path) -> str | None:
    try:
        result = subprocess.run(
            args,
            check=True,
            capture_output=True,
            text=True,
            cwd=cwd,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    value = result.stdout.strip()
    return value or None


def get_git_metadata(repo_root: Path) -> dict[str, str | None]:
    return {
        "branch": run_git_command(["git", "rev-parse", "--abbrev-ref", "HEAD"], repo_root),
        "commit": run_git_command(["git", "rev-parse", "HEAD"], repo_root),
    }


def sanitize_label(label: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in label)
    return cleaned.strip("._") or "unknown"


def try_load_new_cleaner(mode: str):
    try:
        return importlib.import_module("dtcc_core.builder.cleaning")
    except ModuleNotFoundError:
        if mode == "new":
            raise RuntimeError(
                "--mode new requires dtcc_core.builder.cleaning, but that module is not available on this branch."
            ) from None
        return None


def extract_polygon_parts(geometry: BaseGeometry) -> list[Polygon]:
    polygons: list[Polygon] = []
    for part in get_parts(geometry):
        if isinstance(part, Polygon):
            polygons.append(part)
        elif hasattr(part, "geom_type") and "Polygon" in part.geom_type:
            polygons.extend(extract_polygon_parts(part))
    return [poly for poly in polygons if not poly.is_empty and poly.area > 0]


def extract_raw_footprints(buildings: list[Building]) -> tuple[list[Polygon], list[list[int]]]:
    polygons: list[Polygon] = []
    source_map: list[list[int]] = []
    for index, building in enumerate(buildings):
        geom = building.flatten_geometry(GeometryType.LOD0)
        if geom is None:
            continue
        polygon = geom.to_polygon()
        if polygon.is_empty:
            continue
        for part in extract_polygon_parts(make_valid(polygon)):
            polygons.append(part)
            source_map.append([index])
    return polygons, source_map


def compose_index_map(
    parent_map: list[list[int]],
    child_map: list[list[int]],
) -> list[list[int]]:
    composed: list[list[int]] = []
    for child_indices in child_map:
        combined: list[int] = []
        for idx in child_indices:
            if idx < 0 or idx >= len(parent_map):
                continue
            combined.extend(parent_map[idx])
        if combined:
            composed.append(sorted(set(combined)))
    return composed


def extract_lod0_polygons(buildings: list[Building]) -> list[Polygon]:
    polygons: list[Polygon] = []
    for building in buildings:
        geom = building.flatten_geometry(GeometryType.LOD0)
        if geom is None:
            continue
        polygon = geom.to_polygon(simplify=0.0)
        if polygon.is_empty:
            continue
        polygons.extend(extract_polygon_parts(make_valid(polygon)))
    return polygons


def _legacy_merge_building_footprints(
    buildings: list[Building],
    *,
    lod: GeometryType,
    max_distance: float,
    min_area: float,
) -> tuple[list[Building], list[list[int]]]:
    if len(buildings) <= 1:
        return buildings, [[i] for i in range(len(buildings))]

    source_indices: list[int] = []
    footprints: list[Polygon] = []
    building_heights: list[float] = []

    for idx, building in enumerate(buildings):
        flattened_geom = building.get_footprint(lod)
        if flattened_geom is None:
            continue
        footprint = flattened_geom.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        source_indices.append(idx)
        building_heights.append(flattened_geom.zmax)
        footprints.append(footprint)

    merged_footprints, merged_indices = merge_close_polygons(
        footprints,
        max_distance,
        merge_strategy=MergeStrategy.BOUNDARY_EXTENSION,
        preserve_holes=True,
        insert_vertices=True,
        return_mapping=True,
        buffer_cleaning=True,
    )

    merged_buildings: list[Building] = []
    merged_indices_global: list[list[int]] = []

    for footprint, local_indices in zip(merged_footprints, merged_indices):
        if footprint.geom_type == "MultiPolygon" or footprint.is_empty or footprint.area < min_area:
            continue
        global_indices = [source_indices[i] for i in local_indices]
        num = sum(building_heights[i] * footprints[i].area for i in local_indices)
        den = sum(footprints[i].area for i in local_indices)
        height = num / den if den > 0 else 0.0

        surface = Surface()
        surface.from_polygon(footprint, height)
        merged_building = Building()
        merged_building.add_geometry(surface, GeometryType.LOD0)
        merged_building.attributes["height"] = height
        merged_buildings.append(merged_building)
        merged_indices_global.append(global_indices)

    return merged_buildings, merged_indices_global


def _legacy_fix_building_footprint_clearance(
    buildings: list[Building],
    *,
    clearance: float,
    lod: GeometryType,
) -> tuple[list[Building], list[list[int]]]:
    fixed_buildings: list[Building] = []
    index_map: list[list[int]] = []

    for idx, building in enumerate(buildings):
        lod_geom = building.flatten_geometry(lod)
        if lod_geom is None:
            continue
        footprint = lod_geom.to_polygon()
        if footprint is None or footprint.is_empty or footprint.geom_type == "MultiPolygon":
            continue
        footprint = fix_clearance(footprint, clearance)
        if footprint.geom_type == "MultiPolygon":
            footprints = [geom for geom in footprint.geoms if isinstance(geom, Polygon)]
            footprint = unary_union(footprints)
            if footprint.geom_type == "MultiPolygon":
                footprint = max(footprint.geoms, key=lambda polygon: polygon.area)
        if footprint.geom_type != "Polygon":
            continue
        footprint = remove_narrow_wedges(footprint, min_depth=clearance)
        surface = Surface()
        surface.from_polygon(footprint, lod_geom.zmax)
        fixed_building = building.copy()
        fixed_building.add_geometry(surface, GeometryType.LOD0)
        fixed_building.calculate_bounds()
        fixed_buildings.append(fixed_building)
        index_map.append([idx])

    return fixed_buildings, index_map


def _legacy_simplify_building_footprints(
    buildings: list[Building],
    *,
    tolerance: float,
    lod: GeometryType,
) -> tuple[list[Building], list[list[int]]]:
    simplified_buildings: list[Building] = []
    index_map: list[list[int]] = []

    for idx, building in enumerate(buildings):
        lod_geom = building.flatten_geometry(lod)
        if lod_geom is None:
            continue
        footprint = lod_geom.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        footprint = simplify_vwp(footprint, tolerance)
        surface = Surface()
        surface.from_polygon(footprint, lod_geom.zmax)
        simplified_building = building.copy()
        simplified_building.add_geometry(surface, GeometryType.LOD0)
        simplified_building.calculate_bounds()
        simplified_buildings.append(simplified_building)
        index_map.append([idx])

    return simplified_buildings, index_map


def run_legacy_conditioning(
    buildings: list[Building],
    *,
    lod: GeometryType = GeometryType.LOD0,
    merge_buildings: bool,
    merge_tolerance: float,
    min_building_area: float,
    min_building_detail: float,
) -> tuple[list[Polygon], list[list[int]], dict[str, Any]]:
    if merge_buildings:
        merged_buildings, merged_index_map = _legacy_merge_building_footprints(
            buildings,
            lod=lod,
            max_distance=merge_tolerance,
            min_area=min_building_area,
        )
        cleared_buildings, cleared_index_map = _legacy_fix_building_footprint_clearance(
            merged_buildings,
            clearance=min_building_detail,
            lod=lod,
        )
        current_index_map = compose_index_map(merged_index_map, cleared_index_map)
        merged_again, merged_again_index_map = _legacy_merge_building_footprints(
            cleared_buildings,
            lod=lod,
            max_distance=merge_tolerance,
            min_area=min_building_area,
        )
        current_index_map = compose_index_map(current_index_map, merged_again_index_map)
        simplified_buildings, simplified_index_map = _legacy_simplify_building_footprints(
            merged_again,
            tolerance=min_building_detail,
            lod=lod,
        )
        source_map = compose_index_map(current_index_map, simplified_index_map)
        conditioned_polygons = extract_lod0_polygons(simplified_buildings)
    else:
        simplified_buildings, source_map = _legacy_simplify_building_footprints(
            buildings,
            tolerance=min_building_detail,
            lod=lod,
        )
        conditioned_polygons = extract_lod0_polygons(simplified_buildings)

    diagnostics = {
        "pipeline": "legacy_polyforge",
        "output_count": len(conditioned_polygons),
    }
    return conditioned_polygons, source_map, diagnostics


def run_new_conditioning(
    cleaner_module,
    buildings: list[Building],
    *,
    merge_buildings: bool,
    merge_tolerance: float,
    min_building_area: float,
    min_building_detail: float,
    disable_cleaning_diagnostics: bool = False,
) -> tuple[list[Polygon], list[list[int]], dict[str, Any]]:
    options = cleaner_module.ConditioningOptions(
        precision_grid=None,
        min_feature_size=min_building_detail,
        merge_distance=merge_tolerance if merge_buildings else 0.0,
        min_area=min_building_area,
        min_hole_area=min_building_detail**2,
        collect_stage_metrics=not disable_cleaning_diagnostics,
        enable_logging=not disable_cleaning_diagnostics,
    )
    result = cleaner_module.condition_building_footprints(
        buildings,
        lod=GeometryType.LOD0,
        options=options,
    )
    return result.polygons, result.source_map, result.diagnostics


def resolve_mode(mode: str):
    cleaner_module = try_load_new_cleaner(mode)
    if mode == "auto":
        return ("new", cleaner_module) if cleaner_module is not None else ("legacy", None)
    if mode == "new":
        return "new", cleaner_module
    return "legacy", None


def height_for_sources(
    buildings: list[Building],
    source_indices: list[int],
) -> float:
    heights: list[float] = []
    weights: list[float] = []
    for index in source_indices:
        if index < 0 or index >= len(buildings):
            continue
        building = buildings[index]
        geom = building.flatten_geometry(GeometryType.LOD0)
        if geom is None:
            continue
        polygon = geom.to_polygon(simplify=0.0)
        area = float(max(polygon.area, 0.0))
        height = getattr(building, "height", None)
        if height is None or height <= 0:
            continue
        heights.append(float(height))
        weights.append(area if area > 0 else 1.0)
    if not heights:
        return DEFAULT_MAX_MESH_SIZE
    return float(np.average(heights, weights=weights))


def make_conditioned_buildings(
    polygons: list[Polygon],
    source_map: list[list[int]],
    original_buildings: list[Building],
) -> list[Building]:
    conditioned_buildings: list[Building] = []
    for polygon, indices in zip(polygons, source_map):
        surface = Surface()
        surface.from_polygon(polygon, 0.0)
        building = Building()
        building.add_geometry(surface, GeometryType.LOD0)
        building.attributes["height"] = height_for_sources(original_buildings, indices)
        building.attributes["source_map"] = list(indices)
        conditioned_buildings.append(building)
    return conditioned_buildings


def build_mesh_from_conditioned_footprints(
    terrain_raster,
    conditioned_polygons: list[Polygon],
    source_map: list[list[int]],
    source_buildings: list[Building],
    *,
    max_mesh_size: float,
    min_mesh_angle: float,
    min_building_detail: float,
    footprint_diagnostics: dict[str, Any] | None = None,
    disable_cleaning_diagnostics: bool = False,
):
    # This benchmark path inspects stage-2 meshing from the output of stage-1
    # conditioning. Re-running build_city_flat_mesh() would send the footprints
    # back through conditioning with different defaults and confound the
    # comparison.
    del source_buildings

    marker_lookup: dict[tuple[int, ...], int] = {}
    building_markers: list[int] = []
    for source_indices in source_map:
        marker_key = tuple(sorted(set(source_indices)))
        if marker_key not in marker_lookup:
            marker_lookup[marker_key] = len(marker_lookup)
        building_markers.append(marker_lookup[marker_key])

    flat_mesh_bounds = (
        terrain_raster.bounds.xmin,
        terrain_raster.bounds.ymin,
        terrain_raster.bounds.xmax,
        terrain_raster.bounds.ymax,
    )
    region_polygons, region_markers = _condition_flat_mesh_coverage_regions(
        bounds=flat_mesh_bounds,
        building_polygons=conditioned_polygons,
        building_markers=building_markers,
        hole_polygons=[],
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        footprint_diagnostics=footprint_diagnostics or {},
        cleaning_diagnostics=not disable_cleaning_diagnostics,
    )
    flat_mesh, _active_mesher = _build_ground_mesh_from_coverage(
        region_polygons=region_polygons,
        region_markers=region_markers,
        bounds=flat_mesh_bounds,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        mesher=None,
        sort_triangles=True,
    )
    return flat_mesh


def color_from_key(key: str) -> str:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    r = int(digest[0:2], 16) / 255.0
    g = int(digest[2:4], 16) / 255.0
    b = int(digest[4:6], 16) / 255.0
    rgb = np.array([r, g, b])
    rgb = 0.25 + 0.65 * rgb
    return "#{:02x}{:02x}{:02x}".format(*(np.round(rgb * 255).astype(int)))


def make_geodataframe(
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    key_mode: str,
) -> gpd.GeoDataFrame:
    records: list[dict[str, Any]] = []
    for idx, (polygon, indices) in enumerate(zip(polygons, source_map)):
        records.append(
            {
                "footprint_id": idx,
                "source_map": ",".join(str(i) for i in indices),
                "source_count": len(indices),
                "color_key": idx if key_mode == "raw" else tuple(indices),
                "color": color_from_key(
                    str(idx) if key_mode == "raw" else ",".join(str(i) for i in indices)
                ),
                "area": float(polygon.area),
                "geometry": polygon,
            }
        )
    if not records:
        return gpd.GeoDataFrame(
            {
                "footprint_id": pd.Series(dtype=int),
                "source_map": pd.Series(dtype=str),
                "source_count": pd.Series(dtype=int),
                "color_key": pd.Series(dtype=object),
                "color": pd.Series(dtype=str),
                "area": pd.Series(dtype=float),
                "geometry": gpd.GeoSeries([], crs=EPSG),
            },
            geometry="geometry",
            crs=EPSG,
        )
    return gpd.GeoDataFrame(records, geometry="geometry", crs=EPSG)


def make_bounds_geodataframe(bounds: Bounds) -> gpd.GeoDataFrame:
    polygon = box(bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
    return gpd.GeoDataFrame(
        [{"name": "bounds", "geometry": polygon}],
        geometry="geometry",
        crs=EPSG,
    )


def save_geopackage(
    path: Path,
    raw_gdf: gpd.GeoDataFrame,
    conditioned_gdf: gpd.GeoDataFrame,
    bounds: Bounds,
) -> None:
    if path.exists():
        path.unlink()
    raw_gdf.to_file(path, layer="raw", driver="GPKG")
    conditioned_gdf.to_file(path, layer="conditioned", driver="GPKG")
    make_bounds_geodataframe(bounds).to_file(path, layer="bounds", driver="GPKG")


def save_source_map_csv(path: Path, source_map: list[list[int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["conditioned_polygon_id", "source_indices", "source_count"],
        )
        writer.writeheader()
        for polygon_id, indices in enumerate(source_map):
            writer.writerow(
                {
                    "conditioned_polygon_id": polygon_id,
                    "source_indices": ",".join(str(i) for i in indices),
                    "source_count": len(indices),
                }
            )


def coverage_metrics(polygons: list[Polygon]) -> dict[str, float | int | None]:
    if not polygons:
        return {
            "polygon_count": 0,
            "total_area": 0.0,
            "union_area": 0.0,
            "overlap_area": 0.0,
            "minimum_clearance": None,
        }

    total_area = float(sum(polygon.area for polygon in polygons))
    union_geom = unary_union(polygons)
    union_area = float(union_geom.area)
    overlap_area = float(max(total_area - union_area, 0.0))

    clearances: list[float] = []
    for polygon in polygons:
        try:
            clearance = polygon.minimum_clearance
        except Exception:
            continue
        if np.isfinite(clearance):
            clearances.append(float(clearance))

    return {
        "polygon_count": len(polygons),
        "total_area": total_area,
        "union_area": union_area,
        "overlap_area": overlap_area,
        "minimum_clearance": min(clearances) if clearances else None,
    }


def coverage_difference_metrics(
    reference_polygons: list[Polygon],
    candidate_polygons: list[Polygon],
) -> dict[str, float]:
    reference = unary_union(reference_polygons) if reference_polygons else GeometryCollection()
    candidate = unary_union(candidate_polygons) if candidate_polygons else GeometryCollection()
    return {
        "union_area_delta": float(candidate.area - reference.area),
        "symmetric_difference_area": float(reference.symmetric_difference(candidate).area),
        "reference_minus_candidate_area": float(reference.difference(candidate).area),
        "candidate_minus_reference_area": float(candidate.difference(reference).area),
    }


def polygon_boundary_metrics(
    polygons: list[Polygon],
    *,
    short_edge_threshold: float,
) -> dict[str, float | int | None]:
    if not polygons:
        return {
            "polygon_count": 0,
            "vertex_count": 0,
            "edge_count": 0,
            "mean_edge_length": None,
            "min_edge_length": None,
            "short_edge_count": 0,
            "mean_clearance": None,
            "min_clearance": None,
        }

    vertex_count = 0
    edge_count = 0
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
                edge_count += 1
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
        "edge_count": edge_count,
        "mean_edge_length": float(np.mean(edge_lengths)) if edge_lengths else None,
        "min_edge_length": float(np.min(edge_lengths)) if edge_lengths else None,
        "short_edge_count": short_edge_count,
        "mean_clearance": float(np.mean(clearances)) if clearances else None,
        "min_clearance": float(np.min(clearances)) if clearances else None,
    }


def mesh_quality_summary(mesh_quality: dict[str, Any]) -> dict[str, float]:
    return {
        "element_quality_mean": float(mesh_quality["element_quality"]["mean"]),
        "element_quality_worst": float(mesh_quality["element_quality"]["min"]),
        "aspect_ratio_mean": float(mesh_quality["aspect_ratio"]["mean"]),
        "aspect_ratio_worst": float(mesh_quality["aspect_ratio"]["max"]),
        "edge_ratio_mean": float(mesh_quality["edge_ratio"]["mean"]),
        "edge_ratio_worst": float(mesh_quality["edge_ratio"]["max"]),
        "skewness_mean": float(mesh_quality["skewness"]["mean"]),
        "skewness_worst": float(mesh_quality["skewness"]["max"]),
    }


def empty_mesh_quality_summary() -> dict[str, float | None]:
    return {
        "element_quality_mean": None,
        "element_quality_worst": None,
        "aspect_ratio_mean": None,
        "aspect_ratio_worst": None,
        "edge_ratio_mean": None,
        "edge_ratio_worst": None,
        "skewness_mean": None,
        "skewness_worst": None,
    }


def timing_summary(timings_seconds: dict[str, float]) -> dict[str, float]:
    conditioning = float(timings_seconds.get("conditioning", 0.0))
    meshing = float(timings_seconds.get("meshing", 0.0))
    return {
        "conditioning_seconds": conditioning,
        "meshing_seconds": meshing,
        "core_seconds": conditioning + meshing,
    }


def error_details(stage: str, exc: Exception) -> dict[str, str]:
    return {
        "stage": stage,
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def _fmt_metric(value: float | int | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def polygon_metrics_text(
    metrics: dict[str, float | int | None],
    *,
    short_edge_threshold: float,
    delta_metrics: dict[str, float] | None = None,
) -> str:
    lines = [
        f"polys {metrics['polygon_count']}  verts {metrics['vertex_count']}",
        f"edge mean {_fmt_metric(metrics['mean_edge_length'])}  min {_fmt_metric(metrics['min_edge_length'])}",
        f"short < {short_edge_threshold:.2f}m: {metrics['short_edge_count']}  clear {_fmt_metric(metrics['min_clearance'])}",
    ]
    if delta_metrics is not None:
        lines.extend(
            [
                f"symdiff {_fmt_metric(delta_metrics['symmetric_difference_area'])} m2",
                f"missing {_fmt_metric(delta_metrics['reference_minus_candidate_area'])}  extra {_fmt_metric(delta_metrics['candidate_minus_reference_area'])}",
            ]
        )
    return "\n".join(lines)


def mesh_metrics_text(metrics: dict[str, float | None]) -> str:
    return "\n".join(
        [
            (
                f"EQ mean {_fmt_metric(metrics['element_quality_mean'], 3)}  "
                f"worst {_fmt_metric(metrics['element_quality_worst'], 3)}"
            ),
            (
                f"AR mean {_fmt_metric(metrics['aspect_ratio_mean'], 3)}  "
                f"worst {_fmt_metric(metrics['aspect_ratio_worst'], 3)}"
            ),
            (
                f"ER mean {_fmt_metric(metrics['edge_ratio_mean'], 3)}  "
                f"worst {_fmt_metric(metrics['edge_ratio_worst'], 3)}"
            ),
            (
                f"Skew mean {_fmt_metric(metrics['skewness_mean'], 3)}  "
                f"worst {_fmt_metric(metrics['skewness_worst'], 3)}"
            ),
        ]
    )


def timing_metrics_text(metrics: dict[str, float]) -> str:
    return "\n".join(
        [
            f"Cond {metrics['conditioning_seconds']:.3f}s",
            f"Mesh {metrics['meshing_seconds']:.3f}s",
            f"Core {metrics['core_seconds']:.3f}s",
        ]
    )


def mesh_edge_segments(mesh) -> list[np.ndarray]:
    segments: list[np.ndarray] = []
    if mesh is None:
        return segments
    if mesh.faces is None or mesh.vertices is None:
        return segments
    xy = mesh.vertices[:, :2]
    for face in mesh.faces:
        points = xy[np.asarray(face, dtype=int)]
        segments.append(points[[0, 1]])
        segments.append(points[[1, 2]])
        segments.append(points[[2, 0]])
    return segments


def set_panel_extent(axes, bounds: Bounds) -> None:
    padding = 0.05 * BOX_SIZE
    xmin = bounds.xmin - padding
    xmax = bounds.xmax + padding
    ymin = bounds.ymin - padding
    ymax = bounds.ymax + padding
    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])


def plot_case(
    output_path: Path,
    bounds: Bounds,
    raw_gdf: gpd.GeoDataFrame,
    conditioned_gdf: gpd.GeoDataFrame,
    flat_mesh,
    raw_polygon_metrics: dict[str, float | int | None],
    conditioned_polygon_metrics: dict[str, float | int | None],
    delta_metrics: dict[str, float],
    mesh_metrics: dict[str, float | None],
    timing_metrics: dict[str, float],
    short_edge_threshold: float,
    title: str,
    error_info: dict[str, str] | None = None,
    show: bool = False,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 7.5), constrained_layout=True)
    mesh_available = flat_mesh is not None

    for ax, gdf, panel_title in [
        (axes[0], raw_gdf, "Raw footprints"),
        (axes[1], conditioned_gdf, "Conditioned footprints"),
        (
            axes[2],
            conditioned_gdf,
            "Conditioned + flat mesh" if mesh_available else "Conditioned footprints",
        ),
    ]:
        if not gdf.empty:
            gdf.plot(
                ax=ax,
                color=gdf["color"].tolist(),
                edgecolor="black",
                linewidth=0.6,
            )
        ax.set_title(panel_title)

    segments = mesh_edge_segments(flat_mesh)
    if segments:
        axes[2].add_collection(
            LineCollection(segments, colors="black", linewidths=0.25, alpha=0.6)
        )

    text_kwargs = {
        "transform": axes[0].transAxes,
        "va": "bottom",
        "ha": "left",
        "fontsize": 9,
        "bbox": {
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "alpha": 0.88,
            "edgecolor": "black",
            "linewidth": 0.4,
        },
    }
    axes[0].text(
        0.02,
        0.02,
        polygon_metrics_text(
            raw_polygon_metrics,
            short_edge_threshold=short_edge_threshold,
        ),
        **text_kwargs,
    )
    axes[1].text(
        0.02,
        0.02,
        polygon_metrics_text(
            conditioned_polygon_metrics,
            short_edge_threshold=short_edge_threshold,
            delta_metrics=delta_metrics,
        ),
        **{**text_kwargs, "transform": axes[1].transAxes},
    )
    panel_text = mesh_metrics_text(mesh_metrics) + "\n" + timing_metrics_text(timing_metrics)
    if error_info is not None:
        panel_text += (
            "\n"
            f"{error_info['stage']} failed: {error_info['type']}\n"
            f"{error_info['message']}"
        )
    axes[2].text(
        0.02,
        0.02,
        panel_text,
        **{**text_kwargs, "transform": axes[2].transAxes},
    )

    set_panel_extent(axes, bounds)
    fig.suptitle(title)
    fig.savefig(output_path, dpi=200)
    if show:
        plt.show()
    plt.close(fig)


def nested_metric(summary: dict[str, Any], path: tuple[str, ...]) -> float | int | None:
    value: Any = summary
    for key in path:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return value
    return None


def build_summary_grid(
    summaries: list[dict[str, Any]],
    path: tuple[str, ...],
) -> np.ndarray:
    grid = np.full((NY, NX), np.nan)
    for summary in summaries:
        value = nested_metric(summary, path)
        if value is None:
            continue
        ix, iy = case_to_grid(int(summary["case"]))
        grid[iy, ix] = float(value)
    return grid


def plot_overview(
    summaries: list[dict[str, Any]],
    output_path: Path,
    *,
    label: str,
    mode_name: str,
) -> None:
    plt_mod, _, _, _ = load_plot_modules()
    fig, axes = plt_mod.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    axes = np.asarray(axes).reshape(2, 3)

    for ax, (path, title, cmap, fmt) in zip(axes.flat, OVERVIEW_METRICS):
        grid = build_summary_grid(summaries, path)
        image = ax.imshow(grid, origin="lower", cmap=cmap, aspect="equal")
        annotate_heatmap(ax, grid, fmt)
        ax.set_title(title)
        ax.set_xticks(range(NX))
        ax.set_yticks(range(NY))
        ax.set_xlabel("Grid X")
        ax.set_ylabel("Grid Y")
        fig.colorbar(image, ax=ax, shrink=0.82, pad=0.03)

    fig.suptitle(f"Footprint benchmark overview\nlabel={label}, mode={mode_name}")
    fig.savefig(output_path, dpi=180)
    plt_mod.close(fig)


def status_emoji(status: str | None) -> str:
    return "✅" if status == "completed" else "❌"


def aggregate_status_emoji(successes: int, failures: int) -> str:
    if failures == 0:
        return "✅"
    if successes == 0:
        return "❌"
    return "⚠️"


def build_summary_report(
    summaries: list[dict[str, Any]],
    *,
    label: str,
    mode_name: str,
    elapsed_seconds: float,
) -> str:
    if not summaries:
        return ""

    case_rows: list[list[Any]] = []
    completed = [summary for summary in summaries if summary.get("status") == "completed"]
    for summary in summaries:
        if summary.get("status") == "completed":
            detail_tail = _fmt_metric(
                nested_metric(summary, ("timing_summary", "core_seconds")),
                2,
            )
        else:
            error = summary.get("error") or {}
            detail_tail = (
                f"{error.get('stage', '-')}: {error.get('type', '-')}"
                f" ({error.get('message', '-')})"
            )
        case_rows.append(
            [
                f"{summary['case']:03d}",
                f"{status_emoji(summary.get('status'))} {summary['status']}",
                summary["raw_polygon_count"],
                summary["conditioned_polygon_count"],
                _fmt_metric(
                    nested_metric(summary, ("conditioned_polygon_boundary_metrics", "min_clearance"))
                ),
                _fmt_metric(
                    nested_metric(
                        summary,
                        ("raw_to_conditioned_difference_metrics", "symmetric_difference_area"),
                    )
                ),
                _fmt_metric(
                    nested_metric(
                        summary,
                        ("raw_to_conditioned_difference_metrics", "candidate_minus_reference_area"),
                    )
                ),
                _fmt_metric(
                    nested_metric(summary, ("flat_mesh_quality_summary", "element_quality_worst")),
                    3,
                ),
                _fmt_metric(
                    nested_metric(summary, ("flat_mesh_quality_summary", "aspect_ratio_worst")),
                    2,
                ),
                detail_tail,
            ]
        )

    if not completed:
        summary_rows = [[aggregate_status_emoji(0, len(summaries)), 0, len(summaries), "-", "-", "-", "-", "-"]]
        return "\n".join(
            [
                "bench_footprints",
                f"Config: mode={mode_name}, label={label}",
                f"Elapsed time: {elapsed_seconds:.1f}s",
                "",
                format_console_table(
                    ["Status", "Completed", "Failed", "Worst clear", "Median symdiff", "Worst EQ", "Worst AR", "Mean core s"],
                    summary_rows,
                    title="Summary",
                ),
                "",
                format_console_table(
                    ["Case", "Status", "Raw", "Cond", "Clear min", "SymDiff", "Extra", "EQ worst", "AR worst", "Core s/Error"],
                    case_rows,
                    title="Detailed results",
                ),
            ]
        )

    symdiff_values = [
        value
        for value in (
            nested_metric(
                summary,
                ("raw_to_conditioned_difference_metrics", "symmetric_difference_area"),
            )
            for summary in completed
        )
        if value is not None
    ]
    core_seconds = [
        value
        for value in (
            nested_metric(summary, ("timing_summary", "core_seconds")) for summary in completed
        )
        if value is not None
    ]

    failures = sum(1 for summary in summaries if summary.get("status") != "completed")
    aggregate_row = [[
        aggregate_status_emoji(len(completed), failures),
        len(completed),
        failures,
        _fmt_metric(
            min(
                (
                    value
                    for value in (
                        nested_metric(summary, ("conditioned_polygon_boundary_metrics", "min_clearance"))
                        for summary in completed
                    )
                    if value is not None
                ),
                default=None,
            )
        ),
        _fmt_metric(float(np.median(symdiff_values)) if symdiff_values else None),
        _fmt_metric(
            min(
                (
                    value
                    for value in (
                        nested_metric(summary, ("flat_mesh_quality_summary", "element_quality_worst"))
                        for summary in completed
                    )
                    if value is not None
                ),
                default=None,
            ),
            3,
        ),
        _fmt_metric(
            max(
                (
                    value
                    for value in (
                        nested_metric(summary, ("flat_mesh_quality_summary", "aspect_ratio_worst"))
                        for summary in completed
                    )
                    if value is not None
                ),
                default=None,
            ),
            2,
        ),
        _fmt_metric(float(np.mean(core_seconds)) if core_seconds else None, 2),
    ]]
    return "\n".join(
        [
            "bench_footprints",
            f"Config: mode={mode_name}, label={label}",
            f"Elapsed time: {elapsed_seconds:.1f}s",
            "",
            format_console_table(
                [
                    "Status",
                    "Completed",
                    "Failed",
                    "Worst clear",
                    "Median symdiff",
                    "Worst EQ",
                    "Worst AR",
                    "Mean core s",
                ],
                aggregate_row,
                title="Summary",
            ),
            "",
            format_console_table(
                [
                    "Case",
                    "Status",
                    "Raw",
                    "Cond",
                    "Clear min",
                    "SymDiff",
                    "Extra",
                    "EQ worst",
                    "AR worst",
                    "Core s/Error",
                ],
                case_rows,
                title="Detailed results",
            ),
        ]
    )


def prepare_case_inputs(
    bounds: Bounds,
    raster_cell_size: float,
    raster_radius: float,
) -> tuple[Any, list[Building], dict[str, float]]:
    timings: dict[str, float] = {}

    t0 = time.perf_counter()
    pointcloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
    cached_footprint_files = sorted(CACHED_FOOTPRINTS_DIR.glob("*.gpkg"))
    buildings = []
    if cached_footprint_files:
        for path in cached_footprint_files:
            buildings = dtcc_core.io.load_footprints(str(path), bounds=bounds)
            if buildings:
                break
    if not buildings:
        buildings = dtcc_core.io.data.download_footprints(bounds=bounds)
    timings["download"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    terrain_raster = build_terrain_raster(
        pointcloud,
        cell_size=raster_cell_size,
        radius=raster_radius,
        ground_only=True,
    )
    buildings = extract_roof_points(buildings, pointcloud)
    buildings = compute_building_heights(buildings, terrain_raster, overwrite=True)
    timings["terrain_and_buildings"] = time.perf_counter() - t0

    return terrain_raster, buildings, timings


def run_case(
    number: int,
    args: argparse.Namespace,
    git_metadata: dict[str, str | None],
    *,
    mode_name: str,
    cleaner_module,
) -> dict[str, Any]:
    ix, iy = case_to_grid(number)
    bounds = make_bounds(ix, iy)
    case_dir = args.output_dir / f"{number:03d}"
    case_dir.mkdir(parents=True, exist_ok=True)

    terrain_raster, buildings, timings = prepare_case_inputs(
        bounds,
        raster_cell_size=args.raster_cell_size,
        raster_radius=args.raster_radius,
    )

    raw_polygons, raw_source_map = extract_raw_footprints(buildings)
    conditioned_polygons: list[Polygon] = []
    source_map: list[list[int]] = []
    diagnostics: dict[str, Any] = {}
    flat_mesh = None
    mesh_quality: dict[str, Any] = {}
    mesh_metrics: dict[str, float | None] = empty_mesh_quality_summary()
    case_status = "completed"
    error_info: dict[str, str] | None = None

    t0 = time.perf_counter()
    try:
        if mode_name == "legacy":
            conditioned_polygons, source_map, diagnostics = run_legacy_conditioning(
                buildings,
                merge_buildings=not args.no_merge_buildings,
                merge_tolerance=args.merge_tolerance,
                min_building_area=args.min_building_area,
                min_building_detail=args.min_building_detail,
            )
        else:
            conditioned_polygons, source_map, diagnostics = run_new_conditioning(
                cleaner_module,
                buildings,
                merge_buildings=not args.no_merge_buildings,
                merge_tolerance=args.merge_tolerance,
                min_building_area=args.min_building_area,
                min_building_detail=args.min_building_detail,
                disable_cleaning_diagnostics=args.disable_cleaning_diagnostics,
            )
    except Exception as exc:
        case_status = "failed"
        error_info = error_details("conditioning", exc)
    timings["conditioning"] = time.perf_counter() - t0

    if case_status == "completed":
        t0 = time.perf_counter()
        try:
            flat_mesh = build_mesh_from_conditioned_footprints(
                terrain_raster,
                conditioned_polygons,
                source_map,
                buildings,
                max_mesh_size=args.max_mesh_size,
                min_mesh_angle=args.min_mesh_angle,
                min_building_detail=args.min_building_detail,
                footprint_diagnostics=diagnostics,
                disable_cleaning_diagnostics=args.disable_cleaning_diagnostics,
            )
            mesh_quality = flat_mesh.quality()
            mesh_metrics = mesh_quality_summary(mesh_quality)
            if not args.summary_only:
                flat_mesh.save(case_dir / "flat_mesh.vtu")
        except Exception as exc:
            case_status = "failed"
            error_info = error_details("meshing", exc)
        timings["meshing"] = time.perf_counter() - t0
    timing_metrics = timing_summary(timings)

    if not args.summary_only:
        raw_gdf = make_geodataframe(raw_polygons, raw_source_map, key_mode="raw")
        conditioned_gdf = make_geodataframe(
            conditioned_polygons,
            [sorted(set(indices)) for indices in source_map],
            key_mode="conditioned",
        )
        save_geopackage(case_dir / "footprints.gpkg", raw_gdf, conditioned_gdf, bounds)
        save_source_map_csv(case_dir / "source_map.csv", source_map)

    raw_metrics = coverage_metrics(raw_polygons)
    conditioned_metrics = coverage_metrics(conditioned_polygons)
    delta_metrics = coverage_difference_metrics(raw_polygons, conditioned_polygons)
    raw_boundary_metrics = polygon_boundary_metrics(
        raw_polygons,
        short_edge_threshold=args.min_building_detail,
    )
    conditioned_boundary_metrics = polygon_boundary_metrics(
        conditioned_polygons,
        short_edge_threshold=args.min_building_detail,
    )

    if not args.summary_only and not args.no_plots:
        t0 = time.perf_counter()
        plot_case(
            case_dir / "comparison.png",
            bounds,
            raw_gdf,
            conditioned_gdf,
            flat_mesh,
            raw_boundary_metrics,
            conditioned_boundary_metrics,
            delta_metrics,
            mesh_metrics,
            timing_metrics,
            args.min_building_detail,
            title=f"Case {number:03d} | {mode_name}",
            error_info=error_info,
            show=bool(args.show_plot and len(args.cases) == 1),
        )
        timings["plotting"] = time.perf_counter() - t0

    artifacts: dict[str, str] = {}
    if not args.summary_only:
        if not args.no_plots:
            artifacts["comparison_png"] = "comparison.png"
        if flat_mesh is not None:
            artifacts["flat_mesh_vtu"] = "flat_mesh.vtu"
        artifacts["footprints_gpkg"] = "footprints.gpkg"
        artifacts["source_map_csv"] = "source_map.csv"

    summary = {
        "git_branch": git_metadata["branch"],
        "git_commit": git_metadata["commit"],
        "mode": mode_name,
        "requested_mode": args.mode,
        "status": case_status,
        "error": error_info,
        "case": number,
        "bounds": bounds_to_dict(bounds),
        "input_building_count": len(buildings),
        "raw_polygon_count": len(raw_polygons),
        "conditioned_polygon_count": len(conditioned_polygons),
        "source_map_sizes": [len(indices) for indices in source_map],
        "raw_coverage_metrics": raw_metrics,
        "raw_polygon_boundary_metrics": raw_boundary_metrics,
        "conditioned_coverage_metrics": conditioned_metrics,
        "conditioned_polygon_boundary_metrics": conditioned_boundary_metrics,
        "raw_to_conditioned_difference_metrics": delta_metrics,
        "conditioning_diagnostics": diagnostics,
        "flat_mesh_quality": mesh_quality,
        "flat_mesh_quality_summary": mesh_metrics,
        "mesh_available": flat_mesh is not None,
        "timing_summary": timing_metrics,
        "timings_seconds": {name: round(value, 3) for name, value in timings.items()},
        "summary_only": args.summary_only,
        "artifacts": artifacts,
    }
    with (case_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(json_ready(summary), handle, indent=2)

    status_bits = [
        f"Case {number:03d}: {mode_name}",
        f"raw={len(raw_polygons)}",
        f"conditioned={len(conditioned_polygons)}",
    ]
    if case_status != "completed" and error_info is not None:
        status_bits.append(f"status={case_status}")
        status_bits.append(f"stage={error_info['stage']}")
        status_bits.append(f"error={error_info['type']}")
    print(" | ".join(status_bits) + f" -> {case_dir}")
    return summary


def main() -> int:
    args = parse_args()
    benchmark_start = time.perf_counter()
    git_root = args.git_root or repo_root()
    git_metadata = get_git_metadata(git_root)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    label = sanitize_label(args.label or git_metadata["branch"] or "unknown")
    mode_name, cleaner_module = resolve_mode(args.mode)

    summaries: list[dict[str, Any]] = []

    for index, number in enumerate(args.cases):
        summaries.append(
            run_case(
                number,
                args,
                git_metadata,
                mode_name=mode_name,
                cleaner_module=cleaner_module,
            )
        )
        if index < len(args.cases) - 1 and args.delay > 0:
            time.sleep(args.delay)

    total_elapsed = time.perf_counter() - benchmark_start
    results_path = args.output_dir / "results.json"
    save_results(
        results_path,
        {summary["case"]: summary for summary in summaries},
        metadata={
            "benchmark": "bench_footprints",
            "config": {
                "label": label,
                "mode": mode_name,
                "requested_mode": args.mode,
                "cases": args.cases,
                "max_mesh_size": args.max_mesh_size,
                "min_mesh_angle": args.min_mesh_angle,
                "min_building_detail": args.min_building_detail,
                "min_building_area": args.min_building_area,
                "merge_tolerance": args.merge_tolerance,
                "merge_buildings": not args.no_merge_buildings,
                "summary_only": args.summary_only,
                "plots_enabled": not args.no_plots,
            },
            "git": git_metadata,
        },
    )
    report = build_summary_report(
        summaries,
        label=label,
        mode_name=mode_name,
        elapsed_seconds=total_elapsed,
    )
    summary_path = args.output_dir / "summary.txt"
    summary_path.write_text(report + "\n", encoding="utf-8")

    overview_path: Path | None = None
    if not args.no_plots and not args.summary_only:
        overview_path = args.output_dir / "overview.png"
        plot_overview(summaries, overview_path, label=label, mode_name=mode_name)

    print()
    print(report)
    print()
    print(f"Results file: {results_path}")
    if overview_path is not None:
        print(f"Overview plot: {overview_path}")
    print(f"Summary text: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise
