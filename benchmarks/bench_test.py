"""
Benchmark the Helsingborg surface-mesh demo across legacy and current pipelines.

This reproduces the same input area and preprocessing flow as
`dtcc/demos/build_city_surface_mesh.py`, then compares:

- legacy Polyforge conditioning + legacy C++ surface builder
- legacy Polyforge conditioning + current surface builder
- current conditioning + legacy C++ surface builder
- current conditioning + current surface builder

The mixed scenarios make it easier to see whether regressions come from the
conditioning stage, the surface-meshing stage, or both.

Typical usage:
    /Users/logg/scratch/dtcc/venv/bin/python benchmarks/bench_test.py
    /Users/logg/scratch/dtcc/venv/bin/python benchmarks/bench_test.py --no-hybrids
    /Users/logg/scratch/dtcc/venv/bin/python benchmarks/bench_test.py --current-mesher triangle
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from polyforge import MergeStrategy, fix_clearance, merge_close_polygons, simplify_vwp
from polyforge.ops.clearance.protrusions import remove_narrow_wedges
from shapely.geometry import Polygon
from shapely.ops import unary_union

import dtcc_core
from dtcc_core.builder import _dtcc_builder
from dtcc_core.builder.geometry_builders.meshes import (
    _LOD_PRIORITY,
    _build_city_surface_mesh_from_ground_mesh,
    _build_ground_mesh_from_coverage,
    _condition_meshing_footprints,
    _normalize_lod_values,
    _prepare_surface_ground_regions,
    _promote_surface_shell_target_lods,
    _require_city_terrain_raster,
    _resolve_conditioned_target_lods,
    _split_ground_mesh_building_components,
    _surface_region_audit,
    _triangle_mesh_audit,
)
from dtcc_core.builder.model_conversion import (
    builder_mesh_to_mesh,
    create_builder_surface,
    raster_to_builder_gridfield,
)
from dtcc_core.model import Bounds, Building, City, GeometryType, Mesh, Surface
from dtcc_core.model.mixins.mesh.quality import triangle_mesh_quality

try:
    from _stockholm_common import format_console_table, json_ready
except ImportError:
    from benchmarks._stockholm_common import format_console_table, json_ready


DEMO_XMIN = 319_891.0
DEMO_YMIN = 6_399_790.0
DEMO_SIZE = 2_000.0
DEMO_ZMIN = 0.0
DEMO_ZMAX = 200.0

DEFAULT_MIN_BUILDING_DETAIL = 0.5
DEFAULT_MIN_BUILDING_AREA = 15.0
DEFAULT_MERGE_TOLERANCE = 0.5
DEFAULT_BUILDING_TRIANGLE_SIZE = 5.0
DEFAULT_MAX_MESH_SIZE = 10.0
DEFAULT_MIN_MESH_ANGLE = 25.0
DEFAULT_CURRENT_MESHER = "auto"
DEFAULT_OUTPUT_JSON = Path("benchmarks/output_surface_demo/bench_test_results.json")


@dataclass
class ConditioningCase:
    key: str
    label: str
    surfaces: list[Surface]
    source_map: list[list[int]]
    conditioned_resolution: list[float]
    target_lods: list[GeometryType]
    diagnostics: dict[str, Any]
    seconds: float


@dataclass
class Scenario:
    name: str
    conditioning_key: str
    conditioning_label: str
    builder_key: str
    builder_label: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-mesh-size",
        type=float,
        default=DEFAULT_MAX_MESH_SIZE,
        help="Maximum terrain triangle size for the demo run.",
    )
    parser.add_argument(
        "--building-mesh-triangle-size",
        type=float,
        default=DEFAULT_BUILDING_TRIANGLE_SIZE,
        help="Requested building-region triangle size.",
    )
    parser.add_argument(
        "--min-mesh-angle",
        type=float,
        default=DEFAULT_MIN_MESH_ANGLE,
        help="Minimum 2D mesh angle for both legacy and current builders.",
    )
    parser.add_argument(
        "--min-building-detail",
        type=float,
        default=DEFAULT_MIN_BUILDING_DETAIL,
        help="Minimum footprint feature size used by both conditioning pipelines.",
    )
    parser.add_argument(
        "--min-building-area",
        type=float,
        default=DEFAULT_MIN_BUILDING_AREA,
        help="Minimum footprint area kept after conditioning.",
    )
    parser.add_argument(
        "--merge-tolerance",
        type=float,
        default=DEFAULT_MERGE_TOLERANCE,
        help="Merge distance used by both conditioning pipelines.",
    )
    parser.add_argument(
        "--current-mesher",
        choices=("auto", "dtcc_mesher", "triangle", "spade"),
        default=DEFAULT_CURRENT_MESHER,
        help="2D backend for the current surface pipeline.",
    )
    parser.add_argument(
        "--no-merge-buildings",
        action="store_true",
        help="Disable building-footprint merging in both conditioning pipelines.",
    )
    parser.add_argument(
        "--no-hybrids",
        action="store_true",
        help="Only run legacy/legacy and current/current, skipping mixed scenarios.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=DEFAULT_OUTPUT_JSON,
        help="Optional JSON output path for the benchmark results.",
    )
    return parser.parse_args()


def demo_bounds() -> Bounds:
    return Bounds(
        DEMO_XMIN,
        DEMO_YMIN,
        DEMO_XMIN + DEMO_SIZE,
        DEMO_YMIN + DEMO_SIZE,
        DEMO_ZMIN,
        DEMO_ZMAX,
    )


def prepare_demo_city() -> tuple[City, dict[str, Any]]:
    bounds = demo_bounds()
    city = City()
    city.bounds = bounds

    start = time.perf_counter()
    city.download_pointcloud(bounds=bounds, filter_on_z_bounds=True)
    point_count = int(len(city.pointcloud.points)) if city.pointcloud is not None else 0
    city.download_footprints(bounds=bounds)
    raw_building_count = int(len(city.buildings))
    city.building_heights_from_pointcloud()
    prepare_seconds = time.perf_counter() - start

    terrain_raster = city.terrain.raster if city.terrain is not None else None
    height_count = 0
    for building in city.buildings:
        try:
            if building.height is not None and float(building.height) > 0.0:
                height_count += 1
        except (AttributeError, TypeError, ValueError):
            continue

    summary = {
        "bounds": {
            "xmin": float(bounds.xmin),
            "ymin": float(bounds.ymin),
            "xmax": float(bounds.xmax),
            "ymax": float(bounds.ymax),
            "zmin": float(bounds.zmin),
            "zmax": float(bounds.zmax),
        },
        "prepare_seconds": float(prepare_seconds),
        "point_count": point_count,
        "building_count": int(len(city.buildings)),
        "buildings_with_height": int(height_count),
        "terrain_raster_width": int(getattr(terrain_raster, "width", 0) or 0),
        "terrain_raster_height": int(getattr(terrain_raster, "height", 0) or 0),
        "raw_building_count": raw_building_count,
    }
    return city, summary


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
        composed.append(sorted(set(combined)))
    return composed


def _surface_from_building_lod(
    building: Building,
    lod_value: GeometryType,
) -> Surface | None:
    geometry = building.flatten_geometry(lod_value)
    if geometry is None:
        return None
    polygon = geometry.to_polygon(simplify=0.0)
    if polygon is None or polygon.is_empty:
        return None
    surface = Surface()
    surface.from_polygon(polygon, float(getattr(geometry.bounds, "zmax", 0.0)))
    return surface


def _legacy_merge_building_footprints(
    buildings: list[Building],
    lod_values: Sequence[GeometryType],
    *,
    max_distance: float,
    min_area: float,
) -> tuple[list[Building], list[list[int]]]:
    if len(buildings) <= 1:
        return list(buildings), [[index] for index in range(len(buildings))]

    source_indices: list[int] = []
    footprints: list[Polygon] = []
    building_heights: list[float] = []

    for index, (building, lod_value) in enumerate(zip(buildings, lod_values)):
        footprint_surface = _surface_from_building_lod(building, lod_value)
        if footprint_surface is None:
            continue
        footprint = footprint_surface.to_polygon(simplify=0.0)
        if footprint is None or footprint.is_empty:
            continue
        source_indices.append(index)
        building_heights.append(float(getattr(footprint_surface.bounds, "zmax", 0.0)))
        footprints.append(footprint)

    if not footprints:
        return [], []

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
        numerator = sum(building_heights[i] * footprints[i].area for i in local_indices)
        denominator = sum(footprints[i].area for i in local_indices)
        roof_z = numerator / denominator if denominator > 0 else 0.0

        surface = Surface()
        surface.from_polygon(footprint, roof_z)
        merged_building = Building()
        merged_building.add_geometry(surface, GeometryType.LOD0)
        merged_building.attributes = {"height": roof_z}
        merged_building.calculate_bounds()
        merged_buildings.append(merged_building)
        merged_indices_global.append(global_indices)

    return merged_buildings, merged_indices_global


def _legacy_fix_building_footprint_clearance(
    buildings: list[Building],
    *,
    clearance: float,
) -> tuple[list[Building], list[list[int]]]:
    fixed_buildings: list[Building] = []
    index_map: list[list[int]] = []

    for index, building in enumerate(buildings):
        geometry = building.flatten_geometry(GeometryType.LOD0)
        if geometry is None:
            continue
        footprint = geometry.to_polygon(simplify=0.0)
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
        surface.from_polygon(footprint, float(getattr(geometry.bounds, "zmax", 0.0)))
        fixed_building = building.copy()
        fixed_building.add_geometry(surface, GeometryType.LOD0)
        fixed_building.calculate_bounds()
        fixed_buildings.append(fixed_building)
        index_map.append([index])

    return fixed_buildings, index_map


def _legacy_simplify_building_footprints(
    buildings: list[Building],
    *,
    tolerance: float,
) -> tuple[list[Building], list[list[int]]]:
    simplified_buildings: list[Building] = []
    index_map: list[list[int]] = []

    for index, building in enumerate(buildings):
        geometry = building.flatten_geometry(GeometryType.LOD0)
        if geometry is None:
            continue
        footprint = geometry.to_polygon(simplify=0.0)
        if footprint is None or footprint.is_empty:
            continue
        footprint = simplify_vwp(footprint, tolerance)
        surface = Surface()
        surface.from_polygon(footprint, float(getattr(geometry.bounds, "zmax", 0.0)))
        simplified_building = building.copy()
        simplified_building.add_geometry(surface, GeometryType.LOD0)
        simplified_building.calculate_bounds()
        simplified_buildings.append(simplified_building)
        index_map.append([index])

    return simplified_buildings, index_map


def _resolution_from_buildings(
    buildings: Sequence[Building],
    *,
    max_mesh_size: float | None,
    min_building_detail: float,
) -> list[float]:
    resolution: list[float] = []
    for building in buildings:
        try:
            height = float(building.height)
        except (AttributeError, TypeError, ValueError):
            height = 0.0
        if height <= 0.0:
            if max_mesh_size is None:
                height = float(min_building_detail)
            else:
                height = float(max_mesh_size)
        if max_mesh_size is None:
            resolution.append(float(height))
        else:
            resolution.append(min(float(height), float(max_mesh_size)))
    return resolution


def run_legacy_conditioning(
    buildings: list[Building],
    *,
    lod_values: Sequence[GeometryType],
    merge_buildings: bool,
    merge_tolerance: float,
    min_building_area: float,
    min_building_detail: float,
    max_mesh_size: float | None,
) -> ConditioningCase:
    start = time.perf_counter()

    if merge_buildings:
        merged_buildings, merged_index_map = _legacy_merge_building_footprints(
            buildings,
            lod_values,
            max_distance=merge_tolerance,
            min_area=min_building_area,
        )
        cleared_buildings, cleared_index_map = _legacy_fix_building_footprint_clearance(
            merged_buildings,
            clearance=min_building_detail,
        )
        current_index_map = compose_index_map(merged_index_map, cleared_index_map)
        merged_again, merged_again_index_map = _legacy_merge_building_footprints(
            cleared_buildings,
            [GeometryType.LOD0] * len(cleared_buildings),
            max_distance=merge_tolerance,
            min_area=min_building_area,
        )
        current_index_map = compose_index_map(current_index_map, merged_again_index_map)
        simplified_buildings, simplified_index_map = _legacy_simplify_building_footprints(
            merged_again,
            tolerance=min_building_detail,
        )
        source_map = compose_index_map(current_index_map, simplified_index_map)
        conditioned_buildings = simplified_buildings
    else:
        footprint_buildings: list[Building] = []
        source_map = []
        for index, (building, lod_value) in enumerate(zip(buildings, lod_values)):
            surface = _surface_from_building_lod(building, lod_value)
            if surface is None:
                continue
            footprint_building = Building()
            footprint_building.add_geometry(surface, GeometryType.LOD0)
            footprint_building.attributes = dict(getattr(building, "attributes", {}) or {})
            try:
                footprint_building.attributes["height"] = float(building.height)
            except (AttributeError, TypeError, ValueError):
                pass
            footprint_building.calculate_bounds()
            footprint_buildings.append(footprint_building)
            source_map.append([index])
        conditioned_buildings, simplified_index_map = _legacy_simplify_building_footprints(
            footprint_buildings,
            tolerance=min_building_detail,
        )
        source_map = compose_index_map(source_map, simplified_index_map)

    surfaces: list[Surface] = []
    for building in conditioned_buildings:
        surface = building.get_footprint(GeometryType.LOD0)
        if surface is not None:
            surfaces.append(surface)

    target_lods = [
        min((lod_values[index] for index in indices), key=lambda value: _LOD_PRIORITY[value])
        if indices
        else GeometryType.LOD0
        for indices in source_map
    ]
    diagnostics = {
        "pipeline": "legacy_polyforge",
        "input_count": int(len(buildings)),
        "output_count": int(len(surfaces)),
        "merged_group_count": int(sum(1 for indices in source_map if len(indices) > 1)),
    }
    return ConditioningCase(
        key="legacy",
        label="legacy_polyforge",
        surfaces=surfaces,
        source_map=source_map,
        conditioned_resolution=_resolution_from_buildings(
            conditioned_buildings,
            max_mesh_size=max_mesh_size,
            min_building_detail=min_building_detail,
        ),
        target_lods=target_lods,
        diagnostics=diagnostics,
        seconds=float(time.perf_counter() - start),
    )


def run_current_conditioning(
    buildings: list[Building],
    *,
    lod: GeometryType | Sequence[GeometryType] | None,
    min_building_detail: float,
    min_building_area: float,
    merge_tolerance: float,
    merge_buildings: bool,
    max_mesh_size: float | None,
) -> ConditioningCase:
    start = time.perf_counter()
    surfaces, source_map, conditioned_resolution, diagnostics = _condition_meshing_footprints(
        buildings,
        lod=lod,
        min_building_detail=min_building_detail,
        min_building_area=min_building_area,
        merge_tolerance=merge_tolerance,
        merge_buildings=merge_buildings,
        max_mesh_size=max_mesh_size,
        cleaning_diagnostics=True,
    )
    seconds = time.perf_counter() - start
    return ConditioningCase(
        key="current",
        label="current_conditioning",
        surfaces=surfaces,
        source_map=source_map,
        conditioned_resolution=conditioned_resolution,
        target_lods=_resolve_conditioned_target_lods(buildings, lod, source_map),
        diagnostics=diagnostics,
        seconds=float(seconds),
    )


def _face_marker_stats(mesh: Mesh) -> dict[str, int]:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    markers_raw = getattr(mesh, "markers", None)
    if markers_raw is None or len(markers_raw) != len(faces):
        markers = np.empty(0, dtype=np.int64)
    else:
        markers = np.asarray(markers_raw, dtype=np.int64)

    if len(markers) == 0:
        return {
            "building_face_count": 0,
            "halo_face_count": 0,
            "ground_face_count": 0,
            "building_region_count": 0,
        }

    building_mask = markers >= 0
    return {
        "building_face_count": int(np.count_nonzero(building_mask)),
        "halo_face_count": int(np.count_nonzero(markers == -1)),
        "ground_face_count": int(np.count_nonzero(markers == -2)),
        "building_region_count": int(np.unique(markers[building_mask]).size),
    }


def _mesh_metrics(mesh: Mesh) -> dict[str, Any]:
    quality = triangle_mesh_quality(mesh.vertices, mesh.faces)
    audit = _triangle_mesh_audit(mesh)
    metrics = {
        "num_vertices": int(len(mesh.vertices)),
        "num_faces": int(len(mesh.faces)),
        "marker_count": int(audit.get("marker_count", 0)),
        "element_quality_min": float(quality["element_quality"]["min"]),
        "element_quality_mean": float(quality["element_quality"]["mean"]),
        "aspect_ratio_max": float(quality["aspect_ratio"]["max"]),
        "edge_ratio_max": float(quality["edge_ratio"]["max"]),
        "skewness_max": float(quality["skewness"]["max"]),
        "edge_length_min": float(audit.get("edge_length_min", float("nan"))),
        "edge_length_p01": float(audit.get("edge_length_p01", float("nan"))),
        "edge_length_p05": float(audit.get("edge_length_p05", float("nan"))),
        "area_min": float(audit.get("area_min", float("nan"))),
        "area_p01": float(audit.get("area_p01", float("nan"))),
        "area_p05": float(audit.get("area_p05", float("nan"))),
    }
    metrics.update(_face_marker_stats(mesh))
    return metrics


def build_surface_with_legacy_cpp(
    city: City,
    conditioned: ConditioningCase,
    *,
    building_mesh_triangle_size: float,
    max_mesh_size: float,
    min_mesh_angle: float,
    smoothing: int = 0,
    merge_meshes: bool = True,
    sort_triangles: bool = False,
    treat_lod0_as_holes: bool = False,
) -> tuple[Mesh, dict[str, Any]]:
    _terrain, terrain_raster = _require_city_terrain_raster(
        city,
        max_mesh_size=max_mesh_size,
    )
    builder_dem = raster_to_builder_gridfield(terrain_raster)
    default_priority = _LOD_PRIORITY[GeometryType.LOD3]
    target_lods = (
        list(conditioned.target_lods)
        if treat_lod0_as_holes
        else _promote_surface_shell_target_lods(conditioned.target_lods)
    )

    building_surfaces = []
    hole_surfaces = []
    building_directives = []
    building_resolution = []

    for surface, lod_value in zip(conditioned.surfaces, target_lods):
        builder_surface = create_builder_surface(surface)
        if treat_lod0_as_holes and lod_value == GeometryType.LOD0:
            hole_surfaces.append(builder_surface)
            continue
        building_surfaces.append(builder_surface)
        building_directives.append(_LOD_PRIORITY.get(lod_value, default_priority))
        building_resolution.append(float(building_mesh_triangle_size))

    start = time.perf_counter()
    builder_meshes = _dtcc_builder.build_city_surface_mesh(
        building_surfaces,
        hole_surfaces,
        building_directives,
        building_resolution,
        builder_dem,
        float(max_mesh_size),
        float(min_mesh_angle),
        int(smoothing),
        bool(merge_meshes),
        bool(sort_triangles),
    )
    build_seconds = time.perf_counter() - start

    if merge_meshes:
        mesh = builder_mesh_to_mesh(builder_meshes[0])
    else:
        raise ValueError("This benchmark expects merge_meshes=True.")

    return mesh, {
        "active_mesher": "legacy_cpp_internal",
        "build_seconds": float(build_seconds),
    }


def build_surface_with_current_pipeline(
    city: City,
    conditioned: ConditioningCase,
    *,
    building_mesh_triangle_size: float,
    max_mesh_size: float,
    min_mesh_angle: float,
    current_mesher: str,
    smoothing: int = 0,
    merge_meshes: bool = True,
    sort_triangles: bool = False,
    treat_lod0_as_holes: bool = False,
    min_building_detail: float,
) -> tuple[Mesh, dict[str, Any]]:
    _terrain, terrain_raster = _require_city_terrain_raster(
        city,
        max_mesh_size=max_mesh_size,
    )
    surface_mesh_bounds = (
        terrain_raster.bounds.xmin,
        terrain_raster.bounds.ymin,
        terrain_raster.bounds.xmax,
        terrain_raster.bounds.ymax,
    )
    base_resolution = [
        min(float(resolution), float(building_mesh_triangle_size))
        if building_mesh_triangle_size > 0.0
        else float(resolution)
        for resolution in conditioned.conditioned_resolution
    ]
    target_lods = (
        list(conditioned.target_lods)
        if treat_lod0_as_holes
        else _promote_surface_shell_target_lods(conditioned.target_lods)
    )

    region_prep_start = time.perf_counter()
    (
        building_surfaces,
        building_lod_switches,
        region_polygons,
        region_markers,
        region_triangle_sizes,
        region_points,
    ) = _prepare_surface_ground_regions(
        conditioned_surfaces=conditioned.surfaces,
        conditioned_resolution=base_resolution,
        target_lods=target_lods,
        bounds=surface_mesh_bounds,
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        footprint_diagnostics=conditioned.diagnostics,
        cleaning_diagnostics=True,
        treat_lod0_as_holes=treat_lod0_as_holes,
    )
    region_prep_seconds = time.perf_counter() - region_prep_start
    region_audit = _surface_region_audit(
        building_surfaces=building_surfaces,
        region_polygons=region_polygons,
        region_markers=region_markers,
        region_triangle_sizes=region_triangle_sizes,
    )

    ground_start = time.perf_counter()
    ground_mesh, active_mesher = _build_ground_mesh_from_coverage(
        region_polygons=region_polygons,
        region_markers=region_markers,
        region_points=region_points,
        bounds=surface_mesh_bounds,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        mesher=current_mesher,
        sort_triangles=sort_triangles,
        region_triangle_sizes=region_triangle_sizes,
        add_halo_markers=False,
    )
    ground_seconds = time.perf_counter() - ground_start
    ground_audit = _triangle_mesh_audit(ground_mesh)

    extrusion_start = time.perf_counter()
    ground_mesh, building_surfaces, building_lod_switches = _split_ground_mesh_building_components(
        ground_mesh=ground_mesh,
        building_surfaces=building_surfaces,
        meshing_directives=building_lod_switches,
    )
    mesh = _build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=terrain_raster,
        building_surfaces=building_surfaces,
        meshing_directives=building_lod_switches,
        smoothing=smoothing,
        merge_meshes=merge_meshes,
    )
    extrusion_seconds = time.perf_counter() - extrusion_start

    if not isinstance(mesh, Mesh):
        raise ValueError("This benchmark expects merge_meshes=True.")

    return mesh, {
        "active_mesher": active_mesher,
        "region_prep_seconds": float(region_prep_seconds),
        "ground_mesh_seconds": float(ground_seconds),
        "surface_extrusion_seconds": float(extrusion_seconds),
        "build_seconds": float(region_prep_seconds + ground_seconds + extrusion_seconds),
        "region_audit": region_audit,
        "ground_mesh_audit": ground_audit,
    }


def build_scenarios(include_hybrids: bool) -> list[Scenario]:
    scenarios = [
        Scenario(
            name="legacy_legacy",
            conditioning_key="legacy",
            conditioning_label="legacy_polyforge",
            builder_key="legacy_cpp",
            builder_label="legacy_cpp_surface",
        ),
        Scenario(
            name="current_current",
            conditioning_key="current",
            conditioning_label="current_conditioning",
            builder_key="current_surface",
            builder_label="current_surface_pipeline",
        ),
    ]
    if include_hybrids:
        scenarios.insert(
            1,
            Scenario(
                name="legacy_current",
                conditioning_key="legacy",
                conditioning_label="legacy_polyforge",
                builder_key="current_surface",
                builder_label="current_surface_pipeline",
            ),
        )
        scenarios.insert(
            2,
            Scenario(
                name="current_legacy",
                conditioning_key="current",
                conditioning_label="current_conditioning",
                builder_key="legacy_cpp",
                builder_label="legacy_cpp_surface",
            ),
        )
    return scenarios


def run_scenario(
    scenario: Scenario,
    city: City,
    conditioning_cases: dict[str, ConditioningCase],
    *,
    building_mesh_triangle_size: float,
    max_mesh_size: float,
    min_mesh_angle: float,
    current_mesher: str,
    min_building_detail: float,
) -> dict[str, Any]:
    conditioned = conditioning_cases[scenario.conditioning_key]
    try:
        if scenario.builder_key == "legacy_cpp":
            mesh, build_info = build_surface_with_legacy_cpp(
                city,
                conditioned,
                building_mesh_triangle_size=building_mesh_triangle_size,
                max_mesh_size=max_mesh_size,
                min_mesh_angle=min_mesh_angle,
            )
        elif scenario.builder_key == "current_surface":
            mesh, build_info = build_surface_with_current_pipeline(
                city,
                conditioned,
                building_mesh_triangle_size=building_mesh_triangle_size,
                max_mesh_size=max_mesh_size,
                min_mesh_angle=min_mesh_angle,
                current_mesher=current_mesher,
                min_building_detail=min_building_detail,
            )
        else:
            raise ValueError(f"Unknown builder key: {scenario.builder_key}")
    except Exception as exc:
        message = str(exc)
        status = "failed"
        if (
            scenario.builder_key == "legacy_cpp"
            and "No triangulation backend is available in this build." in message
        ):
            status = "unavailable"
        return {
            "name": scenario.name,
            "conditioning": scenario.conditioning_label,
            "builder": scenario.builder_label,
            "status": status,
            "conditioning_seconds": float(conditioned.seconds),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }

    metrics = _mesh_metrics(mesh)
    total_seconds = float(conditioned.seconds + build_info["build_seconds"])
    return {
        "name": scenario.name,
        "conditioning": scenario.conditioning_label,
        "builder": scenario.builder_label,
        "status": "success",
        "conditioning_seconds": float(conditioned.seconds),
        "build_seconds": float(build_info["build_seconds"]),
        "total_seconds": total_seconds,
        "active_mesher": build_info.get("active_mesher"),
        "conditioned_footprints": int(len(conditioned.surfaces)),
        "source_group_count": int(len(conditioned.source_map)),
        "conditioning_diagnostics": conditioned.diagnostics,
        "builder_info": build_info,
        "metrics": metrics,
    }


def print_city_summary(city_summary: dict[str, Any]) -> None:
    rows = [[
        f"{city_summary['bounds']['xmin']:.0f}",
        f"{city_summary['bounds']['ymin']:.0f}",
        f"{city_summary['bounds']['xmax']:.0f}",
        f"{city_summary['bounds']['ymax']:.0f}",
        city_summary["point_count"],
        city_summary["building_count"],
        city_summary["buildings_with_height"],
        f"{city_summary['prepare_seconds']:.2f}s",
    ]]
    print(
        format_console_table(
            [
                "xmin",
                "ymin",
                "xmax",
                "ymax",
                "Points",
                "Buildings",
                "Heights",
                "Prep",
            ],
            rows,
            title="Demo Input",
        )
    )


def print_conditioning_summary(conditioning_cases: Sequence[ConditioningCase]) -> None:
    rows = []
    for case in conditioning_cases:
        rows.append([
            case.label,
            len(case.surfaces),
            sum(1 for indices in case.source_map if len(indices) > 1),
            case.diagnostics.get("output_grid", "-"),
            case.diagnostics.get("mesher_regularized_polygon_count", "-"),
            f"{case.seconds:.2f}s",
        ])
    print()
    print(
        format_console_table(
            [
                "Conditioning",
                "Footprints",
                "Merged groups",
                "Output grid",
                "Mesher regularized",
                "Time",
            ],
            rows,
            title="Conditioning Summary",
        )
    )


def print_result_summary(results: Sequence[dict[str, Any]]) -> None:
    rows = []
    for result in results:
        if result["status"] != "success":
            rows.append([
                result["name"],
                result["conditioning"],
                result["builder"],
                result.get("active_mesher", "-"),
                result["status"],
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
            ])
            continue
        metrics = result["metrics"]
        rows.append([
            result["name"],
            result["conditioning"],
            result["builder"],
            result.get("active_mesher", "-"),
            f"{result['conditioning_seconds']:.2f}s",
            f"{result['build_seconds']:.2f}s",
            f"{result['total_seconds']:.2f}s",
            metrics["num_vertices"],
            metrics["num_faces"],
            result["conditioned_footprints"],
            metrics["building_region_count"],
        ])
    print()
    print(
        format_console_table(
            [
                "Scenario",
                "Conditioning",
                "Builder",
                "Mesher",
                "Cond",
                "Build",
                "Total",
                "Vertices",
                "Faces",
                "Footprints",
                "Regions",
            ],
            rows,
            title="Scenario Summary",
        )
    )

    quality_rows = []
    counts_rows = []
    for result in results:
        if result["status"] != "success":
            quality_rows.append([
                result["name"],
                result["status"],
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
            ])
            counts_rows.append([
                result["name"],
                result["status"],
                "-",
                "-",
                "-",
                "-",
                "-",
            ])
            continue
        metrics = result["metrics"]
        quality_rows.append([
            result["name"],
            "ok",
            f"{metrics['element_quality_min']:.4f}",
            f"{metrics['element_quality_mean']:.4f}",
            f"{metrics['aspect_ratio_max']:.2f}",
            f"{metrics['edge_ratio_max']:.2f}",
            f"{metrics['skewness_max']:.2f}",
            f"{metrics['edge_length_p01']:.3f}",
            f"{metrics['area_p01']:.3f}",
        ])
        counts_rows.append([
            result["name"],
            "ok",
            metrics["building_face_count"],
            metrics["ground_face_count"],
            metrics["halo_face_count"],
            metrics["marker_count"],
            result["conditioned_footprints"],
        ])

    print()
    print(
        format_console_table(
            [
                "Scenario",
                "Status",
                "EQ min",
                "EQ mean",
                "AR max",
                "ER max",
                "Skew max",
                "Edge p01",
                "Area p01",
            ],
            quality_rows,
            title="Mesh Quality",
        )
    )
    print()
    print(
        format_console_table(
            [
                "Scenario",
                "Status",
                "Building faces",
                "Ground faces",
                "Halo faces",
                "Markers",
                "Footprints",
            ],
            counts_rows,
            title="Marker / Region Counts",
        )
    )


def save_results(
    path: Path | None,
    *,
    args: argparse.Namespace,
    city_summary: dict[str, Any],
    conditioning_cases: dict[str, ConditioningCase],
    results: list[dict[str, Any]],
) -> None:
    if path is None:
        return

    payload = {
        "config": {
            "max_mesh_size": float(args.max_mesh_size),
            "building_mesh_triangle_size": float(args.building_mesh_triangle_size),
            "min_mesh_angle": float(args.min_mesh_angle),
            "min_building_detail": float(args.min_building_detail),
            "min_building_area": float(args.min_building_area),
            "merge_tolerance": float(args.merge_tolerance),
            "merge_buildings": bool(not args.no_merge_buildings),
            "current_mesher": args.current_mesher,
            "include_hybrids": bool(not args.no_hybrids),
        },
        "city": city_summary,
        "conditioning": {
            key: {
                "label": case.label,
                "seconds": case.seconds,
                "conditioned_footprints": len(case.surfaces),
                "source_group_count": len(case.source_map),
                "target_lods": [lod.name for lod in case.target_lods],
                "diagnostics": case.diagnostics,
            }
            for key, case in conditioning_cases.items()
        },
        "results": results,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2))


def main() -> int:
    args = parse_args()
    merge_buildings = not args.no_merge_buildings

    city, city_summary = prepare_demo_city()
    buildings = list(city.buildings)
    lod = None
    lod_values = _normalize_lod_values(buildings, lod)

    conditioning_cases = {
        "legacy": run_legacy_conditioning(
            buildings,
            lod_values=lod_values,
            merge_buildings=merge_buildings,
            merge_tolerance=args.merge_tolerance,
            min_building_area=args.min_building_area,
            min_building_detail=args.min_building_detail,
            max_mesh_size=args.max_mesh_size,
        ),
        "current": run_current_conditioning(
            buildings,
            lod=lod,
            min_building_detail=args.min_building_detail,
            min_building_area=args.min_building_area,
            merge_tolerance=args.merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=args.max_mesh_size,
        ),
    }

    scenarios = build_scenarios(include_hybrids=not args.no_hybrids)
    results = [
        run_scenario(
            scenario,
            city,
            conditioning_cases,
            building_mesh_triangle_size=args.building_mesh_triangle_size,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            current_mesher=args.current_mesher,
            min_building_detail=args.min_building_detail,
        )
        for scenario in scenarios
    ]

    print_city_summary(city_summary)
    print_conditioning_summary(list(conditioning_cases.values()))
    print_result_summary(results)

    if args.output_json is not None:
        save_results(
            args.output_json,
            args=args,
            city_summary=city_summary,
            conditioning_cases=conditioning_cases,
            results=results,
        )
        print()
        print(f"Results written to {args.output_json.resolve()}")

    failed = any(result["status"] == "failed" for result in results)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
