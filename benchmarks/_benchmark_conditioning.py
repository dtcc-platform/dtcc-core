from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator
from unittest.mock import patch

from dtcc_core.logging import info
from dtcc_core.model import Building, GeometryType, Surface

try:
    import bench_footprints
except ImportError:
    from benchmarks import bench_footprints


CONDITIONING_MODES = ("new", "legacy")


def _resolved_benchmark_lod(
    buildings: list[Building],
    lod: GeometryType | list[GeometryType] | None,
) -> GeometryType:
    import dtcc_core.builder.geometry_builders.meshes as mesh_builders

    lod_values = mesh_builders._normalize_lod_values(buildings, lod)
    unique_lods = {lod_value for lod_value in lod_values}
    if len(unique_lods) != 1:
        raise ValueError(
            "Benchmark legacy conditioning requires a single resolved lod value "
            "for all buildings in the case."
        )
    return lod_values[0]


def _legacy_conditioned_meshing_footprints(
    buildings: list[Building],
    *,
    lod: GeometryType | list[GeometryType] | None,
    min_building_detail: float,
    min_building_area: float,
    merge_tolerance: float,
    merge_buildings: bool,
    max_mesh_size: float | None,
    cleaning_diagnostics: bool = True,
    pipeline_mode: str = "strict",
):
    import dtcc_core.builder.geometry_builders.meshes as mesh_builders

    if not buildings:
        return [], [], [], {
            "conditioning_mode": "legacy",
            "pipeline_mode": pipeline_mode,
            "output_grid": float(min_building_detail),
            "output_count": 0,
        }

    resolved_lod = _resolved_benchmark_lod(buildings, lod)
    normalized_mesh_size = mesh_builders._normalize_max_mesh_size(max_mesh_size)

    source_areas = [0.0] * len(buildings)
    source_heights: list[float | None] = [None] * len(buildings)
    source_roof_z = [0.0] * len(buildings)
    for index, building in enumerate(buildings):
        _polygon, area, height, roof_z = mesh_builders._extract_meshing_polygon(
            building,
            resolved_lod,
        )
        source_areas[index] = area
        source_heights[index] = height
        source_roof_z[index] = roof_z

    conditioned_polygons, source_map, diagnostics = bench_footprints.run_legacy_conditioning(
        buildings,
        lod=resolved_lod,
        merge_buildings=merge_buildings,
        merge_tolerance=merge_tolerance,
        min_building_area=min_building_area,
        min_building_detail=min_building_detail,
    )

    conditioned_surfaces: list[Surface] = []
    conditioned_source_map: list[list[int]] = []
    subdomain_resolution: list[float] = []
    conservative_roof_count = 0
    conservative_roof_max_span = 0.0

    for polygon, source_indices in zip(conditioned_polygons, source_map):
        height_default = normalized_mesh_size or float(min_building_detail)
        roof_z, height, conservative_roof, roof_z_span = (
            mesh_builders._resolve_merged_group_roof_metadata(
                source_indices,
                source_areas=source_areas,
                source_roof_z=source_roof_z,
                source_heights=source_heights,
                default_height=height_default,
            )
        )
        if conservative_roof:
            conservative_roof_count += 1
            conservative_roof_max_span = max(conservative_roof_max_span, roof_z_span)

        surface = Surface()
        surface.from_polygon(polygon, roof_z)
        conditioned_surfaces.append(surface)

        unique_sources = sorted(set(source_indices))
        conditioned_source_map.append(unique_sources)
        if normalized_mesh_size is None:
            subdomain_resolution.append(height)
        else:
            subdomain_resolution.append(min(height, normalized_mesh_size))

    diagnostics = dict(diagnostics)
    diagnostics.setdefault("output_grid", float(min_building_detail))
    diagnostics["conditioning_mode"] = "legacy"
    diagnostics["pipeline_mode"] = pipeline_mode
    diagnostics["conservative_merged_roof_count"] = conservative_roof_count
    diagnostics["conservative_merged_roof_max_span"] = conservative_roof_max_span
    diagnostics["mesher_ready_coverage_revalidation_enabled"] = False
    diagnostics["mesher_ready_coverage_revalidation_attempted"] = False

    if cleaning_diagnostics:
        info(
            "Legacy meshing footprint conditioning complete: "
            f"{len(buildings)} buildings -> {len(conditioned_surfaces)} footprints, "
            f"output_grid={diagnostics.get('output_grid')} m, "
            f"conservative_merged_roofs={conservative_roof_count}."
        )

    return (
        conditioned_surfaces,
        conditioned_source_map,
        subdomain_resolution,
        diagnostics,
    )


@contextmanager
def benchmark_conditioning_override(mode_name: str) -> Iterator[None]:
    if mode_name != "legacy":
        yield
        return

    import dtcc_core.builder.geometry_builders.meshes as mesh_builders

    with patch.object(
        mesh_builders,
        "_condition_meshing_footprints",
        _legacy_conditioned_meshing_footprints,
    ):
        yield
