from __future__ import annotations

from contextlib import nullcontext
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

import dtcc_core
from dtcc_core.model import Bounds, City, GeometryType, Surface


@dataclass(frozen=True)
class CityMeshingFootprints:
    """Meshing-ready conditioned footprints prepared from a city model."""

    terrain: Any
    terrain_raster: Any
    footprints: list[Surface]
    source_map: list[list[int]]
    subdomain_resolution: list[float]
    diagnostics: dict[str, Any]
    contract: dict[str, Any]

    @property
    def polygons(self):
        return [footprint.to_polygon(simplify=0.0) for footprint in self.footprints]


def ground_level_from_raster(raster) -> float:
    valid_mask = np.isfinite(raster.data)
    if not np.isnan(raster.nodata):
        valid_mask &= raster.data != raster.nodata
    return float(raster.data[valid_mask].min())


def prepare_city_from_bounds(
    bounds: Bounds,
    *,
    raster_cell_size: float,
    raster_radius: float,
    remove_outliers: bool,
    outlier_threshold: float,
    flat_ground: bool = False,
    ground_level: float | None = None,
    progress: Any | None = None,
) -> City:
    """Prepare a meshing-ready city model from live data sources.

    This helper centralizes the public mesh-dataset preparation path so the
    end-user datasets and the benchmark scripts can share the same
    bounds-to-City implementation.
    """

    def phase(name: str, message: str):
        if progress is None:
            return nullcontext()
        return progress.phase(name, message)

    with phase("download_pointcloud", "Downloading point cloud data..."):
        pointcloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)

    with phase("download_footprints", "Downloading building footprints..."):
        buildings = dtcc_core.io.data.download_footprints(bounds=bounds)

    with phase(
        "remove_outliers",
        (
            "Removing outliers..."
            if remove_outliers
            else "Skipping outlier removal..."
        ),
    ):
        if remove_outliers:
            pointcloud = pointcloud.remove_global_outliers(outlier_threshold)

    with phase("build_terrain", "Building terrain raster..."):
        raster = dtcc_core.builder.build_terrain_raster(
            pointcloud,
            cell_size=raster_cell_size,
            radius=raster_radius,
            ground_only=True,
        )

    with phase("extract_roof_points", "Extracting roof points..."):
        buildings = dtcc_core.builder.extract_roof_points(buildings, pointcloud)

    with phase("compute_building_heights", "Computing building heights..."):
        buildings = dtcc_core.builder.compute_building_heights(
            buildings,
            raster,
            overwrite=True,
        )

    with phase(
        "build_city",
        (
            "Preparing flat-ground city model..."
            if flat_ground
            else "Assembling city model..."
        ),
    ):
        if flat_ground:
            raster = dtcc_core.builder.flatten_terrain_raster(
                raster,
                height=ground_level,
            )
            buildings = dtcc_core.builder.set_building_heights_from_attribute(
                buildings,
                raster,
                height_attribute="height",
                default_ground_height=ground_level_from_raster(raster),
                always_use_default_ground=True,
            )

        city = City()
        city.add_terrain(raster)
        city.add_buildings(buildings, remove_outside_terrain=True)

    return city


def condition_city_meshing_footprints(
    city: City,
    *,
    lod: GeometryType | Sequence[GeometryType] | None = None,
    min_building_detail: float,
    min_building_area: float,
    merge_tolerance: float,
    merge_buildings: bool,
    max_mesh_size: float | None,
    cleaning_diagnostics: bool,
    show_footprints: bool = False,
    footprint_cleaning_plot_block: bool = True,
    pipeline_mode: str = "strict",
) -> CityMeshingFootprints:
    """Run the same meshing-footprint preparation stage used by mesh builders."""

    import dtcc_core.builder.geometry_builders.meshes as mesh_builders

    (
        terrain,
        terrain_raster,
        building_footprints,
        conditioned_source_map,
        subdomain_resolution,
        diagnostics,
    ) = mesh_builders._prepare_city_meshing_inputs(
        city,
        lod=lod,
        min_building_detail=min_building_detail,
        min_building_area=min_building_area,
        merge_tolerance=merge_tolerance,
        merge_buildings=merge_buildings,
        max_mesh_size=max_mesh_size,
        cleaning_diagnostics=cleaning_diagnostics,
        show_footprints=show_footprints,
        footprint_cleaning_plot_block=footprint_cleaning_plot_block,
        pipeline_mode=pipeline_mode,
    )

    contract = mesh_builders._conditioned_footprint_contract(
        surfaces=building_footprints,
        min_building_detail=min_building_detail,
        diagnostics=diagnostics,
    )
    mesh_builders._raise_stage_contract_errors("Conditioned footprints", contract)

    return CityMeshingFootprints(
        terrain=terrain,
        terrain_raster=terrain_raster,
        footprints=building_footprints,
        source_map=conditioned_source_map,
        subdomain_resolution=subdomain_resolution,
        diagnostics=diagnostics,
        contract=contract,
    )
