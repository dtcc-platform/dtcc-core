from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import numpy as np

from dtcc_core.builder.geometry_builders import meshes as meshes_module
from dtcc_core.datasets._city_mesh_common import (
    CityMeshingFootprints,
    condition_city_meshing_footprints,
    prepare_city_from_bounds,
)
from dtcc_core.model import Bounds


@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets._city_mesh_common.City")
def test_prepare_city_from_bounds_default_path(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    city = Mock(name="city")

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_city_cls.return_value = city

    result = prepare_city_from_bounds(
        Bounds(0.0, 0.0, 1.0, 1.0),
        raster_cell_size=2.0,
        raster_radius=3.0,
        remove_outliers=True,
        outlier_threshold=3.0,
    )

    assert result is city
    mock_download_pointcloud.assert_called_once_with(bounds=ANY)
    mock_download_footprints.assert_called_once_with(bounds=ANY)
    raw_pointcloud.remove_global_outliers.assert_called_once_with(3.0)
    mock_build_terrain_raster.assert_called_once_with(
        filtered_pointcloud,
        cell_size=2.0,
        radius=3.0,
        ground_only=True,
    )
    mock_extract_roof_points.assert_called_once_with(buildings, filtered_pointcloud)
    mock_compute_building_heights.assert_called_once_with(
        roof_buildings,
        raster,
        overwrite=True,
    )
    city.add_terrain.assert_called_once_with(raster)
    city.add_buildings.assert_called_once_with(
        heighted_buildings,
        remove_outside_terrain=True,
    )


@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets._city_mesh_common.City")
def test_prepare_city_from_bounds_skips_outlier_removal_when_disabled(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_city_cls.return_value = Mock(name="city")

    prepare_city_from_bounds(
        Bounds(0.0, 0.0, 1.0, 1.0),
        raster_cell_size=2.0,
        raster_radius=3.0,
        remove_outliers=False,
        outlier_threshold=3.0,
    )

    raw_pointcloud.remove_global_outliers.assert_not_called()
    mock_build_terrain_raster.assert_called_once_with(
        raw_pointcloud,
        cell_size=2.0,
        radius=3.0,
        ground_only=True,
    )


@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.set_building_heights_from_attribute")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.flatten_terrain_raster")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets._city_mesh_common.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets._city_mesh_common.City")
def test_prepare_city_from_bounds_flat_ground_rebuilds_city_on_flat_raster(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_flatten_terrain_raster,
    mock_set_building_heights_from_attribute,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    flat_buildings = Mock(name="flat_buildings")
    flat_raster = Mock(name="flat_raster")
    city = Mock(name="city")

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    flat_raster.data = np.array([[17.5]])
    flat_raster.nodata = np.nan

    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_flatten_terrain_raster.return_value = flat_raster
    mock_set_building_heights_from_attribute.return_value = flat_buildings
    mock_city_cls.return_value = city

    result = prepare_city_from_bounds(
        Bounds(0.0, 0.0, 1.0, 1.0),
        raster_cell_size=2.0,
        raster_radius=3.0,
        remove_outliers=True,
        outlier_threshold=3.0,
        flat_ground=True,
        ground_level=17.5,
    )

    assert result is city
    mock_flatten_terrain_raster.assert_called_once_with(raster, height=17.5)
    mock_set_building_heights_from_attribute.assert_called_once_with(
        heighted_buildings,
        flat_raster,
        height_attribute="height",
        default_ground_height=17.5,
        always_use_default_ground=True,
    )
    city.add_terrain.assert_called_once_with(flat_raster)
    city.add_buildings.assert_called_once_with(
        flat_buildings,
        remove_outside_terrain=True,
    )


@patch(
    "dtcc_core.builder.geometry_builders.meshes._prepare_city_meshing_inputs"
)
@patch("dtcc_core.builder.geometry_builders.meshes._raise_stage_contract_errors")
def test_condition_city_meshing_footprints_uses_shared_meshing_stage(
    mock_raise_stage_contract_errors,
    mock_prepare_city_meshing_inputs,
):
    city = Mock(name="city")
    terrain = Mock(name="terrain")
    terrain_raster = Mock(name="terrain_raster")
    footprint = Mock(name="footprint_surface")
    polygon = Mock(name="footprint_polygon")
    footprint.to_polygon.return_value = polygon
    diagnostics = {"output_grid": 0.03125}
    contract = {"requirements": {"scale_contract_satisfied": True}}
    conditioned = meshes_module.ConditionedFootprints(
        surfaces=[footprint],
        source_map=[[1, 2]],
        subdomain_resolution=[7.5],
        diagnostics=diagnostics,
        declared_scale=0.5,
        contract=contract,
    )

    mock_prepare_city_meshing_inputs.return_value = (
        terrain,
        terrain_raster,
        conditioned,
    )

    result = condition_city_meshing_footprints(
        city,
        min_building_detail=0.5,
        min_building_area=15.0,
        merge_tolerance=0.5,
        merge_buildings=True,
        max_mesh_size=10.0,
        cleaning_diagnostics=True,
    )

    assert isinstance(result, CityMeshingFootprints)
    assert result.terrain is terrain
    assert result.terrain_raster is terrain_raster
    assert result.footprints == [footprint]
    assert result.source_map == [[1, 2]]
    assert result.subdomain_resolution == [7.5]
    assert result.diagnostics is diagnostics
    assert result.declared_scale == 0.5
    assert result.contract is contract
    assert result.polygons == [polygon]

    mock_prepare_city_meshing_inputs.assert_called_once_with(
        city,
        lod=None,
        min_building_detail=0.5,
        min_building_area=15.0,
        merge_tolerance=0.5,
        merge_buildings=True,
        max_mesh_size=10.0,
        cleaning_diagnostics=True,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
        pipeline_mode="strict",
    )
    mock_raise_stage_contract_errors.assert_called_once_with(
        "Conditioned footprints",
        contract,
    )
