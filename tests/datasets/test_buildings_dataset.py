"""Tests for the buildings dataset wrapper."""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.buildings import BuildingArgs, BuildingDataset
from dtcc_core.model import Building, BuildingCollection


def test_buildings_registered_name():
    """The dataset class should expose the correct registration name."""
    assert BuildingDataset().name == "buildings"


def test_buildings_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("buildings")
    assert ds is not None
    assert ds.name == "buildings"


def test_buildings_module_attribute():
    """The dataset should be exposed as a callable module attribute."""
    assert hasattr(datasets, "buildings")
    assert callable(datasets.buildings)


def test_buildings_prepare_result_returns_building_collection():
    """Public no-format calls adapt raw building lists to a semantic collection."""
    buildings = [Building(), Building()]
    dataset = BuildingDataset()
    args = BuildingArgs(bounds=(0.0, 0.0, 1.0, 1.0))

    result = dataset.prepare_result(buildings, args)

    assert isinstance(result, BuildingCollection)
    assert list(result) == buildings
    assert result[0] is buildings[0]


@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_lod1_buildings")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.buildings.City")
def test_buildings_default_build_returns_lod1_buildings(
    mock_city_cls,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_lod1_buildings,
):
    """The default wrapper path should return the built LoD1 buildings list."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    terrain_raster = Mock(name="terrain_raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    lod1_buildings = [Mock(name="lod1_building_1"), Mock(name="lod1_building_2")]
    initial_buildings = [Mock(name="initial_building")]

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_build_terrain_raster.return_value = terrain_raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_build_lod1_buildings.return_value = lod1_buildings

    city = Mock(name="city")
    city.pointcloud = raw_pointcloud
    city.buildings = initial_buildings

    def add_point_cloud_side_effect(pc):
        city.pointcloud = pc

    city.add_point_cloud.side_effect = add_point_cloud_side_effect
    mock_city_cls.return_value = city

    dataset = BuildingDataset()
    result = dataset.build(BuildingArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is lod1_buildings
    city.download_pointcloud.assert_called_once_with()
    city.download_footprints.assert_called_once_with()
    raw_pointcloud.remove_global_outliers.assert_called_once_with(3.0)
    city.add_point_cloud.assert_called_once_with(filtered_pointcloud)
    mock_build_terrain_raster.assert_called_once_with(
        filtered_pointcloud,
        cell_size=2.0,
        ground_only=True,
        _report_progress=False,
    )
    mock_extract_roof_points.assert_called_once_with(initial_buildings, filtered_pointcloud)
    mock_compute_building_heights.assert_called_once_with(
        roof_buildings,
        terrain_raster,
        overwrite=True,
    )
    mock_build_lod1_buildings.assert_called_once_with(heighted_buildings)


@patch("dtcc_core.datasets.buildings.dtcc_core.builder.meshing.merge_meshes")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_lod1_buildings")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.buildings.City")
def test_buildings_obj_export_returns_bytes(
    mock_city_cls,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_lod1_buildings,
    mock_merge_meshes,
):
    """The export branch should merge building meshes and return bytes."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    terrain_raster = Mock(name="terrain_raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    merged_mesh = Mock(name="merged_mesh")

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_build_terrain_raster.return_value = terrain_raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_merge_meshes.return_value = merged_mesh

    building_mesh_1 = Mock(name="building_mesh_1")
    building_mesh_2 = Mock(name="building_mesh_2")
    lod1_building_1 = Mock(name="lod1_building_1")
    lod1_building_2 = Mock(name="lod1_building_2")
    lod1_building_1.lod1.mesh.return_value = building_mesh_1
    lod1_building_2.lod1.mesh.return_value = building_mesh_2
    mock_build_lod1_buildings.return_value = [lod1_building_1, lod1_building_2]

    city = Mock(name="city")
    city.pointcloud = raw_pointcloud
    city.buildings = [Mock(name="initial_building")]
    city.add_point_cloud.side_effect = lambda pc: setattr(city, "pointcloud", pc)
    mock_city_cls.return_value = city

    dataset = BuildingDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"obj-bytes") as mock_export:
        result = dataset.build(
            BuildingArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"obj-bytes"
    lod1_building_1.lod1.mesh.assert_called_once_with(weld=True, snap=0.005)
    lod1_building_2.lod1.mesh.assert_called_once_with(weld=True, snap=0.005)
    mock_merge_meshes.assert_called_once_with([building_mesh_1, building_mesh_2])
    mock_export.assert_called_once_with(merged_mesh, "obj")


@patch("dtcc_core.datasets.buildings.dtcc_core.builder.meshing.merge_meshes")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_lod1_buildings")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.buildings.City")
def test_buildings_place_on_zero_offsets_meshes_before_merge(
    mock_city_cls,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_lod1_buildings,
    mock_merge_meshes,
):
    """place_on_zero=True should offset each building mesh before merge."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_build_terrain_raster.return_value = Mock(name="terrain_raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_merge_meshes.return_value = Mock(name="merged_mesh")

    building_mesh_1 = Mock(name="building_mesh_1")
    building_mesh_1.bounds = SimpleNamespace(zmin=2.5)
    building_mesh_2 = Mock(name="building_mesh_2")
    building_mesh_2.bounds = SimpleNamespace(zmin=-1.0)
    lod1_building_1 = Mock(name="lod1_building_1")
    lod1_building_2 = Mock(name="lod1_building_2")
    lod1_building_1.lod1.mesh.return_value = building_mesh_1
    lod1_building_2.lod1.mesh.return_value = building_mesh_2
    mock_build_lod1_buildings.return_value = [lod1_building_1, lod1_building_2]

    city = Mock(name="city")
    city.pointcloud = raw_pointcloud
    city.buildings = [Mock(name="initial_building")]
    city.add_point_cloud.side_effect = lambda pc: setattr(city, "pointcloud", pc)
    mock_city_cls.return_value = city

    dataset = BuildingDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"stl-bytes"):
        dataset.build(
            BuildingArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="stl",
                place_on_zero=True,
            )
        )

    building_mesh_1.offset.assert_called_once_with([0, 0, -2.5])
    building_mesh_2.offset.assert_called_once_with([0, 0, 1.0])
