"""Tests for the city flat mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_flat_mesh import CityFlatMeshArgs, CityFlatMeshDataset
from dtcc_core.model import GeometryType


def test_city_flat_mesh_registered_name():
    """The dataset class should expose the correct registration name."""
    assert CityFlatMeshDataset().name == "city_flat_mesh"


def test_city_flat_mesh_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("city_flat_mesh")
    assert ds is not None
    assert ds.name == "city_flat_mesh"


def test_city_flat_mesh_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "city_flat_mesh")
    assert callable(datasets.city_flat_mesh)


@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_flat_mesh.City")
def test_city_flat_mesh_default_build_returns_flat_mesh(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_flat_mesh,
):
    """The default wrapper path should return the built flat mesh."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    flat_mesh = Mock(name="flat_mesh")

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_build_city_flat_mesh.return_value = flat_mesh

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityFlatMeshDataset()
    result = dataset.build(CityFlatMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is flat_mesh
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
    mock_build_city_flat_mesh.assert_called_once_with(
        city,
        lod=GeometryType.LOD0,
        max_mesh_size=10.0,
        min_mesh_angle=25.0,
        min_building_detail=0.5,
        min_building_area=15.0,
        merge_buildings=True,
    )


@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_flat_mesh.City")
def test_city_flat_mesh_parameter_plumbing(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_flat_mesh,
):
    """Flat-mesh arguments should be forwarded to the builder call."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_flat_mesh.return_value = Mock(name="flat_mesh")

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityFlatMeshDataset()
    dataset.build(
        CityFlatMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            max_mesh_size=14.0,
            min_mesh_angle=29.0,
            min_building_detail=0.85,
            min_building_area=19.0,
            merge_buildings=False,
        )
    )

    mock_build_city_flat_mesh.assert_called_once_with(
        city,
        lod=GeometryType.LOD0,
        max_mesh_size=14.0,
        min_mesh_angle=29.0,
        min_building_detail=0.85,
        min_building_area=19.0,
        merge_buildings=False,
    )


@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_flat_mesh.City")
def test_city_flat_mesh_obj_export_returns_bytes(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_flat_mesh,
):
    """Export formats should return bytes via export_to_bytes()."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = Mock(name="filtered_pointcloud")
    flat_mesh = Mock(name="flat_mesh")

    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_flat_mesh.return_value = flat_mesh
    mock_city_cls.return_value = Mock(name="city")

    dataset = CityFlatMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"flat-mesh") as mock_export:
        result = dataset.build(
            CityFlatMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"flat-mesh"
    mock_export.assert_called_once_with(flat_mesh, "obj")
