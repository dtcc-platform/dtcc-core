"""Tests for the city surface mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import numpy as np
import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_surface_mesh import (
    CitySurfaceMeshArgs,
    CitySurfaceMeshDataset,
)


def test_city_surface_mesh_registered_name():
    """The dataset class should expose the correct registration name."""
    assert CitySurfaceMeshDataset().name == "city_surface_mesh"


def test_city_surface_mesh_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("city_surface_mesh")
    assert ds is not None
    assert ds.name == "city_surface_mesh"


def test_city_surface_mesh_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "city_surface_mesh")
    assert callable(datasets.city_surface_mesh)


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_surface_mesh.City")
def test_city_surface_mesh_default_build_returns_surface_mesh(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_surface_mesh,
):
    """The default wrapper path should return the built surface mesh."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    surface_mesh = Mock(name="surface_mesh")

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_build_city_surface_mesh.return_value = surface_mesh

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CitySurfaceMeshDataset()
    result = dataset.build(CitySurfaceMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is surface_mesh
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
    mock_build_city_surface_mesh.assert_called_once_with(
        city,
        max_mesh_size=10.0,
        min_mesh_angle=25.0,
        min_building_detail=0.5,
        min_building_area=15.0,
        merge_buildings=True,
        smoothing=0,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
    )


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_surface_mesh.City")
def test_city_surface_mesh_remove_outliers_false_skips_filtering(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_surface_mesh,
):
    """Outlier removal should be skipped when remove_outliers=False."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_surface_mesh.return_value = Mock(name="surface_mesh")
    mock_city_cls.return_value = Mock(name="city")

    dataset = CitySurfaceMeshDataset()
    dataset.build(
        CitySurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            remove_outliers=False,
        )
    )

    raw_pointcloud.remove_global_outliers.assert_not_called()
    mock_build_terrain_raster.assert_called_once_with(
        raw_pointcloud,
        cell_size=2.0,
        radius=3.0,
        ground_only=True,
    )
    mock_extract_roof_points.assert_called_once_with(
        mock_download_footprints.return_value,
        raw_pointcloud,
    )


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_surface_mesh.City")
def test_city_surface_mesh_parameter_plumbing(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_surface_mesh,
):
    """Surface-mesh arguments should be forwarded to the builder call."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_surface_mesh.return_value = Mock(name="surface_mesh")

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CitySurfaceMeshDataset()
    dataset.build(
        CitySurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            max_mesh_size=18.0,
            min_mesh_angle=31.0,
            min_building_detail=1.25,
            min_building_area=27.0,
            merge_buildings=False,
            smoothing=4,
        )
    )

    mock_build_city_surface_mesh.assert_called_once_with(
        city,
        max_mesh_size=18.0,
        min_mesh_angle=31.0,
        min_building_detail=1.25,
        min_building_area=27.0,
        merge_buildings=False,
        smoothing=4,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
    )


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_surface_mesh.City")
def test_city_surface_mesh_obj_export_returns_bytes(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_surface_mesh,
):
    """Export formats should return bytes via export_to_bytes()."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = Mock(name="filtered_pointcloud")
    surface_mesh = Mock(name="surface_mesh")

    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_surface_mesh.return_value = surface_mesh
    mock_city_cls.return_value = Mock(name="city")

    dataset = CitySurfaceMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"surface-mesh") as mock_export:
        result = dataset.build(
            CitySurfaceMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"surface-mesh"
    mock_export.assert_called_once_with(surface_mesh, "obj")


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.set_building_heights_from_attribute")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.flatten_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_surface_mesh.City")
def test_city_surface_mesh_flat_ground_replaces_raster_and_resets_building_heights(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_flatten_terrain_raster,
    mock_set_building_heights_from_attribute,
    mock_build_city_surface_mesh,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    flat_buildings = Mock(name="flat_buildings")
    surface_mesh = Mock(name="surface_mesh")
    flat_raster = Mock(name="flat_raster")

    flat_raster.data = np.array([[17.5]])
    flat_raster.nodata = np.nan

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_flatten_terrain_raster.return_value = flat_raster
    mock_set_building_heights_from_attribute.return_value = flat_buildings
    mock_build_city_surface_mesh.return_value = surface_mesh

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CitySurfaceMeshDataset()
    result = dataset.build(
        CitySurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            flat_ground=True,
            ground_level=17.5,
        )
    )

    assert result is surface_mesh
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


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.set_building_heights_from_attribute")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.flatten_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_surface_mesh.City")
def test_city_surface_mesh_default_path_skips_flat_ground_rebuild(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_flatten_terrain_raster,
    mock_set_building_heights_from_attribute,
    mock_build_city_surface_mesh,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = Mock(name="filtered_pointcloud")
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_surface_mesh.return_value = Mock(name="surface_mesh")
    mock_city_cls.return_value = Mock(name="city")

    dataset = CitySurfaceMeshDataset()
    dataset.build(CitySurfaceMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    mock_flatten_terrain_raster.assert_not_called()
    mock_set_building_heights_from_attribute.assert_not_called()
