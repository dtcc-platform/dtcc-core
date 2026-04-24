"""Tests for the city volume mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import numpy as np
import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_volume_mesh import (
    CityVolumeMeshArgs,
    CityVolumeMeshDataset,
    _regular_tet_volume,
)


def test_city_volume_mesh_registered_name():
    """The dataset class should expose the correct registration name."""
    assert CityVolumeMeshDataset().name == "city_volume_mesh"


def test_city_volume_mesh_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("city_volume_mesh")
    assert ds is not None
    assert ds.name == "city_volume_mesh"


def test_city_volume_mesh_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "city_volume_mesh")
    assert callable(datasets.city_volume_mesh)


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_default_build_returns_volume_mesh_and_uses_regular_tet_default_max_volume(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_volume_mesh,
):
    """The default path should derive TetGen max_volume from max_mesh_size."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    volume_mesh = Mock(name="volume_mesh")

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_build_city_volume_mesh.return_value = volume_mesh

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityVolumeMeshDataset()
    result = dataset.build(CityVolumeMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is volume_mesh
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
    mock_build_city_volume_mesh.assert_called_once_with(
        city,
        max_mesh_size=25.0,
        top_cap_max_mesh_size=None,
        domain_height=80.0,
        min_building_detail=0.5,
        boundary_face_markers=True,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
        mesher=None,
        max_volume=_regular_tet_volume(25.0),
        tetgen_switches={
            "extra": "",
        },
    )


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_explicit_max_volume_and_tetgen_extra_are_forwarded(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_volume_mesh,
):
    """Explicit TetGen settings should be forwarded through tetgen_switches."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_volume_mesh.return_value = Mock(name="volume_mesh")

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityVolumeMeshDataset()
    dataset.build(
        CityVolumeMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            max_mesh_size=30.0,
            domain_height=95.0,
            boundary_face_markers=False,
            min_building_detail=1.25,
            max_volume=12.5,
            tetgen_extra="VV",
        )
    )

    mock_build_city_volume_mesh.assert_called_once_with(
        city,
        max_mesh_size=30.0,
        top_cap_max_mesh_size=None,
        domain_height=95.0,
        min_building_detail=1.25,
        boundary_face_markers=False,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
        mesher=None,
        max_volume=12.5,
        tetgen_switches={
            "extra": "VV",
        },
    )


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_top_cap_max_mesh_size_is_forwarded(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_volume_mesh,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_volume_mesh.return_value = Mock(name="volume_mesh")

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityVolumeMeshDataset()
    dataset.build(
        CityVolumeMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            max_mesh_size=30.0,
            top_cap_max_mesh_size=45.0,
        )
    )

    assert mock_build_city_volume_mesh.call_args.kwargs["top_cap_max_mesh_size"] == 45.0


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_mesher_is_forwarded(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_volume_mesh,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_volume_mesh.return_value = Mock(name="volume_mesh")
    mock_city_cls.return_value = Mock(name="city")

    dataset = CityVolumeMeshDataset()
    dataset.build(CityVolumeMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0), mesher="dtcc_mesher"))

    assert mock_build_city_volume_mesh.call_args.kwargs["mesher"] == "dtcc_mesher"


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_vtu_export_returns_bytes(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_city_volume_mesh,
):
    """Export formats should return bytes via export_to_bytes()."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = Mock(name="filtered_pointcloud")
    volume_mesh = Mock(name="volume_mesh")

    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_volume_mesh.return_value = volume_mesh
    mock_city_cls.return_value = Mock(name="city")

    dataset = CityVolumeMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"volume-mesh") as mock_export:
        result = dataset.build(
            CityVolumeMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="vtu",
            )
        )

    assert result == b"volume-mesh"
    mock_export.assert_called_once_with(volume_mesh, "vtu")


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.set_building_heights_from_attribute")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.flatten_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_flat_ground_replaces_raster_and_resets_building_heights(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_flatten_terrain_raster,
    mock_set_building_heights_from_attribute,
    mock_build_city_volume_mesh,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    filtered_pointcloud = Mock(name="filtered_pointcloud")
    buildings = Mock(name="buildings")
    raster = Mock(name="raster")
    roof_buildings = Mock(name="roof_buildings")
    heighted_buildings = Mock(name="heighted_buildings")
    flat_buildings = Mock(name="flat_buildings")
    volume_mesh = Mock(name="volume_mesh")
    flat_raster = Mock(name="flat_raster")

    flat_raster.data = np.array([[21.0]])
    flat_raster.nodata = np.nan

    raw_pointcloud.remove_global_outliers.return_value = filtered_pointcloud
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = buildings
    mock_build_terrain_raster.return_value = raster
    mock_extract_roof_points.return_value = roof_buildings
    mock_compute_building_heights.return_value = heighted_buildings
    mock_flatten_terrain_raster.return_value = flat_raster
    mock_set_building_heights_from_attribute.return_value = flat_buildings
    mock_build_city_volume_mesh.return_value = volume_mesh

    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityVolumeMeshDataset()
    result = dataset.build(
        CityVolumeMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            flat_ground=True,
            ground_level=21.0,
        )
    )

    assert result is volume_mesh
    mock_flatten_terrain_raster.assert_called_once_with(raster, height=21.0)
    mock_set_building_heights_from_attribute.assert_called_once_with(
        heighted_buildings,
        flat_raster,
        height_attribute="height",
        default_ground_height=21.0,
        always_use_default_ground=True,
    )
    city.add_terrain.assert_called_once_with(flat_raster)
    city.add_buildings.assert_called_once_with(
        flat_buildings,
        remove_outside_terrain=True,
    )


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.set_building_heights_from_attribute")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.flatten_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_footprints")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.io.data.download_pointcloud")
@patch("dtcc_core.datasets.city_volume_mesh.City")
def test_city_volume_mesh_default_path_skips_flat_ground_rebuild(
    mock_city_cls,
    mock_download_pointcloud,
    mock_download_footprints,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_flatten_terrain_raster,
    mock_set_building_heights_from_attribute,
    mock_build_city_volume_mesh,
):
    raw_pointcloud = Mock(name="raw_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = Mock(name="filtered_pointcloud")
    mock_download_pointcloud.return_value = raw_pointcloud
    mock_download_footprints.return_value = Mock(name="buildings")
    mock_build_terrain_raster.return_value = Mock(name="raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")
    mock_build_city_volume_mesh.return_value = Mock(name="volume_mesh")
    mock_city_cls.return_value = Mock(name="city")

    dataset = CityVolumeMeshDataset()
    dataset.build(CityVolumeMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    mock_flatten_terrain_raster.assert_not_called()
    mock_set_building_heights_from_attribute.assert_not_called()
