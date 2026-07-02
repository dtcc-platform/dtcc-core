"""Tests for the building footprints dataset wrapper."""

from __future__ import annotations

import json
from unittest.mock import Mock, patch

import numpy as np

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.footprints import FootprintsArgs, FootprintsDataset
from dtcc_core.model import Building, City, FootprintCollection, GeometryType, Surface


def test_building_footprints_registered_name():
    """The dataset class should expose the correct registration name."""
    assert FootprintsDataset().name == "building_footprints"


def test_building_footprints_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("building_footprints")
    assert ds is not None
    assert ds.name == "building_footprints"


def test_building_footprints_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "building_footprints")
    assert callable(datasets.building_footprints)


def test_building_footprints_prepare_result_returns_footprint_collection():
    """Public no-format calls adapt raw buildings to semantic footprints."""
    city = _city_with_one_footprint()
    dataset = FootprintsDataset()
    args = FootprintsArgs(bounds=(319720.0, 6397660.0, 320220.0, 6398160.0))

    result = dataset.prepare_result(city.buildings, args)

    assert isinstance(result, FootprintCollection)
    assert len(result) == 1
    assert result[0].vertices.shape == (4, 3)


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_default_build_returns_city_buildings(mock_city_cls):
    """Without export, the wrapper should return city.buildings directly."""
    buildings = [Mock(name="building1"), Mock(name="building2")]
    city = Mock(name="city")
    city.buildings = buildings
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    result = dataset.build(FootprintsArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is buildings
    city.download_footprints.assert_called_once_with()
    city.download_pointcloud.assert_not_called()
    city.building_heights_from_pointcloud.assert_not_called()


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_calculate_heights_false_skips_pointcloud(mock_city_cls):
    """Point-cloud download should be skipped when calculate_heights=False."""
    city = Mock(name="city")
    city.buildings = []
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    dataset.build(
        FootprintsArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            calculate_heights=False,
        )
    )

    city.download_footprints.assert_called_once_with()
    city.download_pointcloud.assert_not_called()
    city.building_heights_from_pointcloud.assert_not_called()


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_calculate_heights_true_runs_height_pipeline(mock_city_cls):
    """Height calculation should download the point cloud and compute building heights."""
    city = Mock(name="city")
    city.buildings = []
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    dataset.build(
        FootprintsArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            calculate_heights=True,
        )
    )

    city.download_footprints.assert_called_once_with()
    city.download_pointcloud.assert_called_once_with()
    city.building_heights_from_pointcloud.assert_called_once_with(
        keep_roof_points=False
    )


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_geojson_export_returns_bytes(mock_city_cls):
    """Vector export formats should return bytes via export_to_bytes()."""
    city = Mock(name="city")
    city.buildings = []
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"footprints-bytes") as mock_export:
        result = dataset.build(
            FootprintsArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="geojson",
            )
        )

    assert result == b"footprints-bytes"
    mock_export.assert_called_once()
    args, kwargs = mock_export.call_args
    assert args[:2] == (city, "geojson")
    assert kwargs["save_callable"] is not None


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_crs_forwarded_to_export(mock_city_cls):
    """The crs argument should reach export_to_bytes as the output_crs kwarg."""
    city = Mock(name="city")
    city.buildings = []
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"x") as mock_export:
        dataset.build(
            FootprintsArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="geojson",
                crs="EPSG:3006",
            )
        )

    assert mock_export.call_args.kwargs["output_crs"] == "EPSG:3006"


def _city_with_one_footprint() -> City:
    """A real City with one EPSG:3006 footprint inside the table bounds."""
    surface = Surface()
    surface.vertices = np.array(
        [
            [319900.0, 6397900.0, 0.0],
            [319960.0, 6397900.0, 0.0],
            [319960.0, 6397940.0, 0.0],
            [319900.0, 6397940.0, 0.0],
        ]
    )
    surface.transform.srs = "EPSG:3006"
    building = Building()
    building.id = "bldg-1"
    building.add_geometry(surface, GeometryType.LOD0)
    city = City()
    city.add_buildings([building])
    return city


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_geojson_crs_3006_is_table_compatible(mock_city_cls):
    """crs='EPSG:3006' GeoJSON keeps meter coordinates and declares SWEREF99 TM.

    The DTCC table (Atlas) requires a crs member naming EPSG:3006 and
    rejects WGS84-degree coordinates.
    """
    city = _city_with_one_footprint()
    city.download_footprints = lambda: None
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    payload = dataset.build(
        FootprintsArgs(
            bounds=(319720.0, 6397660.0, 320220.0, 6398160.0),
            format="geojson",
            crs="EPSG:3006",
        )
    )

    data = json.loads(payload)
    assert data["crs"]["properties"]["name"] == "urn:ogc:def:crs:EPSG::3006"
    x, y = data["features"][0]["geometry"]["coordinates"][0][0][:2]
    assert abs(x) > 180 and abs(y) > 90, (x, y)


@patch("dtcc_core.datasets.footprints.City")
def test_building_footprints_geojson_default_remains_wgs84(mock_city_cls):
    """Without crs, GeoJSON export keeps reprojecting to WGS84 degrees."""
    city = _city_with_one_footprint()
    city.download_footprints = lambda: None
    mock_city_cls.return_value = city

    dataset = FootprintsDataset()
    payload = dataset.build(
        FootprintsArgs(
            bounds=(319720.0, 6397660.0, 320220.0, 6398160.0),
            format="geojson",
        )
    )

    data = json.loads(payload)
    x, y = data["features"][0]["geometry"]["coordinates"][0][0][:2]
    assert abs(x) <= 180 and abs(y) <= 90, (x, y)
