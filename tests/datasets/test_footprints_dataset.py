"""Tests for the building footprints dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.footprints import FootprintsArgs, FootprintsDataset


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
