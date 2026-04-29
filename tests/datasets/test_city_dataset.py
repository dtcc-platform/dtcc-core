"""Tests for the city dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city import CityArgs, CityDataset


def test_city_registered_name():
    """The dataset class should expose the correct registration name."""
    assert CityDataset().name == "city"


def test_city_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("city")
    assert ds is not None
    assert ds.name == "city"


def test_city_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "city")
    assert callable(datasets.city)


@patch("dtcc_core.datasets.city.City")
def test_city_default_build_returns_city_object(mock_city_cls):
    """Without export, the wrapper should return the City object itself."""
    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityDataset()
    result = dataset.build(CityArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is city
    city.download_pointcloud.assert_called_once_with()
    city.download_footprints.assert_called_once_with()
    city.build_terrain.assert_called_once_with(build_mesh=True)
    city.build_lod1_buildings.assert_called_once_with()


@patch("dtcc_core.datasets.city.City")
def test_city_json_export_returns_bytes(mock_city_cls):
    """JSON export should return bytes via export_to_bytes()."""
    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"city-json") as mock_export:
        result = dataset.build(
            CityArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="json",
            )
        )

    assert result == b"city-json"
    mock_export.assert_called_once_with(city, "json")


@patch("dtcc_core.datasets.city.City")
def test_city_cityjson_export_currently_uses_json_export_path(mock_city_cls):
    """The wrapper currently exports cityjson requests through the JSON byte path."""
    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"city-cityjson") as mock_export:
        result = dataset.build(
            CityArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="cityjson",
            )
        )

    assert result == b"city-cityjson"
    mock_export.assert_called_once_with(city, "json")
