"""Tests for the city footprints dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import DatasetDescriptor, get_dataset
from dtcc_core.datasets import list as list_datasets
from dtcc_core.datasets.city_footprints import CityFootprintsArgs, CityFootprintsDataset


def test_city_footprints_internal_name():
    assert CityFootprintsDataset().name == "city_footprints"


def test_city_footprints_is_not_public_dataset():
    with pytest.raises(KeyError):
        get_dataset("city_footprints")
    assert "city_footprints" not in list_datasets()
    public_attr = getattr(datasets, "city_footprints", None)
    assert not isinstance(public_attr, DatasetDescriptor)
    assert not callable(public_attr)


@patch("dtcc_core.datasets.city_footprints.prepare_footprint_city_from_bounds")
def test_city_footprints_default_build_returns_cleaned_city(
    mock_prepare_footprint_city_from_bounds,
):
    city = Mock(name="city")
    cleaned = Mock(name="cleaned_city")
    city.clean_footprints.return_value = cleaned
    mock_prepare_footprint_city_from_bounds.return_value = city

    dataset = CityFootprintsDataset()
    result = dataset.build(CityFootprintsArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is cleaned
    mock_prepare_footprint_city_from_bounds.assert_called_once_with(
        ANY,
        progress=ANY,
    )
    city.clean_footprints.assert_called_once_with(
        min_building_detail=0.5,
        min_area=15.0,
        show=False,
    )


def test_city_footprints_build_from_city_forwards_args():
    city = Mock(name="city")
    cleaned = Mock(name="cleaned_city")
    city.clean_footprints.return_value = cleaned

    dataset = CityFootprintsDataset()
    result = dataset.build_from_city(
        city,
        bounds=(0.0, 0.0, 1.0, 1.0),
        min_building_detail=0.85,
        min_building_area=19.0,
        show_footprints=True,
    )

    assert result is cleaned
    city.clean_footprints.assert_called_once_with(
        min_building_detail=0.85,
        min_area=19.0,
        show=True,
    )
