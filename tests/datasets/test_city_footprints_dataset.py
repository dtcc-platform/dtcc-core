"""Tests for the city footprints dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import DatasetDescriptor, get_dataset, list as list_datasets
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
@patch("dtcc_core.datasets.city_footprints.condition_city_meshing_footprints")
def test_city_footprints_default_build_returns_conditioned_result(
    mock_condition_city_meshing_footprints,
    mock_prepare_footprint_city_from_bounds,
):
    city = Mock(name="city")
    conditioned = Mock(name="conditioned_footprints")

    mock_prepare_footprint_city_from_bounds.return_value = city
    mock_condition_city_meshing_footprints.return_value = conditioned

    dataset = CityFootprintsDataset()
    result = dataset.build(CityFootprintsArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is conditioned
    mock_prepare_footprint_city_from_bounds.assert_called_once_with(
        ANY,
        progress=ANY,
    )
    mock_condition_city_meshing_footprints.assert_called_once_with(
        city,
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


@patch("dtcc_core.datasets.city_footprints.condition_city_meshing_footprints")
def test_city_footprints_build_from_city_forwards_advanced_args(
    mock_condition_city_meshing_footprints,
):
    city = Mock(name="city")
    conditioned = Mock(name="conditioned_footprints")
    mock_condition_city_meshing_footprints.return_value = conditioned

    dataset = CityFootprintsDataset()
    result = dataset.build_from_city(
        city,
        bounds=(0.0, 0.0, 1.0, 1.0),
        max_mesh_size=None,
        min_building_detail=0.85,
        min_building_area=19.0,
        merge_buildings=False,
        merge_tolerance=0.75,
        cleaning_diagnostics=False,
        show_footprints=True,
        footprint_cleaning_plot_block=False,
    )

    assert result is conditioned
    mock_condition_city_meshing_footprints.assert_called_once_with(
        city,
        min_building_detail=0.85,
        min_building_area=19.0,
        merge_tolerance=0.75,
        merge_buildings=False,
        max_mesh_size=None,
        cleaning_diagnostics=False,
        show_footprints=True,
        footprint_cleaning_plot_block=False,
        pipeline_mode="strict",
    )
