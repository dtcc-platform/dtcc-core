"""Tests for the city dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import numpy as np

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city import CityArgs, CityDataset
from dtcc_core.model import Building, GeometryType, Surface


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
    result = dataset.build(
        CityArgs(bounds=(0.0, 0.0, 1.0, 1.0), smallest_building_size=0.0)
    )

    assert result is city
    city.download_pointcloud.assert_called_once_with()
    city.download_footprints.assert_called_once_with(provider="dtcc")
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
                smallest_building_size=0.0,
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
                smallest_building_size=0.0,
            )
        )

    assert result == b"city-cityjson"
    mock_export.assert_called_once_with(city, "json")


@patch("dtcc_core.datasets.city.City")
def test_city_osm_source_uses_osm_footprint_provider(mock_city_cls):
    city = Mock(name="city")
    mock_city_cls.return_value = city

    dataset = CityDataset()
    dataset.build(
        CityArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            source="OSM",
            smallest_building_size=0.0,
        )
    )

    city.download_footprints.assert_called_once_with(provider="OSM")


@patch("dtcc_core.datasets.city.City")
def test_city_filters_small_footprints_before_lod1_build(mock_city_cls):
    small = _building_with_footprint("small", 0.0, 0.0, 2.0, 2.0)
    large = _building_with_footprint("large", 0.0, 0.0, 12.0, 10.0)
    city = Mock(name="city")
    city.buildings = [small, large]
    city.replace_buildings.side_effect = lambda buildings: setattr(
        city, "buildings", buildings
    )
    mock_city_cls.return_value = city

    dataset = CityDataset()
    dataset.build(
        CityArgs(bounds=(0.0, 0.0, 1.0, 1.0), smallest_building_size=50.0)
    )

    city.download_footprints.assert_called_once_with(provider="dtcc")
    city.replace_buildings.assert_called_once()
    assert city.buildings == [large]
    city.build_lod1_buildings.assert_called_once_with()


def test_city_context_documents_terrain_lod1_lineage_and_limitations():
    dataset = CityDataset()
    context = dataset.create_context(
        dataset.validate(
            {
                "bounds": (0.0, 0.0, 1.0, 1.0),
                "source": "LM",
                "smallest_building_size": 25.0,
                "format": "cityjson",
            }
        )
    )
    manifest = context.manifest()

    provider_roles = {
        provider["name"]: provider["role"] for provider in manifest.metadata.provider
    }
    assert provider_roles == {
        "Lantmäteriet": "source_provider",
        "OpenStreetMap": "source_provider",
        "DTCC Platform": "processor",
    }
    assert manifest.metadata.lod.startswith("Terrain mesh plus LoD1")
    assert "Requires review" in manifest.metadata.license
    assert "Requires review" in manifest.metadata.collection_period
    assert {item["name"] for item in manifest.provenance.derived_from} == {
        "point_cloud",
        "building_footprints",
    }
    assert any("Build terrain" in step for step in manifest.provenance.processing_steps)
    assert any("LoD1" in step for step in manifest.provenance.processing_steps)
    assert manifest.presentation.headline == "Terrain and LoD1 City Model"
    assert manifest.presentation.legend["title"] == "City model layers"
    assert manifest.presentation.view_hints["building_lod"] == "LoD1"
    assert any("LoD1 block geometry" in warning for warning in manifest.presentation.warnings)
    assert manifest.presentation.limitations
    assert manifest.request.parameters["source"] == "LM"
    assert manifest.request.parameters["format"] == "cityjson"


def _building_with_footprint(
    building_id: str,
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
) -> Building:
    surface = Surface()
    surface.vertices = np.array(
        [
            [xmin, ymin, 0.0],
            [xmax, ymin, 0.0],
            [xmax, ymax, 0.0],
            [xmin, ymax, 0.0],
        ]
    )
    surface.transform.srs = "EPSG:3006"
    building = Building()
    building.id = building_id
    building.add_geometry(surface, GeometryType.LOD0)
    return building
