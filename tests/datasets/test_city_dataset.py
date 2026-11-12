"""Tests for the city dataset wrapper."""

from __future__ import annotations

import json
from unittest.mock import Mock, patch

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset, load_model_package
from dtcc_core.datasets.city import CityArgs, CityDataset
from dtcc_core.model import Building, BuildingPart, City, GeometryType, Surface


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
    with patch.object(
        dataset, "export_to_bytes", return_value=b"city-json"
    ) as mock_export:
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
    with patch.object(
        dataset, "export_to_bytes", return_value=b"city-cityjson"
    ) as mock_export:
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
    dataset.build(CityArgs(bounds=(0.0, 0.0, 1.0, 1.0), smallest_building_size=50.0))

    city.download_footprints.assert_called_once_with(provider="dtcc")
    city.replace_buildings.assert_called_once()
    assert city.buildings == [large]
    city.build_lod1_buildings.assert_called_once_with()


def test_city_context_documents_terrain_lod_lineage_and_limitations():
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
    assert manifest.metadata.lod == "lod1"
    assert "Requires review" in manifest.metadata.license
    assert "Requires review" in manifest.metadata.collection_period
    assert {item["name"] for item in manifest.provenance.derived_from} == {
        "point_cloud",
        "building_footprints",
    }
    assert any("Build terrain" in step for step in manifest.provenance.processing_steps)
    assert any("LOD2" in step for step in manifest.provenance.processing_steps)
    assert manifest.presentation.headline == "Terrain and Buildings"
    assert manifest.presentation.legend["title"] == "City model layers"
    assert any(
        "rejected fits remain LOD1" in warning
        for warning in manifest.presentation.warnings
    )
    assert manifest.presentation.limitations
    assert manifest.request.parameters["source"] == "LM"
    assert manifest.request.parameters["format"] == "cityjson"


@pytest.mark.parametrize("method", ["catalogue", "ransac"])
@patch("dtcc_core.datasets.city.City")
def test_city_lod2_dispatch_and_metadata(mock_city_cls, method):
    args = CityArgs(
        bounds=(0.0, 0.0, 1.0, 1.0),
        smallest_building_size=0.0,
        lod="lod2",
        lod2_method=method,
    )
    dataset = CityDataset()
    city = dataset.build(args)
    city.build_lod2_buildings.assert_called_once_with(lod2_method=method)
    city.build_lod1_buildings.assert_not_called()
    assert dataset.create_context(args).manifest().metadata.lod == "lod2"


def test_city_rejects_ransac_for_lod1():
    with pytest.raises(ValueError, match="requires lod='lod2'"):
        CityArgs(bounds=(0.0, 0.0, 1.0, 1.0), lod2_method="ransac")


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


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.city.City")
def test_city_canonical_geojson_round_trip(mock_city_cls, tmp_path):
    """Canonical GeoJSON is WGS84, reloads from the package and re-exports equal."""
    city = City()
    whole = _building_with_footprint("whole", 319900.0, 6397900.0, 319960.0, 6397940.0)
    whole.attributes["storeys"] = 4
    split = Building(id="split")
    donor = _building_with_footprint("x", 319800.0, 6397800.0, 319830.0, 6397830.0)
    part = BuildingPart()
    part.add_geometry(donor.get_geometry(GeometryType.LOD0), GeometryType.LOD0)
    split.add_child(part)
    city.add_buildings([whole, split, Building(id="empty")])
    calls = []
    city.download_pointcloud = lambda **kwargs: calls.append("pointcloud")
    city.download_footprints = lambda **kwargs: calls.append("footprints")
    city.build_terrain = lambda **kwargs: calls.append("terrain")
    city.build_lod1_buildings = lambda **kwargs: calls.append("lod1")
    mock_city_cls.return_value = city

    result = datasets.city(
        bounds=(319720.0, 6397660.0, 320220.0, 6398160.0),
        smallest_building_size=0.0,
    )
    package = result.export(tmp_path / "p", canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(tmp_path / "q", canonical=True, format="geojson")
    assert len(calls) == 4

    payload, copy = (
        json.loads((p.path / a.path).read_text())
        for p in (package, again)
        for a in p.artifacts
        if a.format == "geojson"
    )
    assert payload == copy
    assert "crs" not in payload
    assert payload["metadata"] == {"skipped": 1}
    first, second = payload["features"]
    assert (first["id"], second["id"]) == ("whole", "split")
    assert first["properties"]["storeys"] == 4
    for feature in (first, second):
        [ring] = feature["geometry"]["coordinates"]
        assert all(-180 <= lon <= 180 and -90 <= lat <= 90 for lon, lat in ring)
    # Oracle: pyproj EPSG:3006 -> EPSG:4326 of the footprint corner
    # (319900, 6397900); rings are counter-clockwise.
    assert [11.97879311854237, 57.68719536715784] in [
        pytest.approx(position, abs=1e-9)
        for position in first["geometry"]["coordinates"][0]
    ]
