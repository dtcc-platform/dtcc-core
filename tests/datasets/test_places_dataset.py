import importlib
import typing

import geopandas as gpd
import pytest
from shapely.geometry import Point

import dtcc_core.datasets as datasets
from dtcc_core.datasets.places import PlacesArgs
from dtcc_core.io.data import osm_places
from dtcc_core.model import Bounds, PlaceCollection, RoadNetwork

places_dataset = importlib.import_module("dtcc_core.datasets.places")


def _gdf():
    row = {column: None for column in osm_places.PLACE_COLUMNS}
    row.update(osm_id=1, osm_type="node", name="Panos", category="food_and_drink",
               osm_key="amenity", osm_value="restaurant")
    return gpd.GeoDataFrame([row], geometry=[Point(10.0, 20.0)], crs="EPSG:3006")


def test_places_dataset_registered():
    assert "places" in datasets.list()
    assert callable(datasets.places)


def test_dataset_categories_match_download_categories():
    category = typing.get_args(PlacesArgs.model_fields["categories"].annotation)[0]
    assert sorted(typing.get_args(typing.get_args(category)[0])) == list(
        osm_places.PLACE_CATEGORY_NAMES
    )


@pytest.fixture
def fake_downloads(monkeypatch):
    calls = {"places": [], "roads": []}

    def fake_places(bounds, categories):
        calls["places"].append((bounds, categories))
        return osm_places.places_from_geodataframe(_gdf())

    def fake_roads(bounds, provider, epsg):
        calls["roads"].append((bounds, provider, epsg))
        return RoadNetwork()

    monkeypatch.setattr("dtcc_core.io.data.download_places", fake_places)
    monkeypatch.setattr("dtcc_core.io.data.download_roadnetwork", fake_roads)
    return calls


def test_places_dataset_returns_place_collection(fake_downloads):
    places = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), categories=["food_and_drink"])

    assert isinstance(places, PlaceCollection)
    assert isinstance(fake_downloads["places"][0][0], Bounds)
    assert fake_downloads["places"][0][1] == ["food_and_drink"]


def test_places_dataset_downloads_context_roads(fake_downloads):
    places = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0))

    assert isinstance(places.context_roads, RoadNetwork)
    bounds, provider, epsg = fake_downloads["roads"][0]
    assert bounds.tuple == (0.0, 0.0, 100.0, 100.0)
    assert (provider, epsg) == ("OSM", "3006")


def test_places_dataset_can_skip_context_roads(fake_downloads):
    places = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), roads=False)

    assert places.context_roads is None
    assert fake_downloads["roads"] == []


def test_places_dataset_keeps_places_when_roads_fail(fake_downloads, monkeypatch):
    warnings = []

    def failing_roads(bounds, provider, epsg):
        raise RuntimeError("All Overpass endpoints failed: HTTP 504")

    monkeypatch.setattr("dtcc_core.io.data.download_roadnetwork", failing_roads)
    monkeypatch.setattr(places_dataset, "warning", warnings.append)

    places = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0))

    assert [p.attributes["name"] for p in places.places()] == ["Panos"]
    assert places.context_roads is None
    assert "HTTP 504" in warnings[0]


def test_places_dataset_exports_geojson(fake_downloads, monkeypatch, tmp_path):
    monkeypatch.setattr(
        osm_places, "download_places_geodataframe", lambda bounds, categories: _gdf()
    )

    payload = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), format="geojson")
    path = tmp_path / "places.geojson"
    path.write_bytes(payload)

    loaded = gpd.read_file(path)
    assert loaded["name"].tolist() == ["Panos"]
    assert loaded["category"].tolist() == ["food_and_drink"]
    assert fake_downloads["roads"] == []


def test_places_rejects_unknown_category():
    with pytest.raises(Exception):
        datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), categories=["nightlife"])
