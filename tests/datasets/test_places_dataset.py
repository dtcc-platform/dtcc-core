import typing

import geopandas as gpd
import pytest
from shapely.geometry import Point

import dtcc_core.datasets as datasets
from dtcc_core.datasets.places import PlacesArgs
from dtcc_core.io.data import osm_places
from dtcc_core.model import Bounds, PlaceCollection


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


def test_places_dataset_returns_place_collection(monkeypatch):
    calls = []

    def fake_download(bounds, categories):
        calls.append((bounds, categories))
        return osm_places.places_from_geodataframe(_gdf())

    monkeypatch.setattr("dtcc_core.io.data.download_places", fake_download)

    places = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), categories=["food_and_drink"])

    assert isinstance(places, PlaceCollection)
    assert isinstance(calls[0][0], Bounds)
    assert calls[0][1] == ["food_and_drink"]


def test_places_dataset_exports_geojson(monkeypatch, tmp_path):
    monkeypatch.setattr(
        osm_places, "download_places_geodataframe", lambda bounds, categories: _gdf()
    )

    payload = datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), format="geojson")
    path = tmp_path / "places.geojson"
    path.write_bytes(payload)

    loaded = gpd.read_file(path)
    assert loaded["name"].tolist() == ["Panos"]
    assert loaded["category"].tolist() == ["food_and_drink"]


def test_places_rejects_unknown_category():
    with pytest.raises(Exception):
        datasets.places(bounds=(0.0, 0.0, 100.0, 100.0), categories=["nightlife"])
