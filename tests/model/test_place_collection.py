import pytest

from dtcc_core.model import Object, PlaceCollection
from dtcc_core.model.geometry import Point


def _place(name, category, x, y):
    place = Object()
    place.attributes = {"name": name, "category": category}
    place.add_geometry(Point(x=x, y=y, z=0.0), "location")
    return place


def _collection():
    collection = PlaceCollection()
    collection.add_place(_place("Panos", "food_and_drink", 1.0, 2.0))
    collection.add_place(_place("Salong Marie", "shopping", 3.0, 4.0))
    collection.add_place(_place("Café Linsen", "food_and_drink", 5.0, 6.0))
    return collection


def test_places_returns_all_children():
    assert [p.attributes["name"] for p in _collection().places()] == [
        "Panos", "Salong Marie", "Café Linsen"
    ]


def test_add_place_rejects_non_objects():
    with pytest.raises(TypeError, match="Place must be an Object"):
        PlaceCollection().add_place("Panos")


def test_to_geodataframe():
    gdf = _collection().to_geodataframe()

    assert gdf["name"].tolist() == ["Panos", "Salong Marie", "Café Linsen"]
    assert gdf.crs.to_epsg() == 3006
    assert [(p.x, p.y) for p in gdf.geometry] == [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]


def test_summary_and_info_sections():
    collection = _collection()

    assert collection._summary_items() == [("num_places", 3)]
    sections = dict((title, rows) for title, _, rows in collection._info_sections())
    assert sections["Places by category"] == [("food_and_drink", 2), ("shopping", 1)]
    assert sections["Sample places (first 3)"][0][0] == "Panos"
