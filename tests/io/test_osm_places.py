import pytest
from pyproj import Transformer

from dtcc_core.io.data import osm_places
from dtcc_core.model import Bounds, PlaceCollection

BOUNDS = Bounds(319000.0, 6397000.0, 320000.0, 6398000.0)
TO_WGS84 = Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)


def _lonlat(x, y):
    lon, lat = TO_WGS84.transform(x, y)
    return {"lon": lon, "lat": lat}


def _node(osm_id, tags, x=319500.0, y=6397500.0):
    return {"type": "node", "id": osm_id, **_lonlat(x, y), "tags": tags}


ELEMENTS = [
    _node(1, {"amenity": "restaurant", "name": "Panos", "cuisine": "greek",
              "opening_hours": "Tu-Su 17:00-22:00", "addr:street": "Teknologgatan"}),
    _node(2, {"shop": "hairdresser", "name": "Salong Marie"}, x=319600.0),
    {"type": "way", "id": 3, "center": _lonlat(319700.0, 6397700.0),
     "tags": {"amenity": "school", "name": "Vasaskolan"}},
    _node(4, {"amenity": "bench", "name": "Minnesbänk"}),               # street furniture
    _node(5, {"amenity": "cafe"}),                                     # unnamed
    _node(6, {"shop": "bakery", "name": "Outside"}, x=321000.0),       # outside bounds
    _node(1, {"amenity": "restaurant", "name": "Panos"}),              # duplicate
    _node(7, {"healthcare": "physiotherapist", "name": "Sjukgymnast"}),
]


@pytest.fixture
def fake_overpass(monkeypatch, tmp_path):
    queries = []

    def fake_query(query, read_timeout=None):
        queries.append(query)
        return {"elements": ELEMENTS}

    monkeypatch.setattr(osm_places, "query_overpass_with_failover", fake_query)
    monkeypatch.setattr(osm_places, "cache_dir", tmp_path)
    fake_query.queries = queries
    return fake_query


@pytest.mark.parametrize(
    "tags, expected",
    [
        ({"amenity": "restaurant"}, ("food_and_drink", "amenity", "restaurant")),
        ({"shop": "bakery"}, ("shopping", "shop", "bakery")),
        ({"tourism": "hotel"}, ("accommodation", "tourism", "hotel")),
        ({"tourism": "attraction"}, ("tourism", "tourism", "attraction")),
        ({"amenity": "theatre"}, ("culture", "amenity", "theatre")),
        ({"healthcare": "physiotherapist"}, ("health", "healthcare", "physiotherapist")),
        ({"amenity": "pharmacy", "shop": "chemist"}, ("health", "amenity", "pharmacy")),
        ({"amenity": "marketplace"}, ("other", "amenity", "marketplace")),
        ({"amenity": "bench"}, None),
        ({"tourism": "information"}, None),
        ({"highway": "bus_stop"}, None),
    ],
)
def test_place_category(tags, expected):
    assert osm_places.place_category(tags) == expected


def test_geodataframe_filters_and_classifies(fake_overpass):
    gdf = osm_places.download_places_geodataframe(BOUNDS)

    assert gdf["name"].tolist() == ["Panos", "Salong Marie", "Vasaskolan", "Sjukgymnast"]
    assert gdf["category"].tolist() == ["food_and_drink", "shopping", "education", "health"]
    assert gdf["osm_type"].tolist() == ["node", "node", "way", "node"]
    assert gdf.crs.to_epsg() == 3006
    school = gdf[gdf["name"] == "Vasaskolan"].geometry.iloc[0]
    assert (school.x, school.y) == pytest.approx((319700.0, 6397700.0), abs=0.01)
    assert gdf.loc[0, "cuisine"] == "greek"
    assert gdf.loc[0, "addr_street"] == "Teknologgatan"


def test_query_asks_for_named_features_of_every_key(fake_overpass):
    osm_places.download_places_geodataframe(BOUNDS)

    query = fake_overpass.queries[0]
    for key in osm_places.PLACE_KEYS:
        assert f'nwr["{key}"]["name"]' in query
    assert "out tags center" in query


def test_download_places_returns_place_collection(fake_overpass):
    places = osm_places.download_places(BOUNDS)

    assert isinstance(places, PlaceCollection)
    assert len(places.places()) == 4
    first = places.places()[0]
    assert first.attributes["name"] == "Panos"
    assert first.attributes["osm_id"] == 1
    assert first.attributes["website"] is None
    point = first.get_geometries()[0]
    assert (point.x, point.y) == pytest.approx((319500.0, 6397500.0), abs=0.01)
    assert "OpenStreetMap" in places.attributes["attribution"]


def test_category_filter(fake_overpass):
    gdf = osm_places.download_places_geodataframe(BOUNDS, categories=["health", "education"])

    assert gdf["name"].tolist() == ["Vasaskolan", "Sjukgymnast"]


def test_uses_cache_on_second_call(fake_overpass):
    osm_places.download_places_geodataframe(BOUNDS)

    cached = osm_places.download_places(BOUNDS, categories="shopping")

    assert len(fake_overpass.queries) == 1
    assert [p.attributes["name"] for p in cached.places()] == ["Salong Marie"]
    assert cached.places()[0].attributes["osm_id"] == 2


def test_empty_area_warns(monkeypatch, tmp_path):
    warnings = []
    monkeypatch.setattr(osm_places, "query_overpass_with_failover",
                        lambda query, read_timeout=None: {"elements": []})
    monkeypatch.setattr(osm_places, "cache_dir", tmp_path)
    monkeypatch.setattr(osm_places, "warning", warnings.append)

    places = osm_places.download_places(BOUNDS)

    assert places.places() == []
    assert "No OSM places" in warnings[0]


def test_unknown_category_raises():
    with pytest.raises(ValueError, match="Unknown place categories"):
        osm_places.download_places(BOUNDS, categories=["nightlife"])
