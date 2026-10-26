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


# Elements as Overpass returns them for "out body center": nodes carry lat/lon,
# ways carry node refs and a centre, relations carry members and a centre.
ELEMENTS = [
    _node(1, {"amenity": "restaurant", "name": "Panos", "cuisine": "greek",
              "opening_hours": "Tu-Su 17:00-22:00", "addr:street": "Teknologgatan"}),
    _node(2, {"shop": "hairdresser", "name": "Salong Marie"}, x=319600.0),
    {"type": "way", "id": 3, "center": _lonlat(319700.0, 6397700.0),
     "nodes": [11, 12, 13, 14, 11], "tags": {"amenity": "school", "name": "Vasaskolan"}},
    _node(4, {"amenity": "bench", "name": "Minnesbänk"}),               # street furniture
    _node(5, {"amenity": "cafe"}),                                     # unnamed
    _node(6, {"shop": "bakery", "name": "Outside"}, x=321000.0),       # outside bounds
    _node(1, {"amenity": "restaurant", "name": "Panos"}),              # duplicate
    _node(7, {"healthcare": "physiotherapist", "name": "Sjukgymnast"}),
    {"type": "relation", "id": 8, "center": _lonlat(319800.0, 6397800.0),
     "members": [{"type": "way", "ref": 21, "role": "outer"}],
     "tags": {"type": "multipolygon", "tourism": "museum", "name": "Röhsska"}},
]


def _overpass_response(query, elements):
    """Return elements as Overpass would for the output mode in ``query``."""
    if "out tags" in query:
        elements = [
            {k: v for k, v in e.items() if e["type"] != "node" or k not in ("lat", "lon")}
            for e in elements
        ]
    return {"version": 0.6, "elements": elements}


@pytest.fixture
def fake_overpass(monkeypatch, tmp_path):
    queries = []

    def fake_query(query, read_timeout=None):
        queries.append(query)
        return _overpass_response(query, ELEMENTS)

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

    assert gdf["name"].tolist() == [
        "Panos", "Salong Marie", "Vasaskolan", "Sjukgymnast", "Röhsska"
    ]
    assert gdf["category"].tolist() == [
        "food_and_drink", "shopping", "education", "health", "culture"
    ]
    assert gdf["osm_type"].tolist() == ["node", "node", "way", "node", "relation"]
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
    assert "out body center" in query
    assert "out tags" not in query


def test_download_places_returns_place_collection(fake_overpass):
    places = osm_places.download_places(BOUNDS)

    assert isinstance(places, PlaceCollection)
    assert len(places.places()) == 5
    assert places.attributes["requested_bounds"] == list(BOUNDS.tuple)
    assert places.bounds.tuple == BOUNDS.tuple
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


def test_nearby_bounds_do_not_share_a_cache_entry(fake_overpass):
    # These bounds round to the same whole metres but select different places.
    near = Bounds(319000.0, 6397000.0, 319599.6, 6398000.0)
    nearer = Bounds(319000.0, 6397000.0, 319600.4, 6398000.0)

    first = osm_places.download_places_geodataframe(near)
    second = osm_places.download_places_geodataframe(nearer)

    assert len(fake_overpass.queries) == 2
    assert osm_places._cache_path(near) != osm_places._cache_path(nearer)
    assert "Salong Marie" not in first["name"].tolist()
    assert "Salong Marie" in second["name"].tolist()


@pytest.mark.parametrize(
    "payload, message",
    [
        ({"elements": [], "remark": "runtime error: Query timed out in \"query\" at line 1 after 91 seconds."},
         "runtime error: Query timed out"),
        ({"remark": "runtime error: Query run out of memory using about 2048 MB of RAM."},
         "out of memory"),
        ({"version": 0.6}, "without an 'elements' list"),
        ({"elements": ["oops"]}, "element 0 is not an object"),
        ({"elements": [{"type": "node", "lat": 57.7, "lon": 11.9,
                        "tags": {"amenity": "cafe", "name": "Linsen"}}]},
         "element 0 has no integer OSM id"),
        ({"elements": [{"type": "area", "id": 1, "tags": {}}]}, "element 0 has no valid OSM type"),
        ({"elements": [{"type": "node", "id": 1, "tags": ["name"]}]}, "tags that are not an object"),
    ],
)
def test_incomplete_response_raises_and_is_not_cached(monkeypatch, tmp_path, payload, message):
    monkeypatch.setattr(osm_places, "query_overpass_with_failover",
                        lambda query, read_timeout=None: payload)
    monkeypatch.setattr(osm_places, "cache_dir", tmp_path)

    with pytest.raises(RuntimeError, match=message) as error:
        osm_places.download_places_geodataframe(BOUNDS)

    assert "DTCC_PLACES_OVERPASS_TIMEOUT" in str(error.value)
    assert not osm_places._cache_path(BOUNDS).exists()
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_failed_endpoints_raise_actionable_error(monkeypatch, tmp_path):
    def failing_query(query, read_timeout=None):
        raise RuntimeError("All Overpass endpoints failed: HTTP 504")

    monkeypatch.setattr(osm_places, "query_overpass_with_failover", failing_query)
    monkeypatch.setattr(osm_places, "cache_dir", tmp_path)

    with pytest.raises(RuntimeError, match="Could not download OSM places.*HTTP 504.*smaller bounds"):
        osm_places.download_places(BOUNDS)
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_places_without_coordinates_are_reported(monkeypatch, tmp_path):
    warnings = []
    elements = [_node(1, {"amenity": "cafe", "name": "Linsen"}),
                {"type": "node", "id": 2, "tags": {"amenity": "cafe", "name": "Ingen plats"}}]
    monkeypatch.setattr(osm_places, "query_overpass_with_failover",
                        lambda query, read_timeout=None: {"elements": elements})
    monkeypatch.setattr(osm_places, "cache_dir", tmp_path)
    monkeypatch.setattr(osm_places, "warning", warnings.append)

    gdf = osm_places.download_places_geodataframe(BOUNDS)

    assert gdf["name"].tolist() == ["Linsen"]
    assert "Skipped 1 OSM places without coordinates" in warnings[0]


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
