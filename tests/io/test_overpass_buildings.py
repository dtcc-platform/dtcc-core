import re

import pyproj
import pytest

from dtcc_core.io.data import overpass

BBOX = (319900.0, 6398900.0, 320100.0, 6399100.0)


def _square(lon0, lat0, size):
    return [
        {"lat": lat0, "lon": lon0},
        {"lat": lat0, "lon": lon0 + size},
        {"lat": lat0 + size, "lon": lon0 + size},
        {"lat": lat0 + size, "lon": lon0},
        {"lat": lat0, "lon": lon0},
    ]


def test_overpass_buildings_keep_osm_id_and_tags(monkeypatch):
    def fake_query(query):
        assert 'relation["building"]["type"="multipolygon"]' in query
        return {
            overpass.OVERPASS_ENDPOINT_METADATA_KEY: "endpoint-1",
            "elements": [
                {
                    "type": "way",
                    "id": 4242,
                    "nodes": [1, 2, 3, 4, 1],
                    "geometry": _square(11.97, 57.70, 0.0002),
                    "tags": {
                        "building": "apartments",
                        "building:levels": "5",
                        "roof:shape": "gabled",
                        "height": "17 m",
                        "addr:street": "Testgatan",
                        "addr:housenumber": "3",
                    },
                },
                {
                    "type": "way",
                    "id": 4243,
                    "nodes": [5, 6, 7, 8, 5],
                    "geometry": _square(11.971, 57.70, 0.0002),
                    "tags": {"building": "yes"},
                },
            ],
        }

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fake_query)

    buildings = overpass.download_overpass_buildings(BBOX)

    assert buildings["osm_id"].tolist() == [4242, 4243]
    assert buildings["osm_type"].tolist() == ["way", "way"]
    assert buildings["building"].tolist() == ["apartments", "yes"]
    assert buildings["building_levels"].iloc[0] == "5"
    assert buildings["roof_shape"].iloc[0] == "gabled"
    assert buildings["osm_height"].iloc[0] == "17 m"
    assert buildings["addr_street"].iloc[0] == "Testgatan"
    assert buildings["addr_housenumber"].iloc[0] == "3"
    assert buildings["building_levels"].isna().iloc[1]
    assert buildings.crs.to_epsg() == 3006
    assert buildings.attrs["source_endpoint"] == "endpoint-1"


def test_overpass_buildings_assemble_multipolygon_relations(monkeypatch):
    outer = _square(11.97, 57.70, 0.001)
    # Outer ring split over two member ways, plus one courtyard.
    outer_a = outer[:3]
    outer_b = outer[2:]
    courtyard = _square(11.9703, 57.7003, 0.0004)

    def fake_query(query):
        return {
            "elements": [
                {
                    "type": "relation",
                    "id": 77,
                    "members": [
                        {"type": "way", "ref": 1, "role": "outer", "geometry": outer_a},
                        {"type": "way", "ref": 2, "role": "outer", "geometry": outer_b},
                        {
                            "type": "way",
                            "ref": 3,
                            "role": "inner",
                            "geometry": courtyard,
                        },
                        {"type": "node", "ref": 9, "role": "entrance"},
                    ],
                    "tags": {
                        "type": "multipolygon",
                        "building": "university",
                        "building:levels": "4",
                    },
                },
                {
                    "type": "relation",
                    "id": 78,
                    "members": [
                        {"type": "way", "ref": 4, "role": "outer", "geometry": outer_a},
                    ],
                    "tags": {"type": "multipolygon", "building": "yes"},
                },
            ],
        }

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fake_query)

    buildings = overpass.download_overpass_buildings(BBOX)

    # Relation 78 has no closed outer ring and is skipped.
    assert buildings["osm_id"].tolist() == [77]
    assert buildings["osm_type"].tolist() == ["relation"]
    assert buildings["building"].tolist() == ["university"]
    polygon = buildings.geometry.iloc[0]
    assert polygon.geom_type == "Polygon"
    assert len(polygon.interiors) == 1


def test_building_cache_ignores_legacy_records(monkeypatch, tmp_path):
    saved_records = []
    seen_records = []

    class FakeBuildings:
        attrs = {}
        columns = []

        def __len__(self):
            return 0

        def to_file(self, filename, layer, driver):
            return None

    monkeypatch.setattr(overpass, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        overpass,
        "load_cache_metadata",
        lambda: [
            {
                "type": "buildings",
                "bbox": [0.0, 0.0, 100.0, 100.0],
                "filepath": "legacy.gpkg",
                "layer": "buildings",
            }
        ],
    )
    monkeypatch.setattr(
        overpass,
        "find_superset_record",
        lambda bbox, records: seen_records.extend(records),
    )
    monkeypatch.setattr(
        overpass, "download_overpass_buildings", lambda bbox: FakeBuildings()
    )
    monkeypatch.setattr(
        overpass,
        "save_cache_metadata",
        lambda records: saved_records.extend(records),
    )

    overpass.get_buildings_for_bbox((20.0, 20.0, 30.0, 30.0))

    assert seen_records == []
    assert saved_records[-1]["version"] == overpass.BUILDING_CACHE_VERSION


@pytest.mark.parametrize(
    "download", ["download_overpass_buildings", "download_overpass_roads"]
)
def test_overpass_query_box_covers_all_bbox_corners(monkeypatch, download):
    queries = []

    def fake_query(query):
        queries.append(query)
        return {"elements": []}

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fake_query)
    bbox = (319720.0, 6397660.0, 320220.0, 6398160.0)

    getattr(overpass, download)(bbox)

    match = re.search(r"\(([-\d.]+),([-\d.]+),([-\d.]+),([-\d.]+)\)", queries[0])
    south, west, north, east = map(float, match.groups())
    to_wgs84 = pyproj.Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)
    for x in (bbox[0], bbox[2]):
        for y in (bbox[1], bbox[3]):
            lon, lat = to_wgs84.transform(x, y)
            assert west <= lon <= east
            assert south <= lat <= north


def test_relation_keeps_outer_ring_inside_courtyard():
    def ring(x0, y0, size):
        return [
            {"lon": x0, "lat": y0},
            {"lon": x0 + size, "lat": y0},
            {"lon": x0 + size, "lat": y0 + size},
            {"lon": x0, "lat": y0 + size},
            {"lon": x0, "lat": y0},
        ]

    relation = {
        "type": "relation",
        "id": 1,
        "members": [
            {"type": "way", "role": "outer", "geometry": ring(0, 0, 10)},
            {"type": "way", "role": "inner", "geometry": ring(2, 2, 6)},
            {"type": "way", "role": "outer", "geometry": ring(4, 4, 2)},
        ],
    }

    shape = overpass._relation_polygon_lonlat(relation)

    assert shape.area == pytest.approx(100 - 36 + 4)
    assert shape.geom_type == "MultiPolygon"


def _rectangle(xmin, ymin, xmax, ymax):
    corners = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]
    return [{"lon": x, "lat": y} for x, y in corners]


@pytest.mark.parametrize(
    "courtyard",
    [(2, 2, 8, 8), (0, 3, 4, 7), (-1, 3, 4, 7)],
    ids=["inside", "touching-wall", "crossing-wall"],
)
def test_relation_cuts_courtyards_touching_or_crossing_the_wall(courtyard):
    relation = {
        "type": "relation",
        "id": 1,
        "members": [
            {"type": "way", "role": "outer", "geometry": _rectangle(0, 0, 10, 10)},
            {"type": "way", "role": "inner", "geometry": _rectangle(*courtyard)},
        ],
    }
    xmin, ymin, xmax, ymax = courtyard
    inside_wall = (min(xmax, 10) - max(xmin, 0)) * (min(ymax, 10) - max(ymin, 0))

    shape = overpass._relation_polygon_lonlat(relation)

    assert shape.area == pytest.approx(100 - inside_wall)


def _outline_relation(relation_id, outline_ref, outline):
    return {
        "type": "relation",
        "id": relation_id,
        "members": [
            {"type": "way", "ref": outline_ref, "role": "outer", "geometry": outline}
        ],
        "tags": {"type": "multipolygon", "building": "yes"},
    }


def _building_way(way_id, outline):
    return {
        "type": "way",
        "id": way_id,
        "nodes": [1, 2, 3, 4, 1],
        "geometry": outline,
        "tags": {"building": "yes"},
    }


def test_building_way_that_outlines_a_relation_is_not_duplicated(monkeypatch):
    outline = _square(11.97, 57.70, 0.0005)
    monkeypatch.setattr(
        overpass,
        "query_overpass_with_failover",
        lambda query: {
            "elements": [_building_way(5, outline), _outline_relation(9, 5, outline)]
        },
    )

    buildings = overpass.download_overpass_buildings(BBOX)

    assert buildings["osm_type"].tolist() == ["relation"]
    assert buildings["osm_id"].tolist() == [9]


def test_building_way_is_kept_when_its_relation_cannot_be_assembled(monkeypatch):
    outline = _square(11.97, 57.70, 0.0005)
    monkeypatch.setattr(
        overpass,
        "query_overpass_with_failover",
        lambda query: {
            "elements": [
                _building_way(5, outline),
                _outline_relation(9, 5, outline[:3]),
            ]
        },
    )

    buildings = overpass.download_overpass_buildings(BBOX)

    assert buildings["osm_type"].tolist() == ["way"]
    assert buildings["osm_id"].tolist() == [5]


def _tagged(element, **tags):
    element["tags"] = {**element["tags"], **tags}
    return element


def test_relation_row_takes_tags_from_its_outline_way(monkeypatch):
    outline = _square(11.97, 57.70, 0.0005)
    way = _tagged(
        _building_way(5, outline),
        **{
            "building": "apartments",
            "building:levels": "5",
            "height": "17",
            "roof:shape": "gabled",
        },
    )
    relation = _tagged(_outline_relation(9, 5, outline), **{"roof:shape": "flat"})
    monkeypatch.setattr(
        overpass,
        "query_overpass_with_failover",
        lambda query: {"elements": [way, relation]},
    )

    buildings = overpass.download_overpass_buildings(BBOX)

    row = buildings.iloc[0]
    assert len(buildings) == 1
    assert row["osm_id"] == 9
    assert row["building"] == "apartments"
    assert row["building_levels"] == "5"
    assert row["osm_height"] == "17"
    assert row["roof_shape"] == "flat"


def test_relation_row_skips_tags_its_outline_ways_disagree_on(monkeypatch):
    left = _square(11.97, 57.70, 0.0005)
    right = _square(11.9705, 57.70, 0.0005)
    relation = {
        "type": "relation",
        "id": 9,
        "members": [
            {"type": "way", "ref": 5, "role": "outer", "geometry": left},
            {"type": "way", "ref": 6, "role": "outer", "geometry": right},
        ],
        "tags": {"type": "multipolygon", "building": "yes"},
    }
    monkeypatch.setattr(
        overpass,
        "query_overpass_with_failover",
        lambda query: {
            "elements": [
                _tagged(_building_way(5, left), **{"building:levels": "3"}),
                _tagged(_building_way(6, right), **{"building:levels": "4"}),
                relation,
            ]
        },
    )

    buildings = overpass.download_overpass_buildings(BBOX)

    assert buildings["osm_id"].tolist() == [9]
    assert buildings["building_levels"].isna().iloc[0]


def test_way_sharing_an_id_with_a_relation_is_kept(monkeypatch):
    outline = _square(11.97, 57.70, 0.0005)
    other = _square(11.971, 57.70, 0.0005)
    monkeypatch.setattr(
        overpass,
        "query_overpass_with_failover",
        lambda query: {
            "elements": [
                _building_way(5, outline),
                _building_way(9, other),
                _outline_relation(9, 5, outline),
            ]
        },
    )

    buildings = overpass.download_overpass_buildings(BBOX)

    assert sorted(zip(buildings["osm_type"], buildings["osm_id"], strict=True)) == [
        ("relation", 9),
        ("way", 9),
    ]


@pytest.mark.parametrize(
    "download", ["download_overpass_buildings", "download_overpass_roads"]
)
def test_overpass_rejects_invalid_bbox_before_querying(monkeypatch, download):
    def fail_query(query):
        raise AssertionError("an invalid bbox must not reach Overpass")

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fail_query)

    with pytest.raises(ValueError, match="xmin < xmax"):
        getattr(overpass, download)((320000.0, 6398000.0, 320000.0, 6398100.0))
