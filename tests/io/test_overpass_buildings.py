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
                        {"type": "way", "ref": 3, "role": "inner", "geometry": courtyard},
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
