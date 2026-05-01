import pytest

from dtcc_core.io.data import overpass
from dtcc_core.io.data import wrapper
from dtcc_core.model import Bounds


def test_overpass_roads_are_segmented_with_traffic_attributes(monkeypatch):
    def fake_query(query):
        return {
            "elements": [
                {"type": "node", "id": 1, "lat": 57.7000, "lon": 11.9700},
                {"type": "node", "id": 2, "lat": 57.7005, "lon": 11.9705},
                {"type": "node", "id": 3, "lat": 57.7010, "lon": 11.9710},
                {
                    "type": "way",
                    "id": 100,
                    "nodes": [1, 2, 3],
                    "tags": {
                        "highway": "residential",
                        "name": "Testgatan",
                        "lanes": "2",
                        "maxspeed": "30 mph",
                        "oneway": "-1",
                        "bridge": "yes",
                        "tunnel": "no",
                        "junction": "roundabout",
                    },
                },
            ]
        }

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fake_query)

    roads = overpass.download_overpass_roads((319900.0, 6398900.0, 320100.0, 6399100.0))

    assert len(roads) == 2
    assert roads["osm_way_id"].tolist() == [100, 100]
    assert roads["osm_start_node_id"].tolist() == [3, 2]
    assert roads["osm_end_node_id"].tolist() == [2, 1]
    assert roads["segment_index"].tolist() == [0, 1]
    assert roads["highway"].tolist() == ["residential", "residential"]
    assert roads["name"].tolist() == ["Testgatan", "Testgatan"]
    assert roads["lanes"].tolist() == ["2", "2"]
    assert roads["maxspeed"].tolist() == ["30 mph", "30 mph"]
    assert roads["maxspeed_kmh"].iloc[0] == pytest.approx(48.28032)
    assert roads["oneway"].tolist() == [True, True]
    assert roads["oneway_raw"].tolist() == ["-1", "-1"]
    assert roads["bridge"].tolist() == ["yes", "yes"]
    assert roads["tunnel"].tolist() == ["no", "no"]
    assert roads["junction"].tolist() == ["roundabout", "roundabout"]


def test_road_cache_ignores_legacy_records(monkeypatch):
    saved_records = []

    class FakeRoads:
        def __len__(self):
            return 0

        def to_file(self, filename, layer, driver):
            return None

    def fake_find_superset_record(bbox, records):
        assert records == [
            {
                "type": "roads",
                "version": overpass.ROAD_CACHE_VERSION,
                "bbox": [0.0, 0.0, 10.0, 10.0],
                "filepath": "current.gpkg",
                "layer": "roads",
            }
        ]
        return None

    monkeypatch.setattr(
        overpass,
        "load_cache_metadata",
        lambda: [
            {
                "type": "roads",
                "bbox": [0.0, 0.0, 10.0, 10.0],
                "filepath": "legacy.gpkg",
                "layer": "roads",
            },
            {
                "type": "roads",
                "version": overpass.ROAD_CACHE_VERSION,
                "bbox": [0.0, 0.0, 10.0, 10.0],
                "filepath": "current.gpkg",
                "layer": "roads",
            },
        ],
    )
    monkeypatch.setattr(overpass, "find_superset_record", fake_find_superset_record)
    monkeypatch.setattr(overpass, "download_overpass_roads", lambda bbox: FakeRoads())
    monkeypatch.setattr(
        overpass,
        "save_cache_metadata",
        lambda records: saved_records.extend(records),
    )

    overpass.get_roads_for_bbox((20.0, 20.0, 30.0, 30.0))

    assert saved_records[-1]["version"] == overpass.ROAD_CACHE_VERSION


def test_download_roadnetwork_defaults_to_osm(monkeypatch):
    calls = []
    bounds = Bounds(0.0, 0.0, 1.0, 1.0)

    def fake_download_data(data_type, provider, bounds, epsg="3006"):
        calls.append((data_type, provider, bounds, epsg))
        return "roads"

    monkeypatch.setattr(wrapper, "download_data", fake_download_data)

    assert wrapper.download_roadnetwork(bounds) == "roads"
    assert calls == [("roads", "OSM", bounds, "3006")]
