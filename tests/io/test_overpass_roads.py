import pytest

from dtcc_core.io.data import overpass, wrapper
from dtcc_core.model import Bounds


def test_overpass_failover_does_not_warn_when_fallback_succeeds(monkeypatch):
    from dtcc_core.io import HttpTransportError, Response

    warnings, calls = [], []
    endpoints = ["https://first.test", "https://second.test"]

    def request(client, req, **options):
        calls.append(req)
        if len(calls) == 1:
            raise HttpTransportError(req, "http_status", status=406)
        return Response(200, req.target, {}, b'{"elements": []}', req)

    monkeypatch.setattr(overpass.HttpClient, "request", request)
    monkeypatch.setattr(overpass, "OVERPASS_ENDPOINTS", endpoints)
    monkeypatch.setattr(overpass, "warning", warnings.append)
    data = overpass.query_overpass_with_failover("[out:json];")
    assert data["elements"] == []
    assert data[overpass.OVERPASS_ENDPOINT_METADATA_KEY] == endpoints[1]
    assert warnings == []
    assert all(req.idempotent for req in calls)
    assert dict(calls[0].headers) == {
        "User-Agent": overpass.OVERPASS_USER_AGENT,
        "Accept": "application/json",
    }
    assert dict(calls[0].form) == {"data": "[out:json];"}


def test_overpass_uses_bounded_connect_and_read_timeouts():
    with overpass.create_retry_session(connect_timeout=2.5, read_timeout=7.5) as client:
        assert client._state.timeout.connect == 2.5
        assert client._state.timeout.read == 7.5
        assert not client._state.retry.connect
        assert not client._state.retry.read


def test_overpass_roads_are_segmented_with_traffic_attributes(monkeypatch):
    def fake_query(query):
        assert f"[timeout:{overpass.OVERPASS_SERVER_TIMEOUT_SECONDS}]" in query
        return {
            overpass.OVERPASS_ENDPOINT_METADATA_KEY: "endpoint-2",
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
            ],
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
    assert roads["source_endpoint"].tolist() == ["endpoint-2", "endpoint-2"]
    assert roads.attrs["source_endpoint"] == "endpoint-2"


def test_overpass_roads_parse_direct_geometry_payload(monkeypatch):
    def fake_query(query):
        assert "out geom;" in query
        return {
            overpass.OVERPASS_ENDPOINT_METADATA_KEY: "endpoint-1",
            "elements": [
                {
                    "type": "way",
                    "id": 101,
                    "nodes": [10, 11],
                    "geometry": [
                        {"lat": 57.7000, "lon": 11.9700},
                        {"lat": 57.7005, "lon": 11.9705},
                    ],
                    "tags": {
                        "highway": "primary",
                        "oneway": "yes",
                    },
                },
            ],
        }

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fake_query)

    roads = overpass.download_overpass_roads((319900.0, 6398900.0, 320100.0, 6399100.0))

    assert len(roads) == 1
    assert roads["osm_start_node_id"].tolist() == [10]
    assert roads["osm_end_node_id"].tolist() == [11]
    assert roads["highway"].tolist() == ["primary"]
    assert roads["oneway"].tolist() == [True]


def test_road_cache_ignores_legacy_records(monkeypatch):
    saved_records = []

    class FakeRoads:
        attrs = {"source_endpoint": "endpoint-2"}

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
    assert saved_records[-1]["source_endpoint"] == "endpoint-2"


def test_download_roadnetwork_defaults_to_osm(monkeypatch):
    calls = []
    bounds = Bounds(0.0, 0.0, 1.0, 1.0)

    def fake_download_data(data_type, provider, bounds, epsg="3006"):
        calls.append((data_type, provider, bounds, epsg))
        return "roads"

    monkeypatch.setattr(wrapper, "download_data", fake_download_data)

    assert wrapper.download_roadnetwork(bounds) == "roads"
    assert calls == [("roads", "OSM", bounds, "3006")]
