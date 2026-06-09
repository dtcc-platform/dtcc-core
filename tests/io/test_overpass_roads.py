import pytest

from dtcc_core.io.data import overpass
from dtcc_core.io.data import wrapper
from dtcc_core.model import Bounds


def test_overpass_failover_does_not_warn_when_fallback_succeeds(monkeypatch):
    warnings = []

    class Response:
        def __init__(self, status_code, payload=None):
            self.status_code = status_code
            self._payload = payload or {}

        def json(self):
            return self._payload

    class Session:
        def __init__(self):
            self.responses = [
                Response(406),
                Response(200, {"elements": []}),
            ]

        def post(self, endpoint, data, timeout):
            return self.responses.pop(0)

    session = Session()
    monkeypatch.setattr(overpass, "OVERPASS_ENDPOINTS", ["endpoint-1", "endpoint-2"])
    monkeypatch.setattr(overpass, "create_retry_session", lambda: session)
    monkeypatch.setattr(overpass, "warning", lambda message: warnings.append(message))

    data = overpass.query_overpass_with_failover("[out:json];")

    assert data["elements"] == []
    assert data[overpass.OVERPASS_ENDPOINT_METADATA_KEY] == "endpoint-2"
    assert warnings == []


def test_overpass_uses_bounded_connect_and_read_timeouts(monkeypatch):
    timeouts = []

    class Response:
        status_code = 200

        def json(self):
            return {"elements": []}

    class Session:
        def post(self, endpoint, data, timeout):
            timeouts.append(timeout)
            return Response()

    monkeypatch.setattr(overpass, "OVERPASS_ENDPOINTS", ["endpoint-1"])
    monkeypatch.setattr(overpass, "create_retry_session", lambda: Session())

    data = overpass.query_overpass_with_failover(
        "[out:json];",
        connect_timeout=2.5,
        read_timeout=7.5,
    )

    assert data["elements"] == []
    assert timeouts == [(2.5, 7.5)]


def test_overpass_retry_session_sets_user_agent():
    session = overpass.create_retry_session()

    assert session.headers["User-Agent"] == overpass.OVERPASS_USER_AGENT
    assert session.headers["Accept"] == "application/json"


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
