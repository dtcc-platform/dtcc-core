from urllib.parse import urlparse

import pytest

from dtcc_core.io.data import overpass, wrapper
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


def _slow_server_session(seconds_to_answer, timeouts):
    """Fake session whose server answers after ``seconds_to_answer``."""
    import requests

    class Response:
        status_code = 200

        def json(self):
            return {"elements": []}

    class Session:
        def post(self, endpoint, data, timeout):
            timeouts.append(timeout)
            if timeout[1] < seconds_to_answer:
                raise requests.exceptions.ReadTimeout("Read timed out.")
            return Response()

    return Session()


def test_overpass_waits_for_a_server_slower_than_30_seconds(monkeypatch):
    timeouts = []
    monkeypatch.setattr(overpass, "OVERPASS_ENDPOINTS", ["endpoint-1"])
    monkeypatch.setattr(overpass, "OVERPASS_READ_TIMEOUT_SECONDS", None)
    monkeypatch.setattr(
        overpass, "create_retry_session", lambda: _slow_server_session(43, timeouts)
    )

    data = overpass.query_overpass_with_failover(
        f"[out:json][timeout:{overpass.OVERPASS_SERVER_TIMEOUT_SECONDS}];node(1);out;"
    )

    assert data["elements"] == []
    assert len(timeouts) == 1


@pytest.mark.parametrize(
    "query, server_timeout",
    [
        ("[out:json][timeout:25];node(1);out;", 25),
        ("[out:json] [ timeout : 120 ];node(1);out;", 120),
        ("[out:json];node(1);out;", 180),  # Overpass server default
    ],
)
def test_overpass_read_timeout_follows_query_timeout(
    monkeypatch, query, server_timeout
):
    monkeypatch.setattr(overpass, "OVERPASS_READ_TIMEOUT_SECONDS", None)
    monkeypatch.setattr(overpass, "OVERPASS_QUEUE_ALLOWANCE_SECONDS", 60.0)

    assert overpass._overpass_read_timeout(query) == server_timeout + 60.0


def test_overpass_read_timeout_environment_override_is_kept(monkeypatch):
    timeouts = []
    monkeypatch.setattr(overpass, "OVERPASS_ENDPOINTS", ["endpoint-1"])
    monkeypatch.setattr(overpass, "OVERPASS_READ_TIMEOUT_SECONDS", 12.0)
    monkeypatch.setattr(
        overpass, "create_retry_session", lambda: _slow_server_session(0, timeouts)
    )

    overpass.query_overpass_with_failover("[out:json][timeout:25];node(1);out;")

    assert timeouts == [(overpass.OVERPASS_CONNECT_TIMEOUT_SECONDS, 12.0)]


def test_overpass_endpoints_are_distinct_servers():
    hosts = [urlparse(endpoint).hostname for endpoint in overpass.OVERPASS_ENDPOINTS]

    assert len(hosts) == len(set(hosts))
    # DNS alias of overpass.private.coffee.
    assert "overpass.kumi.systems" not in hosts


class _FakeClock:
    """Stands in for the ``time`` module so tests never really wait."""

    def __init__(self):
        self.now = 0.0
        self.sleeps = []

    def monotonic(self):
        return self.now

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds


def _scripted_session(clock, answers, posts):
    """Fake session that answers each endpoint from a script.

    ``answers`` maps an endpoint to a list of ``(status, seconds)`` pairs, or
    ``"timeout"`` for a read timeout, used in order.
    """
    import requests

    class Response:
        def __init__(self, status_code):
            self.status_code = status_code

        def json(self):
            return {"elements": []}

    class Session:
        def post(self, endpoint, data, timeout):
            posts.append(endpoint)
            answer = answers[endpoint].pop(0)
            if answer == "timeout":
                clock.now += timeout[1]
                raise requests.exceptions.ReadTimeout("Read timed out.")
            status, seconds = answer
            clock.now += seconds
            return Response(status)

    return Session()


def _install_scripted_session(monkeypatch, answers):
    clock = _FakeClock()
    posts = []
    monkeypatch.setattr(overpass, "OVERPASS_ENDPOINTS", list(answers))
    monkeypatch.setattr(overpass, "OVERPASS_RETRIES", 2)
    monkeypatch.setattr(overpass, "OVERPASS_QUICK_FAILURE_SECONDS", 10.0)
    monkeypatch.setattr(overpass, "time", clock)
    monkeypatch.setattr(
        overpass,
        "create_retry_session",
        lambda: _scripted_session(clock, answers, posts),
    )
    return clock, posts


def test_overpass_retries_a_quick_502_on_the_same_endpoint(monkeypatch):
    clock, posts = _install_scripted_session(
        monkeypatch,
        {"endpoint-1": [(502, 0.3), (200, 42.8)], "endpoint-2": []},
    )

    data = overpass.query_overpass_with_failover("[out:json][timeout:25];")

    assert data[overpass.OVERPASS_ENDPOINT_METADATA_KEY] == "endpoint-1"
    assert posts == ["endpoint-1", "endpoint-1"]
    assert clock.sleeps == [1.0]


@pytest.mark.parametrize(
    "first_answer",
    [(504, 80.0), "timeout"],
    ids=["slow 504", "read timeout"],
)
def test_overpass_moves_on_after_a_slow_failure(monkeypatch, first_answer):
    clock, posts = _install_scripted_session(
        monkeypatch,
        {"endpoint-1": [first_answer], "endpoint-2": [(200, 1.0)]},
    )

    data = overpass.query_overpass_with_failover("[out:json][timeout:25];")

    assert data[overpass.OVERPASS_ENDPOINT_METADATA_KEY] == "endpoint-2"
    assert posts == ["endpoint-1", "endpoint-2"]
    assert clock.sleeps == []


def test_overpass_gives_up_after_the_retry_limit(monkeypatch):
    clock, posts = _install_scripted_session(
        monkeypatch,
        {"endpoint-1": [(502, 0.3)] * 3, "endpoint-2": [(429, 0.2)] * 3},
    )

    with pytest.raises(RuntimeError, match="HTTP 429"):
        overpass.query_overpass_with_failover("[out:json][timeout:25];")

    assert posts == ["endpoint-1"] * 3 + ["endpoint-2"] * 3
    assert clock.sleeps == [1.0, 2.0, 1.0, 2.0]


def test_overpass_session_leaves_retries_to_the_failover_loop():
    session = overpass.create_retry_session()

    assert session.get_adapter("https://example.org").max_retries.total == 0


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
                        "surface": "asphalt",
                        "lit": "yes",
                        "sidewalk:left": "separate",
                        "cycleway:right": "track",
                        "width": "7",
                        "access": "destination",
                    },
                },
                {
                    "type": "node",
                    "id": 2,
                    "lat": 57.7005,
                    "lon": 11.9705,
                    "tags": {"highway": "crossing", "crossing": "marked"},
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
    assert roads["surface"].tolist() == ["asphalt", "asphalt"]
    assert roads["lit"].tolist() == ["yes", "yes"]
    assert roads["sidewalk_left"].tolist() == ["separate", "separate"]
    assert roads["cycleway_right"].tolist() == ["track", "track"]
    assert roads["width"].tolist() == ["7", "7"]
    assert roads["access"].tolist() == ["destination", "destination"]
    assert roads["sidewalk"].isna().all()
    assert roads["service"].isna().all()
    # Segments run 3->2 and 2->1 because of oneway=-1; node 2 is a crossing.
    assert roads["end_node_highway"].fillna("").tolist() == ["crossing", ""]
    assert roads["end_node_crossing"].fillna("").tolist() == ["marked", ""]
    assert roads["start_node_highway"].fillna("").tolist() == ["", "crossing"]
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


def test_overpass_failure_lists_every_endpoint(monkeypatch):
    import requests

    class Response:
        status_code = 400

    class Session:
        def post(self, endpoint, data, timeout):
            if endpoint == "endpoint-3":
                raise requests.Timeout("Read timed out.")
            return Response()

    monkeypatch.setattr(
        overpass, "OVERPASS_ENDPOINTS", ["endpoint-1", "endpoint-2", "endpoint-3"]
    )
    monkeypatch.setattr(overpass, "create_retry_session", lambda: Session())

    with pytest.raises(RuntimeError) as info:
        overpass.query_overpass_with_failover("[out:json];")

    message = str(info.value)
    assert message.startswith(
        "All Overpass endpoints failed. Last error: Read timed out."
    )
    assert "- endpoint-1: HTTP 400" in message
    assert "- endpoint-2: HTTP 400" in message
    assert "- endpoint-3: Read timed out." in message
