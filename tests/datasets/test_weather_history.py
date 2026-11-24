"""Historical weather selection, source semantics and location reprojection."""

import copy
import json
from urllib.parse import parse_qsl, urlsplit, urlunsplit

import httpx2
import numpy as np
import pytest
from pyproj import Geod, Transformer

from dtcc_core import datasets
from dtcc_core.datasets import _weather_history as common
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.weather_history import (
    NoMatchingStation,
    WeatherHistoryDatasetArgs,
    WeatherHistoryLimitError,
)
from dtcc_core.io._http._httpx import NativeHttpClient
from dtcc_core.model import Bounds, Point, TimeSeries, TimeSeriesCollection

SMHI = "https://opendata-download-metobs.smhi.se/api/version/1.0"
OPEN = "https://archive-api.open-meteo.com/v1/archive"
NAMES = ["air_temperature", "wind_speed", "wind_direction"]
START = int(np.datetime64("2020-01-01T00:00", "ms").astype(np.int64))
END = START + 86_400_000


def resource(url, rel):
    return {"link": [{"href": url, "rel": rel, "type": "application/json"}]}


def csv_history(unit, value=10.0, quality="G"):
    rows = [
        "\ufeffStationsnamn;Stationsnummer;Stationsnät;Mäthöjd (meter över marken)",
        "Test station;20;SMHI;10",
        "",
        "Parameternamn;Beskrivning;Enhet",
        f"Measurement;hourly;{unit}",
        "",
        "Datum;Tid (UTC);Measurement;Kvalitet;;Notes",
    ]
    rows.extend(
        f"2020-01-01;{hour:02}:00:00;{value};{quality};;comment" for hour in range(24)
    )
    return "\n".join(rows)


class Response(httpx2.SyncByteStream):
    def __init__(self, body, status=200, headers=None):
        self.body = body
        self.status_code = status
        self.headers = headers or {}
        self.closed = False

    def close(self):
        self.closed = True

    def __iter__(self):
        yield from self.iter_content(65536)

    def iter_content(self, chunk_size):
        body = self.body
        if isinstance(body, dict):
            body = json.dumps(body).encode()
        elif isinstance(body, str):
            body = body.encode("utf-8")
        for offset in range(0, len(body), chunk_size):
            yield body[offset : offset + chunk_size]


@pytest.fixture
def server(monkeypatch):
    responses, calls = {}, []
    common._CACHE.clear()
    variables = {
        "temperature_2m": [10.0] * 24,
        "wind_speed_10m": [3.0] * 24,
        "wind_direction_10m": [270.0] * 24,
    }
    responses[OPEN] = {
        "latitude": 59.25,
        "longitude": 18.0,
        "elevation": 37,
        "utc_offset_seconds": 0,
        "hourly": {"time": list(range(START // 1000, END // 1000, 3600)), **variables},
        "hourly_units": {
            "time": "unixtime",
            "temperature_2m": "°C",
            "wind_speed_10m": "m/s",
            "wind_direction_10m": "°",
        },
    }
    responses["https://customer-archive-api.open-meteo.com/v1/archive"] = responses[
        OPEN
    ]
    version = {"resource": []}
    responses[
        common_url := "https://opendata-download-metobs.smhi.se/api/version/latest.json"
    ] = version
    for pid, unit in ((1, "celsius"), (4, "meter per sekund"), (3, "grader")):
        parameter_url = f"{SMHI}/parameter/{pid}.json"
        version["resource"].append(
            {"key": str(pid), **resource(parameter_url, "parameter")}
        )
        catalog = {
            "unit": unit,
            "valueType": "SAMPLING",
            "summary": "hourly",
            "station": [],
        }
        responses[parameter_url] = catalog
        # Station 10 is closest but only offers temperature. Station 20 offers all.
        for sid, lon, lat in ((10, 18.05, 59.34), (20, 18.06, 59.35), (30, 18.2, 59.5)):
            if sid == 10 and pid != 1:
                continue
            station_url = f"{SMHI}/parameter/{pid}/station/{sid}.json"
            catalog["station"].append(
                {
                    "key": str(sid),
                    "from": START - 86_400_000,
                    "to": END,
                    "active": False,
                    "longitude": lon,
                    "latitude": lat,
                    **resource(station_url, "station"),
                }
            )
            station = {
                "key": str(sid),
                "title": f"Station {sid}",
                "period": [],
                "position": [
                    {
                        "from": START - 86_400_000,
                        "to": END,
                        "longitude": lon,
                        "latitude": lat,
                        "height": 10,
                    }
                ],
            }
            responses[station_url] = station
            for period_name in ("latest-months", "corrected-archive"):
                period_url = (
                    f"{SMHI}/parameter/{pid}/station/{sid}/period/{period_name}.json"
                )
                data_url = period_url.removesuffix(".json") + "/data.csv"
                station["period"].append(
                    {"key": period_name, **resource(period_url, "period")}
                )
                responses[period_url] = {
                    "from": START,
                    "to": END,
                    "updated": END,
                    "data": [
                        {
                            "link": [
                                {"href": data_url, "type": "text/plain", "rel": "data"}
                            ]
                        }
                    ],
                }
                responses[data_url] = csv_history(
                    unit,
                    10 if period_name == "corrected-archive" else 99,
                    "G" if period_name == "corrected-archive" else "Y",
                )

    def send(request):
        assert request.method == "GET"
        parsed = urlsplit(str(request.url))
        url = urlunsplit(parsed._replace(query=""))
        params = dict(parse_qsl(parsed.query)) or None
        if params:
            for name in ("latitude", "longitude"):
                if name in params:
                    params[name] = float(params[name])
        calls.append((url, copy.deepcopy(params)))
        body = responses[url]
        body = body(params) if callable(body) else body
        if isinstance(body, Exception):
            raise body
        body = body if isinstance(body, Response) else Response(body)
        return httpx2.Response(body.status_code, headers=body.headers, stream=body)

    def build_client(self):
        if self._client is None:
            self._client = httpx2.Client(
                transport=httpx2.MockTransport(send), trust_env=False
            )
        return self._client

    monkeypatch.setattr(NativeHttpClient, "_build_client", build_client)
    return responses, calls, common_url


def history(source="open_meteo", **kwargs):
    args = dict(
        source=source,
        location=(18.05, 59.34),
        crs="EPSG:4326",
        start="2020-01-01",
        end="2020-01-02",
        cache_ttl_s=0,
    )
    if source == "smhi":
        args["max_distance_m"] = 50_000
    args.update(kwargs)
    return datasets.weather_history(**args)


@pytest.mark.parametrize("kind", ["tuple", "point", "bounds"])
def test_default_sweref_location_is_reprojected(server, kind):
    x, y = Transformer.from_crs(4326, 3006, always_xy=True).transform(18.05, 59.34)
    location = {
        "tuple": (x, y),
        "point": Point(x=x, y=y),
        "bounds": Bounds(xmin=x - 100, xmax=x + 100, ymin=y - 50, ymax=y + 50),
    }[kind]
    result = datasets.weather_history(
        source="open_meteo",
        location=location,
        start="2020-01-01",
        end="2020-01-02",
        cache_ttl_s=0,
    )
    query = server[1][0][1]
    assert query["latitude"] == pytest.approx(59.34)
    assert query["longitude"] == pytest.approx(18.05)
    assert result.attributes["crs"] == "EPSG:3006"
    requested = result.attributes["requested_point"]
    assert (requested.x, requested.y) == pytest.approx((x, y))
    point = result.locations["open_meteo:location"]["point"]
    assert (point.x, point.y) == pytest.approx(
        Transformer.from_crs(4326, 3006, always_xy=True).transform(18.0, 59.25)
    )
    assert point.transform.srs == "EPSG:3006"


def test_open_meteo_returns_three_series_with_context_and_units(server):
    result = history()
    assert isinstance(result, TimeSeriesCollection) and len(result) == 3
    assert list(result) == [("open_meteo:location", name) for name in NAMES]
    assert all(
        isinstance(ts, TimeSeries) and ts.values.shape == (24, 1)
        for ts in result.values()
    )
    assert {ts.role for ts in result.values()} == {"reanalysis"}
    assert [ts.unit for ts in result.values()] == ["celsius", "m/s", "degrees"]
    assert result.attributes["coverage_complete"]
    assert result.dataset_context.identity.name == "weather_history"
    assert result.dataset_context.metadata.provider[0]["slug"] == "open-meteo"
    assert server[1][0][1]["start_date"] == server[1][0][1]["end_date"] == "2020-01-01"
    table = result.to_dataframe()
    assert len(table) == 72 and str(table.time.dtype) == "datetime64[ns, UTC]"
    assert table.quality.isna().all()
    assert (
        result.series_metadata[("open_meteo:location", "wind_speed")][
            "height_above_ground_m"
        ]
        == 10
    )
    with pytest.raises(NotImplementedError):
        result.to_proto()
    ts = next(iter(result.values()))
    restored = TimeSeries()
    restored.from_proto(ts.to_proto())
    np.testing.assert_array_equal(restored.values, ts.values)


def test_nanosecond_half_open_boundaries_and_aliases(server):
    result = history(
        start="2020-01-01T01:00:00.000000001+01:00",
        end="2020-01-01T03:00+01:00",
        parameters=["temp", 1, "wind"],
    )
    assert len(result) == 2
    assert all(
        len(ts) == 1 and ts.times[0] == np.datetime64("2020-01-01T01:00", "ns")
        for ts in result.values()
    )


def test_smhi_selects_one_common_station_and_preserves_quality(server):
    result = history("smhi")
    selection = result.attributes["selection"]
    assert selection["status"] == "matched" and selection["station_id"] == 20
    distance = Geod(ellps="WGS84").inv(18.05, 59.34, 18.06, 59.35)[2]
    assert selection["distance_m"] == pytest.approx(distance)
    assert selection["station_location"].x == 18.06
    assert selection["station_location"].y == 59.35
    assert all(key[0] == "smhi:20" for key in result)
    assert all(
        ts.role == "observed" and ts.values[0, 0] == 10 for ts in result.values()
    )
    for key in result:
        assert all(
            row["quality"] == "G" and row["source_period"] == "corrected-archive"
            for row in result.observations[key].values()
        )
    assert not any("station/10" in url for url, _ in server[1])
    table = result.to_dataframe()
    assert set(table.station_id) == {20} and set(table.distance_m) == {distance}
    assert result.dataset_context.metadata.provider[0]["slug"] == "smhi"


def test_smhi_radius_exact_boundary_and_no_match(server):
    distance = Geod(ellps="WGS84").inv(18.05, 59.34, 18.06, 59.35)[2]
    assert (
        history("smhi", max_distance_m=distance).attributes["selection"]["station_id"]
        == 20
    )
    empty = history("smhi", max_distance_m=distance - 1)
    assert (
        not empty and empty.attributes["selection"]["status"] == "no_matching_station"
    )
    assert empty.attributes["upstream_error_count"] == 0
    with pytest.raises(NoMatchingStation):
        history("smhi", max_distance_m=1, strict_live=True)


def test_smhi_corrected_missing_reading_is_not_filled_from_provisional(server):
    responses = server[0]
    url = f"{SMHI}/parameter/1/station/20/period/corrected-archive/data.csv"
    responses[url] = responses[url].replace("00:00:00;10;G", "00:00:00;;G")
    result = history("smhi")
    assert np.isnan(result[("smhi:20", "air_temperature")].values[0, 0])
    assert (
        result.attributes["partial_result"]
        and not result.attributes["coverage_complete"]
    )


def test_smhi_default_sweref_point_and_bounds_select_same_station(server):
    x, y = Transformer.from_crs(4326, 3006, always_xy=True).transform(18.05, 59.34)
    for location in (
        Point(x=x, y=y),
        Bounds(xmin=x - 10, xmax=x + 10, ymin=y - 10, ymax=y + 10),
    ):
        result = datasets.weather_history(
            source="smhi",
            location=location,
            start="2020-01-01",
            end="2020-01-02",
            max_distance_m=50_000,
            cache_ttl_s=0,
        )
        assert result.attributes["selection"]["station_id"] == 20
        actual = result.locations["smhi:20"]["point"]
        assert actual.transform.srs == "EPSG:3006"
        assert (actual.x, actual.y) == pytest.approx(
            Transformer.from_crs(4326, 3006, always_xy=True).transform(18.06, 59.35)
        )


@pytest.mark.parametrize("coverage", ["catalog", "period"])
def test_smhi_requires_temporal_coverage_for_every_parameter(server, coverage):
    if coverage == "catalog":
        catalog = server[0][f"{SMHI}/parameter/4.json"]
        next(s for s in catalog["station"] if s["key"] == "20")["to"] = START
    else:
        for name in ("latest-months", "corrected-archive"):
            server[0][f"{SMHI}/parameter/4/station/20/period/{name}.json"]["to"] = START
    assert history("smhi").attributes["selection"]["station_id"] == 30


def test_smhi_equal_distance_tie_uses_station_id(server):
    for pid in (1, 3, 4):
        position = server[0][f"{SMHI}/parameter/{pid}/station/30.json"]["position"][0]
        position.update(longitude=18.06, latitude=59.35)
    assert history("smhi").attributes["selection"]["station_id"] == 20


def test_smhi_moved_station_reports_maximum_distance_and_history(server):
    midpoint = START + 43_200_000
    for pid in (1, 3, 4):
        station = server[0][f"{SMHI}/parameter/{pid}/station/20.json"]
        first = station["position"][0]
        second = {**first, "from": midpoint, "longitude": 18.07, "latitude": 59.36}
        # SMHI's inclusive "to" is the last whole second before the next "from".
        first["to"] = midpoint - 1000
        station["position"].append(second)
    selection = history("smhi").attributes["selection"]
    assert selection["station_id"] == 20 and len(selection["position_history"]) == 2
    assert selection["distance_m"] == pytest.approx(
        Geod(ellps="WGS84").inv(18.05, 59.34, 18.07, 59.36)[2]
    )
    assert selection["station_location"].x == 18.06
    assert (
        history("smhi", max_distance_m=1500).attributes["selection"]["status"]
        == "no_matching_station"
    )


def test_smhi_unknown_position_excludes_station_and_reports_error(server):
    del server[0][f"{SMHI}/parameter/1/station/20.json"]["position"]
    result = history("smhi")
    selection = result.attributes["selection"]
    assert selection["status"] == "matched" and selection["station_id"] == 30
    assert selection["excluded_stations"] == [
        {"station_id": 20, "parameter": "air_temperature", "reason": "invalid_metadata"}
    ]
    assert result.attributes["upstream_errors"][0]["failure_class"] == "invalid_payload"
    with pytest.raises(DatasetUpstreamError):
        history("smhi", strict_live=True)


@pytest.mark.parametrize(
    ("problem", "reason"),
    [("conflicting_positions", "invalid_metadata"), ("http_503", "fetch_failed")],
)
def test_smhi_bad_station_is_excluded_without_failing_selection(
    server, problem, reason
):
    url = f"{SMHI}/parameter/4/station/30.json"
    if problem == "http_503":
        server[0][url] = Response(b"", status=503)
    else:
        # Live station 53430 has overlapping positions at different coordinates.
        position = server[0][url]["position"][0]
        server[0][url]["position"].append(
            {**position, "from": START, "longitude": 18.3, "latitude": 59.6}
        )
    result = history("smhi")
    selection = result.attributes["selection"]
    assert selection["status"] == "matched" and selection["station_id"] == 20
    assert selection["excluded_stations"] == [
        {"station_id": 30, "parameter": "wind_speed", "reason": reason}
    ]
    assert result.attributes["upstream_error_count"] == 1
    with pytest.raises(DatasetUpstreamError):
        history("smhi", strict_live=True)


def test_smhi_failed_nearby_station_does_not_claim_no_station(server):
    server[0][f"{SMHI}/parameter/1/station/20.json"] = Response(b"", status=503)
    result = history("smhi", max_distance_m=1_300)
    assert not result and result.attributes["selection"]["status"] == "discovery_failed"
    assert result.attributes["upstream_error_count"] == 1


def test_smhi_period_failure_excludes_station_and_selects_next(server):
    for name in ("latest-months", "corrected-archive"):
        url = f"{SMHI}/parameter/3/station/20/period/{name}.json"
        server[0][url] = Response(b"", status=503)
    result = history("smhi")
    selection = result.attributes["selection"]
    assert selection["status"] == "matched" and selection["station_id"] == 30
    assert selection["excluded_stations"] == [
        {"station_id": 20, "parameter": "wind_direction", "reason": "fetch_failed"}
    ]
    assert result.attributes["upstream_error_count"] == 1


def test_smhi_skips_station_documents_beyond_relocation_margin(server):
    # Station 30's catalogue position is max_distance_m + 20 km + 1 m away.
    lon, lat, _ = Geod(ellps="WGS84").fwd(18.05, 59.34, 90, 21_301)
    for pid in (1, 3, 4):
        catalog = server[0][f"{SMHI}/parameter/{pid}.json"]
        entry = next(s for s in catalog["station"] if s["key"] == "30")
        entry.update(longitude=lon, latitude=lat)
    result = history("smhi", max_distance_m=1_300)
    assert result.attributes["selection"]["station_id"] == 20
    assert not any("station/30" in url for url, _ in server[1])


def test_smhi_finds_station_that_moved_within_relocation_margin(server):
    # Station 30 used to stand 1 km from the query; its catalogue position is
    # 15 km beyond max_distance_m, like the largest SMHI relocation (15.1 km).
    geod = Geod(ellps="WGS84")
    old_lon, old_lat, _ = geod.fwd(18.05, 59.34, 0, 1_000)
    new_lon, new_lat, _ = geod.fwd(18.05, 59.34, 0, 17_000)
    for pid in (1, 3, 4):
        catalog = server[0][f"{SMHI}/parameter/{pid}.json"]
        entry = next(s for s in catalog["station"] if s["key"] == "30")
        entry.update(longitude=new_lon, latitude=new_lat)
        position = server[0][f"{SMHI}/parameter/{pid}/station/30.json"]["position"][0]
        position.update(longitude=old_lon, latitude=old_lat)
    selection = history("smhi", max_distance_m=2_000).attributes["selection"]
    assert selection["station_id"] == 30
    assert selection["distance_m"] == pytest.approx(1_000)


def test_smhi_catalogue_time_extent_does_not_hide_position_coverage_gaps(server):
    # A station can have older and newer records without operating in between.
    for pid in (1, 3, 4):
        station = server[0][f"{SMHI}/parameter/{pid}/station/20.json"]
        position = station["position"][0]
        station["position"] = [
            {**position, "to": START - 1},
            {**position, "from": END},
        ]
    result = history("smhi")
    assert result.attributes["selection"]["station_id"] == 30
    assert result.attributes["selection"]["excluded_stations"] == [
        {
            "station_id": 20,
            "parameter": "air_temperature",
            "reason": "position_coverage",
        }
    ]


def test_smhi_duplicate_source_rows_fail_without_losing_other_parameters(server):
    url = f"{SMHI}/parameter/1/station/20/period/corrected-archive/data.csv"
    server[0][url] += "\n2020-01-01;00:00:00;11;G"
    result = history("smhi")
    assert result.attributes["partial_result"]
    assert result[("smhi:20", "air_temperature")].values[0, 0] == 99
    assert result.attributes["upstream_errors"][0]["failure_class"] == "invalid_payload"
    with pytest.raises(DatasetUpstreamError):
        history("smhi", strict_live=True)


def test_smhi_discovery_failure_does_not_claim_no_station(server):
    server[0][server[2]] = httpx2.ConnectError("offline")
    result = history("smhi")
    assert not result and result.attributes["selection"]["status"] == "discovery_failed"
    assert result.attributes["upstream_errors"][0]["failure_class"] == "connection"
    with pytest.raises(DatasetUpstreamError):
        history("smhi", strict_live=True)


def test_smhi_fetch_failure_keeps_selected_station_and_other_parameters(server):
    for period in ("corrected-archive", "latest-months"):
        server[0][f"{SMHI}/parameter/1/station/20/period/{period}/data.csv"] = (
            httpx2.ReadTimeout("slow")
        )
    result = history("smhi")
    assert len(result) == 2
    assert result.attributes["selection"]["station_id"] == 20
    assert result.attributes["upstream_error_count"] == 2
    assert not any("station/30/period" in url for url, _ in server[1])


@pytest.mark.parametrize("problem", ["unit", "length", "duplicate", "null", "error"])
def test_open_meteo_bad_payloads_and_nulls(server, problem):
    payload = server[0][OPEN]
    if problem == "unit":
        payload["hourly_units"]["wind_speed_10m"] = "km/h"
    elif problem == "length":
        payload["hourly"]["temperature_2m"].pop()
    elif problem == "duplicate":
        payload["hourly"]["time"][1] = payload["hourly"]["time"][0]
    elif problem == "null":
        payload["hourly"]["temperature_2m"][0] = None
    else:
        payload["error"] = True
    result = history()
    if problem == "null":
        assert np.isnan(result[("open_meteo:location", "air_temperature")].values[0, 0])
    else:
        assert not result and result.attributes["upstream_error_count"] == 1
        with pytest.raises(DatasetUpstreamError):
            history(strict_live=True)


def test_open_meteo_chunks_are_half_open_and_can_be_partial(server):
    original = copy.deepcopy(server[0][OPEN])

    def daily(query):
        payload = copy.deepcopy(original)
        offset = (
            int(
                (np.datetime64(query["start_date"]) - np.datetime64("2020-01-01"))
                / np.timedelta64(1, "D")
            )
            * 86400
        )
        payload["hourly"]["time"] = [
            time + offset for time in payload["hourly"]["time"]
        ]
        return payload

    server[0][OPEN] = daily
    result = history(end="2020-01-03", chunk_days=1)
    assert len(server[1]) == 2
    assert all(len(ts) == 48 for ts in result.values())
    assert result.attributes["coverage_complete"]


def test_transport_limits_and_cache(server):
    history(cache_ttl_s=60)
    history(cache_ttl_s=60)
    assert len(server[1]) == 1
    with pytest.raises(WeatherHistoryLimitError, match="max_download_bytes"):
        history(max_download_bytes=10)
    with pytest.raises(WeatherHistoryLimitError, match="max_observations"):
        history(max_observations=5)
    with pytest.raises(WeatherHistoryLimitError, match="max_requests"):
        history("smhi", max_requests=1)


def test_api_key_is_excluded_from_context_and_failures(server):
    customer = "https://customer-archive-api.open-meteo.com/v1/archive"
    result = history(api_key="test-secret", cache_ttl_s=60)
    assert "api_key" not in result.dataset_context.request.parameters
    assert "test-secret" not in str(result.series_metadata)
    history(api_key="test-secret", cache_ttl_s=60)
    assert len(server[1]) == 2
    server[0][customer] = httpx2.ConnectError("URL?apikey=test-secret")
    result = history(api_key="test-secret")
    assert "test-secret" not in str(result.attributes)
    assert result.attributes["upstream_errors"][0]["failure_class"] == "connection"


@pytest.mark.parametrize(
    "status, failure_class, attempts",
    [(404, "http_4xx", 1), (429, "http_4xx", 3), (503, "http_5xx", 3)],
)
def test_transport_status_errors_keep_dataset_contract(
    server, status, failure_class, attempts
):
    response = Response("provider-secret", status=status)
    server[0][OPEN] = response
    result = history()
    error = result.attributes["upstream_errors"][0]
    assert not result
    assert error["failure_class"] == failure_class
    assert error["status_code"] == status
    assert "provider-secret" not in str(result.attributes)
    assert len(server[1]) == attempts and response.closed
    with pytest.raises(DatasetUpstreamError) as caught:
        history(strict_live=True)
    assert caught.value.failure_class == failure_class
    assert caught.value.status_code == status
    assert "provider-secret" not in str(caught.value)


@pytest.mark.parametrize("body", [b"not JSON: provider-secret", b"\xff"])
def test_transport_decoding_errors_are_invalid_payloads(server, body):
    server[0][OPEN] = Response(body)
    result = history()
    error = result.attributes["upstream_errors"][0]
    assert error["failure_class"] == "invalid_payload"
    assert error["operation"] == "parse"
    assert "provider-secret" not in str(result.attributes)
    assert server[0][OPEN].closed


def test_transport_retries_successfully_within_request_budget(server):
    payload = server[0][OPEN]
    unavailable = Response("busy", status=503)
    server[0][OPEN] = lambda query: unavailable if len(server[1]) == 1 else payload
    result = history(max_requests=2)
    assert result.attributes["coverage_complete"]
    assert len(server[1]) == 2 and unavailable.closed
    assert result.attributes["requests"] == 1
    assert result.attributes["download_bytes"] == len(json.dumps(payload).encode())


def test_transport_request_budget_includes_retries(server):
    response = Response("busy", status=503)
    server[0][OPEN] = response
    with pytest.raises(WeatherHistoryLimitError, match="max_requests"):
        history(max_requests=1)
    assert len(server[1]) == 1 and response.closed


def test_transport_byte_budget_includes_failed_reads(server):
    class BrokenResponse(Response):
        def iter_content(self, chunk_size):
            yield b"12345678"
            raise httpx2.ReadTimeout("provider-secret")

    response = BrokenResponse(b"")
    server[0][OPEN] = response
    with pytest.raises(WeatherHistoryLimitError, match="max_download_bytes"):
        history(max_download_bytes=12)
    assert len(server[1]) == 2 and response.closed


def test_transport_byte_budget_includes_error_bodies(server):
    response = Response("provider-secret", status=404)
    server[0][OPEN] = response
    with pytest.raises(WeatherHistoryLimitError, match="max_download_bytes"):
        history(max_download_bytes=10)
    assert len(server[1]) == 1 and response.closed


@pytest.mark.parametrize(
    "headers",
    [
        {"Cache-Control": "NO-STORE"},
        {"Cache-Control": "private"},
        {"Cache-Control": "max-age=0"},
        {"Cache-Control": "max-age=1", "Age": "2"},
    ],
)
def test_transport_cache_honors_provider_restrictions(server, headers):
    server[0][OPEN] = Response(server[0][OPEN], headers=headers)
    history(cache_ttl_s=60)
    history(cache_ttl_s=60)
    assert len(server[1]) == 2


def test_transport_cache_token_and_cached_delivery_budget(server):
    url = f"{SMHI}/parameter/1/station/20/period/corrected-archive/data.csv"
    body = server[0][url]
    size = len(body.encode("utf-8"))
    args = WeatherHistoryDatasetArgs(
        source="smhi",
        location=(18.05, 59.34),
        crs="EPSG:4326",
        start="2020-01-01",
        end="2020-01-02",
        max_distance_m=50_000,
        cache_ttl_s=60,
        max_download_bytes=3 * size,
    )
    with common._Client(args) as client:
        assert client.text(url, cache_token=1) == body.lstrip("\ufeff")
        client.text(url, cache_token=1)
        assert len(server[1]) == 1
        client.text(url, cache_token=2)
        assert len(server[1]) == 2
        assert client.requests == 3 and client.bytes == 3 * size
        with pytest.raises(WeatherHistoryLimitError, match="max_download_bytes"):
            client.text(url, cache_token=2)
        assert len(server[1]) == 2


@pytest.mark.parametrize(
    "updates",
    [
        {"location": (float("nan"), 1)},
        {"location": (1,)},
        {"crs": "invalid"},
        {"start": "NaT"},
        {"end": "2019-01-01"},
        {"parameters": []},
        {"parameters": [7]},
        {"max_distance_m": 1},
        {"max_days": 1, "end": "2020-01-03"},
        {"source": "smhi"},
        {"source": "smhi", "max_distance_m": -1},
    ],
)
def test_invalid_request_is_rejected_before_network(server, updates):
    args = dict(
        source="open_meteo",
        location=(18, 59),
        crs="EPSG:4326",
        start="2020-01-01",
        end="2020-01-02",
    )
    args.update(updates)
    with pytest.raises(ValueError):
        datasets.weather_history(**args)
    assert not server[1]


def test_collection_validation_and_empty_dataframe():
    assert TimeSeriesCollection().to_dataframe().empty
    with pytest.raises(ValueError, match="Missing location"):
        TimeSeriesCollection(series={("x", "temp"): TimeSeries(name="temp")})
    with pytest.raises(TypeError):
        TimeSeriesCollection(series={("x", "temp"): []})
    args = WeatherHistoryDatasetArgs(
        source="open_meteo",
        location=(1, 2, 3),
        crs=None,
        start="2020-01-01",
        end="2020-01-02",
    )
    assert args.crs == "EPSG:3006" and args.location == (1, 2)


def test_python_only_location_dataset_has_explicit_audit_contract():
    from dtcc_core.datasets.qa import audit_dataset

    report = audit_dataset(datasets.weather_history)
    assert not report.failures()
    findings = {(item.section, item.field): item.status for item in report.findings}
    assert findings[("args_schema", "location")] == "present"
    assert findings[("descriptor", "supported_formats")] == "not_applicable"


@pytest.mark.parametrize(
    "updates",
    [
        {"start": ["2020-01-01", "2020-01-02"]},
        {"parameters": [True]},
    ],
)
def test_scalar_boundaries_and_real_parameter_ids_are_required(server, updates):
    with pytest.raises(ValueError):
        history(**updates)
    assert not server[1]


def test_partial_open_meteo_chunks_report_missing_coverage(server):
    original = server[0][OPEN]

    def fail_second(query):
        return (
            original
            if query["start_date"] == "2020-01-01"
            else httpx2.ReadTimeout("slow")
        )

    server[0][OPEN] = fail_second
    result = history(end="2020-01-03", chunk_days=1)
    assert all(len(ts) == 24 for ts in result.values())
    assert (
        result.attributes["partial_result"]
        and not result.attributes["coverage_complete"]
    )
    assert result.attributes["upstream_error_count"] == 1


def _metadata_value(value):
    if isinstance(value, Point):
        return (value.x, value.y, value.z, value.transform.srs)
    if isinstance(value, dict):
        return {
            key: _metadata_value(item)
            for key, item in value.items()
            if key != "retrieval_time"
        }
    if isinstance(value, (list, tuple)):
        return [_metadata_value(item) for item in value]
    return value


@pytest.mark.parametrize("source", ["open_meteo", "smhi"])
def test_scalar_parameter_returns_series_with_all_source_metadata(server, source):
    collection = history(source, parameters=["air_temperature"])
    result = history(source, parameters="air_temperature")
    assert isinstance(collection, TimeSeriesCollection) and len(collection) == 1
    assert isinstance(result, TimeSeries)
    key = next(iter(collection))
    original = collection[key]
    np.testing.assert_array_equal(result.times, original.times)
    np.testing.assert_array_equal(result.values, original.values)
    assert result.name == original.name and result.unit == original.unit
    assert result.role == original.role and result.association == "sample"
    assert _metadata_value(
        {name: result.attributes[name] for name in collection.attributes}
    ) == _metadata_value(collection.attributes)
    assert result.attributes["location_id"] == key[0]
    assert _metadata_value(result.attributes["location"]) == _metadata_value(
        collection.locations[key[0]]
    )
    assert result.attributes["series_metadata"] == collection.series_metadata[key]
    assert result.attributes["observations"] == collection.observations[key]
    assert result.dataset_context.identity.name == "weather_history"
    assert result.metadata.result_kind == "time_series"
    assert result.metadata.python_return_type == "dtcc_core.model.TimeSeries"
    assert result.metadata.data_types == ["time_series", "weather_history"]
    assert result.dataset_context.request.parameters["parameters"] == "air_temperature"
    assert collection.metadata.result_kind == "time_series_collection"
    assert (
        collection.metadata.python_return_type == "dtcc_core.model.TimeSeriesCollection"
    )
    assert collection.dataset_context.request.parameters["parameters"] == [
        "air_temperature"
    ]
    assert result.provenance == collection.provenance
    if source == "smhi":
        assert all(
            row["quality"] == "G" for row in result.attributes["observations"].values()
        )
        assert (
            result.attributes["selection"]["station_id"]
            == result.attributes["location"]["station_id"]
        )


@pytest.mark.parametrize("parameters", [["air_temperature"], ("temp",), ["temp", 1]])
def test_list_and_tuple_selections_keep_collection_shape(server, parameters):
    result = history(parameters=parameters)
    assert isinstance(result, TimeSeriesCollection) and len(result) == 1
    assert result.dataset_context.request.parameters["parameters"] == [
        "air_temperature"
    ]


@pytest.mark.parametrize(
    "parameter,expected",
    [
        ("temp", "air_temperature"),
        ("wind", "wind_speed"),
        ("wind_direction", "wind_direction"),
    ],
)
def test_scalar_alias_is_normalized_without_losing_request_shape(
    server, parameter, expected
):
    result = history(parameters=parameter)
    assert isinstance(result, TimeSeries) and result.name == expected
    assert result.dataset_context.request.parameters["parameters"] == expected
    assert result.attributes["requested_parameters"] == [expected]


@pytest.mark.parametrize(
    "parameters", ["", "air_temperature,wind_speed", 1, True, [], [True]]
)
def test_bad_scalar_or_list_parameters_fail_before_network(server, parameters):
    with pytest.raises(ValueError):
        history(parameters=parameters)
    assert not server[1]


@pytest.mark.parametrize("source", ["open_meteo", "smhi"])
def test_failed_scalar_fetch_returns_empty_series_with_health_metadata(server, source):
    url = OPEN if source == "open_meteo" else server[2]
    server[0][url] = Response("unavailable", status=404)
    result = history(source, parameters="air_temperature")
    assert isinstance(result, TimeSeries) and len(result) == 0
    assert result.values.shape == (0, 1) and result.unit == "celsius"
    assert result.name == "air_temperature"
    assert result.attributes["location"] is None
    assert result.attributes["observations"] == {}
    assert result.attributes["partial_result"]
    assert not result.attributes["coverage_complete"]
    assert result.attributes["upstream_error_count"] == 1
    assert result.metadata.result_kind == "time_series"
    with pytest.raises(DatasetUpstreamError):
        history(source, parameters="air_temperature", strict_live=True)


def test_scalar_no_matching_station_retains_empty_series_and_selection(server):
    result = history(
        "smhi", parameters="air_temperature", location=(18.0, 59.0), max_distance_m=1
    )
    assert isinstance(result, TimeSeries) and not len(result)
    assert result.role == "observed" and result.unit == "celsius"
    assert result.attributes["selection"]["status"] == "no_matching_station"
    assert result.attributes["location"] is None
    assert result.attributes["upstream_error_count"] == 0
    assert not result.attributes["coverage_complete"]
    with pytest.raises(NoMatchingStation):
        history(
            "smhi",
            parameters="air_temperature",
            location=(18.0, 59.0),
            max_distance_m=1,
            strict_live=True,
        )


def test_scalar_failed_archive_keeps_selected_station_metadata(server):
    for period in ("corrected-archive", "latest-months"):
        server[0][f"{SMHI}/parameter/1/station/10/period/{period}/data.csv"] = Response(
            "unavailable", status=404
        )
    result = history("smhi", parameters="air_temperature")
    assert isinstance(result, TimeSeries) and not len(result)
    assert result.attributes["location_id"] == "smhi:10"
    assert result.attributes["location"]["station_id"] == 10
    assert result.attributes["selection"]["status"] == "matched"
    assert result.attributes["upstream_error_count"] == 2


def test_partial_multiple_parameters_never_unwrap_based_on_returned_count(server):
    for pid in (3, 4):
        for period in ("corrected-archive", "latest-months"):
            server[0][f"{SMHI}/parameter/{pid}/station/20/period/{period}/data.csv"] = (
                Response("unavailable", status=404)
            )
    result = history("smhi")
    assert isinstance(result, TimeSeriesCollection) and len(result) == 1
    assert result.attributes["partial_result"]
    assert result.attributes["upstream_error_count"] == 4


def test_partial_scalar_chunks_keep_series_metadata_and_observations(server):
    original = server[0][OPEN]
    server[0][OPEN] = lambda query: (
        original
        if query["start_date"] == "2020-01-01"
        else Response("unavailable", status=404)
    )
    result = history(parameters="air_temperature", end="2020-01-03", chunk_days=1)
    assert isinstance(result, TimeSeries) and len(result) == 24
    assert (
        result.attributes["partial_result"]
        and not result.attributes["coverage_complete"]
    )
    assert result.attributes["series_metadata"]["expected_hourly_count"] == 48
    assert result.attributes["series_metadata"]["observation_count"] == 24
    assert result.attributes["upstream_error_count"] == 1


def test_scalar_metadata_excludes_api_key_and_survives_copy(server, tmp_path):
    result = history(parameters="air_temperature", api_key="test-secret")
    assert "test-secret" not in str(result.attributes)
    assert "test-secret" not in result.dataset_context.model_dump_json()
    cloned = result.copy()
    assert _metadata_value(cloned.attributes) == _metadata_value(result.attributes)
    assert cloned.dataset_context == result.dataset_context
    cloned.attributes["series_metadata"]["source"] = "edited"
    assert result.attributes["series_metadata"]["source"] == "open_meteo"
    derived = result.between()
    assert type(derived) is TimeSeries
    restored = TimeSeries()
    restored.from_proto(derived.to_proto())
    np.testing.assert_array_equal(restored.values, result.values)
    with pytest.raises(NotImplementedError, match="unsupported"):
        result.to_proto()
    target = tmp_path / "weather.dtccpkg"
    with pytest.raises(NotImplementedError, match="unsupported"):
        result.export(target, canonical=True)
    assert not target.exists()
