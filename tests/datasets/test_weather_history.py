"""Historical weather selection, source semantics and location reprojection."""

import copy
import json

import numpy as np
import pytest
import requests
from pyproj import Geod, Transformer

from dtcc_core import datasets
from dtcc_core.datasets import _weather_history as common
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.weather_history import (
    NoMatchingStation,
    WeatherHistoryDatasetArgs,
    WeatherHistoryLimitError,
)
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


class Response:
    def __init__(self, body):
        self.body = body
        self.headers = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def raise_for_status(self):
        if isinstance(self.body, Exception):
            raise self.body

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

    def get(self, url, params=None, **kwargs):
        calls.append((url, copy.deepcopy(params)))
        body = responses[url]
        return Response(body(params) if callable(body) else body)

    monkeypatch.setattr(requests.Session, "get", get)
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
        first["to"] = midpoint - 1
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


def test_smhi_unknown_position_does_not_silently_select_another_station(server):
    del server[0][f"{SMHI}/parameter/1/station/20.json"]["position"]
    result = history("smhi")
    assert not result and result.attributes["selection"]["status"] == "discovery_failed"


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
    server[0][server[2]] = requests.ConnectionError("offline")
    result = history("smhi")
    assert not result and result.attributes["selection"]["status"] == "discovery_failed"
    assert result.attributes["upstream_errors"][0]["failure_class"] == "connection"
    with pytest.raises(DatasetUpstreamError):
        history("smhi", strict_live=True)


def test_smhi_fetch_failure_keeps_selected_station_and_other_parameters(server):
    for period in ("corrected-archive", "latest-months"):
        server[0][f"{SMHI}/parameter/1/station/20/period/{period}/data.csv"] = (
            requests.Timeout("slow")
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
    result = history(api_key="test-secret")
    assert "api_key" not in result.dataset_context.request.parameters
    assert "test-secret" not in str(result.series_metadata)
    server[0][customer] = requests.ConnectionError("URL?apikey=test-secret")
    result = history(api_key="test-secret")
    assert "test-secret" not in str(result.attributes)
    assert result.attributes["upstream_errors"][0]["failure_class"] == "connection"


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
            else requests.Timeout("slow")
        )

    server[0][OPEN] = fail_second
    result = history(end="2020-01-03", chunk_days=1)
    assert all(len(ts) == 24 for ts in result.values())
    assert (
        result.attributes["partial_result"]
        and not result.attributes["coverage_complete"]
    )
    assert result.attributes["upstream_error_count"] == 1
