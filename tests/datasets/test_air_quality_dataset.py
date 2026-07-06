"""Tests for the air quality dataset."""

from __future__ import annotations

import copy
from datetime import datetime
import json
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest
import requests

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.air_quality import (
    AirQualityDataset,
    AirQualityDatasetArgs,
    _extract_latest_value,
    _get_json,
    _resolve_phenomenon_id,
)
from dtcc_core.datasets.dataset import DatasetUpstreamError


BASE_URL = "http://test-api.example.com"
DEFAULT_BOUNDS = (17.9, 59.2, 18.2, 59.5)
DEFAULT_CRSS = "EPSG:4326"

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "smhi" / "air_quality"


def _fixture_json(name: str):
    path = FIXTURE_DIR / name
    return json.loads(path.read_text(encoding="utf-8"))


MOCK_STATIONS_LIST = _fixture_json("stations.json")
MOCK_STATION1_DETAIL = _fixture_json("station1.json")
MOCK_STATION2_DETAIL = _fixture_json("station2.json")
MOCK_TS1_DATA = _fixture_json("timeseries_ts1.json")
MOCK_TS2_DATA = _fixture_json("timeseries_ts2.json")
MOCK_PHENOMENA = _fixture_json("phenomena.json")


def _upstream_error(
    target: str,
    *,
    operation: str = "fetch_json",
    failure_class: str = "connection",
    status_code: int | None = None,
    message: str | None = None,
) -> DatasetUpstreamError:
    return DatasetUpstreamError(
        dataset="air_quality",
        operation=operation,
        target=target,
        failure_class=failure_class,
        status_code=status_code,
        message=message or f"air_quality {operation} failed for {target}",
    )


def _default_payloads() -> dict[str, object]:
    return {
        f"{BASE_URL}/stations": copy.deepcopy(MOCK_STATIONS_LIST),
        f"{BASE_URL}/stations/station1": copy.deepcopy(MOCK_STATION1_DETAIL),
        f"{BASE_URL}/stations/station2": copy.deepcopy(MOCK_STATION2_DETAIL),
        f"{BASE_URL}/timeseries/ts1": copy.deepcopy(MOCK_TS1_DATA),
        f"{BASE_URL}/timeseries/ts2": copy.deepcopy(MOCK_TS2_DATA),
        f"{BASE_URL}/timeseries/ts1/getData": {"values": [], "uom": "µg/m³"},
        f"{BASE_URL}/timeseries/ts2/getData": {"values": [], "uom": "µg/m³"},
        f"{BASE_URL}/phenomena": copy.deepcopy(MOCK_PHENOMENA),
    }


def _make_get_json_side_effect(
    *,
    payloads: dict[str, object] | None = None,
    errors: dict[str, Exception] | None = None,
    calls: list[tuple[str, dict | None, float]] | None = None,
):
    api_payloads = _default_payloads()
    if payloads:
        api_payloads.update(payloads)
    api_errors = errors or {}

    def side_effect(url: str, params=None, timeout_s: float = 10.0):
        if calls is not None:
            copied_params = dict(params) if isinstance(params, dict) else params
            calls.append((url, copied_params, timeout_s))

        if url in api_errors:
            raise api_errors[url]

        if url not in api_payloads:
            raise AssertionError(f"Unexpected URL in mock: {url}")

        return copy.deepcopy(api_payloads[url])

    return side_effect


def _make_args(**overrides) -> AirQualityDatasetArgs:
    args = {
        "bounds": DEFAULT_BOUNDS,
        "crs": DEFAULT_CRSS,
        "phenomenon": "NO2",
        "base_url": BASE_URL,
    }
    args.update(overrides)
    return AirQualityDatasetArgs(**args)


def _build_dataset(
    *,
    payloads: dict[str, object] | None = None,
    errors: dict[str, Exception] | None = None,
    calls: list[tuple[str, dict | None, float]] | None = None,
    **arg_overrides,
):
    dataset = AirQualityDataset()
    args = _make_args(**arg_overrides)
    with patch(
        "dtcc_core.datasets.air_quality._get_json",
        side_effect=_make_get_json_side_effect(
            payloads=payloads, errors=errors, calls=calls
        ),
    ):
        return dataset.build(args)


def test_air_quality_dataset_with_mocked_api():
    """The happy path should return a populated SensorCollection."""
    result = _build_dataset()

    stations = result.stations()
    assert len(stations) == 2
    assert result.attributes["phenomenon"] == "NO2"
    assert result.attributes["phenomenon_id"] == "8"
    assert result.attributes["total_stations_found"] == 2
    assert result.attributes["stations_used"] == 2
    assert result.attributes["partial_result"] is False
    assert result.attributes["upstream_error_count"] == 0
    assert result.attributes["upstream_errors"] == []
    assert result.attributes["stations_skipped_upstream"] == 0


def test_air_quality_collection_attributes():
    """The collection metadata should expose the result contract consistently."""
    result = _build_dataset()
    attrs = result.attributes

    assert attrs["dtcc_type"] == "sensor_collection"
    assert attrs["source"] == "datavardluft.smhi.se"
    assert attrs["phenomenon"] == "NO2"
    assert attrs["phenomenon_id"] == "8"
    assert attrs["crs"] == "EPSG:4326"
    assert attrs["bounds"] == str(DEFAULT_BOUNDS)
    assert attrs["total_stations_found"] == 2
    assert attrs["stations_used"] == 2
    assert attrs["stations_skipped_upstream"] == 0
    assert attrs["stations_skipped_no_coords"] == 0
    assert attrs["stations_skipped_no_timeseries"] == 0
    assert attrs["stations_skipped_no_value"] == 0
    assert attrs["partial_result"] is False
    assert attrs["upstream_error_count"] == 0
    assert attrs["upstream_errors"] == []
    assert isinstance(attrs["retrieval_time"], str)
    assert attrs["retrieval_time"]


def test_air_quality_station_attributes_and_string_summary():
    """Stations should expose the expected metadata and string rendering should work."""
    result = _build_dataset()

    station = result.stations()[0]
    point = station.geometry["location"]

    assert station.attributes["station_id"] == "station1"
    assert station.attributes["station_name"] == "Test Station 1"
    assert station.attributes["operator"] == "SMHI"
    assert station.attributes["phenomenon"] == "NO2"
    assert station.attributes["phenomenon_id"] == "8"
    assert station.attributes["unit"] == "µg/m³"
    assert station.attributes["value"] == pytest.approx(25.5)
    assert station.attributes["timestamp"]
    assert point.fields[0].name == "NO2"
    assert point.fields[0].unit == "µg/m³"
    assert point.fields[0].values[0] == pytest.approx(25.5)

    text = str(result)
    assert "SensorCollection" in text
    assert "station" in text.lower()


def test_air_quality_dataset_reprojects_and_filters_in_output_crs():
    """Requested output CRS should be reflected in point geometry and filtering."""
    result = _build_dataset(
        bounds=(670300, 6576800, 671300, 6577800),
        crs="EPSG:3006",
    )

    stations = result.stations()
    assert len(stations) == 1

    station = stations[0]
    point = station.geometry["location"]
    assert station.attributes["station_id"] == "station1"
    assert 670300 <= point.x <= 671300
    assert 6576800 <= point.y <= 6577800


def test_air_quality_dataset_format_pb_returns_bytes():
    """format='pb' should exercise the real protobuf export path."""
    result = _build_dataset(format="pb")
    assert isinstance(result, bytes)
    assert len(result) > 0


def test_air_quality_empty_bbox_is_not_marked_partial():
    """A genuinely empty station search should remain a clean empty result."""
    result = _build_dataset(payloads={f"{BASE_URL}/stations": []})

    assert len(result.stations()) == 0
    assert result.attributes["total_stations_found"] == 0
    assert result.attributes["stations_used"] == 0
    assert result.attributes["partial_result"] is False
    assert result.attributes["upstream_error_count"] == 0
    assert result.attributes["upstream_errors"] == []


def test_air_quality_max_stations_caps_processing():
    """The max_stations safety limit should stop processing after the cap."""
    calls: list[tuple[str, dict | None, float]] = []
    result = _build_dataset(max_stations=1, calls=calls)

    assert len(result.stations()) == 1
    assert result.attributes["total_stations_found"] == 2
    assert result.attributes["stations_used"] == 1
    assert f"{BASE_URL}/stations/station2" not in {url for url, _params, _timeout in calls}


def test_air_quality_missing_coordinates_increment_counter():
    """Stations without coordinates should be skipped before detail fetches."""
    stations_list = copy.deepcopy(MOCK_STATIONS_LIST)
    stations_list[1]["geometry"] = {"coordinates": []}
    calls: list[tuple[str, dict | None, float]] = []

    result = _build_dataset(
        payloads={f"{BASE_URL}/stations": stations_list},
        calls=calls,
    )

    assert len(result.stations()) == 1
    assert result.attributes["stations_skipped_no_coords"] == 1
    assert result.attributes["stations_used"] == 1
    assert f"{BASE_URL}/stations/station2" not in {url for url, _params, _timeout in calls}


def test_air_quality_missing_matching_timeseries_increments_counter():
    """Stations with no timeseries for the requested phenomenon should be skipped cleanly."""
    station1_detail = copy.deepcopy(MOCK_STATION1_DETAIL)
    station1_detail["properties"]["timeseries"] = {
        "ts1": {"phenomenon": {"id": "999", "label": "OTHER"}}
    }

    result = _build_dataset(
        payloads={
            f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
            f"{BASE_URL}/stations/station1": station1_detail,
        }
    )

    assert len(result.stations()) == 0
    assert result.attributes["stations_skipped_no_timeseries"] == 1
    assert result.attributes["stations_skipped_upstream"] == 0
    assert result.attributes["partial_result"] is False


def test_air_quality_drop_missing_true_skips_station_without_value():
    """Stations with no current value should be dropped when drop_missing=True."""
    ts1_without_value = {"id": "ts1", "label": "NO2", "uom": "µg/m³"}
    calls: list[tuple[str, dict | None, float]] = []

    result = _build_dataset(
        payloads={
            f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
            f"{BASE_URL}/timeseries/ts1": ts1_without_value,
            f"{BASE_URL}/timeseries/ts1/getData": {"values": [], "uom": "µg/m³"},
        },
        calls=calls,
    )

    assert len(result.stations()) == 0
    assert result.attributes["stations_skipped_no_value"] == 1
    assert result.attributes["partial_result"] is False
    assert any(url.endswith("/timeseries/ts1/getData") for url, _params, _timeout in calls)


def test_air_quality_drop_missing_false_keeps_nan_station_without_value():
    """When drop_missing=False, stations with no value should be kept with NaN."""
    ts1_without_value = {"id": "ts1", "label": "NO2", "uom": "µg/m³"}

    result = _build_dataset(
        payloads={
            f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
            f"{BASE_URL}/timeseries/ts1": ts1_without_value,
            f"{BASE_URL}/timeseries/ts1/getData": {"values": [], "uom": "µg/m³"},
        },
        drop_missing=False,
    )

    assert len(result.stations()) == 1
    station = result.stations()[0]
    point = station.geometry["location"]
    assert np.isnan(station.attributes["value"])
    assert np.isnan(point.fields[0].values[0])
    assert result.attributes["stations_skipped_no_value"] == 0
    assert result.attributes["partial_result"] is False


def test_air_quality_non_strict_station_detail_failures_are_recorded():
    """Station-detail failures should become partial results in default mode."""
    result = _build_dataset(
        errors={
            f"{BASE_URL}/stations/station2": _upstream_error(
                f"{BASE_URL}/stations/station2"
            )
        }
    )

    stations = result.stations()
    assert len(stations) == 1
    assert stations[0].attributes["station_id"] == "station1"
    assert result.attributes["partial_result"] is True
    assert result.attributes["upstream_error_count"] == 1
    assert result.attributes["stations_skipped_upstream"] == 1
    assert result.attributes["stations_used"] == 1
    assert result.attributes["total_stations_found"] == 2
    assert result.attributes["stations_skipped_no_timeseries"] == 0
    assert result.attributes["upstream_errors"][0]["dataset"] == "air_quality"
    assert result.attributes["upstream_errors"][0]["failure_class"] == "connection"


def test_air_quality_non_strict_timeseries_failures_are_recorded():
    """Timeseries-fetch failures should also become partial results in default mode."""
    result = _build_dataset(
        errors={
            f"{BASE_URL}/timeseries/ts2": _upstream_error(
                f"{BASE_URL}/timeseries/ts2"
            )
        }
    )

    assert len(result.stations()) == 1
    assert result.attributes["stations_skipped_upstream"] == 1
    assert result.attributes["upstream_error_count"] == 1
    assert result.attributes["partial_result"] is True


def test_air_quality_non_strict_station_list_failure_returns_partial_empty_result():
    """A station-list failure should degrade to an empty partial result."""
    result = _build_dataset(
        errors={f"{BASE_URL}/stations": _upstream_error(f"{BASE_URL}/stations")}
    )

    assert len(result.stations()) == 0
    assert result.attributes["total_stations_found"] == 0
    assert result.attributes["stations_used"] == 0
    assert result.attributes["partial_result"] is True
    assert result.attributes["upstream_error_count"] == 1
    assert result.attributes["stations_skipped_upstream"] == 0


def test_air_quality_strict_live_raises_on_station_list_failure():
    """strict_live should surface station-list failures immediately."""
    with pytest.raises(DatasetUpstreamError) as exc_info:
        _build_dataset(
            strict_live=True,
            errors={f"{BASE_URL}/stations": _upstream_error(f"{BASE_URL}/stations")},
        )

    assert exc_info.value.failure_class == "connection"


def test_air_quality_strict_live_raises_on_station_detail_failure():
    """strict_live should surface station-detail failures immediately."""
    with pytest.raises(DatasetUpstreamError) as exc_info:
        _build_dataset(
            strict_live=True,
            errors={
                f"{BASE_URL}/stations/station1": _upstream_error(
                    f"{BASE_URL}/stations/station1"
                )
            },
        )

    assert exc_info.value.failure_class == "connection"


def test_air_quality_strict_live_raises_on_timeseries_failure():
    """strict_live should surface timeseries fetch failures immediately."""
    with pytest.raises(DatasetUpstreamError) as exc_info:
        _build_dataset(
            strict_live=True,
            errors={
                f"{BASE_URL}/timeseries/ts1": _upstream_error(
                    f"{BASE_URL}/timeseries/ts1"
                )
            },
        )

    assert exc_info.value.failure_class == "connection"


def test_extract_latest_value_valid_payload():
    """_extract_latest_value should parse the standard lastValue payload."""
    value, timestamp, unit = _extract_latest_value(copy.deepcopy(MOCK_TS1_DATA))
    expected_timestamp = datetime.fromtimestamp(1738569600000 / 1000.0).isoformat()

    assert value == pytest.approx(25.5)
    assert timestamp == expected_timestamp
    assert unit == "µg/m³"


def test_extract_latest_value_missing_last_value_returns_none():
    """Missing lastValue should produce no extracted value."""
    assert _extract_latest_value({"id": "ts1", "uom": "µg/m³"}) is None


def test_extract_latest_value_malformed_last_value_returns_none():
    """Malformed lastValue payloads should be ignored safely."""
    assert (
        _extract_latest_value(
            {"id": "ts1", "lastValue": {"timestamp": "bad", "value": "not-a-number"}}
        )
        is None
    )


def test_air_quality_fallback_getdata_success():
    """If lastValue is missing, the dataset should use /getData as a fallback."""
    ts1_without_value = {"id": "ts1", "label": "NO2", "uom": "µg/m³"}
    calls: list[tuple[str, dict | None, float]] = []

    result = _build_dataset(
        payloads={
            f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
            f"{BASE_URL}/timeseries/ts1": ts1_without_value,
            f"{BASE_URL}/timeseries/ts1/getData": {
                "values": [{"timestamp": "2025-02-03T00:00:00Z", "value": "31.0"}],
                "uom": "µg/m³",
            },
        },
        calls=calls,
    )

    station = result.stations()[0]
    point = station.geometry["location"]
    assert station.attributes["value"] == pytest.approx(31.0)
    assert station.attributes["timestamp"] == "2025-02-03T00:00:00Z"
    assert station.attributes["unit"] == "µg/m³"
    assert point.fields[0].values[0] == pytest.approx(31.0)
    assert (
        f"{BASE_URL}/timeseries/ts1/getData",
        {"timespan": "P7D"},
        10.0,
    ) in calls


def test_air_quality_fallback_getdata_upstream_failure_marks_partial_in_non_strict_mode():
    """Fallback failures should mark the result partial and skip the station."""
    ts1_without_value = {"id": "ts1", "label": "NO2", "uom": "µg/m³"}

    result = _build_dataset(
        payloads={
            f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
            f"{BASE_URL}/timeseries/ts1": ts1_without_value,
        },
        errors={
            f"{BASE_URL}/timeseries/ts1/getData": _upstream_error(
                f"{BASE_URL}/timeseries/ts1/getData"
            )
        },
    )

    assert len(result.stations()) == 0
    assert result.attributes["stations_skipped_upstream"] == 1
    assert result.attributes["stations_skipped_no_value"] == 0
    assert result.attributes["partial_result"] is True
    assert result.attributes["upstream_error_count"] == 1


def test_air_quality_fallback_getdata_upstream_failure_raises_in_strict_mode():
    """strict_live should raise when the fallback /getData request fails."""
    ts1_without_value = {"id": "ts1", "label": "NO2", "uom": "µg/m³"}

    with pytest.raises(DatasetUpstreamError) as exc_info:
        _build_dataset(
            payloads={
                f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
                f"{BASE_URL}/timeseries/ts1": ts1_without_value,
            },
            errors={
                f"{BASE_URL}/timeseries/ts1/getData": _upstream_error(
                    f"{BASE_URL}/timeseries/ts1/getData"
                )
            },
            strict_live=True,
        )

    assert exc_info.value.failure_class == "connection"


def test_air_quality_dataset_phenomenon_resolution_common_names():
    """Common hardcoded phenomenon names should map to the expected IDs."""
    assert _resolve_phenomenon_id("https://test.api", "1", 10.0) == "1"
    assert _resolve_phenomenon_id("https://test.api", "NO2", 10.0) == "8"
    assert _resolve_phenomenon_id("https://test.api", "PM10", 10.0) == "5"
    assert _resolve_phenomenon_id("https://test.api", "O3", 10.0) == "7"


def test_air_quality_api_lookup_resolves_unknown_but_valid_phenomenon():
    """The API lookup path should resolve valid phenomena not in the hardcoded map."""
    with patch(
        "dtcc_core.datasets.air_quality._get_json",
        return_value=[{"id": "42", "label": "CUSTOM"}],
    ):
        assert _resolve_phenomenon_id(BASE_URL, "CUSTOM", 10.0) == "42"


def test_air_quality_non_strict_lookup_failure_falls_back_to_raw_phenomenon_and_records_error():
    """Lookup failures in non-strict mode should record the error and continue."""
    station1_detail = copy.deepcopy(MOCK_STATION1_DETAIL)
    station1_detail["properties"]["timeseries"] = {
        "ts1": {"phenomenon": {"id": "EXOTIC", "label": "EXOTIC"}}
    }

    result = _build_dataset(
        phenomenon="EXOTIC",
        payloads={
            f"{BASE_URL}/stations": [copy.deepcopy(MOCK_STATIONS_LIST[0])],
            f"{BASE_URL}/stations/station1": station1_detail,
        },
        errors={f"{BASE_URL}/phenomena": _upstream_error(f"{BASE_URL}/phenomena")},
    )

    assert len(result.stations()) == 1
    assert result.attributes["phenomenon_id"] == "EXOTIC"
    assert result.attributes["partial_result"] is True
    assert result.attributes["upstream_error_count"] == 1


def test_air_quality_strict_lookup_failure_raises_upstream_error():
    """strict_live should raise when phenomenon lookup fails upstream."""
    with pytest.raises(DatasetUpstreamError) as exc_info:
        _build_dataset(
            phenomenon="EXOTIC",
            strict_live=True,
            errors={f"{BASE_URL}/phenomena": _upstream_error(f"{BASE_URL}/phenomena")},
        )

    assert exc_info.value.target == f"{BASE_URL}/phenomena"


def test_air_quality_unknown_phenomenon_raises_value_error():
    """Unknown phenomena should fail fast when the lookup endpoint succeeds."""
    with patch(
        "dtcc_core.datasets.air_quality._get_json",
        return_value=[{"id": "8", "label": "NO2"}],
    ):
        with pytest.raises(ValueError, match="Unknown air quality phenomenon"):
            _resolve_phenomenon_id(BASE_URL, "NOT_A_REAL_PHENOMENON", 10.0)


def test_get_json_invalid_json_raises_typed_error():
    """_get_json should wrap JSON decoding failures in a typed upstream error."""
    response = Mock()
    response.raise_for_status = Mock()
    response.json.side_effect = json.JSONDecodeError("bad json", "{}", 0)

    with patch("requests.get", return_value=response):
        with pytest.raises(DatasetUpstreamError) as exc_info:
            _get_json(f"{BASE_URL}/stations", timeout_s=3.0)

    assert exc_info.value.failure_class == "invalid_json"
    assert exc_info.value.operation == "decode_json"
    assert exc_info.value.target == f"{BASE_URL}/stations"


def test_air_quality_registered_name():
    """The dataset class should expose the correct registration name."""
    assert AirQualityDataset().name == "air_quality"


def test_air_quality_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("air_quality")
    assert ds is not None
    assert ds.name == "air_quality"


def test_air_quality_module_attribute():
    """The dataset should be exposed as a callable module attribute."""
    assert hasattr(datasets, "air_quality")
    assert callable(datasets.air_quality)
