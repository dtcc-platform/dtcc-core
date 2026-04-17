"""Test AirQuality dataset."""

import pytest
import requests
from unittest.mock import Mock, patch
import numpy as np

from dtcc_core.datasets.dataset import DatasetUpstreamError


MOCK_STATIONS_LIST = [
    {
        "id": "station1",
        "properties": {"label": "Test Station 1"},
        "geometry": {"coordinates": [18.0, 59.3, 0.0]},
    },
    {
        "id": "station2",
        "properties": {"label": "Test Station 2"},
        "geometry": {"coordinates": [18.1, 59.4, 0.0]},
    },
]

MOCK_STATION1_DETAIL = {
    "id": "station1",
    "properties": {
        "label": "Test Station 1",
        "operator": "SMHI",
        "timeseries": {"ts1": {"phenomenon": {"id": "8", "label": "NO2"}}},
    },
    "geometry": {"coordinates": [18.0, 59.3, 0.0]},
}

MOCK_STATION2_DETAIL = {
    "id": "station2",
    "properties": {
        "label": "Test Station 2",
        "operator": "SMHI",
        "timeseries": {"ts2": {"phenomenon": {"id": "8", "label": "NO2"}}},
    },
    "geometry": {"coordinates": [18.1, 59.4, 0.0]},
}

MOCK_TS1_DATA = {
    "id": "ts1",
    "label": "NO2",
    "lastValue": {"timestamp": 1738569600000, "value": 25.5},
    "uom": "µg/m³",
}

MOCK_TS2_DATA = {
    "id": "ts2",
    "label": "NO2",
    "lastValue": {"timestamp": 1738569600000, "value": 30.0},
    "uom": "µg/m³",
}


def _install_mock_api(
    mock_get,
    stations_list=None,
    station1_detail=None,
    station2_detail=None,
    ts1_data=None,
    ts2_data=None,
):
    """Configure requests.get to emulate the air-quality API."""
    if stations_list is None:
        stations_list = MOCK_STATIONS_LIST
    if station1_detail is None:
        station1_detail = MOCK_STATION1_DETAIL
    if station2_detail is None:
        station2_detail = MOCK_STATION2_DETAIL
    if ts1_data is None:
        ts1_data = MOCK_TS1_DATA
    if ts2_data is None:
        ts2_data = MOCK_TS2_DATA

    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.raise_for_status = Mock()

    def json_side_effect():
        url = mock_get.call_args[0][0]
        if "/stations" in url and "/station1" not in url and "/station2" not in url:
            return stations_list
        if "/stations/station1" in url:
            return station1_detail
        if "/stations/station2" in url:
            return station2_detail
        if "/timeseries/ts1" in url:
            return ts1_data
        if "/timeseries/ts2" in url:
            return ts2_data
        return {}

    mock_response.json = json_side_effect
    mock_get.return_value = mock_response


def _json_response(data):
    """Build a mocked requests response with a fixed JSON payload."""
    response = Mock()
    response.status_code = 200
    response.raise_for_status = Mock()
    response.json = Mock(return_value=data)
    return response


def test_air_quality_dataset_with_mocked_api():
    """Test AirQuality dataset with fully mocked HTTP API."""

    with patch("requests.get") as mock_get:
        _install_mock_api(mock_get)

        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        args = AirQualityDatasetArgs(
            bounds=(17.9, 59.2, 18.2, 59.5),
            crs="EPSG:4326",
            phenomenon="NO2",
            base_url="http://test-api.example.com",
        )

        dataset = AirQualityDataset()
        result = dataset.build(args)

        # Verify result
        assert result is not None
        assert hasattr(result, "stations")
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

        # Check station data
        for station in stations:
            assert "station_id" in station.attributes
            assert "value" in station.attributes
            assert len(station.geometry) > 0


def test_air_quality_non_strict_upstream_failures_are_recorded():
    """Default mode should return partial results with explicit error metadata."""

    def side_effect(url, params=None, timeout=None):
        if "/stations" in url and "/station1" not in url and "/station2" not in url:
            return _json_response(MOCK_STATIONS_LIST)
        if "/stations/station1" in url:
            return _json_response(MOCK_STATION1_DETAIL)
        if "/stations/station2" in url:
            raise requests.ConnectionError("station detail down")
        if "/timeseries/ts1" in url:
            return _json_response(MOCK_TS1_DATA)
        raise AssertionError(f"Unexpected URL in mock: {url}")

    with patch("requests.get", side_effect=side_effect):
        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        dataset = AirQualityDataset()
        result = dataset.build(
            AirQualityDatasetArgs(
                bounds=(17.9, 59.2, 18.2, 59.5),
                crs="EPSG:4326",
                phenomenon="NO2",
                base_url="http://test-api.example.com",
            )
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


def test_air_quality_dataset_reprojects_and_filters_in_output_crs():
    """Requested output CRS should be reflected in point geometry and filtering."""

    with patch("requests.get") as mock_get:
        _install_mock_api(mock_get)

        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        dataset = AirQualityDataset()
        result = dataset.build(
            AirQualityDatasetArgs(
                bounds=(670300, 6576800, 671300, 6577800),
                crs="EPSG:3006",
                phenomenon="NO2",
                base_url="http://test-api.example.com",
            )
        )

        stations = result.stations()
        assert len(stations) == 1
        station = stations[0]
        point = station.geometry["location"]
        assert station.attributes["station_id"] == "station1"
        assert 670300 <= point.x <= 671300
        assert 6576800 <= point.y <= 6577800


def test_air_quality_strict_live_raises_typed_upstream_error():
    """strict_live should surface transient upstream failures as typed errors."""

    with patch("requests.get", side_effect=requests.ConnectionError("network down")):
        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        dataset = AirQualityDataset()
        with pytest.raises(DatasetUpstreamError) as exc_info:
            dataset.build(
                AirQualityDatasetArgs(
                    bounds=(17.9, 59.2, 18.2, 59.5),
                    crs="EPSG:4326",
                    phenomenon="NO2",
                    base_url="http://test-api.example.com",
                    strict_live=True,
                )
            )

        assert exc_info.value.failure_class == "connection"


def test_air_quality_empty_bbox_is_not_marked_partial():
    """An empty station search should stay empty without upstream-error metadata."""

    with patch("requests.get") as mock_get:
        _install_mock_api(mock_get, stations_list=[])

        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        dataset = AirQualityDataset()
        result = dataset.build(
            AirQualityDatasetArgs(
                bounds=(17.9, 59.2, 18.2, 59.5),
                crs="EPSG:4326",
                phenomenon="NO2",
                base_url="http://test-api.example.com",
            )
        )

        assert len(result.stations()) == 0
        assert result.attributes["total_stations_found"] == 0
        assert result.attributes["partial_result"] is False
        assert result.attributes["upstream_error_count"] == 0
        assert result.attributes["upstream_errors"] == []


def test_air_quality_dataset_phenomenon_resolution():
    """Test phenomenon name to ID resolution."""
    from dtcc_core.datasets.air_quality import _resolve_phenomenon_id

    # Test numeric ID (should return as-is)
    assert _resolve_phenomenon_id("https://test.api", "1", 10.0) == "1"

    # Test common names (hardcoded mappings - updated to match actual API)
    assert _resolve_phenomenon_id("https://test.api", "NO2", 10.0) == "8"
    assert _resolve_phenomenon_id("https://test.api", "PM10", 10.0) == "5"
    assert _resolve_phenomenon_id("https://test.api", "O3", 10.0) == "7"


def test_air_quality_unknown_phenomenon_raises_value_error():
    """Unknown phenomenon names should fail fast when the lookup succeeds."""
    from dtcc_core.datasets.air_quality import _resolve_phenomenon_id

    with patch("requests.get", return_value=_json_response([{"id": "8", "label": "NO2"}])):
        with pytest.raises(ValueError, match="Unknown air quality phenomenon"):
            _resolve_phenomenon_id("https://test.api", "NOT_A_REAL_PHENOMENON", 10.0)


def test_air_quality_dataset_format_pb():
    """Test that format='pb' is handled (even if it doesn't fully work yet)."""

    mock_stations = [
        {
            "id": "s1",
            "properties": {"label": "Station"},
            "geometry": {"coordinates": [18.0, 59.0, 0.0]},
        }
    ]

    mock_station_detail = {
        "id": "s1",
        "properties": {
            "label": "Station",
            "timeseries": {"ts1": {"phenomenon": {"id": "8", "label": "NO2"}}},
        },
        "geometry": {"coordinates": [18.0, 59.0, 0.0]},
    }

    mock_ts_data = {
        "id": "ts1",
        "lastValue": {"timestamp": 1738569600000, "value": 30.0},
        "uom": "µg/m³",
    }

    with patch("requests.get") as mock_get:
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.raise_for_status = Mock()

        def json_func():
            url = mock_get.call_args[0][0]
            if "/stations" in url and "/s1" not in url:
                return mock_stations
            if "/stations/s1" in url:
                return mock_station_detail
            if "/timeseries/ts1" in url:
                return mock_ts_data
            return {}

        mock_response.json = json_func
        mock_get.return_value = mock_response

        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        args = AirQualityDatasetArgs(
            bounds=(17.9, 58.9, 18.1, 59.1),
            crs="EPSG:4326",
            phenomenon="NO2",
            base_url="http://test.api",
        )

        # Just test that we can build the dataset (format='pb' not yet supported)
        dataset = AirQualityDataset()
        result = dataset.build(args)
        assert result is not None
        assert hasattr(result, "stations")
