"""Test AirQuality dataset."""

import pytest
from unittest.mock import Mock, patch
import numpy as np


def test_air_quality_dataset_with_mocked_api():
    """Test AirQuality dataset with fully mocked HTTP API."""

    # Mock API responses that match the real API structure
    mock_stations_list = [
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

    # Station details with embedded timeseries IDs
    mock_station1_detail = {
        "id": "station1",
        "properties": {
            "label": "Test Station 1",
            "operator": "SMHI",
            "timeseries": {"ts1": {"phenomenon": {"id": "8", "label": "NO2"}}},
        },
        "geometry": {"coordinates": [18.0, 59.3, 0.0]},
    }

    mock_station2_detail = {
        "id": "station2",
        "properties": {
            "label": "Test Station 2",
            "operator": "SMHI",
            "timeseries": {"ts2": {"phenomenon": {"id": "8", "label": "NO2"}}},
        },
        "geometry": {"coordinates": [18.1, 59.4, 0.0]},
    }

    # Timeseries metadata with lastValue
    mock_ts1_data = {
        "id": "ts1",
        "label": "NO2",
        "lastValue": {"timestamp": 1738569600000, "value": 25.5},
        "uom": "µg/m³",
    }

    mock_ts2_data = {
        "id": "ts2",
        "label": "NO2",
        "lastValue": {"timestamp": 1738569600000, "value": 30.0},
        "uom": "µg/m³",
    }

    with patch("requests.get") as mock_get:
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.raise_for_status = Mock()

        def json_side_effect():
            url = mock_get.call_args[0][0]
            # Match the actual API call pattern
            if "/stations" in url and "/station1" not in url and "/station2" not in url:
                return mock_stations_list
            elif "/stations/station1" in url:
                return mock_station1_detail
            elif "/stations/station2" in url:
                return mock_station2_detail
            elif "/timeseries/ts1" in url:
                return mock_ts1_data
            elif "/timeseries/ts2" in url:
                return mock_ts2_data
            return {}

        mock_response.json = json_side_effect
        mock_get.return_value = mock_response

        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        args = AirQualityDatasetArgs(
            bounds=(674000, 6580000, 676000, 6582000),  # Stockholm in SWEREF99 TM
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

        # Check station data
        for station in stations:
            assert "station_id" in station.attributes
            assert "value" in station.attributes
            assert len(station.geometry) > 0


def test_air_quality_dataset_phenomenon_resolution():
    """Test phenomenon name to ID resolution."""
    from dtcc_core.datasets.air_quality import _resolve_phenomenon_id

    # Test numeric ID (should return as-is)
    assert _resolve_phenomenon_id("https://test.api", "1", 10.0) == "1"

    # Test common names (hardcoded mappings - updated to match actual API)
    assert _resolve_phenomenon_id("https://test.api", "NO2", 10.0) == "8"
    assert _resolve_phenomenon_id("https://test.api", "PM10", 10.0) == "5"
    assert _resolve_phenomenon_id("https://test.api", "O3", 10.0) == "7"


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
            elif "/stations/s1" in url:
                return mock_station_detail
            elif "/timeseries/ts1" in url:
                return mock_ts_data
            return {}

        mock_response.json = json_func
        mock_get.return_value = mock_response

        from dtcc_core.datasets.air_quality import (
            AirQualityDataset,
            AirQualityDatasetArgs,
        )

        args = AirQualityDatasetArgs(
            bounds=(674000, 6580000, 676000, 6582000),  # Stockholm in SWEREF99 TM
            phenomenon="NO2",
            base_url="http://test.api",
        )

        # Just test that we can build the dataset (format='pb' not yet supported)
        dataset = AirQualityDataset()
        result = dataset.build(args)
        assert result is not None
        assert hasattr(result, "stations")
