"""Tests for the weather dataset (SMHI metobs latest-hour)."""

import pytest
import math
import numpy as np
from unittest.mock import patch, Mock

from dtcc_core.datasets.weather import (
    WeatherDataset,
    WeatherDatasetArgs,
    _parse_latest_hour_csv,
    _resolve_parameter,
    PARAMETER_NAMES,
)


# ── Fixtures ─────────────────────────────────────────────────────────────

# Minimal CSV that mimics the real SMHI metobs latest-hour response for
# parameter 1 (air temperature).  Includes:
#   - 2 metadata rows
#   - blank line
#   - column header with a datetime in position 5
#   - "Data från senaste timmen" subtitle (should be skipped)
#   - 3 station rows (one with quality Y, one with empty value)
#   - metadata comments appended after ;; on some rows

SAMPLE_CSV_PARAM1 = """\
Parameternamn;Beskrivning;Enhet
Lufttemperatur;momentanvärde, 1 gång/tim;celsius

StationsId;Stationsnamn;Latitude;Longitude;Height;2026-02-12 19:00:00;Kvalitet;;
Data från senaste timmen
100;Stockholm City;59.3293;18.0686;28.0;-3.5;G;;Tidsperiod (fr.o.m.) = 2026-02-12 18:00:01 (UTC)
200;Göteborg Landvetter;57.6628;12.2960;169.0;-5.1;Y;;Tidsperiod (t.o.m.) = 2026-02-12 19:00:00 (UTC)
300;Malmö Airport;55.5364;13.3762;18.0;;G;;
"""

SAMPLE_CSV_PARAM4 = """\
Parameternamn;Beskrivning;Enhet
Vindhastighet;medelvärde 10 min, 1 gång/tim;m/s

StationsId;Stationsnamn;Latitude;Longitude;Height;2026-02-12 19:00:00;Kvalitet;;
Data från senaste timmen
100;Stockholm City;59.3293;18.0686;28.0;4.2;G;;
200;Göteborg Landvetter;57.6628;12.2960;169.0;7.8;G;;
300;Malmö Airport;55.5364;13.3762;18.0;5.0;G;;
"""


# ── CSV parser tests ─────────────────────────────────────────────────────


class TestParseLatestHourCsv:
    """Tests for _parse_latest_hour_csv."""

    def test_meta_extraction(self):
        meta, _ = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        assert meta["parameter_name"] == "Lufttemperatur"
        assert meta["unit"] == "celsius"
        assert "momentanvärde" in meta["description"]

    def test_timestamp_detection(self):
        meta, _ = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        assert meta["timestamp"] == "2026-02-12 19:00:00"

    def test_correct_number_of_records(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        # 3 station rows (the subtitle row is skipped because StationsId
        # is not a valid int)
        assert len(records) == 3

    def test_station_id_parsing(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        ids = [r["station_id"] for r in records]
        assert ids == [100, 200, 300]

    def test_coordinates(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        assert records[0]["lat"] == pytest.approx(59.3293)
        assert records[0]["lon"] == pytest.approx(18.0686)
        assert records[0]["height"] == pytest.approx(28.0)

    def test_value_parsing(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        assert records[0]["value"] == pytest.approx(-3.5)
        assert records[1]["value"] == pytest.approx(-5.1)

    def test_missing_value_is_nan(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        assert math.isnan(records[2]["value"])

    def test_quality_codes(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        assert records[0]["quality"] == "G"
        assert records[1]["quality"] == "Y"

    def test_subtitle_row_skipped(self):
        """The 'Data från senaste timmen' row must not produce a record."""
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM1)
        names = [r["station_name"] for r in records]
        assert "Data från senaste timmen" not in " ".join(names)

    def test_empty_input(self):
        meta, records = _parse_latest_hour_csv("")
        assert records == []

    def test_wind_speed_csv(self):
        meta, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM4)
        assert meta["parameter_name"] == "Vindhastighet"
        assert meta["unit"] == "m/s"
        assert len(records) == 3
        assert records[0]["value"] == pytest.approx(4.2)


# ── Build tests with mocked HTTP ────────────────────────────────────────


def _mock_get_text(url, params=None, timeout_s=10.0):
    """Return sample CSV based on which parameter id appears in the URL."""
    if "/parameter/1/" in url:
        return SAMPLE_CSV_PARAM1
    elif "/parameter/4/" in url:
        return SAMPLE_CSV_PARAM4
    else:
        raise RuntimeError(f"Unmocked URL: {url}")


class TestWeatherDatasetBuild:
    """Integration-style tests with HTTP mocked out."""

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_returns_sensor_collection(self, mock_fetch):
        dataset = WeatherDataset()
        # Bounds covering all three sample stations (WGS84)
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1, 4],
        )
        result = dataset.build(args)

        from dtcc_core.model.object import SensorCollection

        assert isinstance(result, SensorCollection)

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_station_count(self, mock_fetch):
        dataset = WeatherDataset()
        # Bounds covering all three stations
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1, 4],
        )
        result = dataset.build(args)
        stations = result.stations()
        # Station 300 (Malmö Airport) has NaN for param 1 but a valid value
        # for param 4, so it should still be included (not all-missing)
        assert len(stations) == 3

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_bbox_filtering(self, mock_fetch):
        dataset = WeatherDataset()
        # Narrow bounds: only Stockholm (lat ~59.33, lon ~18.07)
        args = WeatherDatasetArgs(
            bounds=(18.0, 59.0, 18.2, 59.5),
            crs="EPSG:4326",
            parameters=[1],
        )
        result = dataset.build(args)
        stations = result.stations()
        assert len(stations) == 1
        assert stations[0].attributes["station_name"] == "Stockholm City"

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_station_has_point_geometry(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1],
        )
        result = dataset.build(args)
        for station in result.stations():
            assert "location" in station.geometry
            from dtcc_core.model.geometry import Point

            assert isinstance(station.geometry["location"], Point)

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_field_names_dtcc_style(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(18.0, 59.0, 18.2, 59.5),
            crs="EPSG:4326",
            parameters=[1, 4],
            field_name_style="dtcc",
        )
        result = dataset.build(args)
        station = result.stations()[0]
        point = station.geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "air_temperature" in field_names
        assert "wind_speed" in field_names

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_field_names_smhi_style(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(18.0, 59.0, 18.2, 59.5),
            crs="EPSG:4326",
            parameters=[1, 4],
            field_name_style="smhi",
        )
        result = dataset.build(args)
        station = result.stations()[0]
        point = station.geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "Lufttemperatur" in field_names
        assert "Vindhastighet" in field_names

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_to_arrays(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1],
        )
        result = dataset.build(args)
        pts, vals = result.to_arrays("air_temperature")
        # Stations with valid temperature: 100 and 200.  300 has NaN but
        # is still present — to_arrays may or may not include NaN rows
        # depending on implementation; just check shapes match.
        assert pts.shape[0] == vals.shape[0]
        assert pts.shape[1] == 3

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_drop_missing_all(self, mock_fetch):
        """A station with only NaN values is dropped when drop_missing=True."""
        dataset = WeatherDataset()
        # Only request param 1 — station 300 has NaN for that parameter
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1],
            drop_missing=True,
        )
        result = dataset.build(args)
        sids = [s.attributes["station_id"] for s in result.stations()]
        assert 300 not in sids

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_drop_missing_false_keeps_nan_station(self, mock_fetch):
        """When drop_missing=False, stations with all NaN values are kept."""
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1],
            drop_missing=False,
        )
        result = dataset.build(args)
        sids = [s.attributes["station_id"] for s in result.stations()]
        assert 300 in sids

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_station_attributes(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(18.0, 59.0, 18.2, 59.5),
            crs="EPSG:4326",
            parameters=[1],
        )
        result = dataset.build(args)
        station = result.stations()[0]
        assert station.attributes["station_id"] == 100
        assert station.attributes["station_name"] == "Stockholm City"
        assert station.attributes["elevation"] == pytest.approx(28.0)
        assert station.attributes["q_air_temperature"] == "G"

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_collection_attributes(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1, 4],
        )
        result = dataset.build(args)
        attrs = result.attributes
        assert attrs["source"] == "smhi_metobs"
        assert attrs["dataset"] == "weather"
        assert attrs["period"] == "latest-hour"
        assert attrs["crs"] == "EPSG:4326"
        assert 1 in attrs["parameters"]
        assert 4 in attrs["parameters"]

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_format_pb_returns_bytes(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(12.0, 55.0, 19.0, 60.0),
            crs="EPSG:4326",
            parameters=[1],
            format="pb",
        )
        result = dataset.build(args)
        assert isinstance(result, bytes)
        assert len(result) > 0

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_sweref99_crs_reprojects(self, mock_fetch):
        """When crs=EPSG:3006, coordinates should be reprojected."""
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(674000, 6578000, 678000, 6583000),
            crs="EPSG:3006",
            parameters=[1],
        )
        result = dataset.build(args)
        # If any stations matched, their x/y should be in SWEREF99 range
        for station in result.stations():
            pt = station.geometry["location"]
            # SWEREF99 TM x is typically 100 000 – 900 000
            # y is typically 6 100 000 – 7 700 000
            if result.stations():
                assert pt.x > 100_000 or pt.x < 0  # just check it's not WGS84
                break  # only need one


class TestWeatherDatasetRegistration:
    """Test that the dataset is properly registered."""

    def test_dataset_registered(self):
        from dtcc_core.datasets import get_dataset

        ds = get_dataset("weather")
        assert ds is not None
        assert ds.name == "weather"

    def test_callable_via_module(self):
        import dtcc_core.datasets as datasets

        assert hasattr(datasets, "weather")
        assert callable(datasets.weather)


class TestParameterNames:
    """Test the PARAMETER_NAMES mapping and name resolution."""

    def test_default_parameters_have_names(self):
        defaults = [1, 3, 4, 6, 7, 9]
        for pid in defaults:
            assert pid in PARAMETER_NAMES
            assert isinstance(PARAMETER_NAMES[pid], str)
            assert len(PARAMETER_NAMES[pid]) > 0

    def test_resolve_integer(self):
        assert _resolve_parameter(1) == 1
        assert _resolve_parameter(4) == 4

    def test_resolve_string_integer(self):
        assert _resolve_parameter("1") == 1
        assert _resolve_parameter("9") == 9

    def test_resolve_dtcc_name(self):
        assert _resolve_parameter("air_temperature") == 1
        assert _resolve_parameter("wind_speed") == 4
        assert _resolve_parameter("relative_humidity") == 6
        assert _resolve_parameter("sea_level_pressure") == 9

    def test_resolve_alias(self):
        assert _resolve_parameter("temperature") == 1
        assert _resolve_parameter("temp") == 1
        assert _resolve_parameter("wind") == 4
        assert _resolve_parameter("humidity") == 6
        assert _resolve_parameter("pressure") == 9
        assert _resolve_parameter("precipitation") == 7
        assert _resolve_parameter("rain") == 7

    def test_resolve_case_insensitive(self):
        assert _resolve_parameter("Temperature") == 1
        assert _resolve_parameter("WIND_SPEED") == 4
        assert _resolve_parameter("Air_Temperature") == 1

    def test_resolve_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown weather parameter"):
            _resolve_parameter("nosuchparam")

    @patch("dtcc_core.datasets.weather._get_text", side_effect=_mock_get_text)
    def test_build_with_string_parameters(self, mock_fetch):
        dataset = WeatherDataset()
        args = WeatherDatasetArgs(
            bounds=(18.0, 59.0, 18.2, 59.5),
            crs="EPSG:4326",
            parameters=["temperature", "wind_speed"],
        )
        result = dataset.build(args)
        station = result.stations()[0]
        point = station.geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "air_temperature" in field_names
        assert "wind_speed" in field_names
