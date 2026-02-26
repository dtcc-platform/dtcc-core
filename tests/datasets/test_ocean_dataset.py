"""Tests for the ocean dataset (SMHI OcObs latest-hour)."""

import pytest
import math
import numpy as np
from unittest.mock import patch, Mock

from dtcc_core.datasets.ocean import (
    OceanDataset,
    OceanDatasetArgs,
    _parse_latest_hour_csv,
    _resolve_parameter,
    PARAMETER_NAMES,
    _NAME_TO_ID,
)


# ── Fixtures ─────────────────────────────────────────────────────────────

# Minimal CSV that mimics the real SMHI OcObs latest-hour response for
# parameter 5 (sea temperature).  Structure:
#   - 2 metadata rows
#   - blank line
#   - column header (StationsId;Stationsnamn;Latitude;Longitude;<Param>;Kvalitet)
#   - station rows with metadata comments appended after ;;

SAMPLE_CSV_PARAM5 = """\
Parameternamn;Beskrivning;Enhet
Havstemperatur;null;°C

StationsId;Stationsnamn;Latitude;Longitude;Havstemperatur;Kvalitet;;Tidsperiod (fr.o.m.) = 2010-06-17 11:00:00 (UTC)
2541;UDDEVALLA;58.3475;11.8948;0.9;O;;Tidsperiod (t.o.m.) = 2026-02-13 12:00:00 (UTC)
33084;ONSALA;57.392;11.919;0.21;G;;
33089;Göteborg-Krossholmen;57.6913;11.7712;-0.52;O;;Kvalitetskoderna:
2507;LANDSORT NORRA;58.7687;17.8589;0.26;O;;
2088;KUNGSHOLMSFORT;56.1052;15.5893;-0.28;O;;
2099;BARSEBÄCK;55.7564;12.9033;4.6;O;;
"""

SAMPLE_CSV_PARAM6 = """\
Parameternamn;Beskrivning;Enhet
Havsvattenstånd;null;cm

StationsId;Stationsnamn;Latitude;Longitude;Havsvattenstånd;Kvalitet;;Tidsperiod (fr.o.m.) = 1886-12-01 00:00:00 (UTC)
33097;Göteborg-Hisingsbron;57.7149;11.9687;-3.8;O;;Tidsperiod (t.o.m.) = 2026-02-13 13:00:00 (UTC)
2541;UDDEVALLA;58.3475;11.8948;-9.0;O;;
33084;ONSALA;57.392;11.919;-7.5;G;;
2507;LANDSORT NORRA;58.7687;17.8589;-35.7;O;;
"""

SAMPLE_CSV_EMPTY = """\
Parameternamn;Beskrivning;Enhet
Våghöjd, signifikant 30 min;null;m

StationsId;Stationsnamn;Latitude;Longitude;Våghöjd, signifikant 30 min;Kvalitet;;
"""

SAMPLE_CSV_MISSING_VALUES = """\
Parameternamn;Beskrivning;Enhet
Havstemperatur;null;°C

StationsId;Stationsnamn;Latitude;Longitude;Havstemperatur;Kvalitet;;
100;TEST STATION;59.0;18.0;;G;;
"""


# ── Mock helper ──────────────────────────────────────────────────────────


def _mock_get_text(url: str, timeout_s: float = 10.0) -> str:
    """Mock _get_text that returns sample data based on URL patterns."""
    if "/parameter/5/" in url:
        return SAMPLE_CSV_PARAM5
    if "/parameter/6/" in url:
        return SAMPLE_CSV_PARAM6
    if "/parameter/1/" in url:
        return SAMPLE_CSV_EMPTY
    raise RuntimeError(f"Unexpected URL in mock: {url}")


# ── CSV parser tests ─────────────────────────────────────────────────────


class TestParseLatestHourCsv:
    """Tests for _parse_latest_hour_csv."""

    def test_meta_extraction(self):
        meta, _ = _parse_latest_hour_csv(SAMPLE_CSV_PARAM5)
        assert meta["parameter_name"] == "Havstemperatur"
        assert meta["unit"] == "°C"

    def test_meta_description(self):
        meta, _ = _parse_latest_hour_csv(SAMPLE_CSV_PARAM5)
        assert meta["description"] == "null"

    def test_station_count(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM5)
        assert len(records) == 6

    def test_station_fields(self):
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM5)
        rec = records[0]
        assert rec["station_id"] == 2541
        assert rec["station_name"] == "UDDEVALLA"
        assert rec["lat"] == pytest.approx(58.3475)
        assert rec["lon"] == pytest.approx(11.8948)
        assert rec["value"] == pytest.approx(0.9)
        assert rec["quality"] == "O"

    def test_negative_value(self):
        """Negative values (e.g. -0.52°C) should parse correctly."""
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM5)
        krossholmen = [
            r for r in records if r["station_name"] == "Göteborg-Krossholmen"
        ][0]
        assert krossholmen["value"] == pytest.approx(-0.52)

    def test_sea_level_csv(self):
        meta, records = _parse_latest_hour_csv(SAMPLE_CSV_PARAM6)
        assert meta["parameter_name"] == "Havsvattenstånd"
        assert meta["unit"] == "cm"
        assert len(records) == 4

    def test_empty_csv(self):
        """CSV with no station rows should return empty list."""
        meta, records = _parse_latest_hour_csv(SAMPLE_CSV_EMPTY)
        assert len(records) == 0
        assert meta["parameter_name"] == "Våghöjd, signifikant 30 min"

    def test_missing_value(self):
        """Empty value cell should parse as NaN."""
        _, records = _parse_latest_hour_csv(SAMPLE_CSV_MISSING_VALUES)
        assert len(records) == 1
        assert math.isnan(records[0]["value"])

    def test_garbage_input(self):
        """Non-CSV input should return empty results."""
        meta, records = _parse_latest_hour_csv("this is not csv data")
        assert len(records) == 0


# ── Parameter resolution tests ───────────────────────────────────────────


class TestResolveParameter:
    """Tests for _resolve_parameter."""

    def test_int_passthrough(self):
        assert _resolve_parameter(5) == 5
        assert _resolve_parameter(6) == 6

    def test_str_numeric(self):
        assert _resolve_parameter("5") == 5
        assert _resolve_parameter("  6  ") == 6

    def test_canonical_name(self):
        assert _resolve_parameter("sea_temperature") == 5
        assert _resolve_parameter("sea_level") == 6
        assert _resolve_parameter("salinity") == 4

    def test_alias(self):
        assert _resolve_parameter("temperature") == 5
        assert _resolve_parameter("level") == 6
        assert _resolve_parameter("waves") == 1
        assert _resolve_parameter("oxygen") == 15

    def test_case_insensitive(self):
        assert _resolve_parameter("SEA_TEMPERATURE") == 5
        assert _resolve_parameter("Sea_Level") == 6

    def test_hyphen_and_space(self):
        assert _resolve_parameter("sea-temperature") == 5
        assert _resolve_parameter("sea level") == 6

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown ocean parameter"):
            _resolve_parameter("nonexistent")

    def test_all_parameter_names_resolve(self):
        """Every entry in PARAMETER_NAMES must be resolvable."""
        for pid, name in PARAMETER_NAMES.items():
            assert _resolve_parameter(name) == pid

    def test_all_aliases_resolve(self):
        """Every entry in _NAME_TO_ID must be resolvable."""
        for name, pid in _NAME_TO_ID.items():
            assert _resolve_parameter(name) == pid


# ── Build tests with mocked HTTP ─────────────────────────────────────────


class TestOceanBuild:
    """Tests for OceanDataset.build with mocked HTTP."""

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_basic_build(self, mock_text):
        """Build with default parameters using WGS84 bounds around west coast."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        stations = sc.stations()
        # Should include stations in the west coast bbox
        assert len(stations) >= 1
        snames = {s.attributes["station_name"] for s in stations}
        assert "UDDEVALLA" in snames or "ONSALA" in snames

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_field_on_station(self, mock_text):
        """Each station should have a field for each parameter."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        st = sc.stations()[0]
        point = st.geometry["location"]
        assert len(point.fields) >= 1
        f = point.fields[0]
        assert f.name == "sea_temperature"
        assert f.unit == "°C"

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_multi_parameter(self, mock_text):
        """Build with two parameters should yield two fields for shared stations."""
        ds = OceanDataset()
        # Wide bbox to capture stations present in both params
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5, 6],
            )
        )
        # Find a station that has both params (UDDEVALLA or ONSALA)
        for st in sc.stations():
            point = st.geometry["location"]
            field_names = {f.name for f in point.fields}
            if "sea_temperature" in field_names and "sea_level" in field_names:
                break
        else:
            # At least one station should have both
            shared = {s.attributes["station_name"] for s in sc.stations()}
            pytest.fail(f"No station has both fields. Stations: {shared}")

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_bbox_filter(self, mock_text):
        """Only stations within bbox should be included."""
        ds = OceanDataset()
        # Tight bbox that only includes ONSALA (57.392, 11.919)
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.8, 57.3, 12.0, 57.5),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "ONSALA" in snames
        assert "UDDEVALLA" not in snames
        assert "LANDSORT NORRA" not in snames

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_drop_missing(self, mock_text):
        """Stations with NaN values should be dropped when drop_missing=True."""

        # Use the missing-values fixture by injecting it for param 5
        def mock_missing(url, timeout_s=10.0):
            if "/parameter/5/" in url:
                return SAMPLE_CSV_MISSING_VALUES
            return _mock_get_text(url, timeout_s)

        mock_text.side_effect = mock_missing
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(17.0, 58.0, 19.0, 60.0),
                crs="EPSG:4326",
                parameters=[5],
                drop_missing=True,
            )
        )
        # TEST STATION has empty value → should be dropped
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "TEST STATION" not in snames

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_keep_missing(self, mock_text):
        """Stations with NaN values kept when drop_missing=False."""

        def mock_missing(url, timeout_s=10.0):
            if "/parameter/5/" in url:
                return SAMPLE_CSV_MISSING_VALUES
            return _mock_get_text(url, timeout_s)

        mock_text.side_effect = mock_missing
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(17.0, 58.0, 19.0, 60.0),
                crs="EPSG:4326",
                parameters=[5],
                drop_missing=False,
            )
        )
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "TEST STATION" in snames

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_smhi_field_name_style(self, mock_text):
        """field_name_style='smhi' should use Swedish parameter names."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
                field_name_style="smhi",
            )
        )
        st = sc.stations()[0]
        point = st.geometry["location"]
        assert point.fields[0].name == "Havstemperatur"

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_protobuf_export(self, mock_text):
        """format='pb' should return bytes."""
        ds = OceanDataset()
        result = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
                format="pb",
            )
        )
        assert isinstance(result, bytes)
        assert len(result) > 0

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_max_stations_limit(self, mock_text):
        """max_stations should cap the number of result stations."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(10.0, 55.0, 20.0, 60.0),
                crs="EPSG:4326",
                parameters=[5],
                max_stations=2,
            )
        )
        assert len(sc.stations()) <= 2

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_quality_attribute(self, mock_text):
        """Quality codes should be stored as station attributes."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        st = sc.stations()[0]
        # At least one quality attribute should exist
        q_attrs = [k for k in st.attributes if k.startswith("q_")]
        assert len(q_attrs) >= 1

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_coordinate_reproject(self, mock_text):
        """EPSG:3006 output should have reprojected coordinates."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(270000, 6350000, 330000, 6500000),
                crs="EPSG:3006",
                parameters=[5],
            )
        )
        stations = sc.stations()
        if len(stations) > 0:
            pt = stations[0].geometry["location"]
            # EPSG:3006 coordinates should be in typical Swedish range
            assert 200000 < pt.x < 900000
            assert 6100000 < pt.y < 7700000

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_collection_attributes(self, mock_text):
        """SensorCollection should have expected metadata attributes."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5, 6],
            )
        )
        assert sc.attributes["source"] == "smhi_ocobs"
        assert sc.attributes["dataset"] == "ocean"
        assert "parameter_fields" in sc.attributes
        pf = sc.attributes["parameter_fields"]
        assert "sea_temperature" in pf
        assert "sea_level" in pf

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_parameter_name_strings(self, mock_text):
        """Parameters can be specified as name strings."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=["temperature", "level"],
            )
        )
        assert len(sc.stations()) >= 1
        st = sc.stations()[0]
        point = st.geometry["location"]
        field_names = {f.name for f in point.fields}
        # Should have at least one of these
        assert "sea_temperature" in field_names or "sea_level" in field_names

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_http_failure_graceful(self, mock_text):
        """Fetch failure for a parameter should be skipped gracefully."""

        def _fail_on_p6(url, timeout_s=10.0):
            if "/parameter/6/" in url:
                raise RuntimeError("Network error")
            return _mock_get_text(url, timeout_s=timeout_s)

        mock_text.side_effect = _fail_on_p6
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5, 6],
            )
        )
        assert len(sc.stations()) >= 1
        st = sc.stations()[0]
        point = st.geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "sea_temperature" in field_names
        assert "sea_level" not in field_names

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_z_is_zero(self, mock_text):
        """Ocean stations should have z=0."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        for st in sc.stations():
            pt = st.geometry["location"]
            assert pt.z == pytest.approx(0.0)


# ── Registration tests ───────────────────────────────────────────────────


class TestOceanRegistration:
    """Tests for dataset registration."""

    def test_registered_name(self):
        ds = OceanDataset()
        assert ds.name == "ocean"

    def test_get_dataset(self):
        from dtcc_core.datasets import get_dataset

        ds = get_dataset("ocean")
        assert ds is not None
        assert ds.name == "ocean"

    def test_module_attribute(self):
        import dtcc_core.datasets as datasets

        assert hasattr(datasets, "ocean")


# ── Str representation tests ─────────────────────────────────────────────


class TestOceanStr:
    """Test that SensorCollection prints nicely."""

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_str_contains_station_info(self, mock_text):
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        text = str(sc)
        assert "SensorCollection" in text


# ── Edge case tests ──────────────────────────────────────────────────────


class TestOceanEdgeCases:
    """Edge case tests."""

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_empty_bbox(self, mock_text):
        """Bbox with no stations should return empty collection."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        assert len(sc.stations()) == 0

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_empty_parameter_response(self, mock_text):
        """Parameter with no stations should produce empty collection."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(10.0, 55.0, 20.0, 70.0),
                crs="EPSG:4326",
                parameters=[1],  # wave height — empty in our mock
            )
        )
        assert len(sc.stations()) == 0

    @patch("dtcc_core.datasets.ocean._get_text", side_effect=_mock_get_text)
    def test_value_attribute_set(self, mock_text):
        """Station 'value' attribute should be set for __str__ display."""
        ds = OceanDataset()
        sc = ds.build(
            OceanDatasetArgs(
                bounds=(11.0, 57.0, 12.5, 59.0),
                crs="EPSG:4326",
                parameters=[5],
            )
        )
        st = sc.stations()[0]
        assert "value" in st.attributes
