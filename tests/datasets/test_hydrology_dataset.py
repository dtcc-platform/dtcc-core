"""Tests for the hydrology dataset (SMHI HydroObs latest-day)."""

import pytest
import math
import json
import numpy as np
from unittest.mock import patch, Mock, MagicMock

from dtcc_core.datasets.hydrology import (
    HydrologyDataset,
    HydrologyDatasetArgs,
    _resolve_parameter,
    _fetch_station_list,
    _fetch_latest_day,
    PARAMETER_NAMES,
    _NAME_TO_ID,
)


# ── Fixtures ─────────────────────────────────────────────────────────────

# Minimal station list JSON mimicking the real SMHI HydroObs response for
# parameter 1 (daily discharge).  Contains stations in different locations
# so bbox filtering can be tested.

SAMPLE_STATION_LIST = {
    "key": "1",
    "title": "Vattenföring (Dygn)",
    "unit": "m³/s",
    "station": [
        {
            "key": "2357",
            "name": "ABISKO",
            "id": 2357,
            "latitude": 68.1936,
            "longitude": 19.9859,
            "active": True,
            "owner": "SMHI",
            "measuringStations": "CORE",
            "catchmentName": "Abiskojokk",
            "catchmentNumber": "123",
            "catchmentSize": "560",
            "from": 946684800000,
            "to": 1739404800000,
            "link": [],
        },
        {
            "key": "2100",
            "name": "STOCKHOLM STN",
            "id": 2100,
            "latitude": 59.33,
            "longitude": 18.07,
            "active": True,
            "owner": "SMHI",
            "measuringStations": "CORE",
            "catchmentName": "Norrström",
            "catchmentNumber": "456",
            "catchmentSize": "22600",
            "from": 946684800000,
            "to": 1739404800000,
            "link": [],
        },
        {
            "key": "2200",
            "name": "MALMÖ STN",
            "id": 2200,
            "latitude": 55.60,
            "longitude": 13.00,
            "active": False,
            "owner": "SMHI",
            "measuringStations": "ADDITIONAL",
            "catchmentName": "Höje å",
            "catchmentNumber": "789",
            "catchmentSize": "310",
            "from": 946684800000,
            "to": 1739404800000,
            "link": [],
        },
    ],
}

# Minimal latest-day JSON response for a single station.
SAMPLE_LATEST_DAY_P1 = {
    "updated": 1739404800000,
    "parameter": {
        "key": "1",
        "name": "Vattenföring (Dygn)",
        "summary": "Dygnsmedelvärde",
        "unit": "m³/s",
    },
    "station": {
        "key": "2100",
        "name": "STOCKHOLM STN",
    },
    "period": {"key": "latest-day"},
    "position": [],
    "link": [],
    "value": [
        {"date": 1739318400000, "value": 42.5, "quality": "G"},
    ],
}

SAMPLE_LATEST_DAY_P3 = {
    "updated": 1739404800000,
    "parameter": {
        "key": "3",
        "name": "Vattenstånd",
        "summary": "Momentanvärde",
        "unit": "cm",
    },
    "station": {
        "key": "2100",
        "name": "STOCKHOLM STN",
    },
    "period": {"key": "latest-day"},
    "position": [],
    "link": [],
    "value": [
        {"date": 1739318400000, "value": 125.0, "quality": "O"},
    ],
}

SAMPLE_LATEST_DAY_EMPTY = {
    "updated": 1739404800000,
    "parameter": {
        "key": "1",
        "name": "Vattenföring (Dygn)",
        "summary": "Dygnsmedelvärde",
        "unit": "m³/s",
    },
    "station": {
        "key": "2200",
        "name": "MALMÖ STN",
    },
    "period": {"key": "latest-day"},
    "position": [],
    "link": [],
    "value": [],
}


# ── Mock helper ──────────────────────────────────────────────────────────


def _mock_get_json(url: str, timeout_s: float = 10.0) -> dict:
    """Mock _get_json that returns sample data based on URL patterns."""
    if "/parameter/1.json" in url and "/station/" not in url:
        return SAMPLE_STATION_LIST
    if "/parameter/3.json" in url and "/station/" not in url:
        # Same stations but with different parameter key
        data = json.loads(json.dumps(SAMPLE_STATION_LIST))
        data["key"] = "3"
        data["title"] = "Vattenstånd"
        data["unit"] = "cm"
        return data
    if "/station/2100/" in url and "/parameter/1/" in url:
        return SAMPLE_LATEST_DAY_P1
    if "/station/2100/" in url and "/parameter/3/" in url:
        return SAMPLE_LATEST_DAY_P3
    if "/station/2357/" in url and "/parameter/1/" in url:
        # Abisko station with discharge
        data = json.loads(json.dumps(SAMPLE_LATEST_DAY_P1))
        data["station"]["key"] = "2357"
        data["station"]["name"] = "ABISKO"
        data["value"] = [{"date": 1739318400000, "value": 22.0, "quality": "O"}]
        return data
    if "/station/2357/" in url and "/parameter/3/" in url:
        data = json.loads(json.dumps(SAMPLE_LATEST_DAY_P3))
        data["station"]["key"] = "2357"
        data["station"]["name"] = "ABISKO"
        data["value"] = [{"date": 1739318400000, "value": 88.0, "quality": "O"}]
        return data
    if "/station/2200/" in url:
        return SAMPLE_LATEST_DAY_EMPTY
    raise RuntimeError(f"Unexpected URL in mock: {url}")


# ── Parameter resolution tests ───────────────────────────────────────────


class TestResolveParameter:
    """Tests for _resolve_parameter."""

    def test_int_passthrough(self):
        assert _resolve_parameter(1) == 1
        assert _resolve_parameter(3) == 3

    def test_str_numeric(self):
        assert _resolve_parameter("1") == 1
        assert _resolve_parameter("  3  ") == 3

    def test_canonical_name(self):
        assert _resolve_parameter("discharge_daily") == 1
        assert _resolve_parameter("water_level") == 3
        assert _resolve_parameter("water_temperature") == 4

    def test_alias(self):
        assert _resolve_parameter("discharge") == 1
        assert _resolve_parameter("flow") == 1
        assert _resolve_parameter("level") == 3
        assert _resolve_parameter("temperature") == 4

    def test_case_insensitive(self):
        assert _resolve_parameter("DISCHARGE_DAILY") == 1
        assert _resolve_parameter("Water_Level") == 3

    def test_hyphen_and_space(self):
        assert _resolve_parameter("discharge-daily") == 1
        assert _resolve_parameter("water level") == 3

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown hydrology parameter"):
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


class TestHydrologyBuild:
    """Tests for HydrologyDataset.build with mocked HTTP."""

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_basic_build(self, mock_json):
        """Build with default parameters using Stockholm-area bounds (WGS84)."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        # Only Stockholm (59.33, 18.07) is in the bbox and active
        stations = sc.stations()
        assert len(stations) >= 1
        # Check station attributes
        st = stations[0]
        assert st.attributes["station_name"] == "STOCKHOLM STN"
        assert st.attributes["station_id"] == "2100"
        assert st.attributes["catchment_name"] == "Norrström"

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_field_on_station(self, mock_json):
        """Each station should have a field for each parameter."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        st = sc.stations()[0]
        point = st.geometry["location"]
        assert len(point.fields) == 1
        f = point.fields[0]
        assert f.name == "discharge_daily"
        assert f.unit == "m³/s"
        assert f.values[0] == pytest.approx(42.5, rel=1e-3)

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_multi_parameter(self, mock_json):
        """Build with two parameters should yield two fields per station."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1, 3],
            )
        )
        st = sc.stations()[0]
        point = st.geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "discharge_daily" in field_names
        assert "water_level" in field_names

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_active_filter(self, mock_json):
        """Inactive stations should be excluded when active_only=True."""
        ds = HydrologyDataset()
        # Use a huge bbox that covers all sample stations
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(10.0, 50.0, 25.0, 70.0),
                crs="EPSG:4326",
                parameters=[1],
                active_only=True,
            )
        )
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "MALMÖ STN" not in snames
        assert "STOCKHOLM STN" in snames
        assert "ABISKO" in snames

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_inactive_included_when_flag_off(self, mock_json):
        """Inactive stations should be included when active_only=False."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(10.0, 50.0, 25.0, 70.0),
                crs="EPSG:4326",
                parameters=[1],
                active_only=False,
                drop_missing=False,
            )
        )
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "MALMÖ STN" in snames

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_bbox_filter(self, mock_json):
        """Only stations within bbox should be included."""
        ds = HydrologyDataset()
        # Tight bbox around Stockholm only
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "STOCKHOLM STN" in snames
        assert "ABISKO" not in snames

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_drop_missing(self, mock_json):
        """Stations with no valid values should be dropped when drop_missing=True."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(10.0, 50.0, 25.0, 70.0),
                crs="EPSG:4326",
                parameters=[1],
                active_only=False,
                drop_missing=True,
            )
        )
        # MALMÖ STN returns empty values → should be dropped
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "MALMÖ STN" not in snames

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_keep_missing(self, mock_json):
        """Stations with no valid values kept when drop_missing=False."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(10.0, 50.0, 25.0, 70.0),
                crs="EPSG:4326",
                parameters=[1],
                active_only=False,
                drop_missing=False,
            )
        )
        snames = {s.attributes["station_name"] for s in sc.stations()}
        assert "MALMÖ STN" in snames

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_smhi_field_name_style(self, mock_json):
        """field_name_style='smhi' should use Swedish parameter names."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
                field_name_style="smhi",
            )
        )
        st = sc.stations()[0]
        point = st.geometry["location"]
        assert point.fields[0].name == "Vattenföring (Dygn)"

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_protobuf_export(self, mock_json):
        """format='pb' should return bytes."""
        ds = HydrologyDataset()
        result = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
                format="pb",
            )
        )
        assert isinstance(result, bytes)
        assert len(result) > 0

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_max_stations_limit(self, mock_json):
        """max_stations should cap the number of result stations."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(10.0, 50.0, 25.0, 70.0),
                crs="EPSG:4326",
                parameters=[1],
                max_stations=1,
            )
        )
        assert len(sc.stations()) <= 1

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_quality_attribute(self, mock_json):
        """Quality codes should be stored as station attributes."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        st = sc.stations()[0]
        assert st.attributes["q_discharge_daily"] == "G"

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_coordinate_reproject(self, mock_json):
        """EPSG:3006 output should have reprojected coordinates."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                # Use EPSG:3006 bounds roughly around Stockholm
                bounds=(670000, 6570000, 690000, 6600000),
                crs="EPSG:3006",
                parameters=[1],
            )
        )
        stations = sc.stations()
        if len(stations) > 0:
            pt = stations[0].geometry["location"]
            # EPSG:3006 easting should be in the ~670k–690k range
            assert 600000 < pt.x < 800000
            assert 6500000 < pt.y < 6700000

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_collection_attributes(self, mock_json):
        """SensorCollection should have expected metadata attributes."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1, 3],
            )
        )
        assert sc.attributes["source"] == "smhi_hydroobs"
        assert sc.attributes["dataset"] == "hydrology"
        assert "parameter_fields" in sc.attributes
        pf = sc.attributes["parameter_fields"]
        assert "discharge_daily" in pf
        assert "water_level" in pf

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_parameter_name_strings(self, mock_json):
        """Parameters can be specified as name strings."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=["discharge", "level"],
            )
        )
        st = sc.stations()[0]
        point = st.geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "discharge_daily" in field_names
        assert "water_level" in field_names

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_http_failure_graceful(self, mock_json):
        """Station list failure for a parameter should be skipped gracefully."""
        calls = [0]
        original = _mock_get_json

        def _fail_on_p3(url, timeout_s=10.0):
            if "/parameter/3.json" in url and "/station/" not in url:
                raise RuntimeError("Network error")
            return original(url, timeout_s=timeout_s)

        mock_json.side_effect = _fail_on_p3
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1, 3],
            )
        )
        # Should still have stations from parameter 1
        assert len(sc.stations()) >= 1
        point = sc.stations()[0].geometry["location"]
        field_names = {f.name for f in point.fields}
        assert "discharge_daily" in field_names
        # Parameter 3 failed → no water_level field
        assert "water_level" not in field_names


# ── Registration tests ───────────────────────────────────────────────────


class TestHydrologyRegistration:
    """Tests for dataset registration."""

    def test_registered_name(self):
        ds = HydrologyDataset()
        assert ds.name == "hydrology"

    def test_get_dataset(self):
        from dtcc_core.datasets import get_dataset

        ds = get_dataset("hydrology")
        assert ds is not None
        assert ds.name == "hydrology"

    def test_module_attribute(self):
        import dtcc_core.datasets as datasets

        assert hasattr(datasets, "hydrology")


# ── Str representation tests ─────────────────────────────────────────────


class TestHydrologyStr:
    """Test that SensorCollection prints nicely."""

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_str_contains_station_info(self, mock_json):
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        text = str(sc)
        assert "SensorCollection" in text
        assert "1 station" in text or "stations" in text.lower()


# ── Edge case tests ──────────────────────────────────────────────────────


class TestHydrologyEdgeCases:
    """Edge case tests."""

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_empty_bbox(self, mock_json):
        """Bbox with no stations should return empty collection."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        assert len(sc.stations()) == 0

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_value_attribute_set(self, mock_json):
        """Station 'value' attribute should be set for __str__ display."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        st = sc.stations()[0]
        assert "value" in st.attributes
        assert st.attributes["value"] == pytest.approx(42.5, rel=1e-3)

    @patch("dtcc_core.datasets.hydrology._get_json", side_effect=_mock_get_json)
    def test_catchment_metadata(self, mock_json):
        """Station should carry catchment metadata."""
        ds = HydrologyDataset()
        sc = ds.build(
            HydrologyDatasetArgs(
                bounds=(17.5, 59.0, 18.5, 59.5),
                crs="EPSG:4326",
                parameters=[1],
            )
        )
        st = sc.stations()[0]
        assert st.attributes["catchment_name"] == "Norrström"
        assert st.attributes["catchment_number"] == "456"
        assert st.attributes["catchment_size"] == "22600"
