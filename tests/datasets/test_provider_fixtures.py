"""Contract tests for committed provider parser fixtures."""

from __future__ import annotations

import json
from pathlib import Path


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "smhi"

TEXT_FIXTURES = [
    "metobs/latest_hour_parameter_1.csv",
    "metobs/latest_hour_parameter_4.csv",
    "ocobs/latest_hour_parameter_5.csv",
    "ocobs/latest_hour_parameter_6.csv",
    "ocobs/empty_parameter_1.csv",
    "ocobs/missing_values_parameter_5.csv",
]

JSON_FIXTURES = [
    "hydroobs/parameter_1_stations.json",
    "hydroobs/station_2100_parameter_1_latest_day.json",
    "hydroobs/station_2100_parameter_3_latest_day.json",
    "hydroobs/station_2200_parameter_1_latest_day_empty.json",
    "air_quality/stations.json",
    "air_quality/station1.json",
    "air_quality/station2.json",
    "air_quality/timeseries_ts1.json",
    "air_quality/timeseries_ts2.json",
    "air_quality/timeseries_ts1_stale.json",
    "air_quality/timeseries_ts1_no_current.json",
    "air_quality/getdata_ts1_empty.json",
    "air_quality/getdata_ts1_recent.json",
    "air_quality/phenomena.json",
]


def test_text_provider_fixtures_exist_and_are_nonempty():
    for fixture in TEXT_FIXTURES:
        path = FIXTURE_ROOT / fixture
        assert path.exists(), f"Missing provider fixture: {path}"
        assert path.read_text(encoding="utf-8").strip(), (
            f"Provider fixture is empty: {path}"
        )


def test_json_provider_fixtures_exist_and_parse():
    for fixture in JSON_FIXTURES:
        path = FIXTURE_ROOT / fixture
        assert path.exists(), f"Missing provider fixture: {path}"
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload not in ({}, []), f"Provider fixture is empty: {path}"
