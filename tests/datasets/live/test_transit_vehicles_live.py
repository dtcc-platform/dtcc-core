"""Credentialed live public-transport dataset checks.

Run with ``DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live``.
Provider credentials are still required for each credentialed provider case.
"""

from __future__ import annotations

import os

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.model.object import VehicleCollection
from tests.datasets.live.live_policy import skip_if_transient

pytestmark = pytest.mark.live


def _require_env(name: str, provider: str) -> str:
    value = os.environ.get(name)
    if value is None or not value.strip():
        pytest.skip(
            f"Set {name} to run the {provider} live vehicle test. "
            "Do not commit credentials."
        )
    return value.strip()


def _require_one_env(names: tuple[str, ...], provider: str) -> str:
    for name in names:
        value = os.environ.get(name)
        if value is not None and value.strip():
            return value.strip()
    pytest.skip(
        f"Set one of {', '.join(names)} to run the {provider} live vehicle test. "
        "Do not commit credentials."
    )


def _live_buses(provider_name: str, **kwargs) -> VehicleCollection:
    """Fetch live buses, skipping on transient upstream failures."""
    try:
        return datasets.buses(max_vehicles=25, strict_live=True, **kwargs)
    except Exception as exc:
        skip_if_transient(exc, f"buses ({provider_name})")
        raise


def _assert_vehicle_collection_contract(vehicles: VehicleCollection) -> None:
    assert isinstance(vehicles, VehicleCollection)
    assert vehicles.attributes["partial_result"] is False
    assert vehicles.attributes["upstream_error_count"] == 0
    assert isinstance(vehicles.attributes["providers"], list)


def test_trafiklab_live_buses_with_credentials():
    api_key = _require_one_env(
        ("TRAFIKLAB_API_KEY", "SAMTRAFIKEN_API_KEY"),
        "Trafiklab",
    )

    vehicles = _live_buses(
        "Trafiklab",
        bounds=(17.85, 59.25, 18.15, 59.45),
        crs="EPSG:4326",
        provider="trafiklab",
        api_key=api_key,
    )

    _assert_vehicle_collection_contract(vehicles)


def test_vasttrafik_live_buses_with_credentials():
    authentication_key = _require_env(
        "VASTTRAFIK_AUTHENTICATION_KEY",
        "Västtrafik",
    )

    vehicles = _live_buses(
        "Västtrafik",
        bounds=(11.85, 57.65, 12.05, 57.78),
        crs="EPSG:4326",
        provider="vasttrafik",
        vasttrafik_authentication_key=authentication_key,
    )

    _assert_vehicle_collection_contract(vehicles)
