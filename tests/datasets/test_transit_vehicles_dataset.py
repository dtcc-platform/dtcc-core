import importlib

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.transport.base import TransportProviderResult, VehicleRecord
from dtcc_core.datasets.transit_vehicles import TransitVehiclesArgs
from dtcc_core.model.object import VehicleCollection


def test_transit_vehicles_registered():
    assert callable(datasets.transit_vehicles)
    assert callable(datasets.buses)
    assert callable(datasets.trams)
    assert callable(datasets.trains)
    assert callable(datasets.metros)
    assert callable(datasets.ferries)


def test_transit_vehicles_metadata():
    metadata = datasets.transit_vehicles.describe()

    assert metadata["name"] == "transit_vehicles"
    assert metadata["data_category"] == "raw"
    assert metadata["result_kind"] == "vehicle_collection"
    assert metadata["python_return_type"] == "dtcc_core.model.VehicleCollection"
    assert "Planera Resa v4" in metadata["description"]
    assert "VASTTRAFIK_AUTHENTICATION_KEY" in metadata["description"]
    assert "pb" in metadata["supported_formats"]
    assert (
        "Planera Resa v4"
        in metadata["args_schema"]["properties"]["vasttrafik_authentication_key"][
            "description"
        ]
    )


def test_args_accept_modes():
    args = TransitVehiclesArgs(
        bounds=(17.9, 59.2, 18.2, 59.4),
        crs="EPSG:4326",
        modes=("bus", "tram"),
    )
    assert args.modes == ("bus", "tram")


def test_dataset_builds_vehicle_collection(monkeypatch):
    def fake_fetch(**kwargs):
        assert kwargs["requested_modes"] == ("bus",)
        return TransportProviderResult(
            records=[
                VehicleRecord(
                    vehicle_id="vehicle-1",
                    lon=18.0,
                    lat=59.3,
                    provider="trafiklab",
                    mode="bus",
                    route_id="1",
                    line="1",
                    speed=7.5,
                    bearing=90.0,
                )
            ],
            metadata={"provider": "trafiklab", "operators": ["sl"]},
        )

    module = importlib.import_module("dtcc_core.datasets.transit_vehicles")
    monkeypatch.setattr(module, "fetch_trafiklab_gtfs_vehicles", fake_fetch)

    vehicles = datasets.transit_vehicles(
        bounds=(17.9, 59.2, 18.2, 59.4),
        crs="EPSG:4326",
        provider="trafiklab",
        modes=("bus",),
    )

    assert isinstance(vehicles, VehicleCollection)
    assert len(vehicles.vehicles()) == 1
    vehicle = vehicles.vehicles()[0]
    assert vehicle.attributes["mode"] == "bus"
    assert vehicle.attributes["line"] == "1"
    assert vehicle.geometry["location"].x == pytest.approx(18.0)
    assert vehicle.geometry["location"].y == pytest.approx(59.3)


def test_shortcut_sets_default_mode(monkeypatch):
    seen = {}

    def fake_fetch(**kwargs):
        seen["requested_modes"] = kwargs["requested_modes"]
        return TransportProviderResult()

    module = importlib.import_module("dtcc_core.datasets.transit_vehicles")
    monkeypatch.setattr(module, "fetch_trafiklab_gtfs_vehicles", fake_fetch)

    datasets.buses(
        bounds=(17.9, 59.2, 18.2, 59.4),
        crs="EPSG:4326",
        provider="trafiklab",
    )

    assert seen["requested_modes"] == ("bus",)


def test_shortcut_rejects_conflicting_modes():
    with pytest.raises(ValueError):
        datasets.buses(
            bounds=(17.9, 59.2, 18.2, 59.4),
            crs="EPSG:4326",
            provider="trafiklab",
            modes=("tram",),
        )


def test_missing_trafiklab_key_degrades(monkeypatch):
    monkeypatch.delenv("TRAFIKLAB_API_KEY", raising=False)
    monkeypatch.delenv("SAMTRAFIKEN_API_KEY", raising=False)

    vehicles = datasets.transit_vehicles(
        bounds=(17.9, 59.2, 18.2, 59.4),
        crs="EPSG:4326",
        provider="trafiklab",
        modes=("bus",),
    )

    assert isinstance(vehicles, VehicleCollection)
    assert len(vehicles.vehicles()) == 0
    assert vehicles.attributes["partial_result"] is True
    assert vehicles.attributes["upstream_error_count"] == 1
    assert "TRAFIKLAB_API_KEY" in str(vehicles)


def test_missing_vasttrafik_credentials_degrades(monkeypatch):
    monkeypatch.delenv("VASTTRAFIK_AUTHENTICATION_KEY", raising=False)

    vehicles = datasets.buses(
        bounds=(11.9, 57.6, 12.1, 57.8),
        crs="EPSG:4326",
        provider="vasttrafik",
    )

    assert isinstance(vehicles, VehicleCollection)
    assert len(vehicles.vehicles()) == 0
    assert vehicles.attributes["partial_result"] is True
    assert vehicles.attributes["upstream_error_count"] == 1
    assert "VASTTRAFIK_AUTHENTICATION_KEY" in str(vehicles)


def test_vasttrafik_portal_credential_names_are_used():
    seen = {}

    def fake_post(url, data=None, headers=None, auth=None, timeout=None):
        seen["headers"] = headers
        seen["auth"] = auth

        class Response:
            def raise_for_status(self):
                pass

            def json(self):
                return {"access_token": "token"}

        return Response()

    fake_requests = type(
        "Requests",
        (),
        {"post": staticmethod(fake_post)},
    )

    module = importlib.import_module("dtcc_core.datasets.transport.vasttrafik")
    module._fetch_token(
        fake_requests,
        authentication_key=module._resolve_authentication_key(
            "Authorization: Basic copied-from-portal"
        ),
        timeout_s=1.0,
    )

    assert seen["headers"] == {"Authorization": "Basic copied-from-portal"}
    assert seen["auth"] is None


def test_vasttrafik_v4_positions_are_parsed():
    module = importlib.import_module("dtcc_core.datasets.transport.vasttrafik")

    records = module._parse_positions(
        [
            {
                "detailsReference": "journey-1",
                "line": {
                    "name": "16",
                    "detailName": "Rosa Express",
                    "transportMode": "bus",
                },
                "name": "16",
                "direction": "Sahlgrenska",
                "latitude": 57.7,
                "longitude": 11.98,
            }
        ]
    )

    assert len(records) == 1
    assert records[0].vehicle_id == "journey-1"
    assert records[0].mode == "bus"
    assert records[0].line == "16"
    assert records[0].destination == "Sahlgrenska"


def test_vasttrafik_positions_403_explains_subscription():
    class RequestException(Exception):
        pass

    class HTTPError(RequestException):
        def __init__(self):
            self.response = type("Response", (), {"status_code": 403})()

    class Response:
        def raise_for_status(self):
            raise HTTPError()

    Requests = type(
        "Requests",
        (),
        {
            "RequestException": RequestException,
            "get": staticmethod(lambda url, params=None, headers=None, timeout=None: Response()),
        },
    )

    module = importlib.import_module("dtcc_core.datasets.transport.vasttrafik")

    with pytest.raises(DatasetUpstreamError) as exc:
        module._fetch_positions(
            Requests,
            token="token",
            bounds_wgs84=(11.9, 57.6, 12.1, 57.8),
            requested_modes=("bus",),
            timeout_s=1.0,
        )

    assert exc.value.status_code == 403
    assert "Planera Resa v4" in exc.value.message
    assert "Subscribe" in exc.value.message


def test_missing_trafiklab_key_strict_raises(monkeypatch):
    monkeypatch.delenv("TRAFIKLAB_API_KEY", raising=False)
    monkeypatch.delenv("SAMTRAFIKEN_API_KEY", raising=False)

    with pytest.raises(DatasetUpstreamError):
        datasets.transit_vehicles(
            bounds=(17.9, 59.2, 18.2, 59.4),
            crs="EPSG:4326",
            provider="trafiklab",
            modes=("bus",),
            strict_live=True,
        )
