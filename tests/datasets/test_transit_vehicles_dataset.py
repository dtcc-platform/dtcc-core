import importlib

import pytest
from pydantic import ValidationError

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.transport.base import TransportProviderResult, VehicleRecord
from dtcc_core.datasets.transit_vehicles import (
    TransitVehiclesArgs,
    TransitVehiclesDataset,
)
from dtcc_core.model.geometry import Point
from dtcc_core.model.object import Object, VehicleCollection


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


def test_transit_vehicles_context_metadata_and_presentation():
    dataset = TransitVehiclesDataset()
    context = dataset.create_context(
        dataset.validate({"bounds": (17.9, 59.2, 18.2, 59.4)})
    )
    manifest = context.manifest()

    assert manifest.identity.title == "Live Transit Vehicles"
    assert {item["name"] for item in manifest.metadata.provider} == {
        "Trafiklab",
        "Västtrafik",
    }
    assert manifest.metadata.source[0]["service"] == "gtfs-rt"
    assert manifest.metadata.source[1]["service"] == "planera-resa-v4"
    assert "TRAFIKLAB_API_KEY" in manifest.metadata.source[0]["credential_env"]
    assert "VASTTRAFIK_AUTHENTICATION_KEY" in (
        manifest.metadata.source[1]["credential_env"]
    )
    assert manifest.metadata.collection_period.startswith("Live snapshot")
    assert "vehicle_collection" in manifest.metadata.data_types
    assert manifest.presentation.headline == "Live Public Transport Vehicles"
    assert manifest.presentation.legend["title"] == "Transit vehicle attributes"
    assert manifest.presentation.view_hints["table_role"] == "live_mobility"
    assert manifest.presentation.warnings
    assert manifest.presentation.limitations


def test_vehicle_collection_plot_uses_presentation_panel_by_default():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)

    collection = VehicleCollection()
    vehicle = Object()
    vehicle.attributes = {"mode": "bus", "bearing": 90.0}
    vehicle.geometry["position"] = Point(x=18.0, y=59.3, z=0.0)
    collection.add_vehicle(vehicle)
    context = TransitVehiclesDataset().create_context(
        TransitVehiclesDataset().validate({"bounds": (17.9, 59.2, 18.2, 59.4)})
    )
    attach_dataset_context(collection, context)

    ax = collection.plot(show=False)
    simple_ax = collection.plot(show=False, presentation=False)

    assert len(ax.figure.axes) >= 2
    assert len(simple_ax.figure.axes) == 1


@pytest.mark.parametrize(
    ("dataset_name", "mode", "title"),
    [
        ("buses", "bus", "Buses"),
        ("trams", "tram", "Trams"),
        ("trains", "train", "Trains"),
        ("metros", "metro", "Metros"),
        ("ferries", "ferry", "Ferries"),
    ],
)
def test_shortcut_context_documents_mode_preset(dataset_name, mode, title):
    dataset = datasets.get_dataset(dataset_name)
    context = dataset.create_context(
        dataset.validate({"bounds": (17.9, 59.2, 18.2, 59.4)})
    )
    manifest = context.manifest()

    assert manifest.identity.name == dataset_name
    assert manifest.identity.title == title
    assert f"modes=('{mode}',)" in manifest.metadata.description
    assert manifest.presentation.headline == f"Live {title}"
    assert manifest.request.parameters["modes"] is None


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
    assert vehicles.attributes["dataset"] == "transit_vehicles"
    assert vehicles.attributes["provider_selection"] == "trafiklab"
    assert vehicles.attributes["providers"] == [
        {"provider": "trafiklab", "operators": ["sl"]}
    ]
    assert vehicles.attributes["modes"] == ["bus"]
    assert vehicles.attributes["fetched_parameters"] == ["bus"]
    vehicle = vehicles.vehicles()[0]
    assert vehicle.attributes["mode"] == "bus"
    assert vehicle.attributes["line"] == "1"
    assert {field.name for field in vehicle.geometry["location"].fields} == {
        "speed",
        "bearing",
    }
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


def test_shortcut_result_records_dataset_and_default_modes(monkeypatch):
    def fake_fetch(**kwargs):
        return TransportProviderResult(
            records=[
                VehicleRecord(
                    vehicle_id="vehicle-1",
                    lon=18.0,
                    lat=59.3,
                    provider="trafiklab",
                    mode="bus",
                )
            ],
            metadata={"provider": "trafiklab", "operators": ["sl"]},
        )

    module = importlib.import_module("dtcc_core.datasets.transit_vehicles")
    monkeypatch.setattr(module, "fetch_trafiklab_gtfs_vehicles", fake_fetch)

    vehicles = datasets.buses(
        bounds=(17.9, 59.2, 18.2, 59.4),
        crs="EPSG:4326",
        provider="trafiklab",
    )

    assert vehicles.attributes["dataset"] == "buses"
    assert vehicles.attributes["default_modes"] == ["bus"]
    assert vehicles.attributes["modes"] == ["bus"]


def test_shortcut_rejects_conflicting_modes():
    with pytest.raises(ValueError):
        datasets.buses(
            bounds=(17.9, 59.2, 18.2, 59.4),
            crs="EPSG:4326",
            provider="trafiklab",
            modes=("tram",),
        )


def test_invalid_provider_fails_validation():
    with pytest.raises(ValidationError, match="provider"):
        TransitVehiclesArgs(
            bounds=(17.9, 59.2, 18.2, 59.4),
            provider="not-a-provider",
        )


def test_auto_provider_unsupported_region_degrades():
    vehicles = datasets.transit_vehicles(
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        provider="auto",
    )

    assert isinstance(vehicles, VehicleCollection)
    assert len(vehicles.vehicles()) == 0
    assert vehicles.attributes["partial_result"] is True
    assert vehicles.attributes["upstream_errors"][0]["failure_class"] == (
        "unsupported_region"
    )


def test_auto_provider_unsupported_region_strict_raises():
    with pytest.raises(DatasetUpstreamError) as exc:
        datasets.transit_vehicles(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            provider="auto",
            strict_live=True,
        )

    assert exc.value.failure_class == "unsupported_region"


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


def test_missing_vasttrafik_credentials_strict_raises(monkeypatch):
    monkeypatch.delenv("VASTTRAFIK_AUTHENTICATION_KEY", raising=False)

    with pytest.raises(DatasetUpstreamError) as exc:
        datasets.buses(
            bounds=(11.9, 57.6, 12.1, 57.8),
            crs="EPSG:4326",
            provider="vasttrafik",
            strict_live=True,
        )

    assert exc.value.failure_class == "configuration"
    assert "VASTTRAFIK_AUTHENTICATION_KEY" in exc.value.message


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


def test_trafiklab_gtfs_rt_vehicle_position_payload_is_parsed():
    from google.transit import gtfs_realtime_pb2

    module = importlib.import_module("dtcc_core.datasets.transport.trafiklab_gtfs")

    feed = gtfs_realtime_pb2.FeedMessage()
    feed.header.gtfs_realtime_version = "2.0"
    feed.header.timestamp = 1_735_689_600
    entity = feed.entity.add()
    entity.id = "entity-1"
    vehicle = entity.vehicle
    vehicle.trip.route_id = "route-1"
    vehicle.trip.trip_id = "trip-1"
    vehicle.vehicle.id = "vehicle-1"
    vehicle.position.latitude = 59.3
    vehicle.position.longitude = 18.0
    vehicle.position.speed = 7.5
    vehicle.position.bearing = 90.0
    vehicle.timestamp = 1_735_689_660

    records = module._parse_vehicle_positions(
        feed.SerializeToString(),
        operator="sl",
        route_metadata={
            "route-1": {
                "route_type": "3",
                "mode": "bus",
                "line": "1",
                "route_short_name": "1",
                "route_long_name": "Centralen",
            }
        },
    )

    assert len(records) == 1
    assert records[0].provider == "trafiklab"
    assert records[0].vehicle_id == "vehicle-1"
    assert records[0].trip_id == "trip-1"
    assert records[0].mode == "bus"
    assert records[0].line == "1"
    assert records[0].timestamp == "2025-01-01T00:01:00+00:00"
    assert records[0].speed == pytest.approx(7.5)
    assert records[0].bearing == pytest.approx(90.0)


def test_trafiklab_static_route_metadata_normalizes_modes():
    module = importlib.import_module("dtcc_core.datasets.transport.trafiklab_gtfs")

    routes = module._read_routes_csv(
        "route_id,route_short_name,route_long_name,route_type,route_color,"
        "route_text_color\n"
        "bus-1,1,Centralen,3,0055aa,ffffff\n"
        "tram-1,7,Spårvagn,0,aa5500,000000\n"
    )

    assert routes["bus-1"]["mode"] == "bus"
    assert routes["bus-1"]["line"] == "1"
    assert routes["tram-1"]["mode"] == "tram"


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
            "get": staticmethod(
                lambda url, params=None, headers=None, timeout=None: Response()
            ),
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
