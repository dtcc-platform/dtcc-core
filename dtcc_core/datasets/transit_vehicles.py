# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

"""Live public transport vehicle positions.

This module exposes a provider-neutral dataset for near-real-time public
transport vehicles. It follows the same snapshot-oriented live-data contract as
the existing weather, air-quality, hydrology, and ocean datasets.

For Västtrafik, create an application in the Västtrafik developer portal,
subscribe that same application to "Planera Resa v4", and copy the
"Autentiseringsnyckel" value to the VASTTRAFIK_AUTHENTICATION_KEY environment
variable.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Literal, Optional, Tuple

import numpy as np
from pydantic import Field, field_validator

from .dataset import DatasetBaseArgs, DatasetDescriptor, DatasetUpstreamError
from .transport.base import SUPPORTED_TRANSIT_MODES, TransitMode
from .transport.trafiklab_gtfs import (
    fetch_trafiklab_gtfs_vehicles,
    select_trafiklab_operators,
)
from .transport.vasttrafik import (
    bounds_overlap_vasttrafik,
    fetch_vasttrafik_vehicles,
)
from ..common import info, warning
from ..common.progress import ProgressTracker, report_progress
from ..model.geometry import Point
from ..model.object import Object, VehicleCollection
from ..model.values import Field as DtccField
from ..reproject.reproject import reproject_array


ProviderName = Literal["auto", "trafiklab", "vasttrafik"]
GTFSFeed = Literal["regional", "sweden"]


class TransitVehiclesArgs(DatasetBaseArgs):
    """Arguments for the live public transport vehicle dataset."""

    modes: Optional[Tuple[TransitMode, ...]] = Field(
        None,
        description=(
            "Vehicle modes to include. If omitted, all modes from the selected "
            "provider(s) are returned."
        ),
    )
    provider: ProviderName = Field(
        "auto",
        description='Provider to use: "auto", "trafiklab", or "vasttrafik".',
    )
    operators: Optional[Tuple[str, ...]] = Field(
        None,
        description=(
            "Provider-specific operators to query. Trafiklab examples: sl, skane, "
            "otraf. If omitted, operators are selected from bounds."
        ),
    )
    crs: str = Field(
        "EPSG:3006",
        description="Coordinate reference system of input bounds and returned points.",
    )
    feed: GTFSFeed = Field(
        "regional",
        description='Trafiklab GTFS feed family: "regional" or "sweden".',
    )
    api_key: Optional[str] = Field(
        None,
        description=(
            "Trafiklab API key. If omitted, TRAFIKLAB_API_KEY or "
            "SAMTRAFIKEN_API_KEY is used."
        ),
    )
    vasttrafik_authentication_key: Optional[str] = Field(
        None,
        description=(
            "Västtrafik Autentiseringsnyckel / Authentication key. Create an "
            "application at https://developer.vasttrafik.se/, subscribe that "
            "same application to Planera Resa v4, then copy its "
            "Autentiseringsnyckel. If omitted, VASTTRAFIK_AUTHENTICATION_KEY "
            "is used."
        ),
    )
    timeout_s: float = Field(
        10.0,
        gt=0,
        description="Request timeout per upstream call in seconds.",
    )
    max_vehicles: Optional[int] = Field(
        None,
        gt=0,
        description="Maximum number of vehicles returned after filtering.",
    )
    include_static_metadata: bool = Field(
        True,
        description=(
            "Fetch static GTFS route metadata when needed for line names and modes."
        ),
    )
    include_unknown_modes: bool = Field(
        False,
        description=(
            "When a mode filter is set, include vehicles whose provider mode "
            "cannot be resolved."
        ),
    )
    format: Optional[Literal["pb"]] = Field(
        None,
        description='Output format ("pb" for protobuf).',
    )

    @field_validator("operators")
    @classmethod
    def normalize_operators(cls, value):
        if value is None:
            return None
        return tuple(operator.strip().lower() for operator in value if operator.strip())


class TransitVehiclesDataset(DatasetDescriptor):
    """Latest public transport vehicle positions for a bounded area.

    Västtrafik setup: create an application in the Västtrafik developer portal,
    subscribe that same application to Planera Resa v4, and set
    VASTTRAFIK_AUTHENTICATION_KEY to the application's Autentiseringsnyckel.
    """

    name = "transit_vehicles"
    title = "Transit vehicles"
    description = (
        "Fetch a live snapshot of public transport vehicle positions from "
        "regional Swedish transport providers. Returns a VehicleCollection "
        "with one point object per vehicle. For Västtrafik, create an "
        "application at https://developer.vasttrafik.se/, subscribe that same "
        "application to Planera Resa v4, and set VASTTRAFIK_AUTHENTICATION_KEY "
        "to the application's Autentiseringsnyckel."
    )
    ArgsModel = TransitVehiclesArgs
    data_category = "raw"
    result_kind = "vehicle_collection"
    python_return_type = "dtcc_core.model.VehicleCollection"
    timeout_hint = 20

    default_modes: tuple[str, ...] | None = None

    def build(self, args: TransitVehiclesArgs):
        return _build_transit_vehicles(args, self.default_modes, self.name)


class BusesDataset(TransitVehiclesDataset):
    name = "buses"
    title = "Buses"
    description = "Shortcut for transit_vehicles(..., modes=('bus',))."
    default_modes = ("bus",)


class TramsDataset(TransitVehiclesDataset):
    name = "trams"
    title = "Trams"
    description = "Shortcut for transit_vehicles(..., modes=('tram',))."
    default_modes = ("tram",)


class TrainsDataset(TransitVehiclesDataset):
    name = "trains"
    title = "Trains"
    description = "Shortcut for transit_vehicles(..., modes=('train',))."
    default_modes = ("train",)


class MetrosDataset(TransitVehiclesDataset):
    name = "metros"
    title = "Metros"
    description = "Shortcut for transit_vehicles(..., modes=('metro',))."
    default_modes = ("metro",)


class FerriesDataset(TransitVehiclesDataset):
    name = "ferries"
    title = "Ferries"
    description = "Shortcut for transit_vehicles(..., modes=('ferry',))."
    default_modes = ("ferry",)


def _build_transit_vehicles(
    args: TransitVehiclesArgs,
    default_modes: tuple[str, ...] | None,
    dataset_name: str,
):
    bounds = DatasetDescriptor.parse_bounds(args.bounds)
    bounds_tuple = (
        bounds.xmin,
        bounds.ymin,
        bounds.xmax,
        bounds.ymax,
    )
    if (
        default_modes is not None
        and args.modes is not None
        and tuple(args.modes) != tuple(default_modes)
    ):
        raise ValueError(
            f"{dataset_name} is a shortcut for modes={default_modes}; "
            "use transit_vehicles for custom mode filters."
        )
    modes = default_modes if default_modes is not None else args.modes
    if modes is not None:
        modes = tuple(modes)

    progress = ProgressTracker(
        phases={"prepare": 0.1, "fetch": 0.5, "convert": 0.3, "export": 0.1}
    )

    with progress:
        with progress.phase("prepare", "Preparing transport request..."):
            bounds_wgs84 = _transform_bounds_to_wgs84(bounds_tuple, args.crs)
            providers = _select_providers(args.provider, bounds_wgs84)
            report_progress(
                percent=100,
                message=f"Selected providers: {', '.join(providers) or 'none'}",
            )

        upstream_errors = []
        if not providers:
            exc = DatasetUpstreamError(
                dataset=dataset_name,
                operation="select_provider",
                target=str(bounds_wgs84),
                failure_class="unsupported_region",
                message=(
                    "No live public transport provider matched the requested "
                    "bounds. Use provider='vasttrafik' for Västra Götaland or "
                    "provider='trafiklab' with operators=(...) for Trafiklab "
                    "GTFS-RT regions."
                ),
            )
            upstream_errors.append(exc)
            if args.strict_live:
                raise exc
        records = []
        provider_metadata = []
        with progress.phase("fetch", "Fetching live vehicle positions..."):
            for provider in providers:
                if provider == "trafiklab":
                    result = fetch_trafiklab_gtfs_vehicles(
                        dataset_name=dataset_name,
                        bounds_wgs84=bounds_wgs84,
                        requested_modes=modes,
                        api_key=args.api_key,
                        operators=args.operators,
                        feed=args.feed,
                        timeout_s=args.timeout_s,
                        include_static_metadata=args.include_static_metadata,
                        include_unknown_modes=args.include_unknown_modes,
                        strict_live=args.strict_live,
                    )
                elif provider == "vasttrafik":
                    result = fetch_vasttrafik_vehicles(
                        dataset_name=dataset_name,
                        bounds_wgs84=bounds_wgs84,
                        requested_modes=modes,
                        authentication_key=args.vasttrafik_authentication_key,
                        timeout_s=args.timeout_s,
                        include_unknown_modes=args.include_unknown_modes,
                        strict_live=args.strict_live,
                    )
                else:
                    continue
                records.extend(result.records)
                upstream_errors.extend(result.upstream_errors)
                provider_metadata.append(result.metadata)
            report_progress(percent=100, message=f"Fetched {len(records)} vehicles")

        with progress.phase("convert", "Building vehicle collection..."):
            collection = VehicleCollection()
            collection.attributes = {
                "dtcc_type": "vehicle_collection",
                "source": "live public transport APIs",
                "provider": ",".join(providers),
                "providers": provider_metadata,
                "configuration_help": _configuration_help(providers),
                "crs": args.crs,
                "bounds": f"{bounds_tuple}",
                "bounds_wgs84": f"{bounds_wgs84}",
                "modes": list(modes) if modes is not None else list(SUPPORTED_TRANSIT_MODES),
                "retrieval_time": datetime.now(timezone.utc).isoformat(),
                "total_vehicles_found": len(records),
            }
            seen = set()
            vehicles_used = 0
            vehicles_skipped_duplicate = 0
            vehicles_skipped_no_coords = 0
            for record in records:
                dedupe_key = (record.provider, record.vehicle_id, record.trip_id)
                if dedupe_key in seen:
                    vehicles_skipped_duplicate += 1
                    continue
                seen.add(dedupe_key)

                point_coords = _transform_point_from_wgs84(
                    record.lon, record.lat, args.crs
                )
                if point_coords is None:
                    vehicles_skipped_no_coords += 1
                    continue
                x, y, z = point_coords
                if not DatasetDescriptor.point_within_bounds(x, y, bounds):
                    continue

                vehicle = _vehicle_object_from_record(record, x, y, z)
                collection.add_vehicle(vehicle)
                vehicles_used += 1
                if args.max_vehicles is not None and vehicles_used >= args.max_vehicles:
                    break

            collection.attributes.update(
                {
                    "vehicles_used": vehicles_used,
                    "vehicles_skipped_duplicate": vehicles_skipped_duplicate,
                    "vehicles_skipped_no_coords": vehicles_skipped_no_coords,
                }
            )
            DatasetDescriptor.apply_result_health_metadata(
                collection.attributes,
                upstream_errors=upstream_errors,
                requested_parameters=modes,
                fetched_parameters=sorted(
                    {
                        vehicle.attributes.get("mode")
                        for vehicle in collection.vehicles()
                        if vehicle.attributes.get("mode")
                    }
                ),
            )
            collection.calculate_bounds()
            _log_live_data_status(dataset_name, collection.attributes)
            info(f"Retrieved {vehicles_used} live transport vehicles")

        with progress.phase(
            "export",
            "Exporting to protobuf..." if args.format == "pb" else "Preparing result...",
        ):
            if args.format == "pb":
                return collection.to_proto().SerializeToString()
            return collection


def _select_providers(
    provider: ProviderName,
    bounds_wgs84: tuple[float, float, float, float],
) -> list[str]:
    if provider != "auto":
        return [provider]

    if bounds_overlap_vasttrafik(bounds_wgs84):
        return ["vasttrafik"]
    providers = []
    if select_trafiklab_operators(bounds_wgs84):
        providers.append("trafiklab")
    return providers


def _configuration_help(providers: list[str]) -> list[str]:
    help_lines = []
    if "vasttrafik" in providers:
        help_lines.append(
            "Västtrafik: set VASTTRAFIK_AUTHENTICATION_KEY from "
            "'Autentiseringsnyckel'. In the Västtrafik developer portal, "
            "create/select an application, subscribe that same application to "
            "Planera Resa v4, then copy its Autentiseringsnyckel. You can also "
            "pass vasttrafik_authentication_key=<key>."
        )
    if "trafiklab" in providers:
        help_lines.append(
            "Trafiklab: set TRAFIKLAB_API_KEY, or pass api_key=<api_key>. "
            "Create an API key at https://www.trafiklab.se/."
        )
    if help_lines:
        help_lines.append(
            "Do not commit keys to source code; use environment variables, "
            "a local .env loader, or a deployment secret store."
        )
    return help_lines


def _log_live_data_status(dataset_name: str, attributes: dict[str, Any]) -> None:
    errors = attributes.get("upstream_errors") or []
    if not errors:
        return

    first = errors[0]
    warning(
        f"{dataset_name}: returned a partial/empty live result "
        f"({len(errors)} upstream issue(s)). {first.get('message', '')}"
    )


def _transform_bounds_to_wgs84(
    bounds: tuple[float, float, float, float], crs: str
) -> tuple[float, float, float, float]:
    if crs.upper() in ("CRS84", "EPSG:4326", "WGS84"):
        return bounds
    xmin, ymin, xmax, ymax = bounds
    corners = np.array(
        [[xmin, ymin, 0.0], [xmin, ymax, 0.0], [xmax, ymin, 0.0], [xmax, ymax, 0.0]],
        dtype=float,
    )
    transformed = reproject_array(corners, crs, "EPSG:4326")
    xs = transformed[:, 0]
    ys = transformed[:, 1]
    return (float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max()))


def _transform_point_from_wgs84(lon: float, lat: float, crs: str):
    if crs.upper() in ("CRS84", "EPSG:4326", "WGS84"):
        return float(lon), float(lat), 0.0
    transformed = reproject_array(np.array([[lon, lat, 0.0]], dtype=float), "EPSG:4326", crs)
    return float(transformed[0, 0]), float(transformed[0, 1]), float(transformed[0, 2])


def _vehicle_object_from_record(record, x: float, y: float, z: float) -> Object:
    vehicle = Object()
    attrs: dict[str, Any] = {
        "vehicle_id": record.vehicle_id,
        "provider": record.provider,
        "mode": record.mode or "unknown",
        "timestamp": record.timestamp,
        "route_id": record.route_id,
        "trip_id": record.trip_id,
        "line": record.line,
        "destination": record.destination,
        "bearing": record.bearing,
        "speed": record.speed,
    }
    attrs.update(record.attributes)
    vehicle.attributes = {k: v for k, v in attrs.items() if v not in (None, "")}

    point = Point(x=x, y=y, z=z)
    fields = []
    for name, value, unit in (
        ("speed", record.speed, "m/s"),
        ("bearing", record.bearing, "degrees"),
    ):
        if value is None:
            continue
        field = DtccField()
        field.name = name
        field.unit = unit
        field.dim = 1
        field.values = np.array([float(value)], dtype=np.float32)
        fields.append(field)
    point.fields = fields
    vehicle.geometry["location"] = point
    return vehicle


__all__ = [
    "TransitVehiclesArgs",
    "TransitVehiclesDataset",
    "BusesDataset",
    "TramsDataset",
    "TrainsDataset",
    "MetrosDataset",
    "FerriesDataset",
]
