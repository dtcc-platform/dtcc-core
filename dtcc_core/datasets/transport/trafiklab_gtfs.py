# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from __future__ import annotations

import csv
import io
import os
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from platformdirs import user_cache_dir

from ..dataset import DatasetDescriptor, DatasetUpstreamError
from .base import (
    TransportProviderResult,
    VehicleRecord,
    make_config_error,
    normalize_gtfs_route_type,
)


REGIONAL_RT_URL = (
    "https://opendata.samtrafiken.se/gtfs-rt/{operator}/VehiclePositions.pb"
)
REGIONAL_STATIC_URL = "https://opendata.samtrafiken.se/gtfs/{operator}/{operator}.zip"
SWEDEN_RT_URL = (
    "https://opendata.samtrafiken.se/gtfs-rt-sweden/"
    "{operator}/VehiclePositionsSweden.pb"
)
SWEDEN_STATIC_URL = "https://opendata.samtrafiken.se/gtfs-sweden/sweden.zip"


@dataclass(frozen=True)
class TrafiklabOperator:
    code: str
    name: str
    bbox: tuple[float, float, float, float]


# Coarse WGS84 service boxes for the Trafiklab operators with vehicle positions.
# The boxes are intentionally generous: they are provider selectors, not filters.
TRAFIKLAB_OPERATORS: tuple[TrafiklabOperator, ...] = (
    TrafiklabOperator("sl", "SL", (17.2, 58.7, 18.9, 60.2)),
    TrafiklabOperator("ul", "UL", (16.6, 59.4, 18.8, 60.8)),
    TrafiklabOperator("otraf", "Östgötatrafiken", (14.4, 57.6, 16.9, 59.1)),
    TrafiklabOperator("jlt", "JLT", (13.0, 56.8, 15.8, 58.4)),
    TrafiklabOperator("krono", "Kronoberg", (13.2, 56.1, 15.6, 57.4)),
    TrafiklabOperator("klt", "KLT", (15.1, 56.0, 17.8, 58.3)),
    TrafiklabOperator("gotland", "Gotland", (18.0, 56.7, 19.4, 58.0)),
    TrafiklabOperator("blekinge", "Blekingetrafiken", (14.1, 55.8, 16.1, 56.6)),
    TrafiklabOperator("skane", "Skånetrafiken", (12.4, 55.2, 14.8, 56.6)),
    TrafiklabOperator("halland", "Hallandstrafiken", (11.7, 56.2, 13.5, 57.8)),
    TrafiklabOperator("varm", "Värmlandstrafik", (11.7, 58.6, 14.7, 61.2)),
    TrafiklabOperator("orebro", "Örebro", (14.2, 58.6, 16.1, 60.2)),
    TrafiklabOperator("vastmanland", "Västmanland", (15.4, 59.2, 17.6, 60.4)),
    TrafiklabOperator("dt", "Dalatrafik", (12.0, 59.6, 16.5, 62.4)),
    TrafiklabOperator("xt", "X-trafik", (15.4, 60.1, 17.8, 62.4)),
    TrafiklabOperator("dintur", "Din Tur", (14.4, 61.8, 19.2, 64.2)),
)


def fetch_trafiklab_gtfs_vehicles(
    *,
    dataset_name: str,
    bounds_wgs84: tuple[float, float, float, float],
    requested_modes: tuple[str, ...] | None,
    api_key: str | None,
    operators: tuple[str, ...] | None,
    feed: str,
    timeout_s: float,
    include_static_metadata: bool,
    include_unknown_modes: bool,
    strict_live: bool,
) -> TransportProviderResult:
    """Fetch vehicle positions from Trafiklab GTFS-RT feeds."""
    result = TransportProviderResult(
        metadata={
            "provider": "trafiklab",
            "feed": feed,
            "operators": [],
        }
    )

    resolved_key = _resolve_api_key(api_key)
    if not resolved_key:
        exc = make_config_error(
            dataset_name,
            "trafiklab",
            (
                "Trafiklab API key missing. Set TRAFIKLAB_API_KEY, or pass "
                "api_key=<api_key>. Create a Trafiklab account/key at "
                "https://www.trafiklab.se/. Do not commit API keys."
            ),
        )
        result.upstream_errors.append(exc)
        if strict_live:
            raise exc
        return result

    selected = (
        [operator.lower() for operator in operators]
        if operators
        else select_trafiklab_operators(bounds_wgs84)
    )
    result.metadata["operators"] = selected

    if not selected:
        return result

    try:
        import requests
    except ImportError as exc:
        err = make_config_error(
            dataset_name,
            "trafiklab",
            "requests library required for Trafiklab transport data.",
        )
        result.upstream_errors.append(err)
        if strict_live:
            raise err from exc
        return result

    for operator in selected:
        route_metadata = {}
        if include_static_metadata or requested_modes is not None:
            try:
                route_metadata = _fetch_route_metadata(
                    requests,
                    operator=operator,
                    feed=feed,
                    api_key=resolved_key,
                    timeout_s=timeout_s,
                )
            except DatasetUpstreamError as exc:
                result.upstream_errors.append(exc)
                if strict_live:
                    raise

        try:
            payload = _fetch_realtime_payload(
                requests,
                operator=operator,
                feed=feed,
                api_key=resolved_key,
                timeout_s=timeout_s,
            )
            records = _parse_vehicle_positions(
                payload,
                operator=operator,
                route_metadata=route_metadata,
            )
        except DatasetUpstreamError as exc:
            result.upstream_errors.append(exc)
            if strict_live:
                raise
            continue

        for record in records:
            if not _point_in_wgs84_bounds(record.lon, record.lat, bounds_wgs84):
                continue
            if requested_modes is not None:
                if record.mode in requested_modes:
                    pass
                elif record.mode is None and include_unknown_modes:
                    pass
                else:
                    continue
            result.records.append(record)

    return result


def select_trafiklab_operators(
    bounds_wgs84: tuple[float, float, float, float]
) -> list[str]:
    """Select Trafiklab operators whose coarse service box overlaps bounds."""
    return [
        operator.code
        for operator in TRAFIKLAB_OPERATORS
        if _bboxes_overlap(bounds_wgs84, operator.bbox)
    ]


def _resolve_api_key(api_key: str | None) -> str | None:
    if api_key:
        return api_key
    return os.environ.get("TRAFIKLAB_API_KEY") or os.environ.get("SAMTRAFIKEN_API_KEY")


def _fetch_realtime_payload(
    requests,
    *,
    operator: str,
    feed: str,
    api_key: str,
    timeout_s: float,
) -> bytes:
    url = _realtime_url(operator, feed)
    try:
        response = requests.get(url, params={"key": api_key}, timeout=timeout_s)
        response.raise_for_status()
        return response.content
    except requests.RequestException as exc:
        raise DatasetDescriptor.build_upstream_error(
            "transit_vehicles", "fetch_gtfs_realtime", url, exc
        ) from exc


def _fetch_route_metadata(
    requests,
    *,
    operator: str,
    feed: str,
    api_key: str,
    timeout_s: float,
) -> dict[str, dict[str, Any]]:
    cached = _cached_routes_path(operator, feed)
    if cached.exists():
        return _read_routes_csv(cached.read_text(encoding="utf-8"))

    url = _static_url(operator, feed)
    try:
        response = requests.get(url, params={"key": api_key}, timeout=timeout_s)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise DatasetDescriptor.build_upstream_error(
            "transit_vehicles", "fetch_gtfs_static", url, exc
        ) from exc

    try:
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            routes_txt = archive.read("routes.txt").decode("utf-8-sig")
    except (KeyError, zipfile.BadZipFile, UnicodeDecodeError) as exc:
        raise DatasetUpstreamError(
            dataset="transit_vehicles",
            operation="decode_gtfs_static",
            target=url,
            failure_class="invalid_gtfs",
            message=f"transit_vehicles decode_gtfs_static failed for {url}: {exc}",
        ) from exc

    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_text(routes_txt, encoding="utf-8")
    return _read_routes_csv(routes_txt)


def _read_routes_csv(text: str) -> dict[str, dict[str, Any]]:
    reader = csv.DictReader(io.StringIO(text))
    routes: dict[str, dict[str, Any]] = {}
    for row in reader:
        route_id = row.get("route_id")
        if not route_id:
            continue
        route_type = row.get("route_type")
        routes[route_id] = {
            "route_type": route_type,
            "mode": normalize_gtfs_route_type(route_type),
            "line": row.get("route_short_name") or row.get("route_long_name"),
            "route_short_name": row.get("route_short_name"),
            "route_long_name": row.get("route_long_name"),
            "route_color": row.get("route_color"),
            "route_text_color": row.get("route_text_color"),
        }
    return routes


def _parse_vehicle_positions(
    payload: bytes,
    *,
    operator: str,
    route_metadata: dict[str, dict[str, Any]],
) -> list[VehicleRecord]:
    try:
        from google.transit import gtfs_realtime_pb2
    except ImportError as exc:
        raise DatasetUpstreamError(
            dataset="transit_vehicles",
            operation="import_gtfs_realtime",
            target="gtfs-realtime-bindings",
            failure_class="missing_dependency",
            message=(
                "gtfs-realtime-bindings is required to parse Trafiklab "
                "GTFS-RT feeds. Install gtfs-realtime-bindings."
            ),
        ) from exc

    feed = gtfs_realtime_pb2.FeedMessage()
    try:
        feed.ParseFromString(payload)
    except Exception as exc:
        raise DatasetUpstreamError(
            dataset="transit_vehicles",
            operation="decode_gtfs_realtime",
            target=operator,
            failure_class="invalid_protobuf",
            message=f"transit_vehicles decode_gtfs_realtime failed for {operator}: {exc}",
        ) from exc

    records: list[VehicleRecord] = []
    feed_timestamp = _posix_to_iso(getattr(feed.header, "timestamp", None))
    for entity in feed.entity:
        if not _has_field(entity, "vehicle"):
            continue
        vehicle_position = entity.vehicle
        if not _has_field(vehicle_position, "position"):
            continue

        position = vehicle_position.position
        route_id = _message_value(vehicle_position.trip, "route_id")
        trip_id = _message_value(vehicle_position.trip, "trip_id")
        vehicle_id = (
            _message_value(vehicle_position.vehicle, "id")
            or _message_value(vehicle_position.vehicle, "label")
            or entity.id
        )
        route_info = route_metadata.get(route_id or "", {})
        timestamp = _posix_to_iso(_message_value(vehicle_position, "timestamp"))
        if timestamp is None:
            timestamp = feed_timestamp

        attrs = {
            "provider": "trafiklab",
            "operator": operator,
            "feed_entity_id": entity.id,
            "gtfs_route_type": route_info.get("route_type"),
            "route_short_name": route_info.get("route_short_name"),
            "route_long_name": route_info.get("route_long_name"),
            "route_color": route_info.get("route_color"),
            "route_text_color": route_info.get("route_text_color"),
            "current_stop_sequence": _message_value(
                vehicle_position, "current_stop_sequence"
            ),
            "stop_id": _message_value(vehicle_position, "stop_id"),
            "current_status": _enum_name(vehicle_position, "current_status"),
            "occupancy_status": _enum_name(vehicle_position, "occupancy_status"),
            "congestion_level": _enum_name(vehicle_position, "congestion_level"),
        }

        records.append(
            VehicleRecord(
                vehicle_id=str(vehicle_id),
                lon=float(position.longitude),
                lat=float(position.latitude),
                provider="trafiklab",
                mode=route_info.get("mode"),
                timestamp=timestamp,
                route_id=route_id,
                trip_id=trip_id,
                line=route_info.get("line") or route_id,
                bearing=_message_float(position, "bearing"),
                speed=_message_float(position, "speed"),
                attributes={k: v for k, v in attrs.items() if v not in (None, "")},
            )
        )
    return records


def _cached_routes_path(operator: str, feed: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    return (
        Path(user_cache_dir("dtcc-core", "dtcc"))
        / "transport"
        / "trafiklab"
        / feed
        / operator
        / f"routes_{stamp}.txt"
    )


def _realtime_url(operator: str, feed: str) -> str:
    if feed == "sweden":
        return SWEDEN_RT_URL.format(operator=operator)
    return REGIONAL_RT_URL.format(operator=operator)


def _static_url(operator: str, feed: str) -> str:
    if feed == "sweden":
        return SWEDEN_STATIC_URL
    return REGIONAL_STATIC_URL.format(operator=operator)


def _point_in_wgs84_bounds(
    lon: float, lat: float, bounds: tuple[float, float, float, float]
) -> bool:
    xmin, ymin, xmax, ymax = bounds
    return xmin <= lon <= xmax and ymin <= lat <= ymax


def _bboxes_overlap(
    a: tuple[float, float, float, float], b: tuple[float, float, float, float]
) -> bool:
    return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])


def _has_field(message, field_name: str) -> bool:
    try:
        return message.HasField(field_name)
    except (AttributeError, ValueError):
        return bool(getattr(message, field_name, None))


def _message_value(message, field_name: str):
    if message is None:
        return None
    try:
        value = getattr(message, field_name)
    except AttributeError:
        return None
    if value == "":
        return None
    return value


def _message_float(message, field_name: str) -> float | None:
    if not _has_field(message, field_name):
        return None
    value = _message_value(message, field_name)
    if value is None:
        return None
    return float(value)


def _enum_name(message, field_name: str) -> str | None:
    if not _has_field(message, field_name):
        return None
    try:
        descriptor = message.DESCRIPTOR.fields_by_name[field_name]
        return descriptor.enum_type.values_by_number[int(getattr(message, field_name))].name
    except Exception:
        value = _message_value(message, field_name)
        return None if value is None else str(value)


def _posix_to_iso(value) -> str | None:
    if value in (None, "", 0):
        return None
    try:
        return datetime.fromtimestamp(float(value), tz=timezone.utc).isoformat()
    except (TypeError, ValueError, OSError):
        return None
