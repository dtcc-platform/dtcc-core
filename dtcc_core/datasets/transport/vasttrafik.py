# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from __future__ import annotations

import os
from datetime import datetime, timezone
from typing import Any

from ..dataset import DatasetDescriptor, DatasetUpstreamError
from .base import (
    TransportProviderResult,
    VehicleRecord,
    make_config_error,
    normalize_mode,
)


VASTTRAFIK_BBOX = (10.7, 57.0, 14.9, 59.5)
TOKEN_URL = "https://ext-api.vasttrafik.se/token"
POSITIONS_URL = "https://ext-api.vasttrafik.se/pr/v4/positions"
VASTTRAFIK_POSITION_MODES = ("tram", "bus", "ferry", "train")


def fetch_vasttrafik_vehicles(
    *,
    dataset_name: str,
    bounds_wgs84: tuple[float, float, float, float],
    requested_modes: tuple[str, ...] | None,
    authentication_key: str | None,
    timeout_s: float,
    include_unknown_modes: bool,
    strict_live: bool,
) -> TransportProviderResult:
    """Fetch vehicle positions from Västtrafik's Planera Resa v4 API.

    The endpoint is kept behind a provider adapter because Västtrafik is not
    currently exposed as a GTFS-RT vehicle-position feed through Trafiklab.
    The application whose Autentiseringsnyckel is used must be subscribed to
    Planera Resa v4 in the Västtrafik developer portal.
    """
    result = TransportProviderResult(
        metadata={
            "provider": "vasttrafik",
            "operators": ["vt"],
        }
    )

    resolved_authentication_key = _resolve_authentication_key(authentication_key)
    if not resolved_authentication_key:
        exc = make_config_error(
            dataset_name,
            "vasttrafik",
            (
                "Västtrafik credentials missing. Set "
                "VASTTRAFIK_AUTHENTICATION_KEY from the portal field "
                "'Autentiseringsnyckel', or pass "
                "vasttrafik_authentication_key=<key>. In the Västtrafik "
                "developer portal, create/select an application, subscribe "
                "that same application to Planera Resa v4, then copy its "
                "Autentiseringsnyckel. "
                "Do not commit credentials."
            ),
        )
        result.upstream_errors.append(exc)
        if strict_live:
            raise exc
        return result

    try:
        import requests
    except ImportError as exc:
        err = make_config_error(
            dataset_name,
            "vasttrafik",
            "requests library required for Västtrafik transport data.",
        )
        result.upstream_errors.append(err)
        if strict_live:
            raise err from exc
        return result

    try:
        token = _fetch_token(
            requests,
            authentication_key=resolved_authentication_key,
            timeout_s=timeout_s,
        )
        data = _fetch_positions(
            requests,
            token=token,
            bounds_wgs84=bounds_wgs84,
            requested_modes=requested_modes,
            timeout_s=timeout_s,
        )
        records = _parse_positions(data)
    except Exception as exc:
        if hasattr(exc, "failure_class"):
            result.upstream_errors.append(exc)
            if strict_live:
                raise
            return result
        err = make_config_error(
            dataset_name,
            "vasttrafik",
            f"Västtrafik positions failed: {exc}",
        )
        result.upstream_errors.append(err)
        if strict_live:
            raise err from exc
        return result

    for record in records:
        if requested_modes is not None:
            if record.mode in requested_modes:
                pass
            elif record.mode is None and include_unknown_modes:
                pass
            else:
                continue
        result.records.append(record)
    return result


def bounds_overlap_vasttrafik(bounds_wgs84: tuple[float, float, float, float]) -> bool:
    return not (
        bounds_wgs84[2] < VASTTRAFIK_BBOX[0]
        or VASTTRAFIK_BBOX[2] < bounds_wgs84[0]
        or bounds_wgs84[3] < VASTTRAFIK_BBOX[1]
        or VASTTRAFIK_BBOX[3] < bounds_wgs84[1]
    )


def _fetch_token(
    requests,
    *,
    authentication_key: str,
    timeout_s: float,
):
    try:
        response = requests.post(
            TOKEN_URL,
            data={"grant_type": "client_credentials"},
            headers={"Authorization": f"Basic {authentication_key}"},
            timeout=timeout_s,
        )
        response.raise_for_status()
        payload = response.json()
        token = payload.get("access_token")
        if not token:
            raise ValueError("token response did not contain access_token")
        return token
    except requests.RequestException as exc:
        raise DatasetDescriptor.build_upstream_error(
            "transit_vehicles", "fetch_vasttrafik_token", TOKEN_URL, exc
        ) from exc


def _resolve_authentication_key(explicit_key: str | None) -> str | None:
    value = explicit_key or os.environ.get("VASTTRAFIK_AUTHENTICATION_KEY")
    if not value:
        return None
    value = value.strip()
    if value.lower().startswith("authorization:"):
        value = value.split(":", 1)[1].strip()
    if value.lower().startswith("basic "):
        return value[6:].strip()
    return value


def _fetch_positions(
    requests,
    *,
    token: str,
    bounds_wgs84: tuple[float, float, float, float],
    requested_modes: tuple[str, ...] | None,
    timeout_s: float,
):
    xmin, ymin, xmax, ymax = bounds_wgs84
    transport_modes = _vasttrafik_position_modes(requested_modes)
    if not transport_modes:
        return []

    params = {
        "lowerLeftLat": ymin,
        "lowerLeftLong": xmin,
        "upperRightLat": ymax,
        "upperRightLong": xmax,
        "transportModes": list(transport_modes),
        "limit": 200,
    }
    try:
        response = requests.get(
            POSITIONS_URL,
            params=params,
            headers={"Authorization": f"Bearer {token}"},
            timeout=timeout_s,
        )
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        response = getattr(exc, "response", None)
        if getattr(response, "status_code", None) == 403:
            raise DatasetUpstreamError(
                dataset="transit_vehicles",
                operation="fetch_vasttrafik_positions",
                target=POSITIONS_URL,
                failure_class="http_4xx",
                status_code=403,
                message=(
                    "Västtrafik returned 403 Forbidden for Planera Resa v4 "
                    "positions. The authentication key is valid, but the "
                    "application is not authorized for this endpoint. Subscribe "
                    "the same application/key to Planera Resa v4 in the "
                    "Västtrafik developer portal."
                ),
            ) from exc
        raise DatasetDescriptor.build_upstream_error(
            "transit_vehicles", "fetch_vasttrafik_positions", POSITIONS_URL, exc
        ) from exc


def _parse_positions(data: Any) -> list[VehicleRecord]:
    records = []
    for item in _find_vehicle_like_dicts(data):
        lon = _coordinate_value(item, ("lon", "longitude", "x"))
        lat = _coordinate_value(item, ("lat", "latitude", "y"))
        if lon is None or lat is None:
            continue
        if abs(lon) > 180:
            lon = lon / 1_000_000.0
        if abs(lat) > 90:
            lat = lat / 1_000_000.0

        line_details = item.get("line") if isinstance(item.get("line"), dict) else {}
        mode = normalize_mode(
            item.get("type")
            or item.get("transportMode")
            or item.get("journeyType")
            or item.get("vehicleType")
            or line_details.get("transportMode")
        )
        timestamp = item.get("timestamp") or item.get("time") or item.get("created")
        if timestamp is None:
            timestamp = datetime.now(timezone.utc).isoformat()

        details_reference = item.get("detailsReference")
        vehicle_id = (
            item.get("vehicle_id")
            or item.get("vehicleId")
            or item.get("id")
            or item.get("gid")
            or item.get("journeyid")
            or item.get("journeyId")
            or details_reference
            or f"vasttrafik:{lon:.6f}:{lat:.6f}"
        )
        route_id = (
            item.get("route_id")
            or item.get("routeId")
            or item.get("lineid")
            or line_details.get("name")
        )
        line = (
            line_details.get("name")
            or line_details.get("detailName")
            or item.get("name")
            or item.get("sname")
            or route_id
        )
        direction_details = item.get("directionDetails")
        if isinstance(direction_details, dict):
            direction = (
                direction_details.get("fullDirection")
                or direction_details.get("shortDirection")
                or item.get("direction")
            )
        else:
            direction = item.get("direction")
        records.append(
            VehicleRecord(
                vehicle_id=str(vehicle_id),
                lon=float(lon),
                lat=float(lat),
                provider="vasttrafik",
                mode=mode,
                timestamp=str(timestamp),
                route_id=None if route_id is None else str(route_id),
                trip_id=item.get("trip_id") or item.get("tripId"),
                line=None if line is None else str(line),
                destination=item.get("destination") or direction,
                bearing=_float_or_none(item.get("bearing") or item.get("directionAngle")),
                speed=_float_or_none(item.get("speed")),
                attributes={
                    "provider": "vasttrafik",
                    "operator": "vt",
                    "raw_type": item.get("type") or item.get("journeyType"),
                    "details_reference": details_reference,
                },
            )
        )
    return records


def _vasttrafik_position_modes(requested_modes: tuple[str, ...] | None) -> tuple[str, ...]:
    if requested_modes is None:
        return VASTTRAFIK_POSITION_MODES
    return tuple(mode for mode in requested_modes if mode in VASTTRAFIK_POSITION_MODES)


def _find_vehicle_like_dicts(data: Any):
    if isinstance(data, dict):
        if any(key in data for key in ("x", "lon", "longitude")) and any(
            key in data for key in ("y", "lat", "latitude")
        ):
            yield data
        for value in data.values():
            yield from _find_vehicle_like_dicts(value)
    elif isinstance(data, list):
        for value in data:
            yield from _find_vehicle_like_dicts(value)


def _coordinate_value(data: dict[str, Any], keys: tuple[str, ...]) -> float | None:
    for key in keys:
        if key in data:
            return _float_or_none(data.get(key))
    return None


def _float_or_none(value) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
