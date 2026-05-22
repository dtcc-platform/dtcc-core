# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from ..dataset import DatasetUpstreamError


TransitMode = Literal["bus", "tram", "train", "metro", "ferry"]
SUPPORTED_TRANSIT_MODES: tuple[str, ...] = ("bus", "tram", "train", "metro", "ferry")


@dataclass(frozen=True)
class VehicleRecord:
    """Provider-neutral live vehicle position."""

    vehicle_id: str
    lon: float
    lat: float
    provider: str
    mode: str | None = None
    timestamp: str | None = None
    route_id: str | None = None
    trip_id: str | None = None
    line: str | None = None
    destination: str | None = None
    bearing: float | None = None
    speed: float | None = None
    attributes: dict[str, Any] = field(default_factory=dict)


@dataclass
class TransportProviderResult:
    """Result from one or more transport providers."""

    records: list[VehicleRecord] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    upstream_errors: list[DatasetUpstreamError] = field(default_factory=list)


def normalize_mode(value: Any) -> str | None:
    """Normalize common provider mode names to DTCC's mode vocabulary."""
    if value is None:
        return None
    text = str(value).strip().lower()
    if not text:
        return None

    aliases = {
        "bus": "bus",
        "buss": "bus",
        "coach": "bus",
        "express_bus": "bus",
        "rail_replacement_bus": "bus",
        "tram": "tram",
        "spårvagn": "tram",
        "sparvagn": "tram",
        "light_rail": "tram",
        "train": "train",
        "rail": "train",
        "regional_train": "train",
        "commuter_train": "train",
        "metro": "metro",
        "subway": "metro",
        "underground": "metro",
        "ferry": "ferry",
        "boat": "ferry",
        "water": "ferry",
        "ship": "ferry",
    }
    if text in aliases:
        return aliases[text]

    if "bus" in text or "buss" in text:
        return "bus"
    if "tram" in text or "spår" in text or "spar" in text:
        return "tram"
    if "metro" in text or "subway" in text:
        return "metro"
    if "train" in text or "rail" in text or "tåg" in text or "tag" in text:
        return "train"
    if "ferry" in text or "boat" in text or "ship" in text:
        return "ferry"
    return None


def normalize_gtfs_route_type(route_type: Any) -> str | None:
    """Normalize standard and extended GTFS route types."""
    if route_type is None or route_type == "":
        return None
    try:
        value = int(route_type)
    except (TypeError, ValueError):
        return normalize_mode(route_type)

    if value == 0:
        return "tram"
    if value == 1:
        return "metro"
    if value == 2:
        return "train"
    if value == 3:
        return "bus"
    if value == 4:
        return "ferry"

    # Extended GTFS route types used by Trafiklab and many European feeds.
    if 100 <= value < 200:
        return "train"
    if 400 <= value < 500:
        return "metro"
    if 700 <= value < 800:
        return "bus"
    if 900 <= value < 1000:
        return "tram"
    if 1000 <= value < 1100:
        return "ferry"
    return None


def make_config_error(dataset: str, provider: str, message: str) -> DatasetUpstreamError:
    """Create a dataset error for missing provider configuration."""
    return DatasetUpstreamError(
        dataset=dataset,
        operation="configure_provider",
        target=provider,
        failure_class="configuration",
        message=message,
    )
