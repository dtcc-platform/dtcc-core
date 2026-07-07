"""Provider naming helpers for Dataset v2 metadata."""

from __future__ import annotations

from dataclasses import dataclass
import re
import unicodedata


@dataclass(frozen=True)
class ProviderInfo:
    """Canonical provider identity used in dataset metadata."""

    slug: str
    display_name: str


_PROVIDERS: dict[str, ProviderInfo] = {
    "dtcc-platform": ProviderInfo("dtcc-platform", "DTCC Platform"),
    "lantmateriet": ProviderInfo("lantmateriet", "Lantmäteriet"),
    "openstreetmap": ProviderInfo("openstreetmap", "OpenStreetMap"),
    "scb": ProviderInfo("scb", "SCB"),
    "smhi": ProviderInfo("smhi", "SMHI"),
    "trafiklab": ProviderInfo("trafiklab", "Trafiklab"),
    "vasttrafik": ProviderInfo("vasttrafik", "Västtrafik"),
}

_ALIASES: dict[str, str] = {
    "dtcc": "dtcc-platform",
    "dtcc-platform": "dtcc-platform",
    "lantmateriet": "lantmateriet",
    "lantmäteriet": "lantmateriet",
    "lm": "lantmateriet",
    "osm": "openstreetmap",
    "openstreetmap": "openstreetmap",
    "open-street-map": "openstreetmap",
    "scb": "scb",
    "smhi": "smhi",
    "trafiklab": "trafiklab",
    "vasttrafik": "vasttrafik",
    "västtrafik": "vasttrafik",
}


def provider_info(provider: str) -> ProviderInfo:
    """Return canonical provider information for ``provider``."""

    slug = provider_slug(provider)
    try:
        return _PROVIDERS[slug]
    except KeyError as exc:
        raise ValueError(f"Unknown dataset provider {provider!r}.") from exc


def provider_slug(provider: str) -> str:
    """Return the canonical ASCII slug for a provider name or alias."""

    normalized = _normalize_provider_key(provider)
    if not normalized:
        raise ValueError("provider must be a non-empty string.")
    return _ALIASES.get(normalized, normalized)


def provider_display_name(provider: str) -> str:
    """Return the canonical human-facing provider display name."""

    return provider_info(provider).display_name


def provider_entry(provider: str, *, role: str = "source_provider") -> dict[str, str]:
    """Return a normalized Dataset v2 provider metadata item."""

    info = provider_info(provider)
    return {"name": info.display_name, "slug": info.slug, "role": role}


def _normalize_provider_key(provider: str) -> str:
    if not isinstance(provider, str):
        raise ValueError("provider must be a non-empty string.")
    text = provider.strip().lower()
    if not text:
        return ""
    ascii_text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    collapsed = re.sub(r"[\s_/-]+", " ", ascii_text).strip()
    return collapsed.replace(" ", "-")


__all__ = [
    "ProviderInfo",
    "provider_display_name",
    "provider_entry",
    "provider_info",
    "provider_slug",
]
