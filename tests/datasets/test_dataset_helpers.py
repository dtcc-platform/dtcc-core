"""Tests for shared dataset helper modules."""

from __future__ import annotations

import pytest

from dtcc_core.datasets.providers import (
    provider_display_name,
    provider_entry,
    provider_info,
    provider_slug,
)


def test_provider_helpers_normalize_aliases_display_names_and_slugs():
    assert provider_slug("Lantmäteriet") == "lantmateriet"
    assert provider_slug("LM") == "lantmateriet"
    assert provider_slug("Open Street Map") == "openstreetmap"
    assert provider_display_name("vasttrafik") == "Västtrafik"

    smhi = provider_info("SMHI")
    assert smhi.slug == "smhi"
    assert smhi.display_name == "SMHI"

    assert provider_entry("västtrafik") == {
        "name": "Västtrafik",
        "slug": "vasttrafik",
        "role": "source_provider",
    }


def test_provider_helpers_fail_loudly_for_unknown_provider():
    with pytest.raises(ValueError, match="Unknown dataset provider"):
        provider_entry("not-a-provider")
