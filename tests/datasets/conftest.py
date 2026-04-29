"""Pytest configuration for live dataset tests."""

from __future__ import annotations

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add an opt-in switch for live dataset tests."""
    parser.addoption(
        "--run-live",
        action="store_true",
        default=False,
        help="Run only tests marked as live (equivalent to -m live).",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register the live marker and make --run-live behave like -m live."""
    config.addinivalue_line(
        "markers",
        "live: opt-in tests that hit live upstream dataset endpoints",
    )

    if config.getoption("--run-live") and not config.option.markexpr:
        config.option.markexpr = "live"


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Deselect live tests unless explicitly requested."""
    markexpr = (config.option.markexpr or "").strip()
    if config.getoption("--run-live") or "live" in markexpr:
        return

    live_items = [item for item in items if item.get_closest_marker("live")]
    if not live_items:
        return

    keep_items = [item for item in items if item not in live_items]
    config.hook.pytest_deselected(items=live_items)
    items[:] = keep_items
