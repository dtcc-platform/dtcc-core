"""Pytest configuration for live dataset tests."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


LIVE_DATASET_TESTS_ENV = "DTCC_LIVE_DATASET_TESTS"


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add an opt-in switch for live dataset tests."""
    parser.addoption(
        "--run-live",
        action="store_true",
        default=False,
        help=(
            "Run tests marked as live. Also requires "
            f"{LIVE_DATASET_TESTS_ENV}=1."
        ),
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
    live_items = [item for item in items if item.get_closest_marker("live")]
    if not live_items:
        return

    live_requested = (
        config.getoption("--run-live")
        or "live" in markexpr
        or _live_path_requested(config.args)
    )
    if live_requested:
        if os.environ.get(LIVE_DATASET_TESTS_ENV) == "1":
            return
        skip_live = pytest.mark.skip(
            reason=(
                f"Set {LIVE_DATASET_TESTS_ENV}=1 to run live provider tests."
            )
        )
        for item in live_items:
            item.add_marker(skip_live)
        return

    keep_items = [item for item in items if item not in live_items]
    config.hook.pytest_deselected(items=live_items)
    items[:] = keep_items


def _live_path_requested(args: list[str]) -> bool:
    for arg in args:
        parts = Path(str(arg)).parts
        if len(parts) >= 3 and parts[-3:] == ("tests", "datasets", "live"):
            return True
        if "tests/datasets/live" in str(arg).replace("\\", "/"):
            return True
    return False
