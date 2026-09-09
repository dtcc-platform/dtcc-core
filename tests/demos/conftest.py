"""Pytest configuration for demo execution tests.

Full demo execution (see ``test_demo_execution.py``) downloads data from the
DTCC data service and runs meshing.  Those tests are marked ``demo`` and are
deselected by default so normal CI stays fast and offline.  Run them with::

    DTCC_RUN_DEMOS=1 pytest tests/demos/test_demo_execution.py

or::

    pytest tests/demos --run-demos
"""

from __future__ import annotations

import os

import pytest


RUN_DEMOS_ENV = "DTCC_RUN_DEMOS"


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add an opt-in switch for full demo execution tests."""
    parser.addoption(
        "--run-demos",
        action="store_true",
        default=False,
        help=(
            "Run tests marked as demo (full demo scripts: network + meshing). "
            f"Also honoured via {RUN_DEMOS_ENV}=1."
        ),
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register the demo marker."""
    config.addinivalue_line(
        "markers",
        "demo: opt-in tests that execute full demo scripts (network + meshing)",
    )


def _demos_requested(config: pytest.Config) -> bool:
    if config.getoption("--run-demos"):
        return True
    if os.environ.get(RUN_DEMOS_ENV) == "1":
        return True
    markexpr = (config.option.markexpr or "").strip()
    return "demo" in markexpr


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Deselect demo-execution tests unless explicitly requested."""
    demo_items = [item for item in items if item.get_closest_marker("demo")]
    if not demo_items or _demos_requested(config):
        return

    keep_items = [item for item in items if item not in demo_items]
    config.hook.pytest_deselected(items=demo_items)
    items[:] = keep_items
