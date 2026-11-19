"""Pytest configuration for live dataset tests."""

from __future__ import annotations

import importlib
import os
from pathlib import Path

import pytest

from dtcc_core.datasets import get_dataset, load_model_package

LIVE_DATASET_TESTS_ENV = "DTCC_LIVE_DATASET_TESTS"


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add an opt-in switch for live dataset tests."""
    parser.addoption(
        "--run-live",
        action="store_true",
        default=False,
        help=(f"Run tests marked as live. Also requires {LIVE_DATASET_TESTS_ENV}=1."),
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
            reason=(f"Set {LIVE_DATASET_TESTS_ENV}=1 to run live provider tests.")
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


@pytest.fixture
def canonical_round_trip(tmp_path):
    """Export a Dataset result as a canonical package and load it back.

    Asserts the restored type is the Dataset's ``python_return_type`` and the
    whole ``dataset_context`` survives, then returns the restored object.
    """

    def run(dataset_name, result):
        module_name, _, class_name = get_dataset(
            dataset_name
        ).python_return_type.rpartition(".")
        expected = getattr(importlib.import_module(module_name), class_name)
        package = result.export(tmp_path / "x.dtccpkg", canonical=True)
        restored = load_model_package(package.path)
        assert type(restored) is expected
        assert result.dataset_context is not None
        assert restored.dataset_context == result.dataset_context
        return restored

    return run
