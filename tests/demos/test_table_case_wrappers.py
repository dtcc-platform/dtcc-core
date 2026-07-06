"""Compatibility checks for legacy table-case demo entry points."""

from __future__ import annotations

import importlib.util
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parents[2]


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_smoke_table_case_demo_wrapper_reexports_script_api():
    wrapper = _load(_REPO_ROOT / "demos" / "smoke_table_cases.py", "smoke_wrapper")
    script = _load(
        _REPO_ROOT / "scripts" / "table_cases" / "smoke_table_cases.py",
        "smoke_script",
    )

    assert wrapper.BOUNDS == script.BOUNDS
    assert wrapper.CASES == script.CASES
    assert wrapper.main is not None


def test_grid_table_case_demo_wrapper_reexports_script_api():
    wrapper = _load(_REPO_ROOT / "demos" / "grid_table_case.py", "grid_wrapper")
    script = _load(
        _REPO_ROOT / "scripts" / "table_cases" / "grid_table_case.py",
        "grid_script",
    )

    assert wrapper.BOUNDS == script.BOUNDS
    assert wrapper.CASE == script.CASE
    assert wrapper.main is not None


def test_footprints_table_case_demo_wrapper_reexports_script_api():
    wrapper = _load(
        _REPO_ROOT / "demos" / "footprints_table_case.py",
        "footprints_wrapper",
    )
    script = _load(
        _REPO_ROOT / "scripts" / "table_cases" / "footprints_table_case.py",
        "footprints_script",
    )

    assert wrapper.BOUNDS == script.BOUNDS
    assert wrapper.CASE == script.CASE
    assert wrapper.main is not None
