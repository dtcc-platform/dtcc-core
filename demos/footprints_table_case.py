"""Compatibility wrapper for ``scripts/table_cases/footprints_table_case.py``."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


_TARGET = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "table_cases"
    / "footprints_table_case.py"
)


def _load_target():
    spec = importlib.util.spec_from_file_location(
        "_dtcc_footprints_table_case",
        _TARGET,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


_module = _load_target()

BOUNDS = _module.BOUNDS
OUTPUT_DIR = _module.OUTPUT_DIR
CASE = _module.CASE
main = _module.main


if __name__ == "__main__":
    sys.exit(main())
