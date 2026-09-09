"""Opt-in execution tests for workflow demos.

These run each workflow demo as a script.  They download footprints and point
clouds from the DTCC data service and run meshing, so they are marked ``demo``
and are deselected by default.  Run them explicitly with::

    DTCC_RUN_DEMOS=1 pytest tests/demos/test_demo_execution.py

or ``pytest tests/demos --run-demos``.
"""

import runpy
from pathlib import Path

import pytest

from tests.demos.demo_inventory import WORKFLOW_DEMOS

REPO_ROOT = Path(__file__).resolve().parents[2]
DEMO_DIR = REPO_ROOT / "demos"


@pytest.mark.demo
@pytest.mark.parametrize("demo_name", WORKFLOW_DEMOS)
def test_workflow_demo_executes_cleanly(demo_name, monkeypatch):
    """
    Run each workflow demo as a script to verify it imports properly,
    parses its bounds, downloads necessary data, and generates output
    without throwing exceptions.

    To avoid opening UI windows, we ensure ``--view`` is not in sys.argv.
    """
    demo_path = DEMO_DIR / demo_name

    monkeypatch.setattr("sys.argv", ["test_demo_execution"])

    runpy.run_path(str(demo_path), run_name="__main__")
