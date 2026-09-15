"""Static (offline, CI-friendly) checks for every demo script.

These never execute a demo.  They parse each file to catch syntax errors,
forbidden imports of the retired ``dtcc`` umbrella package, and unguarded
viewer calls in workflow demos.
"""

import ast
from pathlib import Path

import pytest

from tests.demos.demo_inventory import DATASET_DEMOS, WORKFLOW_DEMOS

REPO_ROOT = Path(__file__).resolve().parents[2]
DEMO_DIR = REPO_ROOT / "demos"

ALL_DEMOS = tuple(sorted(set(DATASET_DEMOS) | set(WORKFLOW_DEMOS)))

# Package names that must never appear in a migrated demo.
RETIRED_IMPORT_ROOTS = {"dtcc", "dtcc_platform", "dtcc_data", "dtcc_io"}


def _demo_source(demo_name: str) -> str:
    return (DEMO_DIR / demo_name).read_text(encoding="utf-8")


@pytest.mark.parametrize("demo_name", ALL_DEMOS)
def test_demo_parses(demo_name):
    """Every demo compiles and parses cleanly."""
    source = _demo_source(demo_name)
    compile(source, demo_name, "exec")
    ast.parse(source)


@pytest.mark.parametrize("demo_name", ALL_DEMOS)
def test_demo_does_not_import_retired_dtcc(demo_name):
    """No demo may depend on the retired ``dtcc`` umbrella package."""
    tree = ast.parse(_demo_source(demo_name))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".")[0]
                assert root not in RETIRED_IMPORT_ROOTS, (
                    f"{demo_name} imports retired package {alias.name!r}"
                )
        elif isinstance(node, ast.ImportFrom):
            root = (node.module or "").split(".")[0]
            assert root not in RETIRED_IMPORT_ROOTS, (
                f"{demo_name} imports from retired package {node.module!r}"
            )


@pytest.mark.parametrize("demo_name", WORKFLOW_DEMOS)
def test_workflow_demo_viewer_is_optional(demo_name):
    """Workflow demos must gate every ``.view(`` call behind ``--view``."""
    source = _demo_source(demo_name)
    if ".view(" in source:
        assert '"--view" in sys.argv' in source or "'--view' in sys.argv" in source, (
            f"{demo_name} calls .view() without a '--view' in sys.argv guard"
        )
