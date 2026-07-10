"""Static hygiene checks for curated public demos."""

from pathlib import Path

from tests.demos.demo_inventory import (
    NO_DEMO_REASONS,
    NORMAL_DEMOS,
    TABLE_DEMO_MAPPING,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEMO_DIR = REPO_ROOT / "demos"

FORBIDDEN_PATTERNS = (
    "from __future__ import annotations",
    "from pathlib import Path",
    "import pathlib",
    "import os",
    "import tempfile",
    "MPLCONFIGDIR",
    "MPLBACKEND",
    "os.environ",
    "mkdir",
    "publish",
    "upload",
    "DTCC_UPLOAD",
    "dataset_key",
    "importlib.util",
    "spec_from_file_location",
)


def test_design_document_exists_and_demo_catalog_references_it():
    design = REPO_ROOT / "DESIGN.md"
    catalog = REPO_ROOT / "docs" / "datasets" / "demo-catalog.md"

    assert design.is_file()
    assert catalog.is_file()
    assert "DESIGN.md" in catalog.read_text(encoding="utf-8")


def test_legacy_table_case_entry_points_are_removed():
    legacy_demos = sorted(DEMO_DIR.glob("*table*case*.py"))
    legacy_scripts = sorted((REPO_ROOT / "scripts").glob("**/table_cases/**/*.py"))

    assert legacy_demos == []
    assert legacy_scripts == []


def test_every_demo_is_intentionally_classified():
    actual = {
        path.name
        for path in DEMO_DIR.glob("*.py")
        if path.name != "__init__.py"
    }

    assert actual == set(NORMAL_DEMOS)


def test_normal_demos_are_small_and_boring():
    for demo_name in NORMAL_DEMOS:
        path = DEMO_DIR / demo_name
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()

        assert "import dtcc_core as dtcc" in lines, demo_name
        extra_imports = [
            line
            for line in lines
            if line.startswith("import ") or line.startswith("from ")
            if line != "import dtcc_core as dtcc"
        ]
        assert extra_imports == [], f"{demo_name} has extra imports: {extra_imports}"
        assert len(lines) <= 40, f"{demo_name} should stay readable in under a minute"

        for pattern in FORBIDDEN_PATTERNS:
            assert pattern not in text, f"{demo_name} contains forbidden pattern {pattern!r}"


def test_table_demo_mapping_is_complete_and_points_to_real_files():
    mapped = set(TABLE_DEMO_MAPPING)
    explained = set(NO_DEMO_REASONS)

    assert mapped.isdisjoint(explained)
    assert mapped | explained

    for dataset_id, demo_path in TABLE_DEMO_MAPPING.items():
        path = REPO_ROOT / demo_path
        assert path.is_file(), f"{dataset_id} maps to missing demo {demo_path}"
        assert path.name in NORMAL_DEMOS

    for dataset_id, reason in NO_DEMO_REASONS.items():
        assert reason.strip(), f"{dataset_id} needs an explicit no-demo reason"
