"""docs/datasets/export-coverage.md contains the generated coverage table."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "docs" / "datasets" / "export-coverage.md"
SCRIPT = ROOT / "scripts" / "export_coverage_table.py"

spec = importlib.util.spec_from_file_location("export_coverage_table", SCRIPT)
table = importlib.util.module_from_spec(spec)
spec.loader.exec_module(table)


def test_doc_contains_the_generated_table_verbatim():
    assert table.coverage_table() in DOC.read_text(encoding="utf-8")


def test_every_builtin_product_appears_once_in_the_doc():
    declared = [
        (name, entry["product"])
        for name, dataset in table.builtin_datasets().items()
        for entry in dataset.describe()["export"]
    ]
    assert len(declared) == 30
    assert len({name for name, _ in declared}) == 25
    lines = DOC.read_text(encoding="utf-8").splitlines()
    for name, product in declared:
        label = f"`{product}`" if product else "-"
        prefix = f"| `{name}` | {label} | "
        assert sum(line.startswith(prefix) for line in lines) == 1, (name, product)


def test_a_product_without_web_format_states_its_exclusion():
    for dataset, product, _, _, formats, _, exclusion in table.rows():
        if formats == "none":
            assert exclusion != "none", (dataset, product)


def test_tiles_native_type_comes_from_the_explicit_mapping():
    from dtcc_core.datasets import registry
    from dtcc_core.model.values.raster_tiles import RasterTileCollection

    assert table.PRODUCT_TYPES == {("orthophoto", "tiles"): "RasterTileCollection"}
    described = registry.list_datasets()["orthophoto"].describe()
    assert described["python_return_type"] == "dtcc_core.model.Raster"
    assert "products" not in described or not described["products"]
    assert table.native_type("orthophoto", described, "tiles") == (
        RasterTileCollection.__name__
    )
    assert table.native_type("orthophoto", described, "raster") == "Raster"


def test_mapping_never_overrides_a_per_product_type():
    from dtcc_core.datasets import registry

    for name, product in table.PRODUCT_TYPES:
        described = registry.list_datasets()[name].describe()
        assert product not in {p["name"] for p in described.get("products", [])}


@pytest.mark.parametrize("name", ["smoke", "terrain_surface_mesh", "trees"])
def test_per_product_types_are_used_where_describe_supplies_them(name):
    from dtcc_core.datasets import registry

    described = registry.list_datasets()[name].describe()
    for product in described["products"]:
        expected = product["python_return_type"].rpartition(".")[2]
        assert table.native_type(name, described, product["name"]) == expected
