"""Print the per-product export coverage table of the built-in Datasets.

The table is generated from the Dataset registry and each Dataset's
``describe()["export"]`` declarations; ``docs/datasets/export-coverage.md``
contains it verbatim.
"""

from dtcc_core.datasets import registry

HEADER = (
    "| Dataset | Product | Native type | Canonical | Web formats | Constraints "
    "| Exclusion |"
)

# Native types for products whose describe() entry carries no per-product type.
# orthophoto reports Raster for the Dataset but its tiles product returns a
# RasterTileCollection.
PRODUCT_TYPES = {("orthophoto", "tiles"): "RasterTileCollection"}


def builtin_datasets():
    """Registered Datasets defined under dtcc_core.datasets, sorted by name."""
    return {
        name: dataset
        for name, dataset in sorted(registry.list_datasets().items())
        if type(dataset).__module__.startswith("dtcc_core.datasets.")
    }


def native_type(name, described, product):
    """Class name of what the product returns."""
    per_product = {
        p["name"]: p["python_return_type"] for p in described.get("products", [])
    }
    qualified = per_product.get(product)
    if qualified is None:
        qualified = PRODUCT_TYPES.get((name, product), described["python_return_type"])
    return qualified.rpartition(".")[2]


def _constraints(entry):
    parts = []
    for fmt, constraint in entry["constraints"].items():
        needs = [f"{k}={v}" for k, v in constraint["arguments"].items()]
        needs += constraint["preconditions"]
        parts.append(f"{fmt}: {', '.join(needs)}")
    return "; ".join(parts) or "none"


def _exclusions(entry):
    return (
        "; ".join(
            f"{', '.join(e['formats'])}: {e['reason']}" for e in entry["exclusions"]
        )
        or "none"
    )


def rows():
    """One row per declared (Dataset, product) pair."""
    result = []
    for name, dataset in builtin_datasets().items():
        described = dataset.describe()
        for entry in described["export"]:
            product = entry["product"]
            result.append(
                (
                    f"`{name}`",
                    f"`{product}`" if product else "-",
                    f"`{native_type(name, described, product)}`",
                    "yes" if entry["canonical"] else "no",
                    ", ".join(f"`{f}`" for f in entry["web_formats"]) or "none",
                    _constraints(entry),
                    _exclusions(entry),
                )
            )
    return result


def coverage_table():
    """The markdown table, without a trailing newline."""
    lines = [HEADER, "|---|---|---|---|---|---|---|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows()]
    return "\n".join(lines)


if __name__ == "__main__":
    print(coverage_table())
