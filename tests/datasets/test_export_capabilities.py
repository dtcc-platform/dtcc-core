"""Each Dataset's declared export capabilities match what its products do."""

from __future__ import annotations

import json

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import load_model_package, registry
from tests.datasets.export_fixtures import PRODUCTS, VIOLATIONS, build

WEB_FORMATS = ("geojson", "png", "tif")
ENTRY_KEYS = {
    "product",
    "canonical",
    "web_formats",
    "constraints",
    "options",
    "exclusions",
}
# Valid values of every export option; field covers a scalar and a vector component.
OPTION_VALUES = {
    "field": ("pressure", "velocity:u"),
    "band": (1,),
    "colormap": ("magma",),
    "value_range": ((0.0, 1.0),),
    "unit": ("K",),
}
# Defaults that describe a computed value rather than name it.
DESCRIBED_DEFAULTS = {"valid min/max", "Field unit"}
VIOLATION_MESSAGES = {
    "slice_axis": "slice_axis",
    "translation_only_root_transform": "translation",
    "flat_in_exported_frame": "flat",
    "crs_resolvable": "source CRS",
}
REASON_3D = "3D-only product: no map-placed or RFC 7946 web form"
REASON_TILES = "web output out of scope by operator decision on #194"

BUILTINS = {
    name: dataset
    for name, dataset in registry.list_datasets().items()
    if type(dataset).__module__.startswith("dtcc_core.datasets.")
}
ENTRIES = [
    pytest.param(name, entry, id=f"{name}-{entry['product']}")
    for name, dataset in sorted(BUILTINS.items())
    for entry in dataset.describe().get("export", [])
]


def _params(each):
    return [
        pytest.param(name, entry, *values, id="-".join(map(str, (param.id, *values))))
        for param in ENTRIES
        for name, entry in [param.values]
        for values in each(name, entry)
    ]


def _export(model, tmp_path, label, **options):
    return model.export(tmp_path / label, canonical=True, **options)


def _derived(package):
    [artifact] = [a for a in package.artifacts if a.role == "derived"]
    return artifact


def test_every_builtin_dataset_declares_each_product():
    assert len(BUILTINS) == 25
    for name, dataset in BUILTINS.items():
        product = dataset.describe()["args_schema"]["properties"].get("product")
        expected = product["enum"] if product else [None]
        assert [e["product"] for e in dataset.describe()["export"]] == expected, name
    assert {(p.values[0], p.values[1]["product"]) for p in ENTRIES} == set(PRODUCTS)
    for name, entry in (param.values for param in ENTRIES):
        declared = {key for (key,) in _violations(entry)}
        assert declared == VIOLATIONS.get((name, entry["product"]), set()), name
        reason = (
            REASON_TILES
            if (name, entry["product"]) == ("orthophoto", "tiles")
            else REASON_3D
        )
        assert [e["reason"] for e in entry["exclusions"]] in ([], [reason]), name


def test_every_registered_product_has_a_web_format_or_an_exclusion():
    """Test registrations are included; a Dataset that declares nothing has no
    product to check."""
    for name, dataset in registry.list_datasets().items():
        for entry in dataset.describe()["export"]:
            assert set(entry) == ENTRY_KEYS, name
            assert set(entry["web_formats"]) <= set(WEB_FORMATS), name
            assert entry["web_formats"] or entry["exclusions"], name
            for exclusion in entry["exclusions"]:
                assert exclusion["reason"].strip(), name
                assert not set(exclusion["formats"]) & set(entry["web_formats"])
            for fmt, constraint in entry["constraints"].items():
                assert fmt in entry["web_formats"], name
                assert set(constraint) == {"arguments", "preconditions"}, name
            for spec in entry["options"].values():
                assert set(spec) == {"formats", "default"}, name
                assert set(spec["formats"]) <= set(entry["web_formats"]), name


def test_describe_is_json_serializable_for_every_registered_dataset():
    for dataset in registry.list_datasets().values():
        described = dataset.describe()
        assert json.loads(json.dumps(described, allow_nan=False)) == described


@pytest.mark.parametrize(("name", "entry"), ENTRIES)
def test_canonical_matches_the_populated_fixture(name, entry, tmp_path):
    model = build(name, entry["product"], "populated", tmp_path)
    described = BUILTINS[name].describe()
    products = {
        p["name"]: p["python_return_type"] for p in described.get("products", [])
    }
    expected = (
        described["python_return_type"]
        if entry["product"] is None
        else products.get(entry["product"])
    )
    if expected is not None:
        assert type(model).__name__ == expected.split(".")[-1]
    if not entry["canonical"]:
        with pytest.raises(ValueError):
            _export(model, tmp_path, "pkg")
        return
    restored = load_model_package(_export(model, tmp_path, "pkg").path)
    assert type(restored) is type(model)


@pytest.mark.parametrize(("name", "entry"), ENTRIES)
def test_empty_instance_exports_canonically_and_refuses_images(name, entry, tmp_path):
    model = build(name, entry["product"], "empty", tmp_path)
    if entry["canonical"]:
        load_model_package(_export(model, tmp_path, "pkg").path)
    for fmt in {"png", "tif"} & set(entry["web_formats"]):
        with pytest.raises(ValueError, match="empty|no samples|no valid pixels"):
            _export(model, tmp_path, fmt, format=fmt)


@pytest.mark.parametrize("name", ["terrain_surface_mesh", "trees"])
def test_all_invalid_raster_exports_canonically_and_refuses_png(name, tmp_path):
    model = build(name, "raster", "all_invalid", tmp_path)
    load_model_package(_export(model, tmp_path, "pkg").path)
    with pytest.raises(ValueError, match="no valid pixels"):
        _export(model, tmp_path, "png", format="png")


@pytest.mark.parametrize(
    ("name", "entry", "fmt"),
    _params(lambda name, entry: [(f,) for f in entry["web_formats"]]),
)
def test_advertised_web_format_succeeds(name, entry, fmt, tmp_path):
    model = build(name, entry["product"], "populated", tmp_path)
    package = _export(model, tmp_path, "pkg", format=fmt)
    load_model_package(package.path)
    artifact = _derived(package)
    assert artifact.format == fmt
    if fmt == "geojson":
        assert artifact.crs == "EPSG:4326"
        return
    assert artifact.placement is not None
    for option, spec in entry["options"].items():
        if fmt in spec["formats"] and spec["default"] not in DESCRIBED_DEFAULTS:
            assert getattr(artifact.rendering, option) == spec["default"], option


@pytest.mark.parametrize(
    ("name", "entry", "fmt"),
    _params(
        lambda name, entry: [(f,) for f in WEB_FORMATS if f not in entry["web_formats"]]
    ),
)
def test_unadvertised_web_format_fails(name, entry, fmt, tmp_path):
    model = build(name, entry["product"], "populated", tmp_path)
    with pytest.raises(ValueError):
        _export(model, tmp_path, "pkg", format=fmt)


def _violations(entry):
    keys = {
        key
        for constraint in entry["constraints"].values()
        for key in (*constraint["arguments"], *constraint["preconditions"])
    }
    return [(key,) for key in sorted(keys)]


@pytest.mark.parametrize(
    ("name", "entry", "key"),
    _params(
        lambda name, entry: [
            (key,) for key in sorted(VIOLATIONS.get((name, entry["product"]), ()))
        ]
    ),
)
def test_constraint_violation_fails_only_the_constrained_formats(
    name, entry, key, tmp_path
):
    model = build(name, entry["product"], key, tmp_path)
    load_model_package(_export(model, tmp_path, "pkg").path)
    for constraint in entry["constraints"].values():
        if key in constraint["arguments"]:
            parameters = model.dataset_context.request.parameters
            assert parameters[key] != constraint["arguments"][key]
    for fmt in entry["web_formats"]:
        constraint = entry["constraints"].get(fmt, {})
        if key in (
            *constraint.get("arguments", {}),
            *constraint.get("preconditions", []),
        ):
            with pytest.raises(ValueError, match=VIOLATION_MESSAGES[key]):
                _export(model, tmp_path, fmt, format=fmt)
        else:
            _export(model, tmp_path, fmt, format=fmt)


@pytest.mark.parametrize(
    ("name", "entry", "option"),
    _params(lambda name, entry: [(option,) for option in OPTION_VALUES]),
)
def test_option_applies_to_exactly_its_declared_formats(name, entry, option, tmp_path):
    model = build(name, entry["product"], "populated", tmp_path)
    declared = entry["options"].get(option, {"formats": []})["formats"]
    for fmt in (None, *entry["web_formats"]):
        for index, value in enumerate(OPTION_VALUES[option]):
            label = f"{fmt}-{index}"
            if fmt not in declared:
                with pytest.raises(ValueError, match=option):
                    _export(model, tmp_path, label, format=fmt, **{option: value})
                continue
            package = _export(model, tmp_path, label, format=fmt, **{option: value})
            recorded = getattr(_derived(package).rendering, option)
            assert recorded == (list(value) if option == "value_range" else value)


def test_listing_includes_the_export_declaration():
    rows = datasets._group_dataset_summary_rows(
        {
            "smoke": datasets.smoke,
            "orthophoto": datasets.orthophoto,
            "terrain_surface_mesh": datasets.terrain_surface_mesh,
        }
    )
    [smoke] = rows["synthetic"]
    [orthophoto] = rows["raw"]
    [terrain] = rows["derived"]
    assert smoke[-1] == (
        "field: none; slice: geojson*, png*, tif*; streamlines: geojson*"
    )
    assert orthophoto[-1] == "raster: png, tif; tiles: none"
    assert terrain[-1] == "mesh: geojson*; raster: png, tif"
    assert "Web export" in datasets.info(print=False)
