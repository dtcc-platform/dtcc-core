"""FieldSlice.to_raster on the stored local grid."""

from __future__ import annotations

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import load_model_package
from dtcc_core.model import Bounds, Field, FieldSlice

# A 5 x 5 grid whose spacing differs per axis (7.5 in x, 20 in y), so a
# transposed or half-cell-shifted raster cannot match.
DOMAIN = Bounds(xmin=100.0, ymin=200.0, zmin=0.0, xmax=130.0, ymax=280.0, zmax=9.0)
RESOLUTION = 5


def _grid_points(z=3.0):
    xs = np.linspace(DOMAIN.xmin, DOMAIN.xmax, RESOLUTION)
    ys = np.linspace(DOMAIN.ymin, DOMAIN.ymax, RESOLUTION)
    aa, bb = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([aa.ravel(), bb.ravel(), np.full(aa.size, z)])


def _slice(points=None, resolution=RESOLUTION, **overrides):
    points = _grid_points() if points is None else points
    model = FieldSlice(
        points=points,
        resolution=resolution,
        **{"crs": "EPSG:3006", "domain_bounds": DOMAIN, **overrides},
    )
    model.add_field(
        Field(
            association="sample",
            name="speed",
            unit="m/s",
            values=(points[:, 0] + 1000.0 * points[:, 1]).reshape(-1, 1),
            dim=1,
        )
    )
    model.add_field(
        Field(
            association="sample",
            name="velocity",
            unit="m/s",
            values=points.copy(),
            dim=3,
        )
    )
    return model


def _shifted(raster):
    return raster.georef.c, raster.georef.f


def test_values_sit_at_the_cell_of_their_sample():
    model = _slice()
    raster, rendering = model.to_raster()

    spacing = (DOMAIN.width / 4, DOMAIN.height / 4)
    assert raster.data.shape == (5, 5)
    assert raster.crs == "EPSG:3006"
    assert np.isnan(raster.nodata)
    for x, y, _ in model.points:
        col = round((x - DOMAIN.xmin) / spacing[0])
        row = 4 - round((y - DOMAIN.ymin) / spacing[1])
        assert raster.data[row, col] == x + 1000.0 * y
        assert raster.get_value(x, y) == x + 1000.0 * y
    assert rendering == {"field": "speed", "unit": "m/s"}


def test_pixel_centres_coincide_with_the_first_and_last_samples():
    model = _slice()
    raster, _ = model.to_raster()
    first, last = model.points[0], model.points[-1]

    # first sample is the south-west corner, the last the north-east one
    assert raster.georef * (0.5, 4 + 0.5) == (first[0], first[1])
    assert raster.georef * (4 + 0.5, 0.5) == (last[0], last[1])
    assert raster.data[4, 0] == first[0] + 1000.0 * first[1]
    assert raster.data[0, 4] == last[0] + 1000.0 * last[1]
    dx, dy = DOMAIN.width / 4, DOMAIN.height / 4
    assert tuple(raster.georef)[:6] == (
        dx,
        0.0,
        DOMAIN.xmin - dx / 2,
        0.0,
        -dy,
        DOMAIN.ymax + dy / 2,
    )
    bounds = raster.bounds
    assert (bounds.xmin, bounds.xmax) == (DOMAIN.xmin - dx / 2, DOMAIN.xmax + dx / 2)
    assert (bounds.ymin, bounds.ymax) == (DOMAIN.ymin - dy / 2, DOMAIN.ymax + dy / 2)


def test_columns_follow_x_and_rows_follow_decreasing_y():
    model = _slice()
    raster, _ = model.to_raster()
    # x varies along columns (7.5 per step), y down the rows (20 per step)
    np.testing.assert_allclose(np.diff(raster.data[0]), 7.5)
    np.testing.assert_allclose(np.diff(raster.data[:, 0]), -20000.0)


def test_sample_order_does_not_matter():
    model = _slice()
    expected, _ = model.to_raster()
    order = np.random.default_rng(0).permutation(len(model.points))
    shuffled = _slice(model.points[order])

    raster, _ = shuffled.to_raster()
    np.testing.assert_array_equal(raster.data, expected.data)


def test_a_missing_outer_row_leaves_nodata_and_keeps_the_rest_in_place():
    points = _grid_points()
    kept = points[points[:, 1] < DOMAIN.ymax]
    full, _ = _slice().to_raster()

    raster, _ = _slice(kept).to_raster()

    assert np.isnan(raster.data[0]).all()
    np.testing.assert_array_equal(raster.data[1:], full.data[1:])
    assert raster.georef == full.georef


def test_a_translation_shifts_the_georef_only():
    plain, _ = _slice().to_raster()
    model = _slice()
    model.transform.set_translation(10.0, 20.0, 5.0)

    raster, _ = model.to_raster()

    assert _shifted(raster) == (plain.georef.c + 10.0, plain.georef.f + 20.0)
    assert (raster.georef.a, raster.georef.e) == (plain.georef.a, plain.georef.e)
    np.testing.assert_array_equal(raster.data, plain.data)
    moved, base = raster.bounds, plain.bounds
    assert (moved.xmin, moved.ymin, moved.xmax, moved.ymax) == (
        base.xmin + 10.0,
        base.ymin + 20.0,
        base.xmax + 10.0,
        base.ymax + 20.0,
    )


@pytest.mark.parametrize(
    "linear",
    [
        [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]],
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ],
    ids=["rotation-x", "rotation-z", "scale", "shear"],
)
def test_a_root_transform_that_is_not_a_translation_raises(linear):
    model = _slice()
    model.transform.set_rotation(np.array(linear))

    with pytest.raises(ValueError, match=r"translation.*GeoJSON"):
        model.to_raster()


def test_the_z_slice_is_the_only_axis():
    model = _slice(slice_axis="x")

    with pytest.raises(ValueError, match=r"slice_axis='z'.*GeoJSON"):
        model.to_raster()


def test_a_vector_component_is_selected_by_name():
    model = _slice()

    raster, rendering = model.to_raster("velocity:v")

    assert rendering == {"field": "velocity:v", "unit": "m/s"}
    for x, y, _ in model.points:
        assert raster.get_value(x, y) == y
    raster_u, _ = model.to_raster(field="velocity:u")
    assert raster_u.get_value(DOMAIN.xmax, DOMAIN.ymin) == DOMAIN.xmax


def test_the_default_field_is_speed():
    model = _slice()
    assert model.to_raster()[1]["field"] == "speed"
    np.testing.assert_array_equal(
        model.to_raster()[0].data, model.to_raster("speed")[0].data
    )


@pytest.mark.parametrize(
    ("field", "message"),
    [
        ("pressure", "no field named 'pressure'"),
        ("velocity:q", "no component 'q'"),
        ("speed:u", "no component 'u'"),
        ("velocity", "needs a component"),
        ("speed:", "no component ''"),
        (3, "field must be a string"),
    ],
)
def test_a_missing_field_or_component_raises(field, message):
    with pytest.raises(ValueError, match=message):
        _slice().to_raster(field)


def test_the_raster_crs_comes_from_the_transform_srs_when_crs_is_absent():
    model = _slice(crs=None)
    model.transform.srs = "EPSG:3006"

    assert model.to_raster()[0].crs == "EPSG:3006"


def test_conflicting_crs_declarations_raise():
    model = _slice()
    model.transform.srs = "EPSG:4326"

    with pytest.raises(ValueError, match="Conflicting CRS declarations"):
        model.to_raster()


def test_a_slice_without_any_crs_raises():
    with pytest.raises(ValueError, match="Cannot resolve the source CRS"):
        _slice(crs=None).to_raster()


def test_a_geometry_field_is_not_a_sample_field():
    model = _slice()
    model.add_field(
        Field(association="geometry", name="mean", unit="K", values=[[1.0]], dim=1)
    )

    with pytest.raises(ValueError, match="one value per sample"):
        model.to_raster("mean")


@pytest.mark.parametrize(
    ("make", "message"),
    [
        (lambda: _slice(np.zeros((0, 3))), "no samples"),
        (lambda: _slice(domain_bounds=None), "domain_bounds"),
        (lambda: _slice(resolution=1), "resolution must be at least 2"),
        (lambda: _slice(resolution=0), "resolution must be at least 2"),
        (
            lambda: _slice(np.vstack([_grid_points(), _grid_points()[:1]])),
            "duplicate",
        ),
        (
            lambda: _slice(_grid_points() + [1.0, 0.0, 0.0]),
            "not on the grid",
        ),
        (
            lambda: _slice(_grid_points() + [0.0, 1e-3, 0.0]),
            "not on the grid",
        ),
        (
            lambda: _slice(_grid_points() + [DOMAIN.width + 7.5, 0.0, 0.0]),
            "not on the grid",
        ),
        (
            lambda: _slice(np.vstack([_grid_points()[:-1], [[np.nan, 0.0, 0.0]]]), 5),
            "not on the grid",
        ),
    ],
    ids=[
        "empty",
        "no-domain-bounds",
        "resolution-1",
        "resolution-0",
        "duplicate",
        "off-grid-x",
        "off-grid-y",
        "outside",
        "nan",
    ],
)
def test_an_incomplete_grid_contract_raises(make, message):
    model = make()
    with pytest.raises(ValueError, match=message):
        model.to_raster()


def test_a_sample_within_the_tolerance_maps_to_its_cell():
    points = _grid_points()
    points[:, 0] += 5e-7
    raster, _ = _slice(points).to_raster()
    expected, _ = _slice().to_raster()

    np.testing.assert_allclose(raster.data, expected.data, rtol=1e-6)


def test_a_degenerate_extent_raises():
    flat = Bounds(xmin=1.0, ymin=0.0, zmin=0.0, xmax=1.0, ymax=5.0, zmax=1.0)
    model = _slice(domain_bounds=flat)

    with pytest.raises(ValueError, match="domain_bounds"):
        model.to_raster()


def _smoke_slice(path, transform=None):
    result = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        zmax=40.0,
        resolution=4,
        product="slice",
        crs="EPSG:3006",
    )
    if transform is not None:
        transform(result)
    package = result.export(path, canonical=True)
    restored = load_model_package(package.path)
    assert type(restored) is FieldSlice
    return result, restored


def test_the_smoke_slice_raster_covers_the_domain_and_survives_a_reload(tmp_path):
    result, restored = _smoke_slice(tmp_path / "pkg")

    for model in (result, restored):
        raster, rendering = model.to_raster()
        assert raster.data.shape == (4, 4)
        assert not np.isnan(raster.data).any()
        assert raster.crs == "EPSG:3006"
        assert rendering == {"field": "speed", "unit": "m/s"}
        speed = model.field("speed").values[:, 0]
        for (x, y, _), value in zip(model.points, speed, strict=True):
            assert raster.get_value(x, y) == value
    np.testing.assert_array_equal(
        restored.to_raster()[0].data, result.to_raster()[0].data
    )


def test_translation_and_rotation_behave_the_same_after_a_reload(tmp_path):
    def translate(model):
        model.transform.set_translation(10.0, 20.0, 5.0)

    def rotate(model):
        model.transform.set_rotation(
            np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
        )

    plain, _ = _smoke_slice(tmp_path / "plain")
    base = plain.to_raster()[0]
    _, moved = _smoke_slice(tmp_path / "moved", translate)
    raster = moved.to_raster()[0]
    assert _shifted(raster) == (base.georef.c + 10.0, base.georef.f + 20.0)
    np.testing.assert_array_equal(raster.data, base.data)
    assert (raster.bounds.xmin, raster.bounds.ymax) == (
        base.bounds.xmin + 10.0,
        base.bounds.ymax + 20.0,
    )

    _, turned = _smoke_slice(tmp_path / "turned", rotate)
    with pytest.raises(ValueError, match="GeoJSON"):
        turned.to_raster()


def _integer_slice(**fields):
    model = _slice()
    model.fields = [
        Field(
            association="sample",
            name=name,
            unit="m/s",
            values=np.asarray(values, dtype=np.int64).reshape(len(model.points), -1),
            dim=np.asarray(values).reshape(len(model.points), -1).shape[1],
        )
        for name, values in fields.items()
    ]
    return model


def _counts(value=2**53 + 1):
    counts = np.arange(25, dtype=np.int64)
    counts[7] = value
    return counts


def test_a_scalar_field_float64_would_round_is_refused():
    model = _integer_slice(count=_counts())

    with pytest.raises(ValueError, match=rf"Field 'count'.*{2**53 + 1}"):
        model.to_raster(field="count")


def test_a_selected_component_float64_would_round_is_refused():
    rows = np.zeros((25, 3), dtype=np.int64)
    rows[4, 0] = 2**53 + 1
    model = _integer_slice(flow=rows)

    with pytest.raises(ValueError, match=rf"Field 'flow'.*{2**53 + 1}"):
        model.to_raster(field="flow:u")


def test_an_integer_field_float64_holds_exactly_fills_the_raster():
    values = 2**60 - 2**10 * np.arange(25, dtype=np.int64)
    model = _integer_slice(count=values)

    raster, rendering = model.to_raster(field="count")

    assert rendering == {"field": "count", "unit": "m/s"}
    for (x, y, _), value in zip(model.points, values.tolist(), strict=True):
        assert int(raster.get_value(x, y)) == value
    assert raster.get_value(*model.points[0][:2]) == float(2**60)


def _published(model, tmp_path, name, **kwargs):
    from dtcc_core.datasets import attach_dataset_context
    from dtcc_core.datasets.schema import (
        DatasetContext,
        DatasetIdentity,
        DatasetMetadata,
        DatasetPresentation,
        DatasetProvenance,
        DatasetRequest,
    )

    attach_dataset_context(
        model,
        DatasetContext(
            identity=DatasetIdentity(name="example", title="Example"),
            metadata=DatasetMetadata(crs=["EPSG:3006"]),
            provenance=DatasetProvenance(sources=["example"]),
            presentation=DatasetPresentation(),
            request=DatasetRequest(dataset_name="example", bounds=[0.0, 0.0, 1.0, 1.0]),
            health={"status": "ok"},
            warnings=[],
        ),
    )
    return model.export(tmp_path / name, canonical=True, **kwargs)


@pytest.mark.usefixtures("no_network")
def test_a_canonical_package_keeps_an_integer_but_a_tif_refuses_it(tmp_path):
    model = _integer_slice(speed=_counts())

    package = _published(model, tmp_path, "plain")
    restored = load_model_package(package.path)
    assert int(restored.field("speed").values[7, 0]) == 2**53 + 1

    for format in ("tif", "geojson"):
        target = tmp_path / f"{format}_pkg"
        with pytest.raises(ValueError, match=rf"Field 'speed'.*{2**53 + 1}"):
            model.export(target, canonical=True, format=format)
        assert not target.exists()


@pytest.mark.usefixtures("no_network")
def test_a_canonical_tif_holds_an_integer_float64_holds_exactly(tmp_path):
    import rasterio

    values = 2**60 - 2**10 * np.arange(25, dtype=np.int64)
    assert not np.array_equal(values.astype(np.float32).astype(np.int64), values)
    model = _integer_slice(speed=values)

    package = _published(model, tmp_path, "pkg", format="tif")

    [artifact] = [a for a in package.artifacts if a.format == "tif"]
    with rasterio.open(package.path / artifact.path) as dataset:
        band = dataset.read(1)
        for point, expected in zip(model.points, values, strict=True):
            row, column = dataset.index(point[0], point[1])
            assert int(band[row, column]) == int(expected)
