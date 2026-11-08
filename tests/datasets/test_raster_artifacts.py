"""Map-placed PNG and GeoTIFF artifacts in canonical packages."""

from __future__ import annotations

import json
import math

import matplotlib
import numpy as np
import pytest
import rasterio
import rasterio.warp
from affine import Affine
from PIL import Image
from pyproj import Transformer
from rasterio.crs import CRS
from rasterio.warp import Resampling, reproject

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context, load_model_package
from dtcc_core.datasets.schema import DatasetArtifact, DatasetManifest
from dtcc_core.io.raster import valid_pixels
from dtcc_core.model import FieldSlice, Raster, StreamlineCollection

ROTATED = (
    Affine.translation(674000.0, 6580000.0)
    @ Affine.rotation(30.0)
    @ Affine.scale(10.0, -10.0)
)
NORTH_UP = Affine(10.0, 0.0, 674000.0, 0.0, -10.0, 6580000.0)
TO_SOURCE = Transformer.from_crs("EPSG:3857", "EPSG:3006", always_xy=True)
TO_WGS84 = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)


def _context(dataset):
    return dataset.create_context(dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)}))


def _raster(data, *, georef=ROTATED, nodata=np.nan, dataset=datasets.city):
    raster = Raster(data=data, georef=georef, crs="EPSG:3006", nodata=nodata)
    attach_dataset_context(raster, _context(dataset))
    return raster


def _smooth(height=12, width=16):
    rows, cols = np.mgrid[0:height, 0:width]
    return 100.0 + 3.0 * cols - 2.0 * rows


def _export(model, path, **options):
    package = model.export(path, canonical=True, **options)
    [derived] = [a for a in package.artifacts if a.role == "derived"]
    return package, derived


def _png(package, artifact):
    with Image.open(package.path / artifact.path) as image:
        assert image.mode == "RGBA"
        return np.asarray(image)


def _colour(values, rendering):
    low, high = rendering.value_range
    if high > low:
        t = np.clip((values - low) / (high - low), 0.0, 1.0)
    else:
        t = np.zeros_like(values, dtype=float)
    return matplotlib.colormaps[rendering.colormap](t, bytes=True)


def _source_positions(placement, georef):
    """Source pixel coordinates of every destination pixel centre."""
    grid = Affine(*placement.transform)
    rows, cols = np.mgrid[0 : placement.height, 0 : placement.width]
    x, y = grid @ (cols + 0.5, rows + 0.5)
    sx, sy = TO_SOURCE.transform(x, y)
    return ~georef @ (sx, sy)


def _nearest(placement, georef, shape):
    """Nearest source row/column and whether the centre lies on the source grid."""
    u, v = _source_positions(placement, georef)
    # GDAL and pyproj agree to far better than this margin; a centre this close
    # to a pixel edge would make the nearest source pixel ambiguous.
    assert np.abs(u - np.round(u)).min() > 1e-6
    assert np.abs(v - np.round(v)).min() > 1e-6
    inside = (u >= 0) & (u < shape[1]) & (v >= 0) & (v < shape[0])
    rows = np.clip(np.floor(v).astype(int), 0, shape[0] - 1)
    cols = np.clip(np.floor(u).astype(int), 0, shape[1] - 1)
    return rows, cols, inside


def _assert_3857_grid(placement):
    assert placement.crs == "EPSG:3857"
    grid = Affine(*placement.transform)
    assert grid.b == grid.d == 0.0 and grid.a > 0.0 and grid.e < 0.0
    xmin, ymax = grid.c, grid.f
    xmax, ymin = grid @ (placement.width, placement.height)
    assert placement.bounds == pytest.approx([xmin, ymin, xmax, ymax], abs=1e-6)
    lon0, lat0 = TO_WGS84.transform(xmin, ymin)
    lon1, lat1 = TO_WGS84.transform(xmax, ymax)
    assert placement.bounds_wgs84 == pytest.approx([lon0, lat0, lon1, lat1], abs=1e-9)


def test_png_places_integer_labels_by_nearest_source_pixel(tmp_path):
    labels = np.arange(1, 25, dtype=np.int32).reshape(4, 6) * 7
    labels[2, 3] = -1
    raster = _raster(labels, nodata=-1)

    package, artifact = _export(raster, tmp_path / "pkg", format="png")

    assert (artifact.format, artifact.media_type, artifact.data_kind) == (
        "png",
        "image/png",
        "raster",
    )
    assert artifact.crs is None and artifact.bounds is None
    placement, rendering = artifact.placement, artifact.rendering
    _assert_3857_grid(placement)
    assert placement.nodata == -1.0
    assert rendering.model_dump() == {
        "bands": None,
        "field": None,
        "band": 0,
        "unit": None,
        "value_range": [7.0, 168.0],
        "colormap": "viridis",
        "transparent_nodata": True,
    }
    rows, cols, inside = _nearest(placement, ROTATED, labels.shape)
    label = labels[rows, cols]
    valid = inside & (label != -1)
    assert valid.any() and (inside & ~valid).any() and (~inside).any()
    expected = np.zeros((placement.height, placement.width, 4), dtype=np.uint8)
    expected[valid] = _colour(label[valid].astype(float), rendering)
    np.testing.assert_array_equal(_png(package, artifact), expected)
    # The canonical model keeps the source grid.
    restored = load_model_package(package.path)
    assert restored.georef == ROTATED and restored.crs == "EPSG:3006"
    np.testing.assert_array_equal(restored.data, labels)


def test_png_interpolates_floating_values_bilinearly(tmp_path):
    values = _smooth()
    values[0, 0] = np.nan
    raster = _raster(values)

    package, artifact = _export(raster, tmp_path / "pkg", format="png")

    placement, rendering = artifact.placement, artifact.rendering
    _assert_3857_grid(placement)
    assert placement.nodata == "NaN"
    assert rendering.value_range == [float(np.nanmin(values)), np.nanmax(values)]
    low, high = rendering.value_range
    grid = Affine(*placement.transform)
    shape = (placement.height, placement.width)
    warp = dict(
        src_transform=ROTATED,
        src_crs="EPSG:3006",
        dst_transform=grid,
        dst_crs="EPSG:3857",
        tolerance=0,
    )
    scalars = np.full(shape, np.nan)
    reproject(
        values,
        scalars,
        src_nodata=np.nan,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
        **warp,
    )
    mask = np.zeros(shape, dtype=np.uint8)
    reproject(
        np.isfinite(values).astype(np.uint8),
        mask,
        resampling=Resampling.nearest,
        **warp,
    )
    # Bilinear interpolation reproduces a linear field away from the grid edge
    # and the missing corner.
    u, v = _source_positions(placement, ROTATED)
    interior = (u > 2.0) & (u < values.shape[1] - 1.5)
    interior &= (v > 2.0) & (v < values.shape[0] - 1.5)
    assert interior.sum() > 20
    exact = 100.0 + 3.0 * (u - 0.5) - 2.0 * (v - 0.5)
    assert np.abs(scalars - exact)[interior].max() <= 1e-6 * (high - low)

    image = _png(package, artifact)
    valid = mask == 1
    assert valid.any() and (~valid).any()
    np.testing.assert_array_equal(image[..., 3], mask * 255)
    expected = _colour(scalars, rendering)
    difference = np.abs(image.astype(int) - expected.astype(int))[valid]
    assert difference.max() <= 1
    assert (image[~valid] == 0).all()


def test_png_keeps_the_alpha_of_uint8_rgba(tmp_path):
    rows, cols = np.mgrid[0:4, 0:6]
    data = np.stack(
        [cols * 40, rows * 60, np.full_like(cols, 7), np.where(cols < 3, 255, 128)],
        axis=-1,
    ).astype(np.uint8)
    raster = _raster(data, georef=NORTH_UP)

    package, artifact = _export(raster, tmp_path / "pkg", format="png")

    rendering = artifact.rendering
    assert rendering.bands == [0, 1, 2, 3]
    assert rendering.band is rendering.colormap is rendering.value_range is None
    rows, cols, inside = _nearest(artifact.placement, NORTH_UP, data.shape)
    expected = np.where(inside[..., None], data[rows, cols], 0).astype(np.uint8)
    image = _png(package, artifact)
    np.testing.assert_array_equal(image, expected)
    assert {128, 255} <= set(np.unique(image[..., 3]).tolist())


def test_png_makes_finite_sentinel_pixels_transparent(tmp_path):
    values = _smooth()
    values[:, :8] = -9999.0
    raster = _raster(values, georef=NORTH_UP, nodata=-9999.0)

    package, artifact = _export(raster, tmp_path / "pkg", format="png")

    assert artifact.placement.nodata == -9999.0
    assert artifact.rendering.value_range == [
        values[:, 8:].min(),
        values[:, 8:].max(),
    ]
    rows, cols, inside = _nearest(artifact.placement, NORTH_UP, values.shape)
    valid = inside & (cols >= 8)
    np.testing.assert_array_equal(_png(package, artifact)[..., 3], valid * 255)


def test_png_maps_a_constant_array_to_the_low_colour(tmp_path):
    raster = _raster(np.full((4, 6), 5.0), georef=NORTH_UP)

    package, artifact = _export(raster, tmp_path / "pkg", format="png")

    assert artifact.rendering.value_range == [5.0, 5.0]
    image = _png(package, artifact)
    opaque = image[..., 3] == 255
    assert opaque.any()
    low = matplotlib.colormaps["viridis"](0.0, bytes=True)
    assert (image[opaque] == low).all()


@pytest.mark.parametrize("format", ["png", "tif"])
def test_single_band_data_of_shape_hw_or_hw1_writes_the_same_artifact(tmp_path, format):
    labels = (np.arange(24, dtype=np.int16).reshape(4, 6) - 3) * 7

    written = []
    for name, data in (("hw", labels), ("hw1", labels[..., np.newaxis])):
        package, artifact = _export(_raster(data), tmp_path / name, format=format)
        if format == "png":
            written.append((artifact.rendering, _png(package, artifact)))
        else:
            with rasterio.open(package.path / artifact.path) as src:
                assert (src.count, src.dtypes) == (1, ("int16",))
                written.append((artifact.rendering, src.read(1)))

    (rendering, pixels), (rendering_hw1, pixels_hw1) = written
    assert rendering_hw1 == rendering
    np.testing.assert_array_equal(pixels_hw1, pixels)
    if format == "tif":
        np.testing.assert_array_equal(pixels, labels)


@pytest.mark.parametrize("options", [{}, {"value_range": (0, 1)}])
@pytest.mark.parametrize("low", [2**53, -(2**53) - 1])
def test_png_refuses_integers_beyond_the_exact_float64_range(tmp_path, low, options):
    data = np.full((4, 6), low, dtype=np.int64)
    data[:, 1::2] += 1

    with pytest.raises(ValueError, match=r"within \[-2\*\*53, 2\*\*53\]"):
        _raster(data).export(tmp_path / "pkg", canonical=True, format="png", **options)
    assert not (tmp_path / "pkg").exists()


def test_png_renders_integers_at_the_exact_float64_limit(tmp_path):
    data = np.full((4, 6), 2**53 - 1, dtype=np.int64)
    data[:, 1::2] += 1

    package, artifact = _export(_raster(data), tmp_path / "pkg", format="png")

    assert artifact.rendering.value_range == [2.0**53 - 1, 2.0**53]
    pixels = _png(package, artifact)
    opaque = {tuple(p) for p in pixels[pixels[..., 3] == 255]}
    ends = _colour(np.array(artifact.rendering.value_range), artifact.rendering)
    assert opaque == {tuple(c) for c in ends}


def test_integer_nodata_is_matched_exactly():
    data = np.array([2**53 - 1, 2**53, 2**53 + 1, 2**53], dtype=np.int64)

    np.testing.assert_array_equal(
        valid_pixels(data, float(2**53)), [True, False, True, False]
    )
    np.testing.assert_array_equal(valid_pixels(data, 0.5), [True] * 4)
    np.testing.assert_array_equal(
        valid_pixels(np.array([0, 1], dtype=np.int64), 0.5), [True, True]
    )


def test_png_refuses_integers_that_round_to_the_nodata_value(tmp_path):
    data = np.zeros((4, 6), dtype=np.int64)
    data[0, 0] = 2**60
    data[0, 1] = 2**60 + 1

    with pytest.raises(ValueError, match=r"within \[-2\*\*53, 2\*\*53\]"):
        _raster(data, nodata=float(2**60)).export(
            tmp_path / "pkg", canonical=True, format="png"
        )
    assert not (tmp_path / "pkg").exists()


@pytest.mark.parametrize(
    "data, format, message",
    [
        (np.zeros((0, 3)), "png", "empty"),
        (np.zeros((0, 3)), "tif", "empty"),
        (np.array(np.nan), "png", "empty"),
        (np.array(np.nan), "tif", "empty"),
        (np.full((4, 6), np.nan), "png", "no valid pixels"),
    ],
)
def test_empty_or_all_invalid_rasters_are_refused(tmp_path, data, format, message):
    raster = _raster(data)

    with pytest.raises(ValueError, match=message):
        raster.export(tmp_path / "pkg", canonical=True, format=format)
    assert not (tmp_path / "pkg").exists()


def test_all_invalid_raster_is_written_as_a_nodata_geotiff(tmp_path):
    raster = _raster(np.full((4, 6), np.nan))

    package, artifact = _export(raster, tmp_path / "pkg", format="tif")

    with rasterio.open(package.path / artifact.path) as src:
        assert math.isnan(src.nodata) and np.isnan(src.read(1)).all()


@pytest.mark.parametrize(
    "values, nodata, dtype",
    [(_smooth(), np.nan, "float64"), (np.arange(24).reshape(4, 6), -1, "int16")],
    ids=["float-nan", "int-sentinel"],
)
def test_tif_keeps_the_source_grid_nodata_and_view_hint_unit(
    tmp_path, values, nodata, dtype
):
    values = values.astype(dtype)
    if dtype == "float64":
        values[0, 0] = np.nan
    else:
        values[0, 0] = nodata
    raster = _raster(values, nodata=nodata, dataset=datasets.trees)

    package, artifact = _export(raster, tmp_path / "pkg", format="tif")

    with rasterio.open(package.path / artifact.path) as src:
        assert src.crs == CRS.from_string("EPSG:3006")
        assert src.transform == ROTATED
        assert (src.count, src.dtypes) == (1, (dtype,))
        if dtype == "float64":
            assert math.isnan(src.nodata)
        else:
            assert src.nodata == nodata
        np.testing.assert_array_equal(src.read(1), values)
        assert src.units == ("m",)
        assert src.descriptions == (None,)
    placement = artifact.placement
    assert placement.crs == "EPSG:3006"
    assert placement.transform == list(ROTATED)[:6]
    assert (placement.width, placement.height) == (values.shape[1], values.shape[0])
    h, w = values.shape
    corners = [ROTATED @ c for c in ((0, 0), (w, 0), (0, h), (w, h))]
    xs, ys = zip(*corners, strict=True)
    assert placement.bounds == pytest.approx([min(xs), min(ys), max(xs), max(ys)])
    lon, lat = 17.9, 59.3  # Stockholm area
    west, south, east, north = placement.bounds_wgs84
    assert west < east and south < north
    assert abs(west - lon) < 0.5 and abs(south - lat) < 0.5
    assert placement.nodata == ("NaN" if dtype == "float64" else -1.0)
    assert artifact.rendering.model_dump() == {
        "bands": [0],
        "field": None,
        "band": None,
        "unit": "m",
        "value_range": None,
        "colormap": None,
        "transparent_nodata": False,
    }
    if dtype == "float64":
        text = package.manifest_path.read_text()
        assert '"nodata": "NaN"' in text
        reloaded = load_model_package(package.path)
        assert math.isnan(reloaded.nodata)


@pytest.mark.parametrize("format", ["png", "tif"])
def test_view_hint_unit_is_recorded_and_the_unit_option_overrides_it(tmp_path, format):
    raster = _raster(_smooth(), dataset=datasets.trees)

    _, plain = _export(raster, tmp_path / "plain", format=format)
    package, overridden = _export(raster, tmp_path / "ft", format=format, unit="ft")

    assert plain.rendering.unit == "m"
    assert overridden.rendering.unit == "ft"
    if format == "tif":
        with rasterio.open(package.path / overridden.path) as src:
            assert src.units == ("ft",)


def test_png_without_unit_source_records_none(tmp_path):
    _, artifact = _export(_raster(_smooth()), tmp_path / "pkg", format="png")
    assert artifact.rendering.unit is None


def test_an_empty_unit_option_is_recorded_as_given(tmp_path):
    raster = _raster(_smooth(), dataset=datasets.trees)

    _, artifact = _export(raster, tmp_path / "pkg", format="png", unit="")

    assert artifact.rendering.unit == ""


@pytest.mark.parametrize("name", ["pkg", "pkg.dtccpkg"])
def test_tiff_is_written_as_tif(tmp_path, name):
    package, artifact = _export(_raster(_smooth()), tmp_path / name, format="tiff")

    assert artifact.format == "tif" and artifact.path.endswith(".tif")
    assert artifact.placement.crs == "EPSG:3006"


@pytest.mark.parametrize("name", ["pkg", "pkg.dtccpkg"])
def test_other_image_formats_are_refused_in_canonical_packages(tmp_path, name):
    raster = _raster(np.arange(24, dtype=np.uint8).reshape(4, 6))

    with pytest.raises(ValueError, match='only as format="png" or format="tif"'):
        raster.export(tmp_path / name, canonical=True, format="jpg")
    assert not (tmp_path / name).exists()


def test_png_placement_requires_an_exact_reproject(tmp_path, monkeypatch):
    def reproject(source, destination=None, **kwargs):
        raise AssertionError("must not warp")

    monkeypatch.setattr(rasterio.warp, "reproject", reproject)

    with pytest.raises(RuntimeError, match="reproject accepts tolerance"):
        _raster(_smooth()).export(tmp_path / "pkg", canonical=True, format="png")
    assert not (tmp_path / "pkg").exists()


# Applicability of the options: positive cases


def test_band_selects_a_band_of_a_multiband_raster(tmp_path):
    data = np.stack([_smooth(), -_smooth()], axis=-1)

    _, artifact = _export(_raster(data), tmp_path / "pkg", format="png", band=1)

    assert artifact.rendering.band == 1
    assert artifact.rendering.value_range == [-data[..., 0].max(), -data[..., 0].min()]


def test_colormap_and_value_range_are_used_and_recorded(tmp_path):
    labels = _smooth().astype(np.int32)
    raster = _raster(labels, georef=NORTH_UP)

    package, artifact = _export(
        raster, tmp_path / "pkg", format="png", colormap="magma", value_range=(0, 200)
    )

    rendering = artifact.rendering
    assert (rendering.colormap, rendering.value_range) == ("magma", [0.0, 200.0])
    rows, cols, inside = _nearest(artifact.placement, NORTH_UP, (12, 16))
    expected = np.zeros((*inside.shape, 4), dtype=np.uint8)
    expected[inside] = _colour(labels[rows, cols][inside], rendering)
    np.testing.assert_array_equal(_png(package, artifact), expected)


# Applicability of the options: rejected cases


def _rgba():
    return _raster(np.full((4, 6, 4), 200, dtype=np.uint8))


def _two_bands():
    return _raster(np.stack([_smooth(), _smooth()], axis=-1))


@pytest.mark.parametrize(
    "make, options, message",
    [
        (_rgba, {"format": "png", "band": 1}, "band applies"),
        (_two_bands, {"format": "png", "band": 2}, "band must be"),
        (_two_bands, {"format": "png", "band": -1}, "band must be"),
        (_two_bands, {"format": "png", "band": 1.0}, "band must be"),
        (lambda: _raster(_smooth()), {"format": "png", "band": 0}, "band applies"),
        (_two_bands, {"format": "tif", "band": 0}, 'format="png"'),
        (lambda: _raster(_smooth()), {"format": "tif", "colormap": "magma"}, "png"),
        (lambda: _raster(_smooth()), {"format": "tif", "value_range": (0, 1)}, "png"),
        (_rgba, {"format": "png", "colormap": "magma"}, "RGB"),
        (_rgba, {"format": "png", "value_range": (0, 1)}, "RGB"),
        (lambda: _raster(_smooth()), {"format": "png", "colormap": "nope"}, "colormap"),
        (
            lambda: _raster(_smooth()),
            {"format": "png", "value_range": (10, 0)},
            "value_range",
        ),
        (
            lambda: _raster(_smooth()),
            {"format": "png", "value_range": (0, math.inf)},
            "value_range",
        ),
        (
            lambda: _raster(_smooth()),
            {"format": "png", "value_range": (0, 1, 2)},
            "value_range",
        ),
        (
            lambda: _raster(_smooth()),
            {"format": "png", "value_range": [0, 1, 2]},
            "value_range",
        ),
        (
            lambda: _raster(_smooth()),
            {"format": "png", "value_range": "12"},
            "value_range",
        ),
        (
            lambda: _raster(_smooth()),
            {"format": "png", "value_range": (0, True)},
            "value_range",
        ),
        (lambda: _raster(_smooth()), {"format": "png", "field": "speed"}, "field"),
        (lambda: _raster(_smooth()), {"unit": "m"}, 'format="png"'),
        (lambda: _raster(_smooth()), {"format": "dtcc", "unit": "m"}, 'format="png"'),
        (lambda: _raster(_smooth()), {"format": "png", "unit": 3}, "unit"),
    ],
)
def test_inapplicable_options_are_rejected_before_writing(
    tmp_path, make, options, message
):
    raster = make()

    with pytest.raises(ValueError, match=message):
        raster.export(tmp_path / "pkg", canonical=True, **options)
    assert not (tmp_path / "pkg").exists()


def test_renderer_options_are_rejected_by_legacy_packages(tmp_path):
    raster = _raster(_smooth())

    with pytest.raises(ValueError, match="canonical=True"):
        raster.export(tmp_path / "pkg", format="png", colormap="magma")
    assert not (tmp_path / "pkg").exists()


# Re-export and reader


@pytest.mark.parametrize("format", ["png", "tif"])
def test_reexport_of_a_loaded_package_records_the_same_placement(tmp_path, format):
    options = {"format": format, "unit": "ft"}
    if format == "png":
        options.update(colormap="magma", value_range=(90.0, 140.0))
    raster = _raster(_smooth(), dataset=datasets.trees)
    _, first = _export(raster, tmp_path / "first", **options)

    restored = load_model_package(tmp_path / "first")
    _, second = _export(restored, tmp_path / "second", **options)

    assert second.placement == first.placement
    assert second.rendering == first.rendering
    assert second.sha256 == first.sha256


def _derived(manifest):
    [artifact] = [a for a in manifest["artifacts"] if a["role"] == "derived"]
    return artifact


@pytest.mark.parametrize(
    "change",
    [
        {"width": 0},
        {"height": -1},
        {"transform": [1.0, 0.0, 0.0, 0.0, -1.0]},
        {"bounds_wgs84": "INFINITE"},
        {"bounds": [1.0, 0.0, 0.0, 1.0]},
        {"nodata": "nan"},
        {"crs": ""},
        {"extra": 1},
    ],
)
def test_reader_rejects_a_malformed_placement(tmp_path, change):
    package, _ = _export(_raster(_smooth()), tmp_path / "pkg", format="png")
    manifest = json.loads(package.manifest_path.read_text())
    _derived(manifest)["placement"].update(change)
    text = json.dumps(manifest).replace('"INFINITE"', "[0.0, 0.0, 1e999, 1.0]")
    package.manifest_path.write_text(text)

    with pytest.raises(ValueError):
        load_model_package(package.path)


def test_reader_loads_a_manifest_without_the_records(tmp_path):
    package, _ = _export(_raster(_smooth()), tmp_path / "pkg", format="png")
    manifest = json.loads(package.manifest_path.read_text())
    del _derived(manifest)["placement"]
    del _derived(manifest)["rendering"]
    package.manifest_path.write_text(json.dumps(manifest))

    assert type(load_model_package(package.path)) is Raster


def test_manifest_round_trips_with_and_without_the_records():
    context = _context(datasets.city)
    plain = DatasetArtifact(
        path="artifacts/a.png",
        role="derived",
        format="png",
        media_type="image/png",
        data_kind="raster",
    )
    placed = DatasetArtifact.model_validate(
        {
            **plain.model_dump(),
            "placement": {
                "crs": "EPSG:3857",
                "transform": [2.0, 0.0, 10.0, 0.0, -2.0, 20.0],
                "width": 3,
                "height": 2,
                "bounds": [10.0, 16.0, 16.0, 20.0],
                "bounds_wgs84": [0.0001, 0.0001, 0.0002, 0.0002],
                "nodata": "NaN",
            },
            "rendering": {
                "field": "velocity",
                "band": 0,
                "unit": "m/s",
                "value_range": [0.0, 4.0],
                "colormap": "viridis",
                "transparent_nodata": True,
                "components": ["x", "y"],
            },
        }
    )
    for artifacts in ([plain], [placed]):
        manifest = context.manifest(artifacts=artifacts)
        text = manifest.model_dump_json()
        assert DatasetManifest.model_validate(json.loads(text)) == manifest
        assert DatasetManifest.model_validate(json.loads(text), strict=True) == manifest
    assert "placement" not in plain.model_dump_json()
    assert "rendering" not in plain.model_dump_json()
    assert '"components":["x","y"]' in placed.model_dump_json()
    assert placed.rendering.bands is None


# The to_raster hook


def _field_slice():
    model = FieldSlice()
    attach_dataset_context(model, _context(datasets.smoke))
    return model


def _hook(calls, unit="K", components=None):
    def to_raster(self, field=None):
        calls.append(field)
        raster = Raster(data=_smooth(), georef=ROTATED, crs="EPSG:3006")
        rendering = {"field": field, "unit": unit}
        if components is not None:
            rendering["components"] = components
        return raster, rendering

    return to_raster


def test_hook_receives_only_field_and_its_unit_is_recorded(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        FieldSlice, "to_raster", _hook(calls, components=["x", "y"]), raising=False
    )
    model = _field_slice()

    _, png = _export(
        model,
        tmp_path / "png",
        format="png",
        field="velocity",
        colormap="magma",
        value_range=(0, 1),
    )
    package, tif = _export(model, tmp_path / "tif", format="tif", field="velocity")

    assert calls == ["velocity", "velocity"]
    assert png.rendering.model_dump() == {
        "bands": None,
        "field": "velocity",
        "band": 0,
        "unit": "K",
        "value_range": [0.0, 1.0],
        "colormap": "magma",
        "transparent_nodata": True,
        "components": ["x", "y"],
    }
    assert (tif.rendering.field, tif.rendering.unit) == ("velocity", "K")
    with rasterio.open(package.path / tif.path) as src:
        assert src.units == ("K",) and src.descriptions == ("velocity",)
    assert type(load_model_package(package.path)) is FieldSlice


def test_hook_unit_yields_to_the_unit_option(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(FieldSlice, "to_raster", _hook(calls), raising=False)

    _, artifact = _export(_field_slice(), tmp_path / "pkg", format="png", unit="degC")

    assert calls == [None]
    assert artifact.rendering.unit == "degC"
    assert artifact.rendering.field is None


def test_band_on_a_hook_model_is_rejected(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(FieldSlice, "to_raster", _hook(calls), raising=False)

    with pytest.raises(ValueError, match="band applies"):
        _field_slice().export(tmp_path / "pkg", canonical=True, format="png", band=0)
    assert calls == [] and not (tmp_path / "pkg").exists()


def test_tif_renderer_options_are_rejected_before_the_hook_runs(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(FieldSlice, "to_raster", _hook(calls), raising=False)

    with pytest.raises(ValueError, match='colormap applies to format="png"'):
        _field_slice().export(
            tmp_path / "pkg", canonical=True, format="tif", colormap="magma"
        )
    assert calls == [] and not (tmp_path / "pkg").exists()


@pytest.mark.parametrize(
    "result",
    [
        Raster(data=_smooth(), georef=ROTATED, crs="EPSG:3006"),
        None,
        (Raster(data=_smooth(), georef=ROTATED, crs="EPSG:3006"), {}, {}),
        (Raster(data=_smooth(), georef=ROTATED, crs="EPSG:3006"), None),
    ],
    ids=["bare-raster", "none", "three-tuple", "no-rendering"],
)
def test_a_hook_must_return_a_raster_rendering_pair(tmp_path, monkeypatch, result):
    monkeypatch.setattr(
        FieldSlice, "to_raster", lambda self, field=None: result, raising=False
    )

    with pytest.raises(
        ValueError, match=r"to_raster must return a \(Raster, rendering\) pair"
    ):
        _field_slice().export(tmp_path / "pkg", canonical=True, format="png")
    assert not (tmp_path / "pkg").exists()


@pytest.mark.parametrize("format", ["png", "tif"])
def test_models_without_a_hook_are_refused(tmp_path, format):
    model = StreamlineCollection()
    attach_dataset_context(model, _context(datasets.smoke))

    with pytest.raises(ValueError, match=r'Dataset-level format="png"'):
        model.export(tmp_path / "pkg", canonical=True, format=format)
    assert not (tmp_path / "pkg").exists()
