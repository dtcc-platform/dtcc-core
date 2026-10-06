import pytest
import numpy as np
import os
import tempfile
from pathlib import Path
from dtcc_core import io


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / ".." / "data" / "rasters").resolve()


@pytest.fixture
def dem_raster_path(data_dir):
    return data_dir / "test_dem.tif"


@pytest.fixture
def rgb_image_path(data_dir):
    return data_dir / "14040.png"


@pytest.fixture
def dem_raster(dem_raster_path):
    return io.load_raster(dem_raster_path)


@pytest.fixture
def rgb_image(rgb_image_path):
    return io.load_raster(rgb_image_path)


@pytest.fixture
def test_raster_paths(data_dir):
    return [
        data_dir / "testraster_0_0.tif",
        data_dir / "testraster_0_1.tif",
        data_dir / "testraster_1_0.tif",
        data_dir / "testraster_1_1.tif",
    ]


def test_load_dem_dimensions(dem_raster):
    assert dem_raster.width == 20
    assert dem_raster.height == 40
    assert dem_raster.channels == 1


def test_load_image_dimensions(rgb_image):
    assert rgb_image.width == 228
    assert rgb_image.height == 230
    assert rgb_image.channels == 3


@pytest.mark.parametrize(
    "raster, expected_cell_size",
    [("dem_raster", (2.0, -2.0)), ("rgb_image", (0.08, -0.08))],
)
def test_get_cell_size(request, raster, expected_cell_size):
    raster_obj = request.getfixturevalue(raster)
    assert raster_obj.cell_size == expected_cell_size


def test_write_elevation_model(dem_raster):
    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_file:
        outfile = tmp_file.name

    try:
        dem_raster.save(outfile)
        loaded_em = io.load_raster(outfile)

        assert loaded_em.height == 40
        assert loaded_em.width == 20
        assert loaded_em.cell_size == (2.0, -2.0)
    finally:
        os.unlink(outfile)


def test_load_multiple(data_dir, test_raster_paths):
    # Load and check multiple rasters
    combined_raster = io.load_raster(test_raster_paths)
    assert combined_raster.width == 50
    assert combined_raster.height == 50
    assert combined_raster.cell_size == (0.5, -0.5)

    # Compare with reference raster
    reference_raster = io.load_raster(data_dir / "testraster.tif")
    assert np.all(combined_raster.data == reference_raster.data)


# GeoTIFF writer

import math
import tracemalloc

import rasterio
from affine import Affine
from rasterio.enums import ColorInterp, MaskFlags

from dtcc_core.io import raster as raster_io
from dtcc_core.model import Raster


def samples(shape, seed=0):
    return np.random.default_rng(seed).integers(1, 255, shape, dtype=np.uint8)


def georeferenced(data, **kwargs):
    kwargs.setdefault("georef", Affine(0.5, 0.0, 672500.0, 0.0, -0.5, 6577500.0))
    return Raster(data=data, **kwargs)


@pytest.mark.parametrize("channels", [3, 4])
def test_multiband_raster_is_written_band_first(tmp_path, channels):
    data = samples((5, 7, channels))
    path = tmp_path / "out.tif"
    georeferenced(data).save(path)
    with rasterio.open(path) as src:
        assert (src.count, src.height, src.width) == (channels, 5, 7)
        np.testing.assert_array_equal(src.read(), data.transpose(2, 0, 1))
        assert src.transform == Affine(0.5, 0.0, 672500.0, 0.0, -0.5, 6577500.0)
    np.testing.assert_array_equal(io.load_raster(path).data, data)


@pytest.mark.parametrize("shape", [(5, 7), (5, 7, 3)])
def test_crs_and_nodata_are_written(tmp_path, shape):
    path = tmp_path / "out.tif"
    georeferenced(samples(shape), crs="EPSG:3006", nodata=0.0).save(path)
    with rasterio.open(path) as src:
        assert src.crs.to_epsg() == 3006
        assert set(src.nodatavals) == {0.0}


def test_nan_nodata_and_empty_crs_are_not_written(tmp_path):
    path = tmp_path / "out.tif"
    georeferenced(samples((5, 7, 3))).save(path)
    with rasterio.open(path) as src:
        assert src.crs is None
        assert set(src.nodatavals) == {None}


@pytest.mark.parametrize("channels", [2, 4])
def test_alpha_is_written_only_when_requested(tmp_path, channels):
    data = samples((5, 7, channels))
    with_alpha, without = tmp_path / "alpha.tif", tmp_path / "plain.tif"
    georeferenced(data).save(with_alpha, alpha=True)
    georeferenced(data).save(without)
    with rasterio.open(with_alpha) as src:
        assert src.colorinterp[-1] == ColorInterp.alpha
        assert ColorInterp.alpha not in src.colorinterp[:-1]
        np.testing.assert_array_equal(src.read(), data.transpose(2, 0, 1))
    with rasterio.open(without) as src:
        assert ColorInterp.alpha not in src.colorinterp
        assert all(flags == [MaskFlags.all_valid] for flags in src.mask_flag_enums)


@pytest.mark.parametrize("channels", [1, 3, 5])
def test_alpha_needs_two_or_four_bands(tmp_path, channels):
    shape = (5, 7) if channels == 1 else (5, 7, channels)
    path = tmp_path / "out.tif"
    with pytest.raises(ValueError, match="alpha"):
        georeferenced(samples(shape)).save(path, alpha=True)
    assert not path.exists()


def test_missing_crs_as_loaded_is_not_written(tmp_path):
    # The loader stores a file without a CRS as the text "None".
    path = tmp_path / "out.tif"
    georeferenced(samples((5, 7)), crs="None").save(path)
    with rasterio.open(path) as src:
        assert src.crs is None


@pytest.fixture
def written(monkeypatch):
    """Record the size of every array handed to a GeoTIFF write."""
    sizes = []
    real_write = rasterio.io.DatasetWriter.write

    def recording_write(self, arr, *args, **kwargs):
        sizes.append(np.asarray(arr).nbytes)
        return real_write(self, arr, *args, **kwargs)

    monkeypatch.setattr(rasterio.io.DatasetWriter, "write", recording_write)
    return sizes


def test_strips_respect_the_cap_even_for_wide_rows(tmp_path, monkeypatch, written):
    data = samples((6, 20, 4))
    whole, split = tmp_path / "whole.tif", tmp_path / "split.tif"
    georeferenced(data).save(whole, alpha=True)
    written.clear()
    # Smaller than one 20-pixel row of one band.
    monkeypatch.setattr(raster_io, "WRITE_STRIP_BYTES", 7)
    georeferenced(data).save(split, alpha=True)
    assert written and max(written) <= 7
    with rasterio.open(whole) as a, rasterio.open(split) as b:
        np.testing.assert_array_equal(a.read(), b.read())
        assert a.colorinterp == b.colorinterp


@pytest.mark.parametrize("cap", [None, 4 * 20])
def test_each_window_writes_every_band_once(tmp_path, monkeypatch, cap):
    indexes = []
    real_write = rasterio.io.DatasetWriter.write

    def recording_write(self, arr, *args, **kwargs):
        indexes.append(list(args[0] if args else kwargs["indexes"]))
        return real_write(self, arr, *args, **kwargs)

    monkeypatch.setattr(rasterio.io.DatasetWriter, "write", recording_write)
    if cap is not None:
        # One 20-pixel row of all four bands per window.
        monkeypatch.setattr(raster_io, "WRITE_STRIP_BYTES", cap)
    georeferenced(samples((6, 20, 4))).save(tmp_path / "out.tif", alpha=True)
    assert indexes == [[1, 2, 3, 4]] * (1 if cap is None else 6)


def test_small_gdal_cache_does_not_inflate_the_file(tmp_path):
    raster = georeferenced(samples((256, 256, 4)), crs="EPSG:3006")
    small, large = tmp_path / "small.tif", tmp_path / "large.tif"
    # GDAL reads cache sizes below 100000 as megabytes.
    with rasterio.Env(GDAL_CACHEMAX=200_000):
        raster.save(small, alpha=True)
    with rasterio.Env(GDAL_CACHEMAX=512 * 1024**2):
        raster.save(large, alpha=True)
    # Strips rewritten once per band would make the file about twice as large.
    assert small.stat().st_size <= large.stat().st_size * 1.01
    with rasterio.open(small) as a, rasterio.open(large) as b:
        np.testing.assert_array_equal(a.read(), b.read())


def test_writing_allocates_at_most_one_strip(tmp_path, monkeypatch):
    data = samples((1024, 1024, 4))
    cap = 1024**2
    monkeypatch.setattr(raster_io, "WRITE_STRIP_BYTES", cap)
    raster = georeferenced(data, crs="EPSG:3006")
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        raster.save(tmp_path / "out.tif", alpha=True)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    # A second strip alive at the same time would exceed the small allowance.
    assert peak <= cap + 128 * 1024
