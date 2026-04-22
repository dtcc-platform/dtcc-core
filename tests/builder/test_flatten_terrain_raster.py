import numpy as np
from affine import Affine

import dtcc_core.builder as builder
from dtcc_core.builder.geometry_builders.meshes import _is_effectively_flat_raster
from dtcc_core.model import Raster


def test_flatten_terrain_raster_uses_min_valid_value_and_preserves_metadata():
    raster = Raster()
    raster.data = np.array([[5.0, -9999.0], [3.0, 7.0]])
    raster.nodata = -9999.0
    raster.georef = Affine.translation(10.0, 20.0) * Affine.scale(2.0, -2.0)
    raster.crs = "EPSG:3006"

    flat_raster = builder.flatten_terrain_raster(raster)

    assert flat_raster is not raster
    assert flat_raster.georef == raster.georef
    assert flat_raster.crs == raster.crs
    assert np.array_equal(flat_raster.data, np.array([[3.0, -9999.0], [3.0, 3.0]]))
    assert np.array_equal(raster.data, np.array([[5.0, -9999.0], [3.0, 7.0]]))


def test_flatten_terrain_raster_accepts_explicit_height_and_preserves_nan_nodata():
    raster = Raster()
    raster.data = np.array([[np.nan, 4.0], [5.0, 6.0]])
    raster.nodata = np.nan

    flat_raster = builder.flatten_terrain_raster(raster, height=11.5)

    assert np.isnan(flat_raster.data[0, 0])
    assert np.array_equal(flat_raster.data[0, 1:], np.array([11.5]))
    assert np.array_equal(flat_raster.data[1], np.array([11.5, 11.5]))


def test_is_effectively_flat_raster_detects_uniform_valid_values():
    flat_raster = Raster()
    flat_raster.data = np.array([[np.nan, 3.0], [3.0, 3.0]])
    flat_raster.nodata = np.nan

    varying_raster = Raster()
    varying_raster.data = np.array([[3.0, 3.0], [3.0, 3.01]])
    varying_raster.nodata = np.nan

    assert _is_effectively_flat_raster(flat_raster)
    assert not _is_effectively_flat_raster(varying_raster)
