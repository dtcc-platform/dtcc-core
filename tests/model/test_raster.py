import pytest

import numpy as np
from affine import Affine
from dtcc_core.model import Raster


def test_create_empty_raster():
    raster = Raster()
    assert raster.crs == ""
    assert raster.bounds.tuple == (0, 0, 0, 0)
    assert raster.georef.to_gdal() == (0, 1, 0, 0, 0, 1)


def test_copy_raster():
    raster = Raster()
    raster.data = np.ones((10, 10), dtype=np.uint8)
    raster.crs = "EPSG:3857"

    copy_raster = raster.copy()

    assert copy_raster.data.tolist() == raster.data.tolist()
    assert copy_raster.crs == raster.crs


def test_copy_raster_nodata():
    raster = Raster()
    raster.data = np.ones((10, 10), dtype=np.uint8)
    raster.crs = "EPSG:3857"

    copy_raster = raster.copy(no_data=True)

    assert copy_raster.data.shape == ()
    assert copy_raster.crs == raster.crs


def test_raster_protobuf_roundtrip_2d():
    raster = Raster()
    raster.data = np.arange(6, dtype=np.float32).reshape((2, 3))
    raster.nodata = -9999.0
    raster.crs = "EPSG:3006"
    raster.georef = Affine(2.0, 0.0, 10.0, 0.0, -2.0, 20.0)

    raster2 = Raster()
    raster2.from_proto(raster.to_proto())

    assert np.allclose(raster2.data, raster.data)
    assert raster2.data.dtype == raster.data.dtype
    assert raster2.nodata == raster.nodata
    assert raster2.crs == raster.crs
    assert raster2.georef == raster.georef


def test_raster_protobuf_roundtrip_3d():
    raster = Raster()
    raster.data = np.arange(12, dtype=np.uint8).reshape((2, 3, 2))
    raster.nodata = 255.0
    raster.crs = "EPSG:3857"
    raster.georef = Affine(1.0, 0.1, 5.0, 0.2, -1.0, 8.0)

    raster2 = Raster()
    raster2.from_proto(raster.to_proto().SerializeToString())

    assert np.array_equal(raster2.data, raster.data)
    assert raster2.data.dtype == raster.data.dtype
    assert raster2.nodata == raster.nodata
    assert raster2.crs == raster.crs
    assert raster2.georef == raster.georef


if __name__ == "__main__":
    pytest.main()
