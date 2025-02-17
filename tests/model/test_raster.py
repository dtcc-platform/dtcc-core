import pytest

import numpy as np
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


if __name__ == "__main__":
    pytest.main()
