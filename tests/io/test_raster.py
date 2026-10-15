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
