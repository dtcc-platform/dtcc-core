import pytest
from pathlib import Path
from dtcc_core import io
from dtcc_core.model import Raster
from dtcc_core.builder.model_conversion import raster_to_builder_gridfield


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / ".." / "data" / "rasters").resolve()


@pytest.fixture
def dem_raster(data_dir):
    return io.load_raster(data_dir / "test_dem.tif")


@pytest.fixture
def dem_gridfield(dem_raster):
    return raster_to_builder_gridfield(dem_raster)


def test_gridfield_dimensions(dem_gridfield):
    assert dem_gridfield.grid.xsize == 20
    assert dem_gridfield.grid.ysize == 40
    assert dem_gridfield.grid.xstep == 2.0
    assert dem_gridfield.grid.ystep == 2.0


def test_gridfield_values(dem_gridfield):
    values = dem_gridfield.values
    assert len(values) == 20 * 40

    # Check specific points
    assert pytest.approx(values[0], rel=1e-5) == 76  # bottom left corner
    assert pytest.approx(values[1], rel=1e-5) == 76.2  # next to bottom left
    assert pytest.approx(values[(20 * 40) - 20], rel=1e-5) == 0  # top left corner
    assert pytest.approx(values[-1], rel=1e-5) == 3.8  # top right corner


if __name__ == "__main__":
    pytest.main()
