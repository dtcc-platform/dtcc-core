import pytest

from dtcc_core import io
from dtcc_core.model import City

from pathlib import Path


@pytest.fixture
def data_dir():
    return Path(__file__).parent / ".." / "data"


@pytest.fixture
def building_shp_path(data_dir):
    return str((data_dir / "MinimalCase" / "PropertyMap.shp").resolve())


@pytest.fixture
def las_path(data_dir):
    return str((data_dir / "MinimalCase" / "pointcloud.las").resolve())


@pytest.fixture
def raster_path(data_dir):
    return str((data_dir / "rasters" / "test_dem.tif").resolve())


def test_load_footprints(building_shp_path):
    city = City()
    city.load_footprints(building_shp_path)
    assert len(city.buildings) == 5


def test_load_pointcloud(las_path):
    city = City()
    city.load_pointcloud(las_path)
    assert len(city.pointcloud) == 8148  # Assuming the point cloud is not empty


def test_load_terrain_raster(raster_path):
    city = City()
    city.load_terrain_raster(raster_path)
    assert city.terrain is not None
    assert city.terrain.bounds.area > 0
    assert city.terrain.raster is not None
