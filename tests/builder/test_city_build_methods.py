import pytest

import dtcc_core
from dtcc_core.model import City, Bounds, PointCloud


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


def test_build_terrain(las_path):
    city = City()
    city.load_pointcloud(las_path)
    city.build_terrain()

    assert city.terrain is not None
    assert city.terrain.raster is not None
    assert city.terrain.mesh is not None
