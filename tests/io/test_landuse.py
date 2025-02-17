import pytest
from pathlib import Path
from dtcc_core.model import Landuse, LanduseClasses
from dtcc_core import io


@pytest.fixture
def test_data_path():
    path = Path(__file__).parent / ".." / "data" / "landuse" / "landuse_testdata.shp"
    assert path.is_file()
    return path


@pytest.fixture
def landuse_data(test_data_path):
    return io.load_landuse(test_data_path)


def test_load_landuse_type(landuse_data):
    assert isinstance(landuse_data, Landuse)


def test_landuse_surfaces(landuse_data):
    assert len(landuse_data.surfaces) == 79


def test_landuse_classes(landuse_data):
    assert len(landuse_data.landuses) == 79
    assert all(isinstance(landuse, LanduseClasses) for landuse in landuse_data.landuses)


if __name__ == "__main__":
    pytest.main()
