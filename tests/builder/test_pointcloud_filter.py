import pytest
from pathlib import Path

import dtcc_core
from dtcc_core import io
from dtcc_core.model import Bounds


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / ".." / "data" / "MinimalCase").resolve()


@pytest.fixture
def las_file(data_dir):
    return str((data_dir / "pointcloud.las").resolve())


@pytest.fixture
def point_cloud(las_file):
    return io.load_pointcloud(las_file)


def test_crop_nothing(point_cloud):
    """Test that cropping with original bounds preserves all points."""
    bounds = point_cloud.bounds
    cropped_pc = point_cloud.crop(bounds)

    assert len(cropped_pc) == 8148
    assert len(cropped_pc.classification) == 8148


def test_crop(point_cloud):
    """Test that cropping with smaller bounds reduces point count appropriately."""
    bounds = Bounds(-2, -2, 0, 0)
    cropped_pc = point_cloud.crop(bounds)

    assert len(point_cloud) == 8148
    assert len(cropped_pc) == 64
    assert len(cropped_pc.classification) == 64


if __name__ == "__main__":
    pytest.main()
