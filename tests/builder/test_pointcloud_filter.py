import pytest
from pathlib import Path

import dtcc_core
from dtcc_core import io
from dtcc_core.model import Bounds

from dtcc_core.builder.pointcloud.filter import crop


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


def test_crop_no_mut(point_cloud):
    """Test that cropping does not mutate the original point cloud."""
    bounds = Bounds(-2, -2, 0, 0)
    c = point_cloud.crop(bounds)

    assert len(point_cloud) == 8148
    assert len(point_cloud.classification) == 8148
    assert len(c) == 64


def test_crop_copy(point_cloud):
    bounds = Bounds(-2, -2, 0, 0)
    cropped_pc = crop(point_cloud, bounds)

    assert len(point_cloud) == 8148
    assert len(cropped_pc) == 64
    assert len(cropped_pc.classification) == 64


if __name__ == "__main__":
    pytest.main()
