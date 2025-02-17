import pytest
from pathlib import Path
from dtcc_core import io
from dtcc_core.model import Bounds


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / ".." / "data").resolve()


@pytest.fixture
def pointcloud_dir(data_dir):
    return data_dir / "pointclouds"


@pytest.fixture
def pointcloud_directory(pointcloud_dir):
    return io.load_pointcloud_directory(pointcloud_dir)


@pytest.mark.parametrize(
    "bounds, expected_count",
    [
        (None, 20),  # Get all points
        (Bounds(xmin=0.1, xmax=3.1, ymin=0.1, ymax=1.1), 3),  # Points from 1 file
        (Bounds(xmin=0.1, xmax=3.1, ymin=-1, ymax=3), 6),  # Points from 2 files
    ],
)
def test_get_points_with_bounds(pointcloud_directory, bounds, expected_count):
    if bounds is None:
        bounds = pointcloud_directory.bounds
    pc = pointcloud_directory.pointcloud(bounds)
    assert len(pc) == expected_count


def test_load_pointcloud_directory(pointcloud_directory):
    assert len(pointcloud_directory) == 2


def test_pointcloud_directory_bounds(pointcloud_directory):
    bounds = pointcloud_directory.bounds
    assert bounds.ymin == 0
    assert bounds.ymax == 1


if __name__ == "__main__":
    pytest.main()
