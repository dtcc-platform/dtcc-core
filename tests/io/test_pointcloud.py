import pytest
from pathlib import Path
import json
import tempfile
from dtcc_core import io
from dtcc_core.model import Bounds, PointCloud


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / ".." / "data" / "MinimalCase").resolve()


@pytest.fixture
def las_file(data_dir):
    return str((data_dir / "pointcloud.las").resolve())


@pytest.fixture
def point_cloud(las_file):
    return io.load_pointcloud(las_file)


@pytest.fixture
def expected_bounds():
    return {"xmin": -8.01747, "ymin": -18.850332, "xmax": 15.92373, "ymax": 1.83826}


def test_load_pointcloud(point_cloud):
    assert isinstance(point_cloud, PointCloud)
    assert len(point_cloud.points) == 8148
    assert len(point_cloud.classification) == 8148
    assert len(point_cloud.used_classifications()) == 2


def test_load_pointcloud_from_dir(data_dir):
    pc = io.load_pointcloud(data_dir)
    assert len(pc.points) == 8148


def test_load_pointcloud_bounded(las_file):
    pc = io.load_pointcloud(las_file, bounds=Bounds(-2, -2, 0, 0))
    assert len(pc.points) == 64
    assert len(pc.classification) == 64


def test_point_cloud_bounds(las_file, expected_bounds):
    bounds = io.pointcloud.calc_las_bounds(las_file)
    assert pytest.approx(bounds.xmin, rel=1e-3) == expected_bounds["xmin"]
    assert pytest.approx(bounds.ymin, rel=1e-3) == expected_bounds["ymin"]
    assert pytest.approx(bounds.xmax, rel=1e-3) == expected_bounds["xmax"]
    assert pytest.approx(bounds.ymax, rel=1e-3) == expected_bounds["ymax"]


def test_save_pointcloud(point_cloud):
    with tempfile.NamedTemporaryFile(suffix=".las", delete=False) as outfile:
        outpath = Path(outfile.name)

    try:
        point_cloud.save(outpath)
        loaded_pc = io.load_pointcloud(outpath)
        assert len(point_cloud.points) == len(loaded_pc.points)
        assert len(point_cloud.classification) == len(loaded_pc.classification)
    finally:
        outpath.unlink()


def test_load_pointcloud_list(las_file):
    pc = io.load_pointcloud([las_file, las_file])
    assert len(pc.points) == 8148 * 2
    assert len(pc.classification) == 8148 * 2


if __name__ == "__main__":
    pytest.main()
