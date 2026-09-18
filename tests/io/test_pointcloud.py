import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import builder, io
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


def test_load_pointcloud_from_dir_uses_stable_filename_order(monkeypatch, tmp_path):
    (tmp_path / "tile_b.laz").write_bytes(b"")
    (tmp_path / "tile_a.las").write_bytes(b"")
    captured = {}

    def fake_load_list(path, **kwargs):
        captured["path"] = path
        return PointCloud()

    monkeypatch.setattr(io.pointcloud, "load_list", fake_load_list)

    io.pointcloud.load(tmp_path)

    assert [path.name for path in captured["path"]] == [
        "tile_a.las",
        "tile_b.laz",
    ]


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


@pytest.mark.parametrize("file_count", [1, 2])
def test_load_pointcloud_list_preserves_attributes_for_terrain(las_file, file_count):
    original = io.load_pointcloud(las_file)
    pc = io.load_pointcloud([las_file] * file_count)
    assert len(pc.points) == len(original.points) * file_count
    for name in ("classification", "intensity", "return_number", "num_returns"):
        expected = getattr(original, name)
        actual = getattr(pc, name)
        assert actual.dtype == expected.dtype
        np.testing.assert_array_equal(actual, np.tile(expected, file_count))

    raster = builder.build_terrain_raster(
        pc.remove_global_outliers(3.0),
        cell_size=2.0,
        ground_only=True,
        _report_progress=False,
    )
    assert raster.data.size > 0
    assert np.isfinite(raster.data).all()


if __name__ == "__main__":
    pytest.main()
