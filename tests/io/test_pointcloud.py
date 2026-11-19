import tempfile
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import builder, io
from dtcc_core.io.info import info_pointcloud
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


def test_save_pointcloud_csv_round_trips_points_and_classification(tmp_path):
    rng = np.random.default_rng(0)
    points = rng.random((10, 3)) * 1000.0
    classification = np.arange(10, dtype=np.uint8) % 3 + 1
    pc = PointCloud(points=points, classification=classification)
    outpath = tmp_path / "cloud.csv"

    io.save_pointcloud(pc, outpath)
    loaded = io.load_pointcloud(outpath)

    np.testing.assert_array_equal(loaded.points, points)
    np.testing.assert_array_equal(loaded.classification, classification)


def test_save_pointcloud_csv_without_classification_writes_xyz(tmp_path):
    points = np.array([[1.5, 2.5, 3.5]])
    outpath = tmp_path / "cloud.csv"

    io.save_pointcloud(PointCloud(points=points), outpath)
    loaded = io.load_pointcloud(outpath)

    assert np.loadtxt(outpath, delimiter=",", ndmin=2).shape == (1, 3)
    np.testing.assert_array_equal(loaded.points, points)


def test_load_pointcloud_csv_points_only_ignores_classification_column(tmp_path):
    outpath = tmp_path / "cloud.csv"
    np.savetxt(outpath, [[0.0, 0.0, 0.0, 6.0], [1.0, 1.0, 1.0, 2.0]], delimiter=",")

    assert list(io.load_pointcloud(outpath).classification) == [6, 2]
    assert list(io.load_pointcloud(outpath, points_only=True).classification) == [
        1,
        1,
    ]


@pytest.mark.parametrize(
    ("text", "expected"),
    [("1,2,3\n", (1, 1.0, 1.0)), ("1,2,3\n4,5,6\n", (2, 1.0, 4.0))],
)
def test_info_pointcloud_csv(tmp_path, text, expected):
    path = tmp_path / "points.csv"
    path.write_text(text)
    info = info_pointcloud(path)
    assert (info["points"], info["x_min"], info["x_max"]) == expected


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
