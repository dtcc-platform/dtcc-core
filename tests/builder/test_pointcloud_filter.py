import pytest
import numpy as np
from pathlib import Path

import dtcc_core
from dtcc_core import io
from dtcc_core.model import Bounds, PointCloud, Surface

from dtcc_core.builder.pointcloud.filter import (
    crop, find_points_in_polygons, remove_points_in_polygons, points_in_polygons,
)


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


@pytest.fixture
def polygon_selection():
    pc = PointCloud(
        points=np.array([
            [20, 20, 0], [13, 11, 5], [13, 11, 5],
            [12, 12, 100], [15, 12, -100], [11, 11, 0],
        ], dtype=float),
        classification=np.arange(6, dtype=np.uint8),
        intensity=np.arange(10, 16, dtype=np.uint16),
        return_number=np.array([1, 2, 1, 2, 1, 2], dtype=np.uint8),
        num_returns=np.full(6, 2, dtype=np.uint8),
    )
    first = Surface(
        vertices=np.array([[10, 10, 0], [14, 10, 0], [14, 14, 0], [10, 14, 0]]),
        holes=[np.array([[10.5, 10.5, 0], [11.5, 10.5, 0],
                         [11.5, 11.5, 0], [10.5, 11.5, 0]])],
    )
    second = Surface(vertices=np.array([
        [12.5, 10, 0], [16, 10, 0], [16, 14, 0], [12.5, 14, 0],
    ]))
    return pc, [first, None, second]


def test_polygon_filters_preserve_records_attributes_and_input(polygon_selection):
    pc, polygons = polygon_selection
    before = pc.copy()
    for operation, expected in (
        (find_points_in_polygons, [1, 2, 3, 4]),
        (remove_points_in_polygons, [0, 5]),
    ):
        result = operation(pc, polygons)
        for attribute in ("points", "classification", "intensity", "return_number", "num_returns"):
            original = getattr(before, attribute)
            selected = getattr(result, attribute)
            np.testing.assert_array_equal(selected, original[expected])
            np.testing.assert_array_equal(getattr(pc, attribute), original)
            assert not np.shares_memory(selected, getattr(pc, attribute))


def test_polygon_coordinate_helper_keeps_its_coordinate_api(polygon_selection):
    pc, polygons = polygon_selection
    groups = points_in_polygons(pc, polygons, flatten=False)
    assert len(groups) == 2
    for group, expected in zip(groups, ([1, 2, 3], [1, 2, 4])):
        assert group.shape == (3, 3)
        assert sorted(map(tuple, group)) == sorted(map(tuple, pc.points[expected]))
    # Flattened coordinate output retains its existing coordinate deduplication.
    np.testing.assert_array_equal(
        points_in_polygons(pc, polygons), np.unique(pc.points[[1, 2, 3, 4]], axis=0)
    )


@pytest.mark.parametrize("selection", ["no_polygons", "no_matches", "no_points"])
def test_polygon_selection_empty_results(polygon_selection, selection):
    pc, polygons = polygon_selection
    if selection == "no_polygons":
        polygons = [None]
    elif selection == "no_matches":
        pc.points[:, :2] += 100
    else:
        pc = PointCloud()
    kept = find_points_in_polygons(pc, polygons)
    removed = remove_points_in_polygons(pc, polygons)
    assert kept.points.shape == (0, 3)
    np.testing.assert_array_equal(removed.points, pc.points)
    assert points_in_polygons(pc, polygons).shape == (0, 3)
    groups = points_in_polygons(pc, polygons, flatten=False)
    assert len(groups) == (0 if selection == "no_polygons" else 2)
    assert all(group.shape == (0, 3) for group in groups)


@pytest.mark.parametrize("points", [np.zeros((2, 2)), np.array([[0, 0, np.nan]])])
def test_native_polygon_selection_rejects_invalid_points(points):
    from dtcc_core.builder import _dtcc_builder

    with pytest.raises(ValueError, match="points"):
        _dtcc_builder.points_in_polygons(points, [])


if __name__ == "__main__":
    pytest.main()
