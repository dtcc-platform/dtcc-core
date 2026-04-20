import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.builder.geometry_builders.roof_config import RoofDetectionConfig
from dtcc_core.builder.geometry_builders.roof_detection import (
    detect_roof_planes,
    filter_roof_points,
)


def _flat_pc(n=200):
    np.random.seed(0)
    pts = np.column_stack([
        np.random.rand(n) * 10,
        np.random.rand(n) * 10,
        np.full(n, 5.0) + np.random.randn(n) * 0.1,
    ])
    return PointCloud(points=pts)


def test_filter_roof_points_records_normals_and_filter_timings():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    timings = {}
    filtered, normals = filter_roof_points(pc, cfg, stage_timings=timings)
    assert filtered is not None
    assert "normals" in timings
    assert "filter" in timings
    assert timings["normals"] >= 0.0
    assert timings["filter"] >= 0.0


def test_filter_roof_points_none_stage_timings_is_noop():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    filtered_a, _ = filter_roof_points(pc, cfg)
    filtered_b, _ = filter_roof_points(pc, cfg, stage_timings=None)
    assert len(filtered_a.points) == len(filtered_b.points)


def test_detect_roof_planes_records_ransac_and_region_growing():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    filtered, normals = filter_roof_points(pc, cfg)
    timings = {}
    planes = detect_roof_planes(filtered, normals, cfg, stage_timings=timings)
    assert len(planes) >= 1
    assert "ransac" in timings
    assert "region_growing" in timings
    assert timings["ransac"] >= 0.0
    assert timings["region_growing"] >= 0.0


def test_detect_roof_planes_none_stage_timings_is_noop():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    filtered, normals = filter_roof_points(pc, cfg)
    planes_a = detect_roof_planes(filtered, normals, cfg)
    planes_b = detect_roof_planes(filtered, normals, cfg, stage_timings=None)
    assert len(planes_a) == len(planes_b)


def test_filter_roof_points_records_both_timings_even_when_post_filter_insufficient():
    # Dense points on a vertical YZ plane (x~0). Normal estimation will yield
    # horizontal normals (pointing ~±X), so every point will be dropped by the
    # facade filter (angle_to_z ≈ 90° > facade_angle_threshold=75°). The
    # post-filter count goes to 0, triggering the (None, None) branch. Timings
    # for both 'normals' and 'filter' must be recorded before returning.
    np.random.seed(1)
    n = 30  # above min_roof_points (20)
    pts = np.column_stack([
        np.full(n, 0.0) + np.random.randn(n) * 0.001,
        np.random.rand(n) * 5.0,
        np.random.rand(n) * 5.0,
    ])
    pc = PointCloud(points=pts)
    cfg = RoofDetectionConfig()
    timings = {}
    filtered, normals = filter_roof_points(pc, cfg, stage_timings=timings)
    assert filtered is None
    assert "normals" in timings
    assert "filter" in timings


def test_detect_roof_planes_seeds_zero_timings_when_too_few_points():
    pc = PointCloud(points=np.zeros((5, 3)))  # below min_plane_points (20)
    normals = np.tile(np.array([[0.0, 0.0, 1.0]]), (5, 1))
    cfg = RoofDetectionConfig()
    timings = {}
    planes = detect_roof_planes(pc, normals, cfg, stage_timings=timings)
    assert planes == []
    assert timings.get("ransac") == 0.0
    assert timings.get("region_growing") == 0.0
