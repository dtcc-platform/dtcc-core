import numpy as np
import open3d as o3d
import pytest

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.builder.geometry_builders.roof_config import RoofDetectionConfig
from dtcc_core.builder.geometry_builders.roof_detection import (
    filter_roof_points,
    detect_roof_planes,
)


def _make_flat_roof_points(n=50, z=10.0, noise=0.1):
    xy = np.random.rand(n, 2) * 10
    z_vals = np.full(n, z) + np.random.randn(n) * noise
    return np.column_stack([xy, z_vals])


def _make_facade_points(n=20, x=0.0):
    yz = np.random.rand(n, 2) * 10
    x_vals = np.full(n, x) + np.random.randn(n) * 0.05
    return np.column_stack([x_vals, yz])


def _make_gabled_roof_points(n_per_side=100, ridge_height=12.0, eave_height=10.0):
    pts = []
    for _ in range(n_per_side):
        x = np.random.rand() * 10
        y = np.random.rand() * 5
        z = eave_height + (ridge_height - eave_height) * (y / 5.0)
        pts.append([x, y, z])
    for _ in range(n_per_side):
        x = np.random.rand() * 10
        y = 5 + np.random.rand() * 5
        z = ridge_height - (ridge_height - eave_height) * ((y - 5) / 5.0)
        pts.append([x, y, z])
    return np.array(pts) + np.random.randn(n_per_side * 2, 3) * 0.02


# --- filter_roof_points tests ---


def test_filter_removes_facade_points():
    roof_pts = _make_flat_roof_points(50, z=10.0)
    facade_pts = _make_facade_points(20, x=0.0)
    all_pts = np.vstack([roof_pts, facade_pts])
    pc = PointCloud(points=all_pts)
    config = RoofDetectionConfig()

    filtered_pc, normals = filter_roof_points(pc, config)

    assert filtered_pc is not None
    assert len(filtered_pc.points) < len(all_pts)
    assert len(filtered_pc.points) >= 30


def test_filter_uses_classification_when_available():
    pts = _make_flat_roof_points(50, z=10.0)
    classification = np.array([6] * 25 + [2] * 25)
    pc = PointCloud(points=pts, classification=classification)
    config = RoofDetectionConfig()

    filtered_pc, _ = filter_roof_points(pc, config)

    assert filtered_pc is not None
    assert len(filtered_pc.points) <= 30


def test_filter_returns_normals():
    pts = _make_flat_roof_points(50, z=10.0)
    pc = PointCloud(points=pts)
    config = RoofDetectionConfig()

    filtered_pc, normals = filter_roof_points(pc, config)

    assert normals is not None
    assert normals.shape[0] == len(filtered_pc.points)
    assert normals.shape[1] == 3


def test_filter_returns_none_when_too_few_points():
    pts = np.array([[0, 0, 10], [1, 1, 10]], dtype=float)
    pc = PointCloud(points=pts)
    config = RoofDetectionConfig(min_roof_points=20)

    filtered_pc, normals = filter_roof_points(pc, config)

    assert filtered_pc is None
    assert normals is None


# --- detect_roof_planes tests ---


def _estimate_normals(pts):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamKNN(knn=20))
    normals = np.asarray(pcd.normals)
    flip_mask = normals[:, 2] < 0
    normals[flip_mask] *= -1
    return normals


def test_detect_flat_roof():
    pts = _make_flat_roof_points(100, z=10.0, noise=0.05)
    normals = _estimate_normals(pts)
    pc = PointCloud(points=pts)
    config = RoofDetectionConfig()

    planes = detect_roof_planes(pc, normals, config)

    assert len(planes) >= 1
    assert planes[0].slope_deg < 15


def test_detect_gabled_roof():
    pts = _make_gabled_roof_points(100)
    normals = _estimate_normals(pts)
    pc = PointCloud(points=pts)
    config = RoofDetectionConfig()

    planes = detect_roof_planes(pc, normals, config)

    assert len(planes) >= 2


def test_detect_too_few_points_returns_empty():
    pts = np.array([[0, 0, 10], [1, 1, 10]], dtype=float)
    normals = np.array([[0, 0, 1], [0, 0, 1]], dtype=float)
    pc = PointCloud(points=pts)
    config = RoofDetectionConfig(min_plane_points=20)

    planes = detect_roof_planes(pc, normals, config)

    assert len(planes) == 0
