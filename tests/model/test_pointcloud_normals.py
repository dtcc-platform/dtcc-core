import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud


def test_pointcloud_normals_default_none():
    pc = PointCloud(points=np.array([[0, 0, 0], [1, 1, 1]], dtype=float))
    assert pc.normals is None


def test_pointcloud_normals_assigned():
    pts = np.array([[0, 0, 0], [1, 1, 1]], dtype=float)
    norms = np.array([[0, 0, 1], [0, 0, 1]], dtype=float)
    pc = PointCloud(points=pts, normals=norms)
    assert pc.normals is not None
    assert pc.normals.shape == (2, 3)
