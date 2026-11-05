"""Canopy extraction from a classified point cloud on known flat terrain."""

import numpy as np
import pytest
from affine import Affine
from shapely.geometry import Point

from dtcc_core import builder
from dtcc_core.model import PointCloud, Raster


def test_tree_workflow_preserves_ground_height_and_world_coordinates():
    rows, cols = np.indices((32, 32))
    heights = 10 * np.exp(-((rows - 10) ** 2 + (cols - 10) ** 2) / 12) + 8 * np.exp(
        -((rows - 22) ** 2 + (cols - 22) ** 2) / 12
    )
    terrain = Raster(
        data=np.full((32, 32), 15.0),
        georef=Affine.translation(100, 200) * Affine.scale(0.5, -0.5),
    )
    x, y = terrain.georef * (cols.ravel() + 0.5, rows.ravel() + 0.5)
    pc = PointCloud(
        points=np.column_stack((x, y, 15 + heights.ravel())),
        classification=np.full(rows.size, 5, dtype=np.uint8),
    )
    original = pc.points.copy()
    canopy = builder.tree_raster_from_pointcloud(pc, terrain)
    coords, tops, radii = builder.find_tree_tops(canopy)
    assert len(coords) == 2
    assert np.all((tops > 5) & (tops < 11))
    assert np.all(radii > 0)
    crowns, valid = builder.tree_crown_polygons(canopy, coords)
    assert valid.tolist() == [True, True]
    for polygon, (row, col) in zip(crowns, coords):
        assert polygon.is_valid and polygon.area > 1
        assert polygon.covers(Point(*(canopy.georef * (col + 0.5, row + 0.5))))
    trees = builder.trees_from_pointcloud(pc, terrain)
    assert len(trees) == 2
    for tree, coord, height in zip(trees, coords, tops):
        np.testing.assert_allclose(tree.position[:2], terrain.pixel_to_georef(*coord))
        assert tree.position[2] == 15
        assert tree.height == pytest.approx(height)
    np.testing.assert_array_equal(pc.points, original)
    np.testing.assert_array_equal(terrain.data, np.full((32, 32), 15.0))
    with pytest.raises(ValueError, match="Tree type"):
        builder.find_tree_tops(canopy, tree_type="unknown")
