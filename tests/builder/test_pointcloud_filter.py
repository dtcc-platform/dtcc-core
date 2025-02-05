import unittest

import dtcc_core
from dtcc_core import io
from dtcc_core.model import Bounds
from pathlib import Path

data_dir = (Path(__file__).parent / ".." / "data" / "MinimalCase").resolve()
las_file = str((data_dir / "pointcloud.las").resolve())


class TestCropPointCloud(unittest.TestCase):
    def test_crop_nothing(self, xy_only=False):
        pc = io.load_pointcloud(las_file)
        bounds = pc.bounds
        cropped_pc = pc.crop(bounds)
        self.assertEqual(len(cropped_pc), 8148)
        self.assertEqual(len(cropped_pc.classification), 8148)

    def test_crop(self, xy_only=False):
        pc = io.load_pointcloud(las_file)
        bounds = Bounds(-2, -2, 0, 0)
        cropped_pc = pc.crop(bounds)

        self.assertEqual(len(pc), 8148)
        self.assertEqual(len(cropped_pc), 64)
        self.assertEqual(len(cropped_pc.classification), 64)
