import unittest
from pathlib import Path

from dtcc_core import io
from dtcc_core.model import Bounds

data_dir = (Path(__file__).parent / ".." / "data").resolve()


class TestPointCloudDirectory(unittest.TestCase):
    def test_load_pc_directory(self):
        pcd = io.load_pointcloud_directory(data_dir / "pointclouds")
        self.assertEqual(len(pcd), 2)

    def test_bounds(self):
        pcd = io.load_pointcloud_directory(data_dir / "pointclouds")
        bounds = pcd.bounds
        self.assertEqual(bounds.ymin, 0)
        self.assertEqual(bounds.ymax, 1)

    def test_get_all_points(self):
        pcd = io.load_pointcloud_directory(data_dir / "pointclouds")
        pc = pcd.pointcloud(pcd.bounds)
        self.assertEqual(len(pc), 20)

    def test_get_sub_points_1file(self):
        pcd = io.load_pointcloud_directory(data_dir / "pointclouds")
        bounds = Bounds(xmin=0.1, xmax=3.1, ymin=0.1, ymax=1.1)
        pc = pcd.pointcloud(bounds)
        self.assertEqual(len(pc), 3)

    def test_get_sub_points_2files(self):
        pcd = io.load_pointcloud_directory(data_dir / "pointclouds")
        bounds = Bounds(xmin=0.1, xmax=3.1, ymin=-1, ymax=3)
        pc = pcd.pointcloud(bounds)
        self.assertEqual(len(pc), 6)


if __name__ == "__main__":
    unittest.main()
