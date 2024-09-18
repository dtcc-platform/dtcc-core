import unittest
from pathlib import Path
from dtcc_core.model import Raster, PointCloud

from dtcc_core import io
import dtcc_core.builder

test_data_dir = project_dir = (Path(__file__).parent / ".." / "data" ).resolve()


class TestConvert(unittest.TestCase):
    def test_raster_to_pointcloud(self):
        raster = io.load_raster(test_data_dir/"test_dem.tif")
        raster_bounds = raster.bounds
        cell_size = raster.cell_size
        pointcloud = raster.to_pointcloud()
        self.assertIsInstance(pointcloud, PointCloud)
        self.assertEqual(pointcloud.points[:, 2].max(), raster.data.max())
        self.assertEqual(pointcloud.points[:, 2].min(), raster.data.min())

        self.assertEqual(
            pointcloud.points[:, 0].min(), raster_bounds.xmin + cell_size[0] / 2
        )
        self.assertEqual(
            pointcloud.points[:, 1].min(), raster_bounds.ymin - (cell_size[1] / 2)
        )

        self.assertEqual(
            pointcloud.points[:, 0].max(), raster.bounds.xmax - cell_size[0] / 2
        )
        self.assertEqual(
            pointcloud.points[:, 1].max(), raster.bounds.ymax + cell_size[1] / 2
        )


if __name__ == "__main__":
    unittest.main()
