import unittest

from pathlib import Path

import numpy as np
import os, tempfile
from dtcc_core import io


class TestGridField(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = (Path(__file__).parent / ".." / "data").resolve()
        cls.dem_raster = cls.data_dir / "test_dem.tif"

        cls.rgb_img = cls.data_dir / "14040.png"

    def test_load_dem(self):
        dem = io.load_raster(self.dem_raster)
        self.assertEqual(dem.width, 20)
        self.assertEqual(dem.height, 40)
        self.assertEqual(dem.channels, 1)

    def test_load_image(self):
        image = io.load_raster(self.rgb_img)
        self.assertEqual(image.width, 228)
        self.assertEqual(image.height, 230)
        self.assertEqual(image.channels, 3)

    def test_get_cell_size(self):
        em = io.load_raster(self.dem_raster)
        self.assertEqual(em.cell_size, (2.0, -2.0))

    def test_get_cell_size_wld_file(self):
        em = io.load_raster(self.rgb_img)
        self.assertEqual(em.cell_size, (0.08, -0.08))

    def test_write_elevation_model(self):
        em = io.load_raster(self.dem_raster)
        outfile = tempfile.NamedTemporaryFile(suffix=".tif", delete=False).name
        em.save(outfile)
        em = io.load_raster(outfile)
        self.assertEqual(em.height, 40)
        self.assertEqual(em.width, 20)
        self.assertEqual(em.cell_size, (2.0, -2.0))
        os.unlink(outfile)

    def test_load_multiple(self):
        in_rasters = [
            self.data_dir / "testraster_0_0.tif",
            self.data_dir / "testraster_0_1.tif",
            self.data_dir / "testraster_1_0.tif",
            self.data_dir / "testraster_1_1.tif",
        ]
        raster = io.load_raster(in_rasters)
        self.assertEqual(raster.width, 50)
        self.assertEqual(raster.height, 50)
        self.assertEqual(raster.cell_size, (0.5, -0.5))

        test_raster = io.load_raster(self.data_dir / "testraster.tif")
        self.assertEqual(np.all(raster.data == test_raster.data), True)


if __name__ == "__main__":
    unittest.main()
