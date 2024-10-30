import unittest
from dtcc_core.model import Landuse, LanduseClasses
from dtcc_core import io

from pathlib import Path

testdata = Path(__file__).parent / ".." / "data" / "landuse_testdata.shp"
assert testdata.is_file()


class TestLanduse(unittest.TestCase):
    def test_load(self):
        landuse = io.load_landuse(testdata)
        self.assertIsInstance(landuse, Landuse)

    def test_geometry(self):
        landuse = io.load_landuse(testdata)
        self.assertEqual(len(landuse.surfaces), 79)

    def test_classes(self):
        landuse = io.load_landuse(testdata)
        self.assertEqual(len(landuse.landuses), 79)
        for l in landuse.landuses:
            self.assertIsInstance(l, LanduseClasses)


if __name__ == "__main__":
    unittest.main()
