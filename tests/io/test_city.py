import unittest
from dtcc_core import io
from dtcc_core.model import City
from pathlib import Path

testdata = Path(__file__).parent / ".." / "data" / "cityjson" / "DenHaag_01.city.json"


class MyTestCase(unittest.TestCase):
    def test_load_cityjson(self):
        city = io.load_city(testdata)
        self.assertIsInstance(city, City)
        self.assertEqual(len(city.buildings), 844)

    def test_to_df(self):
        city =io.load_city(testdata)
        df = io.city.buildings_to_df(city)
        self.assertEqual(len(df), 844)

        try:
            import dtcc_core.builder

            # these tests will fail without dtcc-builder installed
            self.assertTrue("geometry" in df.columns)
        except ImportError:
            print("dtcc-builder not installed, skipping tests")
            pass


if __name__ == "__main__":
    unittest.main()
