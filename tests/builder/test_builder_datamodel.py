import unittest

import sys, os
from pathlib import Path
from dtcc_core import builder
from dtcc_core import io

data_dir = (Path(__file__).parent / "../data").resolve()


# class TestBuilderCity(unittest.TestCase):
#     def test_convert_city(self):
#         city = io.load_city(data_dir / "MinimalCase" / "PropertyMap.shp")
#         builder_city = builder.model.create_builder_city(city)
#         self.assertEqual(len(builder_city.buildings), len(city.buildings))
#         self.assertEqual(len(builder_city.buildings), 5)
#
#     def test_convert_city_with_holes(self):
#         city = io.load_city(data_dir / "MinimalCase" / "PropertyMap.shp")
#         builder_city = builder.model.create_builder_city(city)
#         self.assertEqual(len(builder_city.buildings[0].footprint.holes), 0)
#         self.assertEqual(len(builder_city.buildings[4].footprint.holes), 1)

# class TestBuilderPolygon(unittest.TestCase):
#     def test_convert_polygon(self):
#         city = io.load_city(data_dir / "MinimalCase" / "PropertyMap.shp")
#         footprint = city.buildings[0].footprint
#         builder_polygon = builder.model.create_builder_polygon(footprint)
#         self.assertEqual(len(builder_polygon.vertices), 4)
#         self.assertEqual(footprint.exterior.coords[0][0], builder_polygon.vertices[0].x)


if __name__ == "__main__":
    unittest.main()
