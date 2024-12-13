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


class TestBuilderPolygon(unittest.TestCase):
    def test_convert_surface(self):
        footprints = io.load_footprints(data_dir / "MinimalCase" / "PropertyMap.shp")
        footprint = footprints[0].get_footprint()
        builder_surface = builder.model_conversion.create_builder_surface(footprint)
        self.assertEqual(len(builder_surface.vertices), 4)
        self.assertEqual(footprint.vertices[0][0], builder_surface.vertices[0].x)


if __name__ == "__main__":
    unittest.main()
