import unittest

import numpy as np

from dtcc_core.model.object import RoadNetwork, RoadType
from dtcc_core.model.logging import info
from dtcc_core.io import roadnetwork
from dtcc_core.model import Bounds
from pathlib import Path

datapath = Path(__file__).parent / ".." / "data" / "road_network" / "test_road.shp"

try:
    import geopandas as gpd
except ImportError:
    gpd = None


class LoadRoadnetwork(unittest.TestCase):

    def test_load(self):
        rn = roadnetwork.load(datapath)
        self.assertIsInstance(rn, RoadNetwork)
        self.assertEqual(len(rn.vertices), 185)
        self.assertEqual(len(rn.edges), 244)
        self.assertEqual(len(rn.length), 244)
        self.assertEqual(len(rn.linestrings), len(rn.edges))

    def test_load_bounds(self):

        bounds = Bounds(xmin=58150, ymin=6398226, xmax=59150, ymax=6398984)
        rn = roadnetwork.load(datapath, bounds=bounds)
        self.assertIsInstance(rn, RoadNetwork)

        for k, v in rn.attributes.items():
            self.assertEqual(len(v), len(rn.edges))
        self.assertEqual(len(rn.edges), 37)

    def test_load_no_geometry(self):
        rn = roadnetwork.load(datapath, load_geometry=False)
        self.assertIsInstance(rn, RoadNetwork)
        self.assertEqual(len(rn.vertices), 185)
        self.assertEqual(len(rn.edges), 244)
        self.assertEqual(len(rn.length), 244)
        self.assertEqual(len(rn.linestrings), 0)

    def test_attributes(self):
        rn = roadnetwork.load(datapath)
        self.assertEqual(len(rn.attributes), 7)
        for k, v in rn.attributes.items():
            self.assertEqual(len(v), 244)


class TestRoadnetworkProtobuf(unittest.TestCase):

    def test_to_from_protobuf(self):
        rn = roadnetwork.load(datapath)

        pb = rn.to_proto()
        rn2 = RoadNetwork()
        rn2.from_proto(pb)

        self.assertEqual(len(rn.vertices), len(rn2.vertices))
        self.assertTrue(np.allclose(rn.vertices, rn2.vertices))
        self.assertEqual(len(rn.edges), len(rn2.edges))
        self.assertTrue(np.allclose(rn.edges, rn2.edges))
        self.assertEqual(len(rn.length), len(rn2.length))
        self.assertTrue(np.allclose(rn.length, rn2.length))


class TestDataFrame(unittest.TestCase):
    def test_to_df(self):
        if gpd is not None:
            rn = roadnetwork.load(datapath)
            df = rn.to_df()
            self.assertIsInstance(df, gpd.GeoDataFrame)
            self.assertEqual(len(df), len(rn.edges))
        else:
            info("Geopandas not found, skipping test")


if __name__ == "__main__":
    unittest.main()
