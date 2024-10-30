import unittest
from dtcc_core.model.object import RoadNetwork, RoadType
import numpy as np


class TestRoadnetwork(unittest.TestCase):
    def test_create(self):
        rn = RoadNetwork()
        self.assertEqual(len(rn.vertices), 0)
        self.assertEqual(len(rn.edges), 0)

    def test_bounds(self):
        rn = RoadNetwork()
        rn.vertices = np.array([(0, 0), (1, 0), (1, 1), (0, 2)])
        rn.edges = np.array(
            [
                (0, 1),
                (1, 2),
                (2, 3),
            ]
        )
        bounds = rn.bounds
        self.assertEqual(bounds.xmin, 0)
        self.assertEqual(bounds.xmax, 1)
        self.assertEqual(bounds.ymin, 0)
        self.assertEqual(bounds.ymax, 2)


if __name__ == "__main__":
    unittest.main()
