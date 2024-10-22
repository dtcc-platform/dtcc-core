import unittest
from dtcc_core.model import LineString, MultiLineString
import numpy as np


class TestLineString(unittest.TestCase):
    def test_length_2D(self):
        ls = LineString(vertices=np.array([(0, 0), (1, 0), (1, 1), (0, 1)]))
        self.assertEqual(ls.length, 3)

    def test_length_3D(self):
        ls = LineString(vertices=np.array([(0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2)]))
        self.assertEqual(ls.length, 3)

    def test_bounds_2D(self):
        ls = LineString(vertices=np.array([(0, 0), (1, 0), (1, 1), (0, 1)]))
        bounds = ls.bounds
        self.assertEqual(bounds.xmin, 0)
        self.assertEqual(bounds.xmax, 1)
        self.assertEqual(bounds.ymin, 0)
        self.assertEqual(bounds.ymax, 1)

    def test_bounds_3D(self):
        ls = LineString(vertices=np.array([(0, 0, 0), (1, 0, 1), (1, 1, 2), (0, 3, 2)]))
        bounds = ls.bounds
        self.assertEqual(bounds.xmin, 0)
        self.assertEqual(bounds.xmax, 1)
        self.assertEqual(bounds.ymin, 0)
        self.assertEqual(bounds.ymax, 3)
        self.assertEqual(bounds.zmin, 0)
        self.assertEqual(bounds.zmax, 2)

    def test_to_from_pb(self):
        ls = LineString(vertices=np.array([(0, 0, 0), (1, 0, 1), (1, 1, 2), (0, 3, 2)]))
        pb = ls.to_proto()
        ls2 = LineString()
        ls2.from_proto(pb)
        self.assertTrue(np.array_equal(ls.vertices, ls2.vertices))


class TestMultiLineString(unittest.TestCase):
    def test_length_2D(self):
        mls = MultiLineString(
            linestrings=[
                LineString(vertices=np.array([(0, 0), (1, 0), (1, 1), (0, 1)])),
                LineString(vertices=np.array([(1, 1), (2, 1), (2, 2), (1, 2)])),
            ]
        )
        self.assertEqual(mls.length, 6)

    def test_length_3D(self):
        mls = MultiLineString(
            linestrings=[
                LineString(
                    vertices=np.array([(0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2)])
                ),
                LineString(
                    vertices=np.array([(0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2)])
                ),
            ]
        )
        self.assertEqual(mls.length, 6)

    def test_bounds_2D(self):
        mls = MultiLineString(
            linestrings=[
                LineString(vertices=np.array([(0, 0), (1, 0), (1, 1), (0, 1)])),
                LineString(vertices=np.array([(1, 1), (2, 1), (2, 2), (0, 2)])),
            ]
        )
        bounds = mls.bounds
        self.assertEqual(bounds.xmin, 0)
        self.assertEqual(bounds.xmax, 2)
        self.assertEqual(bounds.ymin, 0)
        self.assertEqual(bounds.ymax, 2)

    def test_bounds_3D(self):
        mls = MultiLineString(
            linestrings=[
                LineString(
                    vertices=np.array([(0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2)])
                ),
                LineString(
                    vertices=np.array([(0, 0, 2), (1, 3, 2), (1, 1, 2), (0, 1, 2)])
                ),
            ]
        )
        bounds = mls.bounds
        self.assertEqual(bounds.xmin, 0)
        self.assertEqual(bounds.xmax, 1)
        self.assertEqual(bounds.ymin, 0)
        self.assertEqual(bounds.ymax, 3)
        self.assertEqual(bounds.zmin, 2)
        self.assertEqual(bounds.zmax, 2)


if __name__ == "__main__":
    unittest.main()
