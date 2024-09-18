import unittest
from dtcc_core import builder, io
from dtcc_core.model import RoadNetwork, GeometryType, Surface

from pathlib import Path

test_data = (
    Path(__file__).parent / ".." / "data" / "road_network" / "test_road.shp"
).resolve()

assert test_data.is_file()


class TestRoadnetwork(unittest.TestCase):
    def test_bidirectional_to_matrix(self):
        rn = io.load_roadnetwork(test_data)
        largest_idx = max([e[0] for e in rn.edges] + [e[1] for e in rn.edges])
        num_roads = len(rn.edges)
        rn_matrix = rn.to_matrix()
        self.assertEqual(rn_matrix.shape, (largest_idx + 1, largest_idx + 1))
        self.assertEqual(
            rn_matrix.nnz, (num_roads * 2) - 1
        )  # bidirectional so each edge is counted twice, minus 1
        # because we have one loop that starts that shouldn't be counted twice

        r12 = rn.edges[12]
        l12 = rn.length[12]
        self.assertEqual(rn_matrix[r12[0], r12[1]], l12)
        self.assertEqual(rn_matrix[r12[1], r12[0]], l12)

    def test_unidirectional_to_matrix(self):
        rn = io.load_roadnetwork(test_data)
        largest_idx = max([e[0] for e in rn.edges] + [e[1] for e in rn.edges])
        num_roads = len(rn.edges)
        rn_matrix = rn.to_matrix(bidirectional=False)
        self.assertEqual(rn_matrix.nnz, num_roads)

    def test_to_polygon(self):
        rn = io.load_roadnetwork(test_data)
        rn_poly = rn.to_surfaces(as_shapely=True, widths=4)
        self.assertEqual(len(rn_poly), len(rn.linestrings))
        self.assertEqual(rn_poly[0].geom_type, "Polygon")
        self.assertLessEqual(abs(1 - (rn_poly[0].area / (rn.length[0] * 4))), 0.02)

    def test_to_surfaces(self):
        rn = io.load_roadnetwork(test_data)
        rn_surfaces = rn.to_surfaces(widths=4, as_shapely=False)
        self.assertEqual(len(rn_surfaces), len(rn.linestrings))
        self.assertTrue(isinstance(rn_surfaces[0], Surface))


if __name__ == "__main__":
    unittest.main()
