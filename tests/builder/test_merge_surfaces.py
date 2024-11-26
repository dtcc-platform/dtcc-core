import unittest
from dtcc_core import io, builder
from dtcc_core.model import MultiSurface, Surface
import numpy as np
from pathlib import Path

class TestCoplanar(unittest.TestCase):

    def test_are_simple_coplanar_surfaces(self):
        surface1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        surface2 = Surface(vertices=np.array([[1, 1, 0], [0, 1, 0], [0, 0, 0]]))
        self.assertTrue(builder.geometry.surface.are_coplanar(surface1, surface2))

    def test_are_not_coplanar_surfaces(self):
        surface1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        surface2 = Surface(vertices=np.array([[1, 1, 0], [0, 1, 1], [0, 0, 0]]))
        self.assertFalse(builder.geometry.surface.are_coplanar(surface1, surface2))

    def test_are_parllel_not_coplanar(self):
        surface1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        surface2 = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [1, 1, 1]]))
        self.assertFalse(builder.geometry.surface.are_coplanar(surface1, surface2))


class TestMergeCoplanar(unittest.TestCase):
    def test_merge_simple_coplanar(self):
        surface1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        surface2 = Surface(vertices=np.array([[1, 1, 0], [0, 3, 0], [0, 0, 0]]))
        ms = MultiSurface(surfaces=[surface1, surface2])
        merged = builder.geometry.multisurface.merge_coplanar(ms)
        self.assertEqual(len(merged.surfaces), 1)
        self.assertEqual(merged.bounds.ymax, 3)

    def test_merge_disjoint(self):
        surface1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        surface2 = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [1, 1, 1]]))
        ms = MultiSurface(surfaces=[surface1, surface2])
        merged = builder.geometry.multisurface.merge_coplanar(ms)
        self.assertEqual(len(merged.surfaces), 2)

    def test_merge_complex_coplanar(self):
        surface1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        surface2 = Surface(vertices=np.array([[1, 1, 0], [1, 0, 3], [0, 0, 0]]))
        surface3 = Surface(vertices=np.array([[1, 1, 0], [0, 1, 0], [0, 0, 0]]))

        ms = MultiSurface(surfaces=[surface1, surface2, surface3])
        merged = builder.geometry.multisurface.merge_coplanar(ms)
        self.assertEqual(len(merged.surfaces), 2)

class TestMergeMesh(unittest.TestCase):
    def test_merge_large_cube(self):
        mesh = io.load_mesh((Path(__file__).parent / ".."/ "data" / "large_cube.obj"))
        ms_cube = mesh.to_multisurface()
        merged = builder.geometry.multisurface.merge_coplanar(ms_cube)
        self.assertEqual(len(merged.surfaces), 6)