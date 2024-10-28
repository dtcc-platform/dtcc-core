import unittest

import numpy as np

import dtcc_core.builder
from dtcc_core.model import Surface, MultiSurface, Mesh


class TestSurfaceIntersection(unittest.TestCase):
    def test_ray_hit(self):
        surface = Surface(vertices=np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
        origin = [0.5, 0.5, 1]
        direction = [0, 0, -1]
        intersection = surface.ray_intersection(origin, direction)

        self.assertTrue(isinstance(intersection, np.ndarray))
        self.assertTrue((intersection == [0.5, 0.5, 0]).all())
        # self.assertTrue((intersection == [0, 0, 0]).all())

    def test_ray_miss(self):
        surface = Surface(vertices=np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
        origin = [0.5, 0.5, 1]
        direction = [0, 0, 1]
        intersection = surface.ray_intersection(origin, direction)
        self.assertTrue(isinstance(intersection, np.ndarray))

        self.assertTrue(np.isnan(intersection).all())


class TestMultiSurfaceIntersection(unittest.TestCase):
    def test_ray_hit(self):
        surface1 = Surface(vertices=np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
        surface2 = Surface(vertices=np.array([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]))

        multisurface = MultiSurface(surfaces=[surface1, surface2])
        origin = [0.5, 0.5, 2]
        direction = [0, 0, -1]
        intersection = multisurface.ray_intersection(origin, direction)

        self.assertTrue(isinstance(intersection, np.ndarray))
        self.assertTrue((intersection == [0.5, 0.5, 1]).all())

    def test_ray_miss(self):
        surface1 = Surface(vertices=np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
        surface2 = Surface(vertices=np.array([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]))
        multisurface = MultiSurface(surfaces=[surface1, surface2])

        origin = [0.5, 0.5, 2]
        direction = [0, 0, 1]
        intersection = multisurface.ray_intersection(origin, direction)
        self.assertTrue(isinstance(intersection, np.ndarray))

        self.assertTrue(np.isnan(intersection).all())


class TestMeshIntersection(unittest.TestCase):
    def test_mesh_ray_intersection(self):
        vertices = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 1, 1],
            ]
        )
        faces = np.array([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]])
        mesh = Mesh(vertices=vertices, faces=faces)
        origin = [0.5, 0.5, 2]
        direction = [0, 0, -1]
        intersection = mesh.ray_intersection(origin, direction)
        self.assertTrue(isinstance(intersection, np.ndarray))
        self.assertTrue((intersection == [0.5, 0.5, 1]).all())
