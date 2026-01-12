import pytest
import numpy as np
from pathlib import Path
from shapely.geometry import Polygon
from dtcc_core.builder.meshing import (
    mesh_multisurface,
    mesh_surface,
    mesh_multisurfaces,
    disjoint_meshes,
)

from dtcc_core.model import Mesh, Surface, MultiSurface


@pytest.fixture
def simple_surface():
    return Surface(
        vertices=np.array(
            [
                [0, 0, 5],
                [0, 10, 5],
                [10, 10, 8],
                [10, 0, 8],
            ]
        )
    )


@pytest.fixture
def complex_surface():
    return Surface(
        vertices=np.array(
            [
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [2, 10, 0],
                [2, 12, 0],
                [0, 12, 0],
            ]
        )
    )


@pytest.fixture
def multi_surface(simple_surface, complex_surface):
    ms = MultiSurface()
    ms.surfaces = [simple_surface, complex_surface]
    return ms


@pytest.fixture
def second_surface_pair():
    surface3 = Surface(
        vertices=np.array(
            [
                [0, 0, 4],
                [0, 10, 4],
                [10, 10, 7],
                [10, 0, 7],
            ]
        )
    )
    surface4 = Surface(
        vertices=np.array(
            [
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [2, 10, 0],
                [2, 12, 0],
                [0, 12, 0],
            ]
        )
    )
    ms = MultiSurface()
    ms.surfaces = [surface3, surface4]
    return ms


@pytest.fixture
def disjoint_cubes_mesh():
    cube1_vertices = np.array(
        [
            [0, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
            [1, 0, 1],
        ]
    )
    cube1_faces = [
        [0, 1, 2],
        [0, 2, 3],
        [0, 4, 5],
        [0, 5, 1],
        [1, 5, 6],
        [1, 6, 2],
        [2, 6, 7],
        [2, 7, 3],
        [3, 7, 4],
        [3, 4, 0],
        [4, 7, 6],
        [4, 6, 5],
    ]
    cube2_vertices = cube1_vertices + np.array([2, 2, 2])
    cube2_faces = [[i + 8 for i in face] for face in cube1_faces]
    cubes_vertices = np.vstack([cube1_vertices, cube2_vertices])
    cubes_faces = np.array(cube1_faces + cube2_faces)
    return Mesh(vertices=cubes_vertices, faces=cubes_faces)


def test_mesh_simple_surface(simple_surface):
    mesh = mesh_surface(simple_surface)
    assert len(mesh.vertices) == 4
    assert len(mesh.faces) == 2
    assert pytest.approx(mesh.vertices[:, 2].min()) == 5
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_triangle_surface(simple_surface):
    mesh = mesh_surface(simple_surface, triangle_size=5)
    assert len(mesh.vertices) == 11
    assert len(mesh.faces) == 9
    assert pytest.approx(mesh.vertices[:, 2].min()) == 5
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_multisurface(multi_surface):
    mesh = mesh_multisurface(multi_surface)
    assert len(mesh.vertices) == 10
    assert len(mesh.faces) == 6
    assert pytest.approx(mesh.vertices[:, 2].min()) == 0
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_multisurfaces(multi_surface, second_surface_pair):
    ms = multi_surface.translate(0, 0, 1)
    meshes = mesh_multisurfaces([multi_surface, second_surface_pair])

    assert len(meshes) == 2

    # First mesh checks
    assert len(meshes[0].vertices) == 10
    assert len(meshes[0].faces) == 6
    assert pytest.approx(meshes[0].vertices[:, 2].min()) == 1
    assert pytest.approx(meshes[0].vertices[:, 2].max()) == 9

    # Second mesh checks
    assert len(meshes[1].vertices) == 10
    assert len(meshes[1].faces) == 6
    assert pytest.approx(meshes[1].vertices[:, 2].min()) == 0
    assert pytest.approx(meshes[1].vertices[:, 2].max()) == 7


def test_disjoint_mesh(disjoint_cubes_mesh):
    disjointed_meshes = disjoint_meshes(disjoint_cubes_mesh)

    assert len(disjointed_meshes) == 2
    assert disjointed_meshes[0].faces.max() == 7
    assert disjointed_meshes[1].faces.max() == 7
    assert len(disjointed_meshes[0].vertices) == 8
    assert len(disjointed_meshes[1].vertices) == 8


def test_snap_mesh_vertices():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 0, 0], [1, 1.01, 0]])
    mesh = Mesh(vertices=vertices, faces=np.array([[0, 1, 2], [1, 3, 4]]))
    mesh = mesh.snap_vertices(0.1)
    assert mesh.faces[1][2] == 2  # The last vertex should be snapped to the first one
