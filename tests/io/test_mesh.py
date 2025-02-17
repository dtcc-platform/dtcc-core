import pytest
import json
import os
import tempfile
import pathlib
import numpy as np
from dtcc_core import io
from dtcc_core.model import Mesh, VolumeMesh, City


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).parent / ".." / "data" / "meshes"


@pytest.fixture
def mesh_paths(data_dir):
    return {
        "stl": str((data_dir / "cube.stl").resolve()),
        "vtk": str((data_dir / "cube.vtk").resolve()),
        "fbx": str((data_dir / "cube.fbx").resolve()),
        "quad_cube": str((data_dir / "quad_cube.obj").resolve()),
        "quad_tri_cube": str((data_dir / "quad_and_tri_cube.obj").resolve()),
        "nine_cubes": str((data_dir / "9_cubes.obj").resolve()),
    }


@pytest.fixture
def basic_mesh(mesh_paths):
    return io.load_mesh(mesh_paths["stl"])


@pytest.fixture
def volume_mesh_cube():
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
    cells = np.array(
        [
            [0, 1, 2, 4],
            [1, 2, 3, 5],
            [2, 3, 0, 6],
            [3, 0, 1, 7],
            [0, 4, 5, 1],
            [2, 6, 7, 3],
        ]
    )
    return VolumeMesh(vertices=vertices, cells=cells)


# Basic mesh loading tests
def test_load_mesh_stl(mesh_paths):
    mesh = io.load_mesh(mesh_paths["stl"])
    assert len(mesh.vertices) == 24
    assert len(mesh.faces) == 44


def test_load_mesh_vtk(mesh_paths):
    mesh = io.load_mesh(mesh_paths["vtk"])
    assert len(mesh.vertices) == 24
    assert len(mesh.faces) == 44


# @pytest.mark.skip(reason="This test segfaults")
def test_load_mesh_fbx(mesh_paths):
    if not io.meshes.has_assimp():
        pytest.skip("pyassimp not found, skipping test")
    mesh = io.load_mesh(mesh_paths["fbx"])
    assert len(mesh.vertices) == 24
    assert len(mesh.faces) == 44


# Save and load tests for different formats
@pytest.mark.parametrize("file_extension", [".stl", ".vtk"])
def test_save_load_mesh(basic_mesh, file_extension):
    with tempfile.NamedTemporaryFile(suffix=file_extension, delete=False) as tmp_file:
        path = tmp_file.name

    try:
        io.save_mesh(basic_mesh, path)
        loaded_mesh = io.load_mesh(path)
        assert len(loaded_mesh.vertices) == 24
        assert len(loaded_mesh.faces) == 44
    finally:
        os.unlink(path)


@pytest.mark.parametrize("file_extension", [".glb", ".gltf"])
def test_save_mesh_only(basic_mesh, file_extension):
    with tempfile.NamedTemporaryFile(suffix=file_extension, delete=False) as tmp_file:
        path = tmp_file.name

    try:
        basic_mesh.save(path)
        assert os.path.exists(path)
    finally:
        os.unlink(path)


# Specialized mesh tests
def test_load_quad_mesh(mesh_paths):
    mesh = io.load_mesh(mesh_paths["quad_cube"])
    assert len(mesh.vertices) == 8
    assert len(mesh.faces) == 12
    assert mesh.faces.shape[1] == 3


def test_load_quad_tri_mesh(mesh_paths):
    mesh = io.load_mesh(mesh_paths["quad_tri_cube"])
    assert len(mesh.vertices) == 8
    assert len(mesh.faces) == 12
    assert mesh.faces.shape[1] == 3


def test_load_disjoint_mesh(mesh_paths):
    mesh = io.load_mesh(mesh_paths["nine_cubes"])
    assert len(mesh.vertices) == 9 * 8
    assert len(mesh.faces) == 9 * 12


@pytest.mark.parametrize("merge_coplanar", [False, True])
def test_load_city_mesh(mesh_paths, merge_coplanar):
    city = io.load_mesh_as_city(
        mesh_paths["nine_cubes"], merge_coplanar_surfaces=merge_coplanar
    )
    assert isinstance(city, City)
    assert len(city.buildings) == 9

    if merge_coplanar:
        assert len(city.buildings[0].lod1.surfaces) == 6


# Volume mesh tests
def test_write_read_volume_mesh(volume_mesh_cube):
    with tempfile.NamedTemporaryFile(suffix=".vtk", delete=False) as tmp_file:
        path = tmp_file.name

    try:
        volume_mesh_cube.save(path)
        assert os.path.exists(path)
        mesh = io.load_volume_mesh(path)
        assert len(mesh.vertices) == 8
        assert len(mesh.cells) == 6
        assert mesh.cells.dtype == np.int64
    finally:
        os.unlink(path)
