import pytest
import numpy as np
from pathlib import Path
from dtcc_core import io, builder
from dtcc_core.model import MultiSurface, Surface


@pytest.fixture
def base_surface():
    return Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))


@pytest.fixture
def coplanar_surface():
    return Surface(vertices=np.array([[1, 1, 0], [0, 1, 0], [0, 0, 0]]))


@pytest.fixture
def noncoplanar_surface():
    return Surface(vertices=np.array([[1, 1, 0], [0, 1, 1], [0, 0, 0]]))


@pytest.fixture
def parallel_surface():
    return Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [1, 1, 1]]))


@pytest.fixture
def extended_coplanar_surface():
    return Surface(vertices=np.array([[1, 1, 0], [0, 3, 0], [0, 0, 0]]))


@pytest.fixture
def complex_coplanar_surface():
    return Surface(vertices=np.array([[1, 1, 0], [1, 0, 3], [0, 0, 0]]))


def test_coplanar_surfaces(base_surface, coplanar_surface):
    assert builder.geometry.surface.are_coplanar(base_surface, coplanar_surface)


def test_noncoplanar_surfaces(base_surface, noncoplanar_surface):
    assert not builder.geometry.surface.are_coplanar(base_surface, noncoplanar_surface)


def test_parallel_noncoplanar_surfaces(base_surface, parallel_surface):
    assert not builder.geometry.surface.are_coplanar(base_surface, parallel_surface)


def test_merge_simple_coplanar_surfaces(base_surface, extended_coplanar_surface):
    ms = MultiSurface(surfaces=[base_surface, extended_coplanar_surface])
    merged = builder.geometry.multisurface.merge_coplanar(ms)
    assert len(merged.surfaces) == 1
    assert merged.bounds.ymax == 3


def test_merge_disjoint_surfaces(base_surface, parallel_surface):
    ms = MultiSurface(surfaces=[base_surface, parallel_surface])
    merged = builder.geometry.multisurface.merge_coplanar(ms)
    assert len(merged.surfaces) == 2


def test_merge_complex_coplanar_surfaces(
    base_surface, complex_coplanar_surface, coplanar_surface
):
    ms = MultiSurface(
        surfaces=[base_surface, complex_coplanar_surface, coplanar_surface]
    )
    merged = builder.geometry.multisurface.merge_coplanar(ms)
    assert len(merged.surfaces) == 2


def test_merge_large_cube():
    mesh_path = Path(__file__).parent / ".." / "data" / "meshes" / "large_cube.obj"
    mesh = io.load_mesh(mesh_path)
    ms_cube = mesh.to_multisurface()
    merged = builder.geometry.multisurface.merge_coplanar(ms_cube)
    assert len(merged.surfaces) == 6
