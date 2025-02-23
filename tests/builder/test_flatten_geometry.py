import pytest
from dtcc_core.model import (
    Building,
    BuildingPart,
    Surface,
    MultiSurface,
    Mesh,
    GeometryType,
)
import numpy as np


@pytest.fixture
def unit_cube():
    floor = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]))
    ceiling = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]))
    wall1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]))
    wall2 = Surface(vertices=np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]]))
    wall3 = Surface(vertices=np.array([[1, 1, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1]]))
    wall4 = Surface(vertices=np.array([[0, 1, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]]))
    ceiling = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]))
    return (floor, ceiling, wall1, wall2, wall3, wall4)


@pytest.fixture
def floor_roof(unit_cube):
    return MultiSurface(surfaces=[unit_cube[0], unit_cube[1]])


@pytest.fixture
def walls(unit_cube):
    return MultiSurface(
        surfaces=[unit_cube[2], unit_cube[3], unit_cube[4], unit_cube[5]]
    )


@pytest.fixture
def building1(unit_cube):
    ms1 = MultiSurface(surfaces=unit_cube)
    b = Building()
    b.add_geometry(ms1, GeometryType.LOD1)
    return b


@pytest.fixture
def building2(walls, floor_roof):
    b = Building()
    p1 = BuildingPart()
    p1.add_geometry(walls, GeometryType.LOD1)
    p2 = BuildingPart()
    p2.add_geometry(floor_roof, GeometryType.LOD1)
    b.add_children([p1, p2])
    return b


@pytest.fixture
def building3(walls, floor_roof):
    b = Building()
    p1 = BuildingPart()
    b.add_geometry(walls, GeometryType.LOD1)
    p1.add_geometry(floor_roof, GeometryType.LOD1)
    b.add_child(p1)
    return b


def test_flatten_no_child(building1):
    f = building1.flatten_geometry(GeometryType.LOD1)
    assert f is not None
    assert isinstance(f, MultiSurface)
    assert len(f.surfaces) == 6


def test_flatten_children(building2):
    f = building2.flatten_geometry(GeometryType.LOD1)
    assert f is not None
    assert isinstance(f, MultiSurface)
    assert len(f.surfaces) == 6
    # don't mutate original building
    assert len(building2.geometry) == 0


def test_flatten_root_and_child(building3):
    f = building3.flatten_geometry(GeometryType.LOD1)
    assert f is not None
    assert isinstance(f, MultiSurface)
    assert len(f.surfaces) == 6
    # don't mutate original building
    assert len(building3.geometry) == 1
    orig_geom = building3.geometry.get(GeometryType.LOD1, None)
    assert orig_geom is not None
    assert len(orig_geom.surfaces) == 4
