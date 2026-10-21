import numpy as np
import pytest

from dtcc_core.model import (
    Building,
    BuildingPart,
    GeometryType,
    MultiSurface,
    Surface,
)


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
    orig_geom = building3.get_geometry(GeometryType.LOD1)
    assert orig_geom is not None
    assert len(orig_geom.surfaces) == 4


def _square(x0, y0, size, z=0.0):
    return Surface(
        vertices=np.array(
            [[x0, y0, z], [x0 + size, y0, z], [x0 + size, y0 + size, z], [x0, y0 + size, z]],
            dtype=float,
        )
    )


def test_a_footprint_covers_the_footprints_of_children():
    # The model contract for LOD0: a footprint of an object covers those of
    # its children (the parts of a building lie inside its footprint), so a
    # parent with a footprint gives its own, and the children's are not
    # added. Without one, the children's footprints are collected.
    building = Building()
    building.add_geometry(_square(0, 0, 10), GeometryType.LOD0)
    for x in (0, 5):
        part = BuildingPart()
        part.add_geometry(_square(x, 0, 5), GeometryType.LOD0)
        building.add_child(part)
    footprint = building.flatten_geometry(GeometryType.LOD0)
    assert isinstance(footprint, Surface)
    assert np.allclose(footprint.vertices, _square(0, 0, 10).vertices)

    parts_only = Building()
    for x in (0, 5):
        part = BuildingPart()
        part.add_geometry(_square(x, 0, 5), GeometryType.LOD0)
        parts_only.add_child(part)
    footprints = parts_only.flatten_geometry(GeometryType.LOD0)
    assert isinstance(footprints, MultiSurface)
    assert len(footprints.surfaces) == 2


def test_other_geometry_of_a_parent_and_its_children_is_combined(walls, floor_roof):
    # LOD1 and higher are solids or surfaces of their own: a parent's and its
    # children's are combined, at every depth.
    building = Building()
    building.add_geometry(walls, GeometryType.LOD1)
    part = BuildingPart()
    part.add_geometry(floor_roof, GeometryType.LOD1)
    inner = BuildingPart()
    inner.add_geometry(MultiSurface(surfaces=[_square(0, 0, 1, 2.0)]), GeometryType.LOD1)
    part.add_child(inner)
    building.add_child(part)
    assert len(building.flatten_geometry(GeometryType.LOD1).surfaces) == 4 + 2 + 1
