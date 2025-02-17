import pytest
import numpy as np
from shapely.geometry import Polygon
from dtcc_core.model import Surface  # Assuming this is the correct import


@pytest.fixture
def simple_polygon():
    return Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])


@pytest.fixture
def basic_surface(simple_polygon):
    surface = Surface()
    surface.from_polygon(simple_polygon, 10)
    return surface


@pytest.fixture
def complex_surface():
    verts = np.array([[0, 0, 10], [1, 0, 10], [1, 1, 12], [0, 1, 12]])
    hole = np.array([[0.1, 0.1, 10], [0.9, 0.1, 10], [0.9, 0.9, 12], [0.1, 0.9, 12]])
    return Surface(vertices=verts, holes=[hole])


def test_convert_polygon(basic_surface):
    v1 = basic_surface.vertices[1]
    assert len(basic_surface.vertices) == 4
    assert v1[0] == 1
    assert v1[1] == 0
    assert v1[2] == 10
    assert basic_surface.zmax == 10


def test_to_polygon(basic_surface):
    p = basic_surface.to_polygon()
    coords = list(p.exterior.coords)

    assert len(coords) == 4 + 1  # +1 for closing point
    assert coords[0] == (0, 0)
    assert coords[1] == (1, 0)
    assert coords[2] == (1, 1)
    assert coords[3] == (0, 1)


def test_to_proto(complex_surface):
    pb = complex_surface.to_proto()

    assert len(pb.surface.vertices) == 4 * 3
    assert pb.surface.vertices[0] == 0
    assert pb.surface.vertices[-1] == 12
    assert len(pb.surface.holes) == 1


def test_from_proto(complex_surface):
    # Convert to proto and back
    pb = complex_surface.to_proto()
    pb_string = pb.SerializeToString()

    s2 = Surface()
    s2.from_proto(pb_string)

    assert len(s2.vertices) == 4
    assert len(s2.holes) == 1
    assert list(s2.vertices[0]) == [0, 0, 10]
    assert list(np.round(s2.holes[0][2], 6)) == [0.9, 0.9, 12]


if __name__ == "__main__":
    pytest.main()
