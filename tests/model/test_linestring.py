import pytest
from dtcc_core.model import LineString
import numpy as np


@pytest.fixture
def linstring2D():
    return LineString(vertices=np.array([(0, 0), (1, 0), (1, 1), (0, 1)]))


@pytest.fixture
def linstring3D():
    return LineString(vertices=np.array([(0, 0, 0), (1, 0, 2), (1, 3, 2), (0, 1, 2)]))


def test_length_2D(linstring2D):
    assert linstring2D.length == 3


def test_length_3D(linstring3D):
    line_len = 0
    for i in range(1, len(linstring3D.vertices)):
        line_len += np.linalg.norm(
            linstring3D.vertices[i] - linstring3D.vertices[i - 1]
        )
    assert linstring3D.length == line_len


def test_bounds_2D(linstring2D):
    bounds = linstring2D.bounds
    assert bounds.xmin == 0
    assert bounds.xmax == 1
    assert bounds.ymin == 0
    assert bounds.ymax == 1


def test_bounds_3D(linstring3D):
    bounds = linstring3D.bounds
    assert bounds.xmin == 0
    assert bounds.xmax == 1
    assert bounds.ymin == 0
    assert bounds.ymax == 3
    assert bounds.zmin == 0
    assert bounds.zmax == 2


def test_to_from_pb(linstring3D):
    pb = linstring3D.to_proto()
    ls2 = LineString()
    ls2.from_proto(pb)
    assert np.array_equal(linstring3D.vertices, ls2.vertices)


if __name__ == "__main__":
    pytest.main()
