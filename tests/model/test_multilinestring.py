import pytest
from dtcc_core.model import LineString, MultiLineString
import numpy as np


@pytest.fixture
def mls_2d():
    return MultiLineString(
        linestrings=[
            LineString(vertices=np.array([(0, 0), (1, 0), (1, 1), (0, 1)])),
            LineString(vertices=np.array([(1, 1), (2, 1), (2, 2), (1, 2)])),
        ]
    )


@pytest.fixture
def mls_3d():
    return MultiLineString(
        linestrings=[
            LineString(vertices=np.array([(0, 0, 3), (1, 0, 3), (1, 1, 3), (0, 1, 3)])),
            LineString(vertices=np.array([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])),
        ]
    )


def test_length_2D(mls_2d):
    assert mls_2d.length == 6


def test_length_3D(mls_3d):
    assert mls_3d.length == 6


def test_bounds_2D(mls_2d):
    bounds = mls_2d.bounds
    assert bounds.xmin == 0
    assert bounds.xmax == 2
    assert bounds.ymin == 0
    assert bounds.ymax == 2


def test_bounds_3D(mls_3d):
    bounds = mls_3d.bounds
    assert bounds.xmin == 0
    assert bounds.xmax == 1
    assert bounds.ymin == 0
    assert bounds.ymax == 1
    assert bounds.zmin == 1
    assert bounds.zmax == 3
