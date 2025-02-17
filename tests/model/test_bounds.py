import pytest
from dtcc_core.model import Bounds


def test_create():
    bounds = Bounds(1, 2, 3, 4)
    assert bounds.xmin == 1
    assert bounds.ymin == 2
    assert bounds.xmax == 3
    assert bounds.ymax == 4


def test_area():
    bounds = Bounds(0, 0, 10, 10)
    assert bounds.area == 100


def test_buffer():
    bounds = Bounds(0, 0, 10, 10)
    bounds.buffer(1)
    assert bounds.xmin == -1
    assert bounds.ymin == -1
    assert bounds.xmax == 11
    assert bounds.ymax == 11


def test_union():
    bounds = Bounds(0, 0, 10, 10)
    bounds.union(Bounds(5, 5, 15, 15))
    assert bounds.xmin == 0
    assert bounds.ymin == 0
    assert bounds.xmax == 15
    assert bounds.ymax == 15


def test_intersection():
    bounds = Bounds(0, 0, 10, 10)
    bounds.intersect(Bounds(5, 5, 15, 15))
    assert bounds.xmin == 5
    assert bounds.ymin == 5
    assert bounds.xmax == 10
    assert bounds.ymax == 10


if __name__ == "__main__":
    pytest.main()
