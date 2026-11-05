import pytest

from dtcc_core.io.data.bbox import as_bounds
from dtcc_core.model import Bounds


def test_as_bounds_keeps_a_bounds_object():
    bounds = Bounds(1.0, 2.0, 3.0, 4.0)
    assert as_bounds(bounds) is bounds


def test_as_bounds_reads_four_values_as_min_and_max_xy():
    bounds = as_bounds((317000, 6396000, 321000, 6400000))
    assert (bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax) == (
        317000.0,
        6396000.0,
        321000.0,
        6400000.0,
    )


def test_as_bounds_reads_six_values_in_dataset_order():
    bounds = as_bounds([317000, 6396000, 0, 321000, 6400000, 100])
    assert (bounds.xmin, bounds.ymin, bounds.zmin) == (317000.0, 6396000.0, 0.0)
    assert (bounds.xmax, bounds.ymax, bounds.zmax) == (321000.0, 6400000.0, 100.0)


def test_as_bounds_accepts_a_numpy_array():
    import numpy as np

    bounds = as_bounds(np.array([317000.0, 6396000.0, 321000.0, 6400000.0]))
    assert (bounds.xmax, bounds.ymax) == (321000.0, 6400000.0)


@pytest.mark.parametrize(
    "bounds", [(1, 2, 3), "1,2,3,4", None, 5.0, {"xmin": 1}, ("a", 2, 3, 4)]
)
def test_as_bounds_rejects_other_input(bounds):
    with pytest.raises(TypeError, match="4 values"):
        as_bounds(bounds)


@pytest.mark.parametrize(
    "bounds",
    [
        (321000, 6396000, 317000, 6400000),
        # A six-value box in Bounds order: x0, y0, x1, y1, z0, z1.
        (317000, 6396000, 321000, 6400000, 0, 50),
    ],
)
def test_as_bounds_rejects_a_minimum_above_its_maximum(bounds):
    with pytest.raises(ValueError, match="minimum below its maximum"):
        as_bounds(bounds)
