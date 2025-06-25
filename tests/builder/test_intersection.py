import pytest
import numpy as np
import dtcc_core.builder
from dtcc_core.model import Surface, MultiSurface, Mesh

from dtcc_core.builder.geometry.multisurface import (
    ray_intersection as multisurface_ray_intersection,
)
from dtcc_core.builder.geometry.surface import (
    ray_intersection as surface_ray_intersection,
)


@pytest.fixture
def square_surface():
    return Surface(vertices=np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))


@pytest.fixture
def square_surface_elevated():
    return Surface(vertices=np.array([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]))


@pytest.fixture
def two_layer_surface(square_surface, square_surface_elevated):
    return MultiSurface(surfaces=[square_surface, square_surface_elevated])


@pytest.fixture
def test_origin():
    return [0.5, 0.5, 2]


@pytest.mark.parametrize(
    "direction, expected_intersection, test_id",
    [([0, 0, -1], [0.5, 0.5, 0], "hit"), ([0, 0, 1], None, "miss")],
)
def test_surface_ray_intersection(
    square_surface, test_origin, direction, expected_intersection, test_id
):
    intersection = surface_ray_intersection(square_surface, test_origin, direction)
    assert isinstance(intersection, np.ndarray)

    if test_id == "hit":
        assert np.array_equal(intersection, expected_intersection)
    else:
        assert np.isnan(intersection).all()


@pytest.mark.parametrize(
    "direction, expected_intersection, test_id",
    [([0, 0, -1], [0.5, 0.5, 1], "hit"), ([0, 0, 1], None, "miss")],
)
def test_multisurface_ray_intersection(
    two_layer_surface, test_origin, direction, expected_intersection, test_id
):
    intersection = multisurface_ray_intersection(
        two_layer_surface, test_origin, direction
    )
    assert isinstance(intersection, np.ndarray)

    if test_id == "hit":
        assert np.array_equal(intersection, expected_intersection)
    else:
        assert np.isnan(intersection).all()


if __name__ == "__main__":
    pytest.main()
