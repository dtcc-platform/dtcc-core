import pytest
import numpy as np
from dtcc_core.model import Surface, MultiSurface


@pytest.fixture
def surface_pair():
    s1 = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
    s2 = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [1, 1, 2]]))
    return s1, s2


@pytest.fixture
def multi_surface(surface_pair):
    s1, s2 = surface_pair
    return MultiSurface(surfaces=[s1, s2])


def test_zmax(multi_surface):
    assert multi_surface.zmax == 2


def test_to_proto(multi_surface):
    pb = multi_surface.to_proto()
    assert len(pb.multi_surface.surfaces) == 2


def test_from_proto(multi_surface):
    # Convert to proto and back
    pb = multi_surface.to_proto()
    pb_string = pb.SerializeToString()

    ms2 = MultiSurface()
    ms2.from_proto(pb_string)

    assert len(ms2.surfaces) == 2
    assert ms2.zmax == 2


if __name__ == "__main__":
    pytest.main()
