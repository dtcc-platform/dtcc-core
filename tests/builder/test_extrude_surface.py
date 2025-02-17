import pytest
import numpy as np
from dtcc_core.builder.geometry_builders.surface import extrude_surface
from dtcc_core.model import Building, MultiSurface, Surface


@pytest.fixture
def simple_surface():
    return Surface(
        vertices=np.array([[0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]])
    )


@pytest.fixture
def extruded_surface(simple_surface):
    return extrude_surface(simple_surface, 2)


def test_extrude_surface_count(extruded_surface):
    assert len(extruded_surface.surfaces) == 6


def test_extrude_surface_height_range(extruded_surface):
    assert pytest.approx(extruded_surface.surfaces[2].vertices[:, 2].max()) == 10
    assert pytest.approx(extruded_surface.surfaces[2].vertices[:, 2].min()) == 2
