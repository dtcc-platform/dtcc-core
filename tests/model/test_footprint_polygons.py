import pytest
import numpy as np
from shapely.geometry import Polygon

from dtcc_core.model import City, Building, Surface
from dtcc_core.model.object import GeometryType


def make_building(polygon: Polygon, height: float = 10.0) -> Building:
    """Create a Building with a LOD0 surface from a shapely Polygon."""
    surface = Surface()
    surface.from_polygon(polygon, height)
    b = Building()
    b.add_geometry(surface, GeometryType.LOD0)
    return b


def make_square(x=0.0, y=0.0, size=10.0) -> Polygon:
    """Create a simple square polygon at (x, y) with given size."""
    return Polygon([
        (x, y),
        (x + size, y),
        (x + size, y + size),
        (x, y + size),
    ])


def test_empty_city_returns_empty_list():
    city = City()
    assert city.footprint_polygons() == []


def test_single_building():
    city = City()
    city.add_building(make_building(make_square()))
    result = city.footprint_polygons()
    assert len(result) == 1
    assert isinstance(result[0], Polygon)


def test_multiple_buildings():
    city = City()
    city.add_building(make_building(make_square(0, 0)))
    city.add_building(make_building(make_square(20, 20)))
    city.add_building(make_building(make_square(40, 40)))
    result = city.footprint_polygons()
    assert len(result) == 3
    assert all(isinstance(p, Polygon) for p in result)


def test_building_without_geometry_is_skipped():
    city = City()
    city.add_building(make_building(make_square()))
    city.add_building(Building())  # no geometry
    result = city.footprint_polygons()
    assert len(result) == 1


def test_all_buildings_without_geometry():
    city = City()
    city.add_building(Building())
    city.add_building(Building())
    result = city.footprint_polygons()
    assert len(result) == 0


def test_returned_polygons_are_2d():
    city = City()
    city.add_building(make_building(make_square(), height=50.0))
    result = city.footprint_polygons()
    coords = np.array(result[0].exterior.coords)
    # Shapely 2D polygons have 2 coordinate dimensions
    assert coords.shape[1] == 2


def test_returned_polygon_area_is_correct():
    size = 10.0
    city = City()
    city.add_building(make_building(make_square(size=size)))
    result = city.footprint_polygons()
    assert result[0].area == pytest.approx(size * size, rel=0.1)


if __name__ == "__main__":
    pytest.main([__file__])
