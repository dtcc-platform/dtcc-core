"""Small geometry workflows through the exported builder API."""

import numpy as np
import pytest
from shapely.geometry import Polygon, box

from dtcc_core import builder
from dtcc_core.builder.meshing import backends
from dtcc_core.model import Bounds, Building, City, GeometryType, Surface


def building(polygon, identifier="building"):
    surface = Surface()
    surface.from_polygon(polygon, 7.0)
    result = Building(id=identifier)
    result.add_geometry(surface, GeometryType.LOD0)
    return result


def test_split_footprint_walls_preserves_shape_height_and_source():
    original = building(box(0, 0, 12, 6))
    vertices = original.lod0.vertices.copy()
    result, = builder.split_footprint_walls([original], max_wall_length=3)
    assert result.lod0.to_polygon().equals(original.lod0.to_polygon())
    edges = np.roll(result.lod0.vertices, -1, axis=0) - result.lod0.vertices
    assert np.linalg.norm(edges, axis=1).max() <= 3 + 1e-9
    np.testing.assert_allclose(result.lod0.vertices[:, 2], 7)
    np.testing.assert_array_equal(original.lod0.vertices, vertices)
    with pytest.raises(RuntimeError, match="max_wall_length"):
        builder.split_footprint_walls([original], max_wall_length=[])


def test_merge_buildings_joins_nearby_footprints_without_mutating_city():
    city = City()
    city.add_building(building(box(0, 0, 4, 4), "left"))
    city.add_building(building(box(4.1, 0, 8.1, 4), "right"))
    result = builder.merge_buildings(city, max_distance=0.2, min_area=0, simplify=False)
    assert len(result.buildings) == 1
    merged = result.buildings[0].lod0.to_polygon()
    assert merged.is_valid and merged.area >= 32
    assert len(city.buildings) == 2
    assert [b.lod0.to_polygon().area for b in city.buildings] == [16, 16]


def test_fix_building_clearance_removes_small_footprint_defect():
    polygon = Polygon([(0, 0), (4, 0), (4, 4), (2.01, 4), (2, 3.99), (0, 4)])
    city = City()
    city.add_building(building(polygon))
    result = builder.fix_building_clearance(city, target_clearance=0.2, min_angle=20)
    fixed = result.buildings[0].lod0.to_polygon(simplify=0)
    assert fixed.is_valid and fixed.minimum_clearance >= 0.18
    assert city.buildings[0].lod0.to_polygon(simplify=0).equals(polygon)


def test_clean_building_surfaces_removes_redundant_vertices():
    city = City()
    city.add_building(building(Polygon([(0, 0), (2, 0), (4, 0), (4, 4), (0, 4)])))
    result = builder.clean_building_surfaces(city, GeometryType.LOD0)
    surface = result.buildings[0].lod0
    assert len(surface.vertices) == 4
    assert surface.to_polygon().symmetric_difference(box(0, 0, 4, 4)).area < 1e-12
    np.testing.assert_allclose(surface.vertices[:, 2], 7)


def test_flat_terrain_uses_requested_height_and_extent():
    bounds = Bounds(100, 200, 130, 240)
    terrain = builder.flat_terrain(-2.5, bounds)
    np.testing.assert_array_equal(terrain.raster.data, [[-2.5]])
    assert terrain.raster.bounds.tuple == bounds.tuple
    assert terrain.raster.georef * (0.5, 0.5) == (115, 220)


def test_default_mesher_override_reset_and_invalid_selection(monkeypatch):
    monkeypatch.setattr(backends, "_default_2d_mesher_override", None)
    monkeypatch.setattr(backends, "available_2d_meshers", lambda: ["dtcc_mesher", "triangle"])
    monkeypatch.setenv("DTCC_2D_MESHER", "dtcc_mesher")
    assert builder.get_default_2d_mesher() == "dtcc_mesher"
    assert builder.set_default_2d_mesher(" TRIANGLE ") == "triangle"
    assert builder.get_default_2d_mesher() == "triangle"
    with pytest.raises(ValueError, match="Unsupported 2D mesher"):
        builder.set_default_2d_mesher("invalid")
    assert builder.get_default_2d_mesher() == "triangle"
    assert builder.set_default_2d_mesher(None) == "dtcc_mesher"
