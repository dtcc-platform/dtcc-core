"""Roof levels of merged blocks at the vertical resolution, and meshing per part."""

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core.builder.cleaning.levels import (
    group_levels,
    level_regions,
    pinched_levels,
)
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.model import (
    Bounds,
    Building,
    City,
    GeometryType,
    Raster,
    Surface,
    Terrain,
)


def test_grouped_levels_move_at_most_half_the_resolution_and_keep_their_order():
    rng = np.random.default_rng(0)
    for _ in range(2000):
        n = int(rng.integers(1, 9))
        delta = float(rng.choice([0.25, 0.5, 1.0]))
        steps = rng.uniform(0, 1.5 * delta, n)
        if rng.random() < 0.5:
            # Levels on a grid of delta / 4 hit differences of exactly delta.
            steps = np.floor(steps / (delta / 4)) * (delta / 4)
        levels = np.cumsum(steps) + 5.0
        rng.shuffle(levels)
        areas = rng.uniform(1, 500, n)
        grouped = np.asarray(group_levels(levels, areas, delta))
        assert np.max(np.abs(grouped - levels)) <= delta / 2 + 1e-12
        order = np.argsort(levels, kind="stable")
        assert np.all(np.diff(grouped[order]) >= -1e-12)
        distinct = np.unique(np.round(grouped, 9))
        assert np.all(np.diff(distinct) >= delta - 1e-9)
        # Levels at least delta apart are never unified.
        for i in range(n):
            for j in range(n):
                if abs(levels[i] - levels[j]) >= delta:
                    assert grouped[i] != grouped[j]


def test_grouped_levels_example_and_exact_values():
    levels = [9.0, 9.3, 12.0, 12.2, 12.6, 15.0]
    assert group_levels(levels, [100] * 6, 0.5) == pytest.approx(
        [9.15, 9.15, 12.1, 12.1, 12.6, 15.0]
    )
    assert group_levels([12.0, 12.0, 30.0], [1, 2, 3], 0.5) == [12.0, 12.0, 30.0]
    assert group_levels([12.0, 12.3], [1, 2], 0.0) == [12.0, 12.3]
    # Levels exactly delta apart stay apart.
    assert group_levels([9.0, 9.5, 10.0], [1, 1, 1], 0.5) == [9.0, 9.5, 10.0]


def test_cells_of_one_level_are_dissolved():
    cells = [(box(5 + 10 * i, 10, 15 + 10 * i, 25), i) for i in range(4)]
    regions = level_regions(cells, {0: 9.0, 1: 9.2, 2: 15.0, 3: 6.0}, 0.5)
    assert [(round(level, 9), labels) for _, level, labels in regions] == [
        (6.0, [3]),
        (9.1, [0, 1]),
        (15.0, [2]),
    ]
    assert regions[1][0].equals(box(5, 10, 25, 25))


def test_point_contacts_of_a_level_set_are_reported():
    # Two tall corners of a 2 x 2 block touch at its centre.
    cells = [
        (box(0, 0, 10, 10), 0),
        (box(10, 10, 20, 20), 1),
        (box(10, 0, 20, 10), 2),
        (box(0, 10, 10, 20), 3),
    ]
    regions = level_regions(cells, {0: 30.0, 1: 30.0, 2: 9.0, 3: 9.0}, 0.5)
    pinched = pinched_levels([(region, level) for region, level, _ in regions])
    assert pinched == [{"level": 30.0, "points": [[10.0, 10.0]]}]
    assert pinched_levels([(box(0, 0, 10, 10), 9.0), (box(10, 0, 20, 10), 12.0)]) == []


def city_of(polygons, roofs):
    city = City()
    raster = Raster()
    raster.data = np.zeros((20, 20))
    raster.set_bounds(Bounds(-15, -15, 35, 35))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    for polygon, roof in zip(polygons, roofs):
        surface = Surface()
        surface.from_polygon(polygon, roof)
        building = Building()
        building.add_geometry(surface, GeometryType.LOD0)
        building.estimated_height = roof
        city.add_building(building)
    return city


def test_repairable_point_contacts_keep_the_two_roofs_in_meshing_footprints():
    city = city_of(
        [
            box(0, 0, 10, 10),
            box(10, 10, 20, 20),
            box(10, 0, 20, 10),
            box(0, 10, 10, 20),
        ],
        [30.0, 30.0, 9.0, 9.0],
    )
    footprints = meshes._meshing_footprints(
        city,
        lod=None,
        min_building_detail=0.5,
        min_building_area=1.0,
        max_mesh_size=5.0,
    )
    assert {surface.zmax for surface in footprints.surfaces} == {9, 30}
    assert footprints.diagnostics["levels"]["pinched_level_sets"] == []
    assert footprints.diagnostics["levels"]["one_roof"] == []
    assert footprints.contract["errors"] == []
    assert not any(
        "vertical line contacts" in w for w in footprints.contract["warnings"]
    )


def test_reused_cleaning_meshes_parts_exactly_as_a_new_cleaning():
    pytest.importorskip("dtcc_mesher")
    city = city_of(
        [box(0, 0, 10, 10), box(10.3, 0, 20, 10), box(20, 0, 30, 8), box(0, 15, 8, 25)],
        [9.0, 21.0, 21.2, 12.0],
    )
    cleaned = city.clean_footprints()
    fresh = meshes.build_city_surface_mesh(
        city, report_mesh_quality=False, mesher="dtcc_mesher"
    )
    reused = meshes.build_city_surface_mesh(
        cleaned, report_mesh_quality=False, mesher="dtcc_mesher"
    )
    assert np.array_equal(fresh.vertices, reused.vertices)
    assert np.array_equal(fresh.faces, reused.faces)
    assert np.array_equal(fresh.markers, reused.markers)
    # Roofs: 9 and 12 m, and 21 and 21.2 m grouped at their area-weighted mean.
    block = cleaned.buildings[0]
    areas = [
        p.get_geometry(GeometryType.LOD0).to_polygon(simplify=0.0).area
        for p in block.building_parts
    ]
    grouped = (21.0 * areas[1] + 21.2 * areas[2]) / (areas[1] + areas[2])
    z = fresh.vertices[fresh.faces][:, :, 2]
    flat = (z.max(axis=1) == z.min(axis=1)) & (z[:, 0] > 0)
    roofs = sorted(set(np.round(z[flat, 0], 9)))
    assert roofs == pytest.approx([9.0, 12.0, grouped])
