"""LOD1 per part and closed shells of cleaned blocks."""

from collections import Counter
import math

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core.builder import build_building_shell, build_lod1_buildings, build_lod2_buildings
from dtcc_core.model import Building, City, GeometryType, MultiSurface, Surface


def building(id, polygon, height, ground=0.0):
    surface = Surface()
    surface.from_polygon(polygon, ground + height)
    value = Building(id=id)
    value.add_geometry(surface, GeometryType.LOD0)
    value.estimated_height = float(height)
    value.attributes["ground_height"] = float(ground)
    return value


def cleaned(*buildings):
    city = City()
    city.add_buildings(list(buildings))
    return city.clean_footprints().buildings


def volume(geometry):
    mesh = geometry.mesh(weld=True) if isinstance(geometry, MultiSurface) else geometry
    points = np.asarray(mesh.vertices)[np.asarray(mesh.faces)]
    return float(np.einsum("ij,ij->i", points[:, 0], np.cross(points[:, 1], points[:, 2])).sum() / 6)


def closed_and_oriented(mesh):
    edges, directed = Counter(), Counter()
    for face in np.asarray(mesh.faces):
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edges[(min(a, b), max(a, b))] += 1
            directed[(a, b)] += 1
    return all(count == 2 for count in edges.values()) and all(
        count == 1 for count in directed.values()
    )


def test_lod1_of_a_block_is_one_prism_per_part():
    block, solo = cleaned(
        building("long", box(0, 0, 40, 10), 6),
        building("tower", box(15, -5, 25, 15), 30, ground=1.0),
        building("solo", box(60, 0, 70, 10), 9),
    )
    build_lod1_buildings([block, solo])

    assert block.lod1 is None
    parts = {part.id: part for part in block.building_parts}
    # Each part from its own ground to its own roof; the long building is
    # split in two by the tower.
    assert volume(parts["long"].lod1) == pytest.approx(300 * 6)
    assert volume(parts["tower"].lod1) == pytest.approx(200 * 30)
    assert len(block.flatten_geometry(GeometryType.LOD1).surfaces) == len(
        parts["long"].lod1.surfaces
    ) + len(parts["tower"].lod1.surfaces)
    assert volume(solo.lod1) == pytest.approx(900)

    # Without rebuild, existing prisms are kept; LOD2 falls back to them.
    prism = parts["tower"].lod1
    build_lod2_buildings([block], rebuild=True, build_lod1_fallback=True)
    assert block.lod1 is None and parts["tower"].lod1 is prism


def test_shell_of_a_block_is_closed_with_step_walls():
    block, solo = cleaned(
        building("long", box(0, 0, 40, 10), 6),
        building("tower", box(15, -5, 25, 15), 30),
        building("solo", box(60, 0, 70, 10), 9),
    )
    shell = build_building_shell(block)
    assert closed_and_oriented(shell)
    assert volume(shell) == pytest.approx(300 * 6 + 200 * 30)
    # Regions: the two pieces of the long building and the tower.
    assert sorted(set(shell.markers.tolist())) == list(range(9))

    assert volume(build_building_shell(solo)) == pytest.approx(900)

    # Without a floor, the walls reach down to the base and the shell is open
    # only along the bottom.
    open_shell = build_building_shell(block, base_z=-2.0, floor=False)
    assert open_shell.vertices[:, 2].min() == -2.0
    assert not closed_and_oriented(open_shell)
    assert set(open_shell.markers.tolist()) == set(range(6))


def test_shell_groups_roofs_closer_than_the_resolution():
    terrace = cleaned(
        *[
            building(str(i), box(5 + 10 * i, 10, 15 + 10 * i, 25), height)
            for i, height in enumerate([9.0, 9.2, 15.0, 6.0])
        ]
    )[0]
    shell = build_building_shell(terrace)
    assert closed_and_oriented(shell)
    assert volume(shell) == pytest.approx(150 * (9.1 * 2 + 15 + 6))
    z = shell.vertices[shell.faces][:, :, 2]
    flat = (z.max(axis=1) == z.min(axis=1)) & (z[:, 0] > 0)
    assert sorted(set(np.round(z[flat, 0], 9))) == pytest.approx([6.0, 9.1, 15.0])


def test_shell_needs_a_roof_above_its_base():
    value = building("low", box(0, 0, 10, 10), 3)
    with pytest.raises(ValueError, match="at or below its base"):
        build_building_shell(value, base_z=5.0)
    with pytest.raises(ValueError, match="no LOD0 footprint"):
        build_building_shell(Building(id="empty"))


def wall_levels(mesh):
    """The z levels of the shell's vertices on each vertical line, by (x, y)."""
    levels = {}
    for x, y, z in np.round(np.asarray(mesh.vertices), 9):
        levels.setdefault((x, y), set()).add(z)
    return levels


def quality(mesh):
    from dtcc_core.model.mixins.mesh.quality import tri_element_quality

    return tri_element_quality(np.asarray(mesh.vertices), np.asarray(mesh.faces))


def test_gap_closing_keeps_compact_cells_through_reuse_and_reload(tmp_path):
    from dtcc_core import io
    from dtcc_core.builder.cleaning.blocks import block_region, set_part_roof
    from dtcc_core.builder.cleaning.levels import level_regions
    from dtcc_core.builder.geometry_builders.meshes import footprint_cells

    city = City()
    city.add_buildings([
        building(str(i), box(x, 0, x + 10, 10), height)
        for i, (x, height) in enumerate([(0, 8), (10.2, 16), (20.4, 12)])
    ])
    once = city.clean_footprints(min_building_detail=0.5, isolated=False)
    io.save_model(once, tmp_path / "compact.dtcc")
    loaded = io.load_model(tmp_path / "compact.dtcc")
    original_cells, _ = footprint_cells(once.buildings[0])
    for candidate in (once, once.clean_footprints(), loaded, loaded.clean_footprints()):
        block = candidate.buildings[0]
        assert len(candidate.buildings) == 1 and len(block.building_parts) == 3
        assert block_region(block).equals(box(0, 0, 30.4, 10))
        cells, roofs = footprint_cells(block)
        assert sorted(roofs.values()) == [8, 12, 16]
        assert sum(len(p.exterior.coords) - 1 for p, _ in cells) <= 16
        for (p, i), (expected, j) in zip(cells, original_cells):
            assert i == j and p.equals_exact(expected, 0)
        for p, _, _ in level_regions(cells, roofs, 0.5):
            assert len(p.exterior.coords) <= 6
        shell = build_building_shell(block, max_mesh_size=5.0)
        # A handful of rectangular panels needs hundreds, not tens of thousands,
        # of triangles at this mesh size; leave ample triangulator variation.
        assert len(shell.faces) < 1500
        assert quality(shell).min() > 0.25
        assert closed_and_oriented(shell)
        assert volume(shell) == pytest.approx(3652)

    # Exposing different roof groups rechecks contacts without restoring samples.
    block = loaded.buildings[0]
    set_part_roof(block.building_parts[0], 12)
    cells, roofs = footprint_cells(block)
    assert sum(len(p.exterior.coords) - 1 for p, _ in cells) <= 16
    assert sorted(roofs.values()) == [12, 12, 16]


@pytest.mark.parametrize("isolated", [False, True])
def test_gap_closing_city_surface_mesh_has_compact_building_constraints(isolated):
    from dtcc_core.builder import build_city_surface_mesh
    from dtcc_core.model import Bounds, Raster, Terrain

    city = City()
    raster = Raster()
    raster.data = np.zeros((4, 8))
    raster.set_bounds(Bounds(-5, -5, 35, 15))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings([
        building(str(i), box(x, 0, x + 10, 10), height)
        for i, (x, height) in enumerate([(0, 8), (10.2, 16), (20.4, 12)])
    ])
    # A jointly cleaned block can also be consumed by the isolated shell workflow.
    city = city.clean_footprints(min_building_detail=0.5)
    mesh = build_city_surface_mesh(city, isolated=isolated, max_mesh_size=5, report_mesh_quality=False)
    assert len(mesh.faces) < 2000
    assert quality(mesh).min() > 0.1
    points = mesh.vertices[mesh.faces]
    flat = np.ptp(points[:, :, 2], axis=1) == 0
    assert sorted(set(np.round(points[flat, 0, 2], 9))) == [0, 8, 12, 16]


def panel_widths(mesh):
    """The widths of the wall panels at each vertical line, by (x, y)."""
    points = np.round(np.asarray(mesh.vertices), 9)[np.asarray(mesh.faces)]
    normals = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0])
    widths = {}
    for triangle, normal in zip(points, normals):
        if abs(normal[2]) > 1e-9 * np.linalg.norm(normal):
            continue
        lines = {(x, y) for x, y, _ in triangle}
        if len(lines) == 2:
            a, b = lines
            for line in lines:
                widths.setdefault(line, set()).add(math.dist(a, b))
    return widths


@pytest.mark.parametrize("notch", [0.5, 0.05, 0.005])
def test_a_short_wall_edge_refines_only_the_walls_at_its_ends(notch):
    from shapely.geometry import Polygon

    footprint = Polygon(
        [(0, 0), (20, 0), (20, 10), (10.5, 10), (10.5, 10 + notch), (0, 10 + notch)]
    )
    shell = build_building_shell(building("notched", footprint, 30), min_building_detail=0.5)
    assert closed_and_oriented(shell)
    # Each vertical line is split into the fewest strips no taller than
    # max(5 times its narrowest panel, 0.5 m); panels narrower than a
    # fiftieth of 0.5 m are left whole.
    levels = wall_levels(shell)
    for line, widths in panel_widths(shell).items():
        refined = [width for width in widths if 50 * width >= 0.5]
        strips = math.ceil(30 / max(5 * min(refined), 0.5)) if refined else 1
        assert len(levels[line]) == strips + 1, (line, sorted(widths))
    if notch >= 0.01:
        assert quality(shell).min() > 0.04
        assert len(levels[(10.5, 10.0)]) >= math.ceil(30 / max(5 * notch, 0.5)) + 1


def test_walls_of_a_block_conform_at_steps_with_unequal_edges():
    block, _ = cleaned(
        building("long", box(0, 0, 40, 10), 6),
        building("tower", box(15, -5, 25.4, 15), 30),
        building("solo", box(60, 0, 70, 10), 9),
    )
    shell = build_building_shell(block, min_building_detail=0.5)
    assert closed_and_oriented(shell)
    assert volume(shell) == pytest.approx(29.6 * 10 * 6 + 10.4 * 20 * 30)
    assert quality(shell).min() > 0.04


def test_a_step_near_the_top_of_one_line_gives_no_sliver():
    # The vertical line at (30, 10) has a step at 29.4 m, 0.6 m below the top
    # of the high part's outline wall, whose other line has no level between
    # its bottom and 30 m (wide panels: the wall edges are 5 m). Joining the
    # short interval at the step to the far bottom of the other line gave a
    # triangle of quality 0.011.
    from dtcc_core.builder.geometry_builders.meshes import build_city_surface_mesh
    from dtcc_core.model import Bounds, Raster, Terrain

    city = City()
    raster = Raster()
    raster.data = np.zeros((4, 4))
    raster.set_bounds(Bounds(0, -20, 60, 40))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(
        [building("high", box(20, 0, 30, 10), 30), building("mid", box(30, 5, 40, 10), 29.4)]
    )
    mesh = build_city_surface_mesh(
        city.clean_footprints(), max_mesh_size=10.0, report_mesh_quality=False
    )
    assert sorted(wall_levels(mesh)[(30.0, 10.0)])[-2:] == [29.4, 30.0]
    assert quality(mesh).min() > 0.05
