"""LOD1 per part and closed shells of cleaned blocks."""

from collections import Counter

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
