"""Pinched roof levels: reported by every shell builder, one roof for volumes.

A block whose regions at one roof level touch at a point only above a lower
roof (a checkerboard) has a closed shell with edges of four faces. Surface
meshes keep the step walls and warn; a volume mesh meshes such a block with
one roof, so its boundary is a manifold.
"""

from collections import Counter

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core.builder import build_building_shell
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.model import Bounds, Building, City, GeometryType, Raster, Surface, Terrain

CHECKERBOARD = [
    (box(0, 0, 10, 10), 30.0),
    (box(10, 10, 20, 20), 30.0),
    (box(10, 0, 20, 10), 9.0),
    (box(0, 10, 10, 20), 9.0),
]


def checkerboard_city():
    raster = Raster()
    raster.data = np.zeros((25, 25))
    raster.set_bounds(Bounds(-15, -15, 35, 35))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city = City()
    city.add_terrain(terrain)
    for i, (polygon, height) in enumerate(CHECKERBOARD):
        surface = Surface()
        surface.from_polygon(polygon, height)
        building = Building(id=str(i))
        building.add_geometry(surface, GeometryType.LOD0)
        building.estimated_height = height
        building.attributes["ground_height"] = 0.0
        city.add_building(building)
    return city


def edges_with_four_faces(mesh):
    counts = Counter()
    for face in np.asarray(mesh.faces):
        for i in range(3):
            counts[tuple(sorted((int(face[i]), int(face[(i + 1) % 3]))))] += 1
    return sum(count > 2 for count in counts.values())


def test_shells_of_a_pinched_block_are_reported_by_both_builders(caplog):
    cleaned = checkerboard_city().clean_footprints()
    (block,) = cleaned.buildings
    shell = build_building_shell(block)
    assert edges_with_four_faces(shell) > 0
    assert "touch at a point only" in caplog.text
    caplog.clear()
    footprints = meshes._meshing_footprints(
        cleaned, lod=None, min_building_detail=0.5, min_building_area=0.0, max_mesh_size=5.0
    )
    assert footprints.diagnostics["levels"]["pinched_level_sets"]
    assert any("touch at a point only" in message for message in footprints.contract["warnings"])
    assert len({surface.zmax for surface in footprints.surfaces}) == 2


def test_a_volume_mesh_gives_a_pinched_block_one_roof(caplog):
    cleaned = checkerboard_city().clean_footprints()
    footprints = meshes._meshing_footprints(
        cleaned,
        lod=None,
        min_building_detail=0.5,
        min_building_area=0.0,
        max_mesh_size=5.0,
        volume=True,
    )
    levels = footprints.diagnostics["levels"]
    assert levels["pinched_level_sets"] == []
    (block,) = cleaned.buildings
    assert levels["one_roof"] == [
        {"building": 0, "ids": [block.id], "roof": 30.0, "levels": [30.0]}
    ]
    # Roofs spanning more than 2 m: the highest, so no source is cut down.
    assert [surface.zmax for surface in footprints.surfaces] == [30.0]
    assert "meshed with one roof" in caplog.text
    warnings = footprints.contract["warnings"]
    assert any("meshed with one roof each" in message for message in warnings)
    assert not any("vertical line contacts" in message for message in warnings)


def test_volume_mesh_of_a_pinched_block():
    from dtcc_core.builder import build_city_volume_mesh
    from dtcc_core.builder.meshing.tetgen import is_tetgen_available

    if not is_tetgen_available():
        pytest.skip("TetGen is not available")
    pytest.importorskip("dtcc_mesher")
    mesh = build_city_volume_mesh(
        checkerboard_city(),
        max_mesh_size=5.0,
        domain_height=60.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )
    points = np.asarray(mesh.vertices)[np.asarray(mesh.cells)]
    a, b, c = (points[:, i] - points[:, 0] for i in (1, 2, 3))
    volume = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))).sum() / 6.0
    assert volume == pytest.approx(50 * 50 * 60 - 400 * 30, rel=1e-9)
