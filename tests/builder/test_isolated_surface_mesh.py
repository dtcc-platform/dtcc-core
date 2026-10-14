"""Surface meshes of buildings cleaned in isolation: terrain plus separate shells."""

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import io
from dtcc_core.builder import (
    build_city_flat_mesh,
    build_city_surface_mesh,
    build_city_volume_mesh,
)
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.meshing.tetgen import is_watertight_surface
from dtcc_core.model import Bounds, Building, City, GeometryType, Raster, Surface, Terrain

pytest.importorskip("dtcc_mesher")


def building(id, polygon, height, ground=0.0):
    surface = Surface()
    surface.from_polygon(polygon, ground + height)
    value = Building(id=id)
    value.add_geometry(surface, GeometryType.LOD0)
    value.estimated_height = float(height)
    value.attributes["ground_height"] = float(ground)
    return value


def sloped_city(*buildings):
    city = City()
    raster = Raster()
    raster.data = np.tile(0.15 * np.linspace(0, 100, 51), (31, 1))
    raster.set_bounds(Bounds(-10, -10, 90, 50))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(list(buildings))
    return city


def field(mesh, name):
    return next(f for f in mesh.fields if f.name == name)


def overlapping_city():
    return sloped_city(
        building("a", box(0, 0, 10, 10), 10, ground=1.5),
        building("b", box(9.8, 0, 20, 10), 20, ground=3.0),  # overlaps a
        building("c", box(40, 0, 60, 10), 8, ground=7.5),
    )


def test_buildings_cleaned_in_isolation_get_separate_shells_on_the_terrain():
    city = overlapping_city()
    mesh = build_city_surface_mesh(city, isolated=True, report_mesh_quality=False)

    owner = field(mesh, "building").values
    markers = np.asarray(mesh.markers)
    # The terrain is not cut by the footprints; a and b are not merged.
    assert set(np.unique(owner)) == {-1.0, 0.0, 1.0, 2.0}
    assert set(markers[owner == -1]) == {-2}
    # Three regions: walls 0-2, roofs 3-5.
    assert sorted(set(markers[owner >= 0])) == list(range(6))
    assert field(mesh, "watertight").association == "geometry"
    assert not is_watertight_surface(mesh)

    # Walls reach down to the lowest terrain under each footprint.
    vertices = np.asarray(mesh.vertices)
    faces = np.asarray(mesh.faces)
    raster = city.terrain.raster
    for index, polygon in enumerate([box(0, 0, 10, 10), box(9.8, 0, 20, 10), box(40, 0, 60, 10)]):
        bottom = vertices[faces[owner == index]][:, :, 2].min()
        corners = [raster.get_value(x, y) for x, y in polygon.exterior.coords]
        assert bottom <= min(corners) + 1e-9


def test_unmerged_soup_is_the_terrain_and_one_shell_per_building(tmp_path):
    meshes_list = build_city_surface_mesh(
        overlapping_city(), isolated=True, merge_meshes=False, report_mesh_quality=False
    )
    assert len(meshes_list) == 4
    assert all(not is_watertight_surface(mesh) for mesh in meshes_list)
    assert [np.unique(field(mesh, "building").values).tolist() for mesh in meshes_list] == [
        [-1.0],
        [0.0],
        [1.0],
        [2.0],
    ]

    io.save_model(meshes_list[1], tmp_path / "shell.dtcc")
    assert not is_watertight_surface(io.load_model(tmp_path / "shell.dtcc"))


def test_a_block_gets_one_shell_with_step_walls_and_a_recorded_cleaning_is_reused(monkeypatch):
    city = sloped_city(
        building("low", box(10, 10, 30, 25), 6),
        building("high", box(30, 10, 40, 25), 30),
    )
    block_city = city.clean_footprints()  # one block with two parts
    cleaned = block_city.clean_footprints(isolated=True)

    def forbidden(*args, **kwargs):
        raise AssertionError("cleaned again")

    monkeypatch.setattr(meshes, "clean_footprints", forbidden)
    mesh = build_city_surface_mesh(cleaned, isolated=True, report_mesh_quality=False)
    owner = field(mesh, "building").values
    markers = np.asarray(mesh.markers)
    assert set(np.unique(owner)) == {-1.0, 0.0}
    # Two regions in one shell: walls 0-1, roofs 2-3.
    assert sorted(set(markers[owner == 0])) == [0, 1, 2, 3]
    roofs = mesh.vertices[mesh.faces[markers >= 2]][:, :, 2]
    assert sorted(set(np.round(roofs.ravel(), 9))) == [6.0, 30.0]


def test_flat_and_volume_meshes_need_buildings_cleaned_in_unison():
    city = overlapping_city()
    with pytest.raises(ValueError, match="cleaned in unison"):
        build_city_flat_mesh(city, isolated=True)
    with pytest.raises(ValueError, match="cleaned in unison"):
        build_city_volume_mesh(city, isolated=True)
    with pytest.raises(TypeError):
        build_city_surface_mesh(city, isolated=1)
    assert is_watertight_surface(
        build_city_surface_mesh(city, report_mesh_quality=False, max_mesh_size=10.0)
    )


def test_volume_meshing_refuses_a_soup():
    from dtcc_core.builder.meshing.tetgen import build_volume_mesh, is_tetgen_available

    if not is_tetgen_available():
        pytest.skip("TetGen is not available")
    mesh = build_city_surface_mesh(overlapping_city(), isolated=True, report_mesh_quality=False)
    with pytest.raises(ValueError, match="not watertight"):
        build_volume_mesh(mesh)
