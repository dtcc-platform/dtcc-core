"""Building shells: step walls between touching buildings and closure checks."""

import numpy as np
import pytest
from shapely.geometry import box
from shapely.ops import unary_union

from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.model import Bounds, Mesh, Raster, Surface

BOUNDS = (0.0, 0.0, 50.0, 35.0)


def shell(polygons, heights, *, terrain=None, merge_meshes=True):
    """Surface mesh of buildings given as touching ground-mesh regions.

    The cleaner merges touching buildings, so the shell builder is called
    directly with one region per building.
    """
    pytest.importorskip("dtcc_mesher")
    ground = box(*BOUNDS).difference(unary_union(polygons))
    grounds = list(getattr(ground, "geoms", [ground]))
    regions = [*grounds, *polygons]
    ground_mesh, _ = meshes._build_ground_mesh_from_coverage(
        region_polygons=regions,
        region_markers=[-2] * len(grounds) + list(range(len(polygons))),
        region_points=[
            np.asarray(region.representative_point().coords[0]) for region in regions
        ],
        bounds=BOUNDS,
        max_mesh_size=5.0,
        min_mesh_angle=25.0,
        mesher="dtcc_mesher",
        sort_triangles=False,
        region_triangle_sizes={i: 5.0 for i in range(len(polygons))},
        add_halo_markers=False,
    )
    raster = Raster()
    raster.data = np.zeros((8, 8)) if terrain is None else terrain
    raster.set_bounds(Bounds(*BOUNDS))
    surfaces = []
    for polygon, height in zip(polygons, heights):
        surface = Surface()
        surface.from_polygon(polygon, height)
        surfaces.append(surface)
    return meshes._build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=raster,
        building_surfaces=surfaces,
        meshing_directives=[1] * len(polygons),
        smoothing=0,
        merge_meshes=merge_meshes,
    )


def wall_area_on_plane(mesh, x):
    """Area of faces in the plane x = const, by marker."""
    vertices = np.asarray(mesh.vertices)
    areas = {}
    for face, marker in zip(np.asarray(mesh.faces), np.asarray(mesh.markers)):
        points = vertices[face]
        if np.allclose(points[:, 0], x):
            area = 0.5 * np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0]))
            areas[int(marker)] = areas.get(int(marker), 0.0) + area
    return areas


def enclosed_volume(mesh):
    """Volume between the surface and z = 0 (walls are vertical)."""
    points = np.asarray(mesh.vertices)[np.asarray(mesh.faces)]
    a = points[:, 1, :2] - points[:, 0, :2]
    b = points[:, 2, :2] - points[:, 0, :2]
    projected = 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])
    return float(np.sum(projected * points[:, :, 2].mean(axis=1)))


def closed(mesh, num_buildings):
    report = meshes._shell_closure_audit(mesh, num_buildings=num_buildings)
    return report["closed"], report


def test_touching_buildings_get_one_step_wall_owned_by_the_taller():
    low, high = box(10, 10, 30, 25), box(30, 10, 40, 25)
    mesh = shell([low, high], [6.0, 30.0])

    is_closed, report = closed(mesh, 2)
    assert is_closed, report
    # Only the taller building (wall marker 1) has a wall on x = 30, from the
    # lower roof to its own: 15 m x 24 m.
    assert wall_area_on_plane(mesh, 30.0) == {1: pytest.approx(360.0)}
    assert enclosed_volume(mesh) == pytest.approx(300 * 6 + 150 * 30)


def test_equal_roofs_have_no_wall_between_them():
    mesh = shell([box(10, 10, 30, 25), box(30, 10, 40, 25)], [12.0, 12.0])
    is_closed, report = closed(mesh, 2)
    assert is_closed, report
    assert wall_area_on_plane(mesh, 30.0) == {}


def test_terrace_and_tower_in_podium_close():
    terrace = [box(5 + 10 * i, 10, 15 + 10 * i, 25) for i in range(4)]
    mesh = shell(terrace, [9.0, 9.0, 15.0, 6.0])
    is_closed, report = closed(mesh, 4)
    assert is_closed, report
    assert wall_area_on_plane(mesh, 15.0) == {}  # equal neighbours
    assert wall_area_on_plane(mesh, 25.0) == {2: pytest.approx(15 * 6)}
    assert wall_area_on_plane(mesh, 35.0) == {2: pytest.approx(15 * 9)}

    tower = box(15, 10, 27, 22)
    podium = box(5, 5, 45, 30).difference(tower)
    mesh = shell([podium, tower], [10.0, 45.0])
    is_closed, report = closed(mesh, 2)
    assert is_closed, report
    assert enclosed_volume(mesh) == pytest.approx(podium.area * 10 + tower.area * 45)


def test_touching_buildings_share_one_base_on_sloped_terrain():
    xs = np.linspace(0, 50, 8)
    terrain = np.tile(0.2 * xs, (8, 1)) + np.linspace(0, 3, 8)[:, None]
    mesh = shell([box(10, 10, 30, 25), box(30, 10, 40, 25)], [20.0, 40.0], terrain=terrain)
    is_closed, report = closed(mesh, 2)
    assert is_closed, report
    markers = np.asarray(mesh.markers)
    faces = np.asarray(mesh.faces)
    vertices = np.asarray(mesh.vertices)
    bases = [vertices[np.unique(faces[markers == m])][:, 2].min() for m in (0, 1)]
    assert bases[0] == bases[1]


def test_unmerged_building_meshes_weld_into_a_closed_shell():
    meshes_list = shell(
        [box(10, 10, 30, 25), box(30, 10, 40, 25)], [6.0, 30.0], merge_meshes=False
    )
    assert len(meshes_list) == 3
    is_closed, report = closed(meshes._welded_mesh(meshes_list), 2)
    assert is_closed, report


def test_closure_audit_flags_overlapping_walls_and_cracks():
    # A roof of building 1 over x in [30, 40] at 30 m, and a 6 m wall of
    # building 0 on x = 30 facing +x, i.e. into building 1 below its roof.
    vertices = np.array(
        [
            [30, 10, 30], [40, 10, 30], [40, 25, 30], [30, 25, 30],  # roof of 1
            [30, 10, 0], [30, 25, 0], [30, 10, 6], [30, 25, 6],  # wall of 0
        ],
        dtype=float,
    )
    faces = np.array([[0, 1, 2], [0, 2, 3], [4, 5, 7], [4, 7, 6]])
    markers = np.array([3, 3, 0, 0])  # roofs are 2 + marker
    report = meshes._shell_closure_audit(
        Mesh(vertices=vertices, faces=faces, markers=markers), num_buildings=2
    )
    assert report["buried_wall_faces"] == 2
    assert not report["closed"]

    mesh = shell([box(10, 10, 30, 25)], [10.0])
    faces = np.asarray(mesh.faces)
    markers = np.asarray(mesh.markers)
    keep = np.arange(len(faces)) != int(np.flatnonzero(markers == 0)[0])
    report = meshes._shell_closure_audit(
        Mesh(vertices=mesh.vertices, faces=faces[keep], markers=markers[keep]),
        num_buildings=1,
    )
    # The removed face's two edges shared with other wall faces are open; its
    # base edge now belongs to a terrain face only and is not counted.
    assert report["open_interior_edges"] == 2
    errors, _ = meshes._shell_closure_messages(report)
    assert any("open edge" in message for message in errors)


def test_surface_shell_contract_rejects_an_open_or_buried_shell():
    audit = {
        "num_vertices": 3,
        "num_faces": 1,
        "num_markers": 1,
        "marker_count": 1,
        "closure": {
            "open_interior_edges": 0,
            "nonmanifold_edges": 1,
            "inconsistent_edges": 0,
            "buried_wall_faces": 4,
            "closed": False,
        },
    }
    contract = meshes._surface_shell_contract_from_audit(
        audit, reference_length=0.5, preserve_surface_requested=False
    )
    assert contract["requirements"]["shell_closed"] is False
    assert any("inside another building" in message for message in contract["errors"])
    assert any("more than two faces" in message for message in contract["warnings"])


TOUCHING_VOLUME_CASES = {
    "step": ([box(10, 10, 30, 25), box(30, 10, 40, 25)], [6.0, 30.0]),
    "equal": ([box(10, 10, 30, 25), box(30, 10, 40, 25)], [12.0, 12.0]),
    "terrace": ([box(5 + 10 * i, 10, 15 + 10 * i, 25) for i in range(4)], [9.0, 9.0, 15.0, 6.0]),
    "tower_in_podium": (
        [box(5, 5, 45, 30).difference(box(15, 10, 27, 22)), box(15, 10, 27, 22)],
        [10.0, 45.0],
    ),
}


def city_of_touching(polygons, heights, *, cell=2.0):
    """A city of the buildings on flat ground, with a margin of 15 m."""
    from dtcc_core.model import Building, City, GeometryType, Terrain

    xmin, ymin, xmax, ymax = unary_union(polygons).buffer(15.0, join_style=2).bounds
    raster = Raster()
    raster.data = np.zeros((int(round((ymax - ymin) / cell)), int(round((xmax - xmin) / cell))))
    raster.set_bounds(Bounds(xmin, ymin, xmax, ymax))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city = City()
    city.add_terrain(terrain)
    for polygon, height in zip(polygons, heights):
        surface = Surface()
        surface.from_polygon(polygon, height)
        building = Building()
        building.add_geometry(surface, GeometryType.LOD0)
        building.estimated_height = height
        city.add_building(building)
    return city, (xmin, ymin, xmax, ymax)


def test_merged_buildings_get_a_step_wall_through_the_cleaner():
    from dtcc_core.builder import build_city_surface_mesh

    pytest.importorskip("dtcc_mesher")
    low, high = box(10, 10, 30, 25), box(30, 10, 40, 25)
    city, _ = city_of_touching([low, high], [6.0, 30.0])
    mesh = build_city_surface_mesh(
        city, max_mesh_size=5.0, report_mesh_quality=False, mesher="dtcc_mesher"
    )
    is_closed, report = closed(mesh, 2)
    assert is_closed, report
    walls = wall_area_on_plane(mesh, 30.0)
    assert list(walls.values()) == [pytest.approx(360.0)]
    assert enclosed_volume(mesh) == pytest.approx(300 * 6 + 150 * 30)


@pytest.mark.parametrize("name", sorted(TOUCHING_VOLUME_CASES))
def test_volume_mesh_of_touching_buildings(name):
    """TetGen meshes the air around merged buildings with step walls.

    The buildings are merged into one block by the cleaner, and the block is
    meshed per part at the parts' own roofs. With overlapping full-height walls
    (before step walls) TetGen rejected these shells as self-intersecting.
    """
    from dtcc_core.builder import build_city_volume_mesh
    from dtcc_core.builder.meshing.tetgen import is_tetgen_available

    if not is_tetgen_available():
        pytest.skip("TetGen is not available")
    pytest.importorskip("dtcc_mesher")
    polygons, heights = TOUCHING_VOLUME_CASES[name]
    city, (xmin, ymin, xmax, ymax) = city_of_touching(polygons, heights)

    domain_height = 60.0
    volume_mesh = build_city_volume_mesh(
        city,
        max_mesh_size=5.0,
        domain_height=domain_height,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    points = np.asarray(volume_mesh.vertices)[np.asarray(volume_mesh.cells)]
    a, b, c = (points[:, i] - points[:, 0] for i in (1, 2, 3))
    volume = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))).sum() / 6.0
    air = (xmax - xmin) * (ymax - ymin) * domain_height - sum(
        polygon.area * height for polygon, height in zip(polygons, heights)
    )
    assert volume == pytest.approx(air, rel=1e-9)
