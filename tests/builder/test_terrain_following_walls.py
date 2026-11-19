import numpy as np
import pytest
from shapely.affinity import rotate, translate
from shapely.geometry import Polygon, box

from dtcc_core import builder
from dtcc_core.builder.geometry_builders.meshes import (
    _shell_closure_audit,
    _welded_mesh,
)
from dtcc_core.builder.meshing.roof_validation import validate_shell
from dtcc_core.builder.meshing.terrain_attachment import TerrainContact
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.model import (
    Bounds,
    Building,
    City,
    GeometryType,
    Mesh,
    PointCloud,
    Raster,
    Surface,
)


@pytest.fixture(
    params=[("lod1", "catalogue"), ("lod2", "catalogue"), ("lod2", "ransac")]
)
def lod(request):
    return dict(zip(("lod", "lod2_method"), request.param, strict=True))


def ramp(x, y):
    return 5 + 0.2 * x + 0.1 * y


def make_building(outline, roof, name="building", *, fitted_roof=None):
    building = Building(id=name)
    surface = Surface()
    surface.from_polygon(outline, roof)
    building.add_geometry(surface, GeometryType.LOD0)
    building.estimated_height = roof - 5
    building.attributes["ground_height"] = 5.0
    if fitted_roof is not None:
        xmin, ymin, xmax, ymax = outline.bounds
        x, y = np.meshgrid(
            np.linspace(xmin + 0.1, xmax - 0.1, 18),
            np.linspace(ymin + 0.1, ymax - 0.1, 18),
        )
        points = np.c_[x.ravel(), y.ravel(), np.full(x.size, fitted_roof)]
        building.add_geometry(PointCloud(points=points), GeometryType.POINT_CLOUD)
    return building


def make_city(grade, *buildings):
    nodes = np.arange(61) + 0.5
    raster = Raster(data=grade(nodes[None, :], nodes[::-1, None]))
    raster.set_bounds(Bounds(0, 0, 61, 61))
    city = City()
    city.add_terrain(raster)
    city.add_buildings(list(buildings))
    return city


def mesh_city(city, lod, **kwargs):
    return builder.build_city_surface_mesh(
        city, max_mesh_size=4.0, report_mesh_quality=False, **lod, **kwargs
    )


def seam_vertices(mesh):
    terrain = np.unique(mesh.faces[mesh.markers < 0])
    buildings = np.unique(mesh.faces[mesh.markers >= 0])
    return mesh.vertices[np.intersect1d(terrain, buildings)]


def check_shell(mesh):
    count = (int(mesh.markers.max()) + 1) // 2
    assert validate_shell(mesh, num_patches=count, closure_audit=_shell_closure_audit)[
        "valid"
    ]


@pytest.mark.parametrize("slope", [0, 1])
def test_city_walls_follow_grade_with_welded_contacts(lod, slope):
    def grade(x, y):
        return 5 + slope * (ramp(x, y) - 5)

    city = make_city(grade, make_building(box(15, 20, 45, 40), 25))
    mesh = mesh_city(city, lod)
    seam = seam_vertices(mesh)
    assert len(seam) >= 4
    np.testing.assert_allclose(
        seam[:, 2], grade(seam[:, 0], seam[:, 1]), atol=1e-8, rtol=0
    )
    assert np.ptp(seam[:, 2]) == pytest.approx(8 * slope)
    check_shell(mesh)
    separate = _welded_mesh(mesh_city(city, lod, merge_meshes=False))
    check_shell(separate)
    assert sorted(map(tuple, seam_vertices(separate))) == sorted(map(tuple, seam))


def test_grade_bend_in_middle_of_wall(lod):
    def grade(x, y):
        return ramp(x, y) + 0.3 * np.maximum(x - 30.5, 0)

    city = make_city(grade, make_building(box(15, 20, 45, 40), 30))
    mesh = mesh_city(city, lod)
    seam = seam_vertices(mesh)
    along_wall = seam[np.isclose(seam[:, 1], 20)]
    assert len(along_wall) > 2
    np.testing.assert_allclose(
        seam[:, 2], grade(seam[:, 0], seam[:, 1]), atol=1e-8, rtol=0
    )
    steps = np.diff(along_wall[np.argsort(along_wall[:, 0])], axis=0)
    assert np.ptp(steps[:, 2] / steps[:, 0]) > 0.2
    check_shell(mesh)


def test_courtyard_and_roof_step_over_ramp(lod):
    courtyard = box(19, 24, 28, 40)
    outline = Polygon(box(16, 16, 34, 44).exterior.coords, [courtyard.exterior.coords])
    city = make_city(
        ramp,
        make_building(outline, 30, "west"),
        make_building(box(34, 16, 49, 44), 35, "east"),
    )
    mesh = mesh_city(city, lod)
    seam = seam_vertices(mesh)
    in_courtyard = (
        (seam[:, 0] >= 19)
        & (seam[:, 0] <= 28)
        & (seam[:, 1] >= 24)
        & (seam[:, 1] <= 40)
    )
    assert np.count_nonzero(in_courtyard) >= 4
    np.testing.assert_allclose(
        seam[:, 2], ramp(seam[:, 0], seam[:, 1]), atol=1e-8, rtol=0
    )
    assert any(np.isclose(z, 30) and np.isclose(x, 34) for x, _, z in mesh.vertices)
    check_shell(mesh)


@pytest.mark.parametrize("roof", [15.0, 18.0])
def test_roof_crossing_or_touch_at_uphill_corner(lod, roof):
    city = make_city(ramp, make_building(box(15, 20, 45, 40), roof))
    with pytest.raises(
        RuntimeError, match=r"Roof meets or crosses the terrain.*clearance"
    ):
        mesh_city(city, lod)


def test_roof_crossing_between_footprint_corners(lod):
    def grade(x, y):
        return 5 + np.maximum(0, 9 - abs(x - 30.5)) + np.zeros_like(y)

    city = make_city(grade, make_building(box(15, 20, 45, 40), 11))
    with pytest.raises(
        RuntimeError, match=r"Roof meets or crosses the terrain.*clearance"
    ):
        mesh_city(city, lod)


def test_lod2_uses_flat_fallback_if_fitted_roof_crosses_grade():
    city = make_city(ramp, make_building(box(15, 20, 45, 40), 25, fitted_roof=15))
    audit = {}
    mesh = mesh_city(city, dict(lod="lod2"), stage_audit=audit)
    [source] = audit["roof_catalogue"]["sources"]
    assert source["final_outcome"] == "lod1_fallback"
    [failure] = source["admission"]["rejections"]
    assert failure["clearance_m"] == pytest.approx(-3)
    assert max(mesh.vertices[:, 2]) == pytest.approx(25)
    check_shell(mesh)


@pytest.mark.parametrize(
    ("x0", "y0", "angle"),
    [(0, 0, 0), (0.135400041, 0.44463753, 0), (0.135399959, 0.44463753, 26)],
)
def test_lod2_handles_occupied_city_boundary(x0, y0, angle):
    outline = translate(rotate(box(-5, 20, 18, 40), angle), xoff=x0, yoff=y0)
    city = make_city(ramp, make_building(outline, 25))
    city.terrain.raster.set_bounds(Bounds(x0, y0, x0 + 61, y0 + 61))
    mesh = mesh_city(city, dict(lod="lod2"))
    seam = seam_vertices(mesh)
    assert np.any(np.isclose(seam[:, 0], x0, atol=1e-9, rtol=0))
    expected = ramp(np.maximum(seam[:, 0] - x0, 0.5), seam[:, 1] - y0)
    np.testing.assert_allclose(seam[:, 2], expected, atol=1e-8, rtol=0)
    check_shell(mesh)


@pytest.mark.parametrize(
    ("x0", "y0", "angle"), [(0.135399959, 0.44463753, 26), (0.1354, 0.444637504, 3)]
)
def test_lod2_handles_clipped_block_with_an_oblique_roof_step(x0, y0, angle):
    parts = [
        make_building(
            translate(rotate(outline, angle, origin=(6.5, 30)), xoff=x0, yoff=y0),
            level,
            name,
        )
        for outline, level, name in (
            (box(-5, 20, 18, 30), 30, "south"),
            (box(-5, 30, 18, 40), 34, "north"),
        )
    ]
    city = make_city(ramp, *parts)
    city.terrain.raster.set_bounds(Bounds(x0, y0, x0 + 61, y0 + 61))
    audit = {}
    mesh = mesh_city(city, dict(lod="lod2"), stage_audit=audit)
    patches = audit["roof_catalogue"]["patch_markers"]
    assert {patch["source_key"] for patch in patches} == {"south", "north"}
    seam = seam_vertices(mesh)
    assert np.any(np.isclose(seam[:, 0], x0, atol=1e-9, rtol=0))
    expected = ramp(np.maximum(seam[:, 0] - x0, 0.5), seam[:, 1] - y0)
    np.testing.assert_allclose(seam[:, 2], expected, atol=1e-8, rtol=0)
    check_shell(mesh)


def test_roof_breakpoints_interpolate_existing_terrain_edges():
    terrain = Mesh(
        vertices=np.array([[2, 0, 2], [2, 2, 6], [2, 4, 4], [0, 2, 0]], dtype=float),
        faces=np.array([[0, 1, 3], [1, 2, 3]], dtype=int),
    )
    raster = Raster(data=np.zeros((5, 5)))
    raster.set_bounds(Bounds(0, 0, 5, 5))
    contact = TerrainContact(terrain, raster, np.zeros(2), (0, 0, 5, 5))
    a, b = np.array([2, 1]), np.array([2, 3])
    expected = [[2, 1, 4], [2, 2, 6], [2, 3, 5]]
    np.testing.assert_array_equal(contact.profile(a, b, 0), expected)
    np.testing.assert_array_equal(contact.profile(b, a, 0), expected[::-1])


def test_clipped_wall_profile_includes_dem_grid_breakpoints():
    terrain = Mesh(
        vertices=np.array([[2, 0, 2], [2, 4, 4], [5, 2, 0]], dtype=float),
        faces=np.array([[0, 1, 2]], dtype=int),
    )
    raster = Raster(data=np.tile([0, 0, 6, 0, 0], (5, 1)).T.astype(float))
    raster.set_bounds(Bounds(0, 0, 5, 5))
    contact = TerrainContact(terrain, raster, np.zeros(2), (0, 0, 5, 5))
    profile = contact.profile(np.array([0, 1]), np.array([0, 4]), 0)
    np.testing.assert_array_equal(profile[:, 1], [1, 1.5, 2.5, 3.5, 4])
    np.testing.assert_array_equal(profile[:, 2], [0, 0, 6, 0, 0])


def test_short_wall_breakpoint_matches_only_the_overlapping_terrain_segment():
    ground = Mesh(
        vertices=np.array([[2, 0, 2], [2, 4, 4], [5, 2, 0]], dtype=float),
        faces=np.array([[0, 1, 2]], dtype=int),
    )
    raster = Raster(data=np.zeros((5, 5)))
    raster.set_bounds(Bounds(0, 0, 5, 5))
    contact = TerrainContact(ground, raster, np.zeros(2), (0, 0, 5, 5))
    a, b = np.array([2, 1]), np.array([2 + 5e-10, 1.003])
    np.testing.assert_allclose(
        contact.profile(a, b, 0), [[*a, 2.5], [*b, 2.5015]], atol=1e-12, rtol=0
    )


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not installed")
def test_volume_boundary_keeps_sloped_contact(lod):
    city = make_city(ramp, make_building(box(15, 20, 45, 40), 25))
    mesh = builder.build_city_volume_mesh(
        city, max_mesh_size=5, domain_height=35, report_mesh_quality=False, **lod
    )
    faces, markers = mesh.boundary_faces, mesh.boundary_markers
    terrain = np.unique(faces[markers == -1])
    buildings = np.unique(faces[markers >= 0])
    seam = mesh.vertices[np.intersect1d(terrain, buildings)]
    assert len(seam) >= 4
    np.testing.assert_allclose(
        seam[:, 2], ramp(seam[:, 0], seam[:, 1]), atol=1e-8, rtol=0
    )
    tetrahedra = mesh.vertices[mesh.cells]
    assert np.all(np.linalg.det(tetrahedra[:, 1:] - tetrahedra[:, :1]) > 0)
