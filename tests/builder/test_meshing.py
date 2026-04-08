import pytest
import numpy as np
from dtcc_core.builder.meshing import (
    backends as backends_module,
    mesh_multisurface,
    mesh_surface,
    mesh_multisurfaces,
    disjoint_meshes,
)
from dtcc_core.builder.meshing import tetgen as tetgen_module
from dtcc_core.builder.meshing import tetgen_utils
from dtcc_core.builder.meshing.tetgen import is_tetgen_available

from dtcc_core.model import Mesh, Surface, MultiSurface


@pytest.fixture
def simple_surface():
    return Surface(
        vertices=np.array(
            [
                [0, 0, 5],
                [0, 10, 5],
                [10, 10, 8],
                [10, 0, 8],
            ]
        )
    )


@pytest.fixture
def complex_surface():
    return Surface(
        vertices=np.array(
            [
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [2, 10, 0],
                [2, 12, 0],
                [0, 12, 0],
            ]
        )
    )


@pytest.fixture
def multi_surface(simple_surface, complex_surface):
    ms = MultiSurface()
    ms.surfaces = [simple_surface, complex_surface]
    return ms


@pytest.fixture
def second_surface_pair():
    surface3 = Surface(
        vertices=np.array(
            [
                [0, 0, 4],
                [0, 10, 4],
                [10, 10, 7],
                [10, 0, 7],
            ]
        )
    )
    surface4 = Surface(
        vertices=np.array(
            [
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [2, 10, 0],
                [2, 12, 0],
                [0, 12, 0],
            ]
        )
    )
    ms = MultiSurface()
    ms.surfaces = [surface3, surface4]
    return ms


@pytest.fixture
def disjoint_cubes_mesh():
    cube1_vertices = np.array(
        [
            [0, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
            [1, 0, 1],
        ]
    )
    cube1_faces = [
        [0, 1, 2],
        [0, 2, 3],
        [0, 4, 5],
        [0, 5, 1],
        [1, 5, 6],
        [1, 6, 2],
        [2, 6, 7],
        [2, 7, 3],
        [3, 7, 4],
        [3, 4, 0],
        [4, 7, 6],
        [4, 6, 5],
    ]
    cube2_vertices = cube1_vertices + np.array([2, 2, 2])
    cube2_faces = [[i + 8 for i in face] for face in cube1_faces]
    cubes_vertices = np.vstack([cube1_vertices, cube2_vertices])
    cubes_faces = np.array(cube1_faces + cube2_faces)
    return Mesh(vertices=cubes_vertices, faces=cubes_faces)


def test_mesh_simple_surface(simple_surface):
    mesh = mesh_surface(simple_surface, mesher="spade")
    assert len(mesh.vertices) == 4
    assert len(mesh.faces) == 2
    assert pytest.approx(mesh.vertices[:, 2].min()) == 5
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_triangle_surface(simple_surface):
    mesh = mesh_surface(simple_surface, triangle_size=5, mesher="spade")
    assert len(mesh.vertices) >= 11
    assert len(mesh.faces) >= 9
    assert pytest.approx(mesh.vertices[:, 2].min()) == 5
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_multisurface(multi_surface):
    mesh = mesh_multisurface(multi_surface, mesher="spade")
    assert len(mesh.vertices) == 10
    assert len(mesh.faces) == 6
    assert pytest.approx(mesh.vertices[:, 2].min()) == 0
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_multisurfaces(multi_surface, second_surface_pair):
    ms = multi_surface.translate(0, 0, 1)
    meshes = mesh_multisurfaces([multi_surface, second_surface_pair], mesher="spade")

    assert len(meshes) == 2

    # First mesh checks
    assert len(meshes[0].vertices) == 10
    assert len(meshes[0].faces) == 6
    assert pytest.approx(meshes[0].vertices[:, 2].min()) == 1
    assert pytest.approx(meshes[0].vertices[:, 2].max()) == 9

    # Second mesh checks
    assert len(meshes[1].vertices) == 10
    assert len(meshes[1].faces) == 6
    assert pytest.approx(meshes[1].vertices[:, 2].min()) == 0
    assert pytest.approx(meshes[1].vertices[:, 2].max()) == 7


def test_disjoint_mesh(disjoint_cubes_mesh):
    disjointed_meshes = disjoint_meshes(disjoint_cubes_mesh)

    assert len(disjointed_meshes) == 2
    assert disjointed_meshes[0].faces.max() == 7
    assert disjointed_meshes[1].faces.max() == 7
    assert len(disjointed_meshes[0].vertices) == 8
    assert len(disjointed_meshes[1].vertices) == 8


def test_mesh_surface_dtcc_mesher_backend(simple_surface):
    pytest.importorskip("dtcc_mesher")

    mesh = mesh_surface(simple_surface, mesher="dtcc_mesher")

    assert mesh.vertices.shape[0] >= 4
    assert mesh.faces.shape[0] >= 2
    assert pytest.approx(mesh.vertices[:, 2].min()) == 5
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_multisurface_dtcc_mesher_backend(multi_surface):
    pytest.importorskip("dtcc_mesher")

    mesh = mesh_multisurface(multi_surface, mesher="dtcc_mesher")

    assert mesh.vertices.shape[0] >= 8
    assert mesh.faces.shape[0] >= 4
    assert pytest.approx(mesh.vertices[:, 2].min()) == 0
    assert pytest.approx(mesh.vertices[:, 2].max()) == 8


def test_mesh_surface_dtcc_mesher_triangle_size_controls_density(simple_surface):
    pytest.importorskip("dtcc_mesher")

    coarse = mesh_surface(simple_surface, mesher="dtcc_mesher")
    refined = mesh_surface(simple_surface, triangle_size=5.0, mesher="dtcc_mesher")

    assert refined.faces.shape[0] > coarse.faces.shape[0]
    assert refined.vertices.shape[0] > coarse.vertices.shape[0]


def test_available_2d_meshers_reports_triangle_and_spade(monkeypatch):
    monkeypatch.setattr(
        backends_module.importlib.util,
        "find_spec",
        lambda name: object() if name == "dtcc_mesher" else None,
    )
    monkeypatch.setattr(
        backends_module,
        "_builder_backend_available",
        lambda name: name in {"triangle", "spade"},
    )

    assert backends_module.available_2d_meshers() == ["dtcc_mesher", "triangle", "spade"]
    assert backends_module.resolve_2d_mesher("triangle") == "triangle"
    assert backends_module.resolve_2d_mesher("spade") == "spade"
    assert backends_module.resolve_2d_mesher("auto") == "dtcc_mesher"


def test_snap_mesh_vertices():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 0, 0], [1, 1.01, 0]])
    mesh = Mesh(vertices=vertices, faces=np.array([[0, 1, 2], [1, 3, 4]]))
    mesh = mesh.snap_vertices(0.1)
    assert mesh.faces[1][2] == 2  # The last vertex should be snapped to the first one


def test_inspect_tetgen_plc_accepts_simple_shell():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )
    vertices, boundary_facets = tetgen_utils.compute_boundary_facets(mesh, top_height=20.0)

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        mesh.faces,
        boundary_facets,
    )

    assert diagnostics.errors == []
    assert diagnostics.warnings == []


def test_inspect_tetgen_plc_flags_duplicate_faces():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 1, 2]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )
    vertices, boundary_facets = tetgen_utils.compute_boundary_facets(mesh, top_height=20.0)

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        mesh.faces,
        boundary_facets,
    )

    assert any("duplicate triangles" in message for message in diagnostics.errors)


def test_inspect_tetgen_plc_warns_on_tiny_shell_features():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
                [1.0e-5, 1.0e-5, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 4],
                [1, 2, 4],
                [2, 3, 4],
                [3, 0, 4],
            ],
            dtype=int,
        ),
        markers=np.array([0, 0, 0, 0], dtype=int),
    )
    vertices, boundary_facets = tetgen_utils.compute_boundary_facets(mesh, top_height=20.0)

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        mesh.faces,
        boundary_facets,
    )

    assert diagnostics.errors == []
    assert any("minimum triangle quality" in message for message in diagnostics.warnings)
    assert any("smaller than the median edge length" in message for message in diagnostics.warnings)


def test_compute_boundary_facets_preserves_intermediate_wall_vertices():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 5.0, 0.0],
                [10.0, 10.0, 0.0],
                [5.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 5.0, 0.0],
                [5.0, 5.0, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 8],
                [1, 2, 8],
                [2, 3, 8],
                [3, 4, 8],
                [4, 5, 8],
                [5, 6, 8],
                [6, 7, 8],
                [7, 0, 8],
            ],
            dtype=int,
        ),
        markers=np.zeros(8, dtype=int),
    )

    _, boundary_facets = tetgen_utils.compute_boundary_facets(mesh, top_height=20.0)

    assert list(boundary_facets["south"][:3]) == [0, 1, 2]
    assert list(boundary_facets["east"][:3]) == [2, 3, 4]
    assert list(boundary_facets["north"][:3]) == [4, 5, 6]
    assert list(boundary_facets["west"][:3]) == [6, 7, 0]


def test_compute_boundary_triangle_facets_retriangulates_top_from_outer_boundary(monkeypatch):
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 1.0],
                [10.0, 0.0, 2.0],
                [10.0, 10.0, 3.0],
                [0.0, 10.0, 2.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )
    closure_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
                [5.0, 5.0, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 4],
                [1, 2, 4],
                [2, 3, 4],
                [3, 0, 4],
            ],
            dtype=int,
        ),
        markers=np.array([0, 0, 0, 0], dtype=int),
    )

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        ),
    )

    vertices, boundary_facets = tetgen_utils.compute_boundary_triangle_facets(
        surface_mesh,
        closure_mesh,
        top_height=20.0,
        top_cap_backend="auto",
    )

    assert vertices.shape == (8, 3)
    assert len(boundary_facets) == 10
    assert all(len(facet) == 3 for facet in boundary_facets)

    top_faces = np.asarray(boundary_facets[-2:], dtype=int)
    assert np.all(top_faces >= 4)
    for tri in top_faces:
        points = vertices[tri]
        assert np.cross(points[1] - points[0], points[2] - points[0])[2] > 0.0

@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_volume_mesh_rejects_invalid_shell_before_tetgen(monkeypatch):
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 1, 2]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )
    called = {"tetgen": False}

    def fail_if_called(*args, **kwargs):
        called["tetgen"] = True
        raise AssertionError("TetGen wrapper should not run for invalid PLC input.")

    monkeypatch.setattr(tetgen_module.tetwrap, "tetrahedralize", fail_if_called)

    with pytest.raises(ValueError, match="TetGen PLC precheck failed"):
        tetgen_module.build_volume_mesh(mesh, top_height=20.0)

    assert called["tetgen"] is False


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_volume_mesh_rejects_empty_tetgen_output(monkeypatch):
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )

    class EmptyTetResult:
        points = np.zeros((0, 3), dtype=float)
        tets = np.zeros((0, 4), dtype=np.int32)
        boundary_tri_faces = None
        boundary_tri_markers = None

    monkeypatch.setattr(
        tetgen_module.tetwrap,
        "tetrahedralize",
        lambda **kwargs: EmptyTetResult(),
    )

    with pytest.raises(RuntimeError, match="TetGen returned no linear tetrahedral cells"):
        tetgen_module.build_volume_mesh(mesh, top_height=20.0)
