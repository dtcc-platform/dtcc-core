from collections import defaultdict

import pytest
import numpy as np
from dtcc_core.builder.meshing import (
    backends as backends_module,
    merge_meshes,
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
    assert mesh.vertices.shape[0] == 4


def test_merge_meshes_compacts_unused_vertices():
    mesh_a = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [5.0, 5.0, 5.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([1], dtype=int),
    )
    mesh_b = Mesh(
        vertices=np.array(
            [
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
                [9.0, 9.0, 9.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([2], dtype=int),
    )

    merged = merge_meshes([mesh_a, mesh_b], weld=True)

    used = np.unique(np.asarray(merged.faces, dtype=np.int64).reshape(-1))
    assert len(used) == len(merged.vertices)


def test_merge_meshes_without_weld_preserves_vertex_array():
    mesh_a = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [5.0, 5.0, 5.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([1], dtype=int),
    )
    mesh_b = Mesh(
        vertices=np.array(
            [
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
                [9.0, 9.0, 9.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([2], dtype=int),
    )

    merged = merge_meshes([mesh_a, mesh_b], weld=False)

    assert len(merged.vertices) == len(mesh_a.vertices) + len(mesh_b.vertices)


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


def test_inspect_tetgen_plc_flags_unreferenced_vertices():
    diagnostics = tetgen_utils.inspect_tetgen_plc(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
                [50.0, 50.0, 50.0],
            ]
        ),
        np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        {
            "top": np.array([0, 1, 2, 3], dtype=int),
        },
    )

    assert any("unreferenced vertices" in message for message in diagnostics.errors)


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


def test_compute_boundary_triangle_facets_subdivides_tall_sidewalls_into_stacked_triangles(monkeypatch):
    shell = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [0.0, 4.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, -2], dtype=int),
    )
    closure = Mesh(
        vertices=np.array(shell.vertices, copy=True),
        faces=np.array(shell.faces, copy=True),
        markers=np.array(shell.markers, copy=True),
    )

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [4.0, 0.0, 0.0],
                    [4.0, 4.0, 0.0],
                    [0.0, 4.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        ),
    )

    vertices, boundary_facets = tetgen_utils.compute_boundary_triangle_facets(
        shell,
        closure,
        top_height=100.0,
        top_cap_backend="triangle",
    )

    diagnostics = tetgen_utils.inspect_tetgen_plc(vertices, shell.faces, boundary_facets)
    assert diagnostics.errors == []

    vertices = np.asarray(vertices, dtype=float)
    top_z = float(np.max(vertices[:, 2]))
    sidewall_facets = [
        facet
        for facet in boundary_facets
        if not np.allclose(vertices[np.asarray(facet, dtype=int), 2], top_z)
    ]

    assert len(sidewall_facets) == 16
    assert all(len(facet) == 3 for facet in sidewall_facets)

    vertical_edge_lengths = []
    for facet in sidewall_facets:
        coords = vertices[np.asarray(facet, dtype=int)]
        for pa, pb in zip(coords, np.roll(coords, -1, axis=0)):
            if np.allclose(pa[:2], pb[:2]) and not np.isclose(pa[2], pb[2]):
                vertical_edge_lengths.append(abs(float(pa[2] - pb[2])))

    assert vertical_edge_lengths
    assert set(np.round(vertical_edge_lengths, 8)) == {50.0}


def test_compute_oriented_boundary_triangle_facets_matches_shell_shared_edge_winding(monkeypatch):
    shell = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [0.0, 4.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, -2], dtype=int),
    )
    closure = Mesh(
        vertices=np.array(shell.vertices, copy=True),
        faces=np.array(shell.faces, copy=True),
        markers=np.array(shell.markers, copy=True),
    )

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [4.0, 0.0, 0.0],
                    [4.0, 4.0, 0.0],
                    [0.0, 4.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        ),
    )

    _, raw_boundary_facets = tetgen_utils.compute_boundary_triangle_facets(
        shell,
        closure,
        top_height=20.0,
        top_cap_backend="triangle",
    )
    raw_stats = tetgen_utils._triangle_edge_orientation_stats(
        shell.faces,
        raw_boundary_facets,
    )
    assert raw_stats["shell_boundary_same_direction_edge_count"] > 0

    vertices, oriented_shell_faces, oriented_boundary_facets = tetgen_utils.compute_oriented_boundary_triangle_plc(
        shell,
        closure,
        top_height=20.0,
        top_cap_backend="triangle",
    )
    oriented_stats = tetgen_utils._triangle_edge_orientation_stats(
        oriented_shell_faces,
        oriented_boundary_facets,
    )

    assert oriented_stats["same_direction_shared_edge_count"] == 0
    assert oriented_stats["shell_boundary_same_direction_edge_count"] == 0

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        oriented_shell_faces,
        oriented_boundary_facets,
    )
    assert diagnostics.errors == []


def test_compute_oriented_boundary_plc_remeshes_top_cap_from_shell_boundary(
    monkeypatch,
):
    shell = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [0.0, 4.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, -2], dtype=int),
    )
    closure = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [0.0, 4.0, 0.0],
                [2.0, 2.0, 0.0],
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
        markers=np.array([-2, -2, -2, -2], dtype=int),
    )
    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(shell.vertices, copy=True),
            faces=np.array(shell.faces, copy=True),
            markers=np.array([0, 0], dtype=int),
        ),
    )

    (
        vertices,
        oriented_shell_faces,
        boundary_facets,
        boundary_facet_markers,
        audit_boundary_triangles,
    ) = (
        tetgen_utils.compute_oriented_boundary_plc(
            shell,
            closure,
            top_height=20.0,
            top_cap_backend="triangle",
        )
    )

    assert boundary_facets
    assert all(len(facet) == 3 for facet in boundary_facets)
    assert boundary_facet_markers.count(-2) == 2
    assert set(boundary_facet_markers) == {-2, -3, -4, -5, -6}
    assert audit_boundary_triangles.ndim == 2
    assert audit_boundary_triangles.shape[1] == 3

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        oriented_shell_faces,
        boundary_facets,
    )
    assert diagnostics.errors == []
    assert all(
        info["vertex_count"] == 3 for info in diagnostics.boundary_facets.values()
    )


def test_compute_oriented_boundary_plc_uses_remeshed_top_cap_when_max_mesh_size_is_set(
    monkeypatch,
):
    shell = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [0.0, 4.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, -2], dtype=int),
    )
    closure = Mesh(
        vertices=np.array(shell.vertices, copy=True),
        faces=np.array(shell.faces, copy=True),
        markers=np.array(shell.markers, copy=True),
    )

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [4.0, 0.0, 0.0],
                    [4.0, 4.0, 0.0],
                    [0.0, 4.0, 0.0],
                    [2.0, 2.0, 0.0],
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
        ),
    )

    (
        vertices,
        oriented_shell_faces,
        boundary_facets,
        boundary_facet_markers,
        audit_boundary_triangles,
    ) = tetgen_utils.compute_oriented_boundary_plc(
        shell,
        closure,
        top_height=20.0,
        top_cap_backend="triangle",
        top_cap_max_mesh_size=2.0,
    )

    assert boundary_facets
    assert all(len(facet) == 3 for facet in boundary_facets)
    assert boundary_facet_markers.count(-2) == 4
    assert set(boundary_facet_markers) == {-2, -3, -4, -5, -6}
    assert audit_boundary_triangles.shape == (len(boundary_facets), 3)

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        oriented_shell_faces,
        boundary_facets,
    )
    assert diagnostics.errors == []
    assert all(
        info["vertex_count"] == 3 for info in diagnostics.boundary_facets.values()
    )


def test_compute_oriented_boundary_plc_triangles_sidewalls_for_unstable_shell(
    monkeypatch,
):
    shell = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.1, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 4], [2, 3, 4]], dtype=int),
        markers=np.array([-2, -2, -2], dtype=int),
    )
    closure = Mesh(
        vertices=np.array(shell.vertices, copy=True),
        faces=np.array(shell.faces, copy=True),
        markers=np.array(shell.markers, copy=True),
    )
    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(shell.vertices, copy=True),
            faces=np.array(shell.faces, copy=True),
            markers=np.array([0, 0, 0], dtype=int),
        ),
    )

    (
        vertices,
        oriented_shell_faces,
        boundary_facets,
        boundary_facet_markers,
        audit_boundary_triangles,
    ) = (
        tetgen_utils.compute_oriented_boundary_plc(
            shell,
            closure,
            top_height=20.0,
            top_cap_backend="triangle",
        )
    )

    assert len(boundary_facets) > 5
    assert all(len(facet) == 3 for facet in boundary_facets)
    assert boundary_facet_markers.count(-2) == len(shell.faces)
    assert set(boundary_facet_markers) == {-2, -3, -4, -5, -6}
    assert audit_boundary_triangles.ndim == 2
    assert audit_boundary_triangles.shape[1] == 3

    diagnostics = tetgen_utils.inspect_tetgen_plc(
        vertices,
        oriented_shell_faces,
        boundary_facets,
    )
    assert diagnostics.errors == []
    assert all(
        info["vertex_count"] == 3 for info in diagnostics.boundary_facets.values()
    )


def test_compute_boundary_triangle_facets_keeps_corner_sidewall_strips_consistent(monkeypatch):
    shell = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [2.0, 4.0, 0.0],
                [0.0, 4.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 3], [1, 2, 3], [0, 3, 4]], dtype=int),
        markers=np.array([-2, -2, -2], dtype=int),
    )
    closure = Mesh(
        vertices=np.array(shell.vertices, copy=True),
        faces=np.array(shell.faces, copy=True),
        markers=np.array(shell.markers, copy=True),
    )

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [4.0, 0.0, 0.0],
                    [4.0, 4.0, 0.0],
                    [2.0, 4.0, 0.0],
                    [0.0, 4.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 3], [1, 2, 3], [0, 3, 4]], dtype=int),
            markers=np.array([0, 0, 0], dtype=int),
        ),
    )

    vertices, boundary_facets = tetgen_utils.compute_boundary_triangle_facets(
        shell,
        closure,
        top_height=100.0,
        top_cap_backend="triangle",
    )

    diagnostics = tetgen_utils.inspect_tetgen_plc(vertices, shell.faces, boundary_facets)
    assert diagnostics.errors == []

    vertices = np.asarray(vertices, dtype=float)
    top_z = float(np.max(vertices[:, 2]))
    sidewall_facets = [
        facet
        for facet in boundary_facets
        if not np.allclose(vertices[np.asarray(facet, dtype=int), 2], top_z)
    ]
    vertical_edge_lengths = []
    for facet in sidewall_facets:
        coords = vertices[np.asarray(facet, dtype=int)]
        for pa, pb in zip(coords, np.roll(coords, -1, axis=0)):
            if np.allclose(pa[:2], pb[:2]) and not np.isclose(pa[2], pb[2]):
                vertical_edge_lengths.append(abs(float(pa[2] - pb[2])))

    assert vertical_edge_lengths
    assert set(np.round(vertical_edge_lengths, 8)) == {25.0}


def test_compute_boundary_triangle_facets_avoids_horizontal_sidewall_strip_edges_on_sloped_boundary(monkeypatch):
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
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
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
        top_height=100.0,
        top_cap_backend="triangle",
    )

    vertices = np.asarray(vertices, dtype=float)
    top_z = float(np.max(vertices[:, 2]))
    sidewall_facets = [
        facet
        for facet in boundary_facets
        if not np.allclose(vertices[np.asarray(facet, dtype=int), 2], top_z)
    ]

    assert sidewall_facets
    assert all(len(facet) == 3 for facet in sidewall_facets)

    for facet in sidewall_facets:
        coords = vertices[np.asarray(facet, dtype=int)]
        for pa, pb in zip(coords, np.roll(coords, -1, axis=0)):
            assert not (
                not np.allclose(pa[:2], pb[:2])
                and np.isclose(pa[2], pb[2])
                and not np.isclose(pa[2], float(np.min(vertices[:, 2])))
                and not np.isclose(pa[2], top_z)
            )


def test_compute_boundary_triangle_facets_uses_authoritative_closure_topology(monkeypatch):
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
        top_cap_backend="triangle",
    )

    assert vertices.shape == (8, 3)
    assert len(boundary_facets) == 10
    assert all(len(facet) == 3 for facet in boundary_facets)

    top_z = float(np.max(vertices[:, 2]))
    top_faces = np.asarray(
        [
            facet
            for facet in boundary_facets
            if np.allclose(vertices[np.asarray(facet, dtype=int), 2], top_z)
        ],
        dtype=int,
    )
    assert top_faces.shape == (2, 3)
    assert np.all(top_faces >= 4)
    for tri in top_faces:
        points = vertices[tri]
        assert np.cross(points[1] - points[0], points[2] - points[0])[2] > 0.0


def test_compute_boundary_triangle_facets_rejects_top_cap_backend_that_changes_boundary_loop(monkeypatch):
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 1.0],
                [10.0, 0.0, 1.0],
                [10.0, 10.0, 1.0],
                [0.0, 10.0, 1.0],
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
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )

    def fake_top_cap(**kwargs):
        backend = kwargs["backend"]
        if backend == "triangle":
            return Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [5.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 4], [2, 3, 4]], dtype=int),
                markers=np.array([0, 0, 0], dtype=int),
            )

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        fake_top_cap,
    )
    monkeypatch.setattr(
        "dtcc_core.builder.meshing.backends.available_2d_meshers",
        lambda: ["dtcc_mesher"],
    )

    with pytest.raises(
        ValueError,
        match="Top-cap remeshing changed the prescribed outer boundary vertex set",
    ):
        tetgen_utils.compute_boundary_triangle_facets(
            surface_mesh,
            closure_mesh,
            top_height=20.0,
            top_cap_backend="triangle",
        )


def test_compute_boundary_triangle_facets_ignores_closure_boundary_loop_mismatch(monkeypatch):
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 1.0],
                [10.0, 0.0, 1.0],
                [10.0, 10.0, 1.0],
                [0.0, 10.0, 1.0],
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
                [5.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 3], [0, 3, 4], [1, 2, 3]], dtype=int),
        markers=np.array([0, 0, 0], dtype=int),
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
        top_cap_backend="triangle",
    )

    assert vertices.shape == (8, 3)
    diagnostics = tetgen_utils.inspect_tetgen_plc(vertices, surface_mesh.faces, boundary_facets)
    assert diagnostics.errors == []


def test_compute_boundary_triangle_facets_rejects_explicit_top_cap_backend_that_changes_boundary_loop(monkeypatch):
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 1.0],
                [10.0, 0.0, 1.0],
                [10.0, 10.0, 1.0],
                [0.0, 10.0, 1.0],
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
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )

    def fake_top_cap(**kwargs):
        backend = kwargs["backend"]
        if backend == "triangle":
            return Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [5.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 4], [2, 3, 4]], dtype=int),
                markers=np.array([0, 0, 0], dtype=int),
            )
        raise AssertionError("No backend fallback should run in the strict PLC path.")

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        fake_top_cap,
    )
    monkeypatch.setattr(
        "dtcc_core.builder.meshing.backends.available_2d_meshers",
        lambda: ["dtcc_mesher", "spade"],
    )

    with pytest.raises(
        ValueError,
        match="Top-cap remeshing changed the prescribed outer boundary vertex set",
    ):
        tetgen_utils.compute_boundary_triangle_facets(
            surface_mesh,
            closure_mesh,
            top_height=20.0,
            top_cap_backend="triangle",
        )

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


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_volume_mesh_uses_triangle_only_closure_top_cap(monkeypatch):
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
    closure_mesh = Mesh(
        vertices=np.array(mesh.vertices, copy=True),
        faces=np.array(mesh.faces, copy=True),
        markers=np.array(mesh.markers, copy=True),
    )
    captured: dict[str, object] = {}

    class TetgenResult:
        points = np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 0.0, 10.0],
            ],
            dtype=float,
        )
        tets = np.array([[0, 1, 2, 3]], dtype=np.int32)
        boundary_tri_faces = None
        boundary_tri_markers = None

    def fake_tetrahedralize(**kwargs):
        captured["boundary_facets"] = kwargs["boundary_facets"]
        captured["boundary_facet_markers"] = kwargs["boundary_facet_markers"]
        return TetgenResult()

    monkeypatch.setattr(tetgen_module.tetwrap, "tetrahedralize", fake_tetrahedralize)

    tetgen_module.build_volume_mesh(
        mesh,
        closure_mesh=closure_mesh,
        top_height=20.0,
    )

    boundary_facets = captured["boundary_facets"]
    assert isinstance(boundary_facets, list)
    assert all(len(facet) == 3 for facet in boundary_facets)
    assert captured["boundary_facet_markers"].count(-2) == len(closure_mesh.faces)
    assert set(captured["boundary_facet_markers"]) == {-2, -3, -4, -5, -6}


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_volume_mesh_uses_closure_top_cap_when_max_mesh_size_is_set(
    monkeypatch,
):
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
    closure_mesh = Mesh(
        vertices=np.array(mesh.vertices, copy=True),
        faces=np.array(mesh.faces, copy=True),
        markers=np.array(mesh.markers, copy=True),
    )
    captured: dict[str, object] = {}

    class TetgenResult:
        points = np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 0.0, 10.0],
            ],
            dtype=float,
        )
        tets = np.array([[0, 1, 2, 3]], dtype=np.int32)
        boundary_tri_faces = None
        boundary_tri_markers = None

    monkeypatch.setattr(
        "dtcc_core.builder.meshing.flat_mesh_backends.build_city_flat_mesh_from_coverage",
        lambda **kwargs: Mesh(
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
        ),
    )

    def fake_tetrahedralize(**kwargs):
        captured["boundary_facets"] = kwargs["boundary_facets"]
        captured["boundary_facet_markers"] = kwargs["boundary_facet_markers"]
        return TetgenResult()

    monkeypatch.setattr(tetgen_module.tetwrap, "tetrahedralize", fake_tetrahedralize)

    tetgen_module.build_volume_mesh(
        mesh,
        closure_mesh=closure_mesh,
        top_height=20.0,
        top_cap_backend="triangle",
        top_cap_max_mesh_size=2.0,
    )

    boundary_facets = captured["boundary_facets"]
    assert isinstance(boundary_facets, list)
    assert all(len(facet) == 3 for facet in boundary_facets)
    assert captured["boundary_facet_markers"].count(-2) == 4
    assert set(captured["boundary_facet_markers"]) == {-2, -3, -4, -5, -6}


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_volume_mesh_forwards_named_boundary_markers_without_closure_mesh(monkeypatch):
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
    captured = {}

    class TetgenResult:
        points = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=float,
        )
        tets = np.array([[0, 1, 2, 3]], dtype=np.int32)
        boundary_tri_faces = None
        boundary_tri_markers = None

    def fake_tetrahedralize(**kwargs):
        captured["boundary_facets"] = kwargs["boundary_facets"]
        captured["boundary_facet_markers"] = kwargs["boundary_facet_markers"]
        return TetgenResult()

    monkeypatch.setattr(tetgen_module.tetwrap, "tetrahedralize", fake_tetrahedralize)

    tetgen_module.build_volume_mesh(mesh, top_height=20.0, closure_mesh=None)

    assert list(captured["boundary_facets"].keys()) == ["south", "east", "north", "west", "top"]
    assert captured["boundary_facet_markers"] == {
        "south": -5,
        "east": -4,
        "north": -6,
        "west": -3,
        "top": -2,
    }
