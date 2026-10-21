"""Shared constraints, geometric surfaces and bounded pilot transactions."""

import json
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon

from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.meshing.roof_validation import validate_shell
from dtcc_core.builder.meshing.shared_edges import SharedEdges
from dtcc_core.builder.meshing.surface_faces import extract_surface_faces
from dtcc_core.builder.meshing.surface_refinement import (
    SurfaceRefinementOptions,
    refine_shared_surface,
)
from dtcc_core.builder.meshing.surface_refinement_backend import mesh_face
from dtcc_core.model import Mesh


def require_pilot_mesher():
    dm = pytest.importorskip("dtcc_mesher")
    if "reuse_acute_corners" not in dm.MeshingOptions.__dataclass_fields__:
        pytest.skip("Shared-edge pilot requires the documented DTCC Mesher extension")


def fixture_shell():
    data = json.loads(
        (Path(__file__).parent / "fixtures/wall_short_eave.json").read_text()
    )
    return Mesh(
        vertices=np.array(data["vertices"]),
        faces=np.array(data["faces"]),
        markers=np.array(data["markers"]),
    ), data["owners"]


def test_unequal_subdivisions_and_opposite_orientation_are_noded():
    v = np.array([[0.0, 0, 0], [4, 0, 0], [1, 0, 0], [4, 0, 3], [0, 0, 3]])
    registry = SharedEdges(v, {0: [(0, 1)], 1: [(1, 2), (2, 0)], 2: [(3, 4)]})
    assert registry.boundary(0) == [(0, 2), (2, 1)]
    assert registry.boundary(1) == [(1, 2), (2, 0)]
    ei = registry.face_edges[0][0][0]
    affected, added = registry.reconcile(
        [(ei, 0.5), (ei, 0.5 + 1e-10)], max_vertices=20
    )
    assert added == 1 and affected == {0, 1}
    assert registry.boundary(0) == [e[::-1] for e in registry.boundary(1)[::-1]]
    assert registry.boundary(2) == [(3, 4)]
    assert np.array_equal(np.array(registry.vertices)[: len(v)], v)


def test_distinct_coincident_xy_edges_and_atomic_budget():
    v = np.array([[0.0, 0, 0], [1, 0, 0], [0, 0, 1e-4], [1, 0, 1e-4]])
    r = SharedEdges(v, {0: [(0, 1)], 1: [(2, 3)]})
    with pytest.raises(ValueError, match="budget"):
        r.reconcile([(0, 0.5), (1, 0.5)], max_vertices=5)
    assert len(r.vertices) == 4 and all(e.parameters == [0.0, 1.0] for e in r.edges)
    r.reconcile([(0, 0.5)], max_vertices=5)
    assert r.boundary(1) == [(2, 3)]


def test_single_edge_trial_copies_only_mutable_split_state():
    vertices = np.array([[0.0, 0, 0], [4, 0, 0], [0, 0, 2], [4, 0, 2]])
    original = SharedEdges(vertices, {0: [(0, 1)], 1: [(1, 0)], 2: [(2, 3)]})
    first = original.copy_for_split(0)
    second = original.copy_for_split(0)
    first.reconcile([(0, 0.25)], max_vertices=5)
    second.reconcile([(0, 0.75)], max_vertices=5)
    assert original.boundary(0) == [(0, 1)]
    assert original.edges[0].parameters == [0.0, 1.0]
    assert first.edges[0].parameters == [0.0, 0.25, 1.0]
    assert second.edges[0].parameters == [0.0, 0.75, 1.0]
    assert first.boundary(0) == [e[::-1] for e in first.boundary(1)[::-1]]
    assert len(original.vertices) == 4
    assert first.edges[1] is original.edges[1]
    assert first.vertices[0] is original.vertices[0]
    assert np.array_equal(original.vertices, vertices)


@pytest.mark.parametrize("backend", ["dtcc_mesher", "triangle"])
def test_maintained_quality_backend_holes_and_short_eave(backend):
    if backend == "triangle":
        from dtcc_core.builder.meshing.backends import available_2d_meshers

        if backend not in available_2d_meshers():
            pytest.skip("Triangle is optional and not built")
    else:
        require_pilot_mesher()
    polygon = Polygon(
        [(0, 0), (0.009, 0), (12, 0), (12, 8), (0, 8)],
        holes=[[(4, 2), (4, 5), (7, 5), (7, 2)]],
    )
    import dtcc_mesher as dm

    d = dm.Domain.from_loops(
        np.array(polygon.exterior.coords)[:-1],
        [np.array(r.coords)[:-1] for r in polygon.interiors],
    )
    p, t, a = mesh_face(
        d.points,
        d.segments,
        polygon,
        SurfaceRefinementOptions(backend=backend),
        spacing=4.0,
    )
    assert a["backend"] == backend and a["refinement_enabled"]
    assert (
        polygon.symmetric_difference(shapely.union_all(shapely.polygons(p[t]))).area
        < 1e-7
    )
    xyz = np.pad(p, ((0, 0), (0, 1)))
    tri = xyz[t]
    lengths = np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2)
    assert lengths.max() < 6.0
    assert len(t) > 20


def test_bounded_failure_returns_original_complete_component():
    old, owners = fixture_shell()
    result, audit = refine_shared_surface(
        old, owners, SurfaceRefinementOptions(max_vertices=len(old.vertices))
    )
    assert result is old and not audit["accepted"]
    assert (
        "budget" in audit["error"] and "original complete shell" in audit["transaction"]
    )


def test_mesher_failure_has_actionable_local_reproducer():
    old, owners = fixture_shell()

    def fail(*args, **kwargs):
        raise RuntimeError("controlled failure")

    result, audit = refine_shared_surface(
        old, owners, SurfaceRefinementOptions(), adapter=fail
    )
    assert result is old and not audit["accepted"]
    assert audit["reproducer"]["points"] and audit["reproducer"]["segments"]
    assert "controlled failure" in audit["error"]


@pytest.mark.parametrize("offset", [(0.0, 0.0, 0.0), (386000.0, 6174000.0, 40.0)])
def test_roof_wall_terrain_preserves_domains_planes_and_shared_ids(offset, tmp_path):
    require_pilot_mesher()
    from dtcc_core import io
    from dtcc_core.builder.meshing.wall_faces import separate_wall_mesh

    old, owners = fixture_shell()
    old.vertices += offset
    result, audit = refine_shared_surface(
        old, owners, SurfaceRefinementOptions(max_rounds=40)
    )
    assert audit["accepted"], audit.get("error")
    assert audit["shared_boundary_mismatches"] == 0
    triangles = result.vertices[result.faces]
    u = np.roll(triangles, -1, axis=1) - triangles
    v = np.roll(triangles, 1, axis=1) - triangles
    angles = np.degrees(
        np.arctan2(np.linalg.norm(np.cross(u, v), axis=2), np.sum(u * v, axis=2))
    ).min(axis=1)
    exempt = next(
        f.values.astype(bool)
        for f in result.fields
        if f.name == "refinement_quality_exempt"
    )
    assert angles[~exempt].min() >= 19.999
    wall = (result.markers >= 0) & (result.markers < len(owners))
    assert angles[wall].min() >= 19.999
    assert all(r["geometry"]["coverage_error_m2"] < 1e-7 for r in audit["faces"])
    assert all(r["geometry"]["plane_error_m"] < 1e-7 for r in audit["faces"])
    assert validate_shell(
        result, num_patches=len(owners), closure_audit=meshes._shell_closure_audit
    )["valid"]
    result.save(tmp_path / "pilot.dtcc")
    restored = io.load_mesh(tmp_path / "pilot.dtcc")
    assert len(restored.regions) == len(audit["faces"])
    parts = separate_wall_mesh(restored, len(owners))
    welded = meshes._welded_mesh(parts)
    assert len(welded.faces) == len(result.faces)
    assert validate_shell(
        welded, num_patches=len(owners), closure_audit=meshes._shell_closure_audit
    )["valid"]
    assert sum(len(m.regions) for m in parts) == len(result.regions)


def test_superseded_pilot_is_not_a_public_assembly_route():
    from dtcc_core.model import City

    with pytest.raises(TypeError, match="surface_refinement"):
        meshes.build_city_surface_mesh(
            City(), surface_refinement=SurfaceRefinementOptions()
        )


@pytest.mark.parametrize("configuration", ["step", "crossing", "acute"])
def test_generic_patch_interfaces_preserve_exposed_geometry(configuration):
    require_pilot_mesher()
    from shapely.geometry import box

    from tests.builder.test_roof_catalogue import assemble, model

    if configuration == "acute":
        foot = Polygon([(0.0, 0.0), (10.0, 0.0), (10.0, 0.4)])
        patches = model("flat", [8.0], foot).patches(foot)
    else:
        a, b = box(-10, -5, 0, 5), box(0, -5, 10, 5)
        left = model(
            "shed" if configuration == "crossing" else "flat",
            [8.0, 0.0, 0.5] if configuration == "crossing" else [8.0],
            a,
            key="a",
        )
        right = model(
            "shed" if configuration == "crossing" else "flat",
            [8.0, 0.0, -0.5] if configuration == "crossing" else [11.0],
            b,
            key="b",
        )
        patches = [*left.patches(a), *right.patches(b)]
    old, patches = assemble(patches, slope=0.3)
    result, audit = refine_shared_surface(
        old, [p.source_key for p in patches], SurfaceRefinementOptions(max_rounds=20)
    )
    assert audit["accepted"], audit.get("error")
    assert audit["topology"]["valid"]
    assert len(audit["rounds"]) > 1
    for r in audit["faces"]:
        assert r["geometry"]["coverage_error_m2"] < 1e-7
    if configuration == "crossing":
        assert any(np.allclose(v, [0, 0, 8]) for v in result.vertices)
    assert set(r["source_key"] for r in audit["faces"]) == {
        None,
        *[p.source_key for p in patches],
    }


def test_corner_spacing_is_shared_across_all_incident_planes():
    from types import SimpleNamespace

    from dtcc_core.builder.meshing.shared_edges import seed_corner_shells

    v = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 4]])
    r = SharedEdges(v, {0: [(0, 1), (1, 2), (2, 0)], 1: [(0, 3), (3, 1), (1, 0)]})
    faces = [SimpleNamespace(rings=[[0, 1, 2]]), SimpleNamespace(rings=[[0, 3, 1]])]
    seed = seed_corner_shells(r, faces, nominal_spacing=4, max_vertices=100)
    assert seed["corners"] > 0 and seed["vertices"] > 0
    # The short horizontal feature also controls the vertical wall edge.
    for edge in r.edges:
        if 0 in edge.endpoints:
            distances = np.linalg.norm(
                np.array(r.vertices)[edge.vertices] - v[0], axis=1
            )
            assert np.any(abs(distances - 0.25) < 1e-10)
    assert np.array_equal(np.array(r.vertices)[: len(v)], v)


def test_geometric_overlay_failure_is_transactional(monkeypatch):
    import dtcc_core.builder.meshing.surface_refinement as sr

    old, owners = fixture_shell()

    def fail(*args, **kwargs):
        raise shapely.errors.GEOSException("controlled geometric overlay failure")

    monkeypatch.setattr(sr, "extract_surface_faces", fail)
    result, audit = refine_shared_surface(old, owners, SurfaceRefinementOptions())
    assert result is old and not audit["accepted"]
    assert "controlled geometric overlay" in audit["error"]


@pytest.mark.parametrize(
    "options",
    [
        {"max_seconds": float("inf")},
        {"nominal_spacing": float("nan")},
        {"max_rounds": 1.5},
        {"max_vertices": True},
    ],
)
def test_pilot_budgets_must_be_finite_and_discrete(options):
    with pytest.raises(ValueError):
        SurfaceRefinementOptions(**options)


def test_synchronization_budget_retains_whole_shell_without_blame_on_last_face():
    require_pilot_mesher()
    old, owners = fixture_shell()
    result, audit = refine_shared_surface(
        old, owners, SurfaceRefinementOptions(max_rounds=1)
    )
    assert result is old and not audit["accepted"]
    assert "synchronization budget" in audit["error"]
    assert audit["rounds"][0]["unpropagated_faces"] > 0
    assert audit["failed_face"] is None


def test_missing_pilot_backend_extension_is_explicit_and_transactional(monkeypatch):
    from types import SimpleNamespace

    import dtcc_core.builder.meshing.surface_refinement_backend as backend

    monkeypatch.setattr(
        backend,
        "_load_dtcc_mesher",
        lambda: SimpleNamespace(
            MeshingOptions=SimpleNamespace(__dataclass_fields__={})
        ),
    )
    old, owners = fixture_shell()
    result, audit = refine_shared_surface(old, owners, SurfaceRefinementOptions())
    assert result is old and not audit["accepted"]
    assert "requires DTCC Mesher corner-reuse extension" in audit["error"]
