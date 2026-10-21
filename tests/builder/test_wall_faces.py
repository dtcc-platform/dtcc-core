"""Wall-only refinement preserves shell geometry and immutable boundaries."""

import json
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon

from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.meshing import wall_faces as wf
from dtcc_core.builder.meshing.roof_validation import vertex_link_audit
from dtcc_core.model import Mesh


def panel_ring(*, angle=0.0, offset=(0.0, 0.0, 0.0), distinct_owner=False):
    # A 9 mm roof breakpoint on a long planar facade, with a sloped eave.
    xy = np.array([[0, 0], [0.009, 0], [1.2, 0], [12, 0], [12, 5], [0, 5]], dtype=float)
    vertices, faces, markers = [], [], []
    for i in range(len(xy)):
        start = len(vertices)
        a, b = xy[i], xy[(i + 1) % len(xy)]
        for point in [a, b]:
            for z in np.linspace(0.0, 8.0 + 0.02 * point[0], 5):
                vertices.append([*point, z])
        for k in range(4):
            faces.extend(
                [
                    [start + k, start + 5 + k, start + 6 + k],
                    [start + k, start + 6 + k, start + k + 1],
                ]
            )
            markers.extend([i, i])
    mesh = meshes._welded_mesh(
        [
            Mesh(
                vertices=np.array(vertices),
                faces=np.array(faces),
                markers=np.array(markers),
            )
        ]
    )
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    mesh.vertices = mesh.vertices @ rotation.T + offset
    owners = ["one"] * len(xy)
    if distinct_owner:
        owners[0] = "two"
    return mesh, owners


@pytest.mark.parametrize(
    "angle,offset", [(0.0, (0.0, 0.0, 0.0)), (0.71, (386000.0, 6174000.0, 40.0))]
)
def test_short_eave_seams_are_removed_with_ownership_and_geometry_preserved(
    angle, offset
):
    old, owners = panel_ring(angle=angle, offset=offset)
    assert np.sum(wf.triangle_quality(old.vertices, old.faces) < 0.02) > 0
    new, audit = wf.refine_wall_faces(old, owners)
    assert audit["replaced_faces"] > 0 and audit["quality_fallbacks"] == 0
    assert wf.triangle_quality(new.vertices, new.faces).min() >= 0.02
    assert vertex_link_audit(new)["nonmanifold_vertices"] == 0
    assert np.array_equal(new.vertices[: len(old.vertices)], old.vertices)
    for row in audit["faces"]:
        if row["method"] == "connected_planar_cdt":
            assert row["coverage_error_m2"] < 1e-10
            assert {owners[m] for m in row["contributing_markers"]} == {
                row["source_key"]
            }
    assert any(len(r["contributing_markers"]) > 1 for r in audit["faces"])
    # Actual footprint corners remain different wall faces.
    assert all(
        set(r["contributing_markers"]) <= {0, 1, 2}
        or r["contributing_markers"] in ([3], [4], [5])
        for r in audit["faces"]
    )


def test_distinct_owners_do_not_merge_and_failed_refinement_keeps_original():
    old, owners = panel_ring(distinct_owner=True)
    new, audit = wf.refine_wall_faces(old, owners, max_rounds=0, max_face_points=3)
    assert audit["quality_fallbacks"] > 0
    assert np.array_equal(new.vertices, old.vertices)
    assert np.array_equal(new.faces, old.faces)
    assert np.array_equal(new.markers, old.markers)
    assert all(
        len({owners[m] for m in r["contributing_markers"]}) == 1 for r in audit["faces"]
    )
    assert len(audit["faces"]) == 5


@pytest.mark.parametrize(
    "polygon",
    [
        Polygon([(0, 0), (12, 0), (12, 5), (7, 5), (7, 3), (4, 3), (4, 6), (0, 6)]),
        Polygon(
            [(0, 0), (12, 0), (12, 8), (0, 8)], holes=[[(4, 2), (4, 5), (7, 5), (7, 2)]]
        ),
        Polygon([(0, 0), (10, 0), (10, 3)]),
    ],
)
def test_concave_holes_and_collapsed_step_domains(polygon):
    # Boundary-only constrained triangulation must preserve all components.
    import dtcc_mesher as dm

    domain = dm.Domain.from_loops(
        np.array(polygon.exterior.coords)[:-1],
        [np.array(r.coords)[:-1] for r in polygon.interiors],
    )
    points, triangles, rounds, immutable = wf._triangulate_face(
        domain.points, domain.segments, polygon
    )
    error, overlap = wf._check_replacement(polygon, points, triangles, domain.segments)
    assert error < 1e-10 and overlap < 1e-10


def test_optional_triangle_is_not_used(monkeypatch):
    from dtcc_core.builder import _dtcc_builder

    def forbidden(*args, **kwargs):
        raise AssertionError("Optional Triangle must not be required")

    monkeypatch.setattr(_dtcc_builder, "build_city_flat_mesh", forbidden)
    old, owners = panel_ring()
    new, audit = wf.refine_wall_faces(old, owners)
    assert audit["replaced_faces"] and not audit["quality_fallbacks"]


def test_separate_output_retains_spanning_wall_mapping_and_welds_conformingly():
    old, owners = panel_ring()
    new, audit = wf.refine_wall_faces(old, owners)
    parts = wf.separate_wall_mesh(new, len(owners))
    assert len(parts) == len(owners) + 1
    welded = meshes._welded_mesh(parts)
    assert len(welded.faces) == len(new.faces)
    assert vertex_link_audit(welded)["nonmanifold_vertices"] == 0
    assert wf.triangle_quality(welded.vertices, welded.faces).min() >= 0.02
    assert all(len(m.fields[0].values) == len(m.faces) for m in parts)


def test_acute_corner_quality_bound_is_explicit():
    polygon = Polygon([(0.0, 0.0), (10.0, 0.0), (10.0, 0.05)])
    constraint = wf._angle_constraint(polygon)
    assert constraint["minimum_input_angle_degrees"] < 0.3
    assert constraint["corner_quality_upper_bound"] < 0.02


def test_short_eave_shell_keeps_roof_terrain_and_persistent_ownership(tmp_path):
    from dtcc_core import io

    data = json.loads(
        (Path(__file__).parent / "fixtures/wall_short_eave.json").read_text()
    )
    old = Mesh(
        vertices=np.array(data["vertices"]),
        faces=np.array(data["faces"]),
        markers=np.array(data["markers"]),
    )
    new, audit = wf.refine_wall_faces(old, data["owners"])
    count = len(data["owners"])
    old_walls = (old.markers >= 0) & (old.markers < count)
    new_walls = (new.markers >= 0) & (new.markers < count)
    assert wf.triangle_quality(old.vertices, old.faces[old_walls]).min() < 0.001
    assert wf.triangle_quality(new.vertices, new.faces[new_walls]).min() >= 0.02
    assert np.array_equal(
        old.vertices[old.faces[~old_walls]], new.vertices[new.faces[~new_walls]]
    )
    from dtcc_core.builder.meshing.roof_validation import validate_shell

    assert validate_shell(
        new, num_patches=count, closure_audit=meshes._shell_closure_audit
    )["valid"]
    assert any(
        r["boundary_vertices_added"]
        for r in audit["faces"]
        if r["method"] == "connected_planar_cdt"
    )
    new.save(tmp_path / "walls.dtcc")
    restored = io.load_mesh(tmp_path / "walls.dtcc")
    assert len(restored.regions) == audit["replaced_faces"]
    for region in restored.regions:
        assert region.attributes["source_key"] == "source"
        assert set(region.attributes["contributing_patch_markers"]) <= set(range(count))
        assert np.all(restored.fields[0].values[region.indices] == int(region.id))
    separate = wf.separate_wall_mesh(new, count)
    assert sum(len(m.regions) for m in separate) == len(new.regions)


def test_global_vertex_budget_rejection_is_bounded_and_geometry_preserving():
    old, owners = panel_ring()
    new, audit = wf.refine_wall_faces(old, owners, max_added_points=0)
    assert len(new.vertices) == len(old.vertices)
    assert np.array_equal(new.faces, old.faces)
    assert any("budget" in (r["fallback_reason"] or "") for r in audit["faces"])
