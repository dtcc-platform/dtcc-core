"""Compact geometry, shared simplification and conservative quality budgets."""

from dataclasses import replace

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon

from dtcc_core.builder.meshing.compact_geometry import (
    building_domains,
    geometric_graph,
    planar_face,
)
from dtcc_core.builder.meshing.compact_surface import CompactSurfaceOptions
from dtcc_core.builder.meshing.compact_triangulation import point_trial, triangulate
from dtcc_core.builder.roofs import RoofPatch
from dtcc_core.model import Mesh


def plane(points, source="a", kind="roof", normal=(0, 0, 1)):
    return planar_face(
        [np.array(points, dtype=float)],
        kind=kind,
        source=source,
        markers=[0],
        normal=normal,
    )


def triangulate_polygon(p, backend="dtcc_mesher"):
    rings = [np.array(r.coords)[:-1] for r in [p.exterior, *p.interiors]]
    points = np.concatenate(rings)
    edges = []
    offset = 0
    for r in rings:
        edges.extend((offset + i, offset + (i + 1) % len(r)) for i in range(len(r)))
        offset += len(r)
    return points, np.array(edges), triangulate(points, edges, p, backend=backend)[0]


@pytest.mark.parametrize("backend", ["dtcc_mesher", "geos", "triangle"])
@pytest.mark.parametrize("width,height", [(4, 3), (100, 3), (0.001, 12)])
def test_rectangles_are_two_triangles_even_long_or_thin(backend, width, height):
    if backend == "triangle":
        from dtcc_core.builder.meshing.backends import available_2d_meshers

        if backend not in available_2d_meshers():
            pytest.skip("Triangle is optional")
    p = Polygon([(0, 0), (width, 0), (width, height), (0, height)])
    points, edges, t = triangulate_polygon(p, backend)
    assert len(points) == 4 and len(t) == 2
    assert shapely.union_all(shapely.polygons(points[t])).equals(p)


@pytest.mark.parametrize("backend", ["dtcc_mesher", "geos", "triangle"])
def test_holes_and_mandatory_collinear_junctions(backend):
    if backend == "triangle":
        from dtcc_core.builder.meshing.backends import available_2d_meshers

        if backend not in available_2d_meshers():
            pytest.skip("Triangle is optional")
    p = Polygon(
        [(0, 0), (1, 0), (8, 0), (8, 8), (0, 8)],
        holes=[[(2, 2), (2, 4), (4, 4), (4, 2)]],
    )
    points, edges, t = triangulate_polygon(p, backend)
    assert len(t) == len(points) + 2 - 2
    assert 1 in t


@pytest.mark.parametrize(
    "angle,offset", [(0, (0, 0, 0)), (0.71, (386000, 6174000, 40))]
)
def test_joint_collinear_removal_is_deterministic(angle, offset):
    a = np.array([[0, 0, 0], [1, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]], float)
    b = np.array([[0, 0, 0], [0, 0, 2], [4, 0, 2], [4, 0, 0], [1, 0, 0]], float)
    R = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    results = []
    for _ in range(2):
        fs = [
            plane(a @ R.T + offset),
            plane(b @ R.T + offset, kind="wall", normal=np.array([0, 1, 0]) @ R.T),
        ]
        registry, audit = geometric_graph(fs)
        assert audit["removed_redundant_vertices"] == 1
        assert all(len(f.rings[0]) == 4 for f in fs)
        assert any(len(e.incidents) == 2 for e in registry.edges)
        results.append(audit)
    assert results[0] == results[1]


def test_t_junction_and_sloped_eave_breakpoint_survive():
    fs = [
        plane([[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]),
        plane([[2, 0, 0], [2, -2, 0], [4, -2, 0], [4, 0, 0]], source="b"),
    ]
    _, audit = geometric_graph(fs)
    assert len(fs[0].rings[0]) == 5
    assert any(r["reason"] == "source interface" for r in audit["retained_vertices"])
    wall = plane(
        [[0, 0, 0], [4, 0, 0], [4, 0, 2], [2, 0, 3], [0, 0, 2]],
        kind="wall",
        normal=(0, -1, 0),
    )
    _, audit = geometric_graph([wall])
    assert len(wall.rings[0]) == 5


def test_coplanar_contact_can_simplify_but_bent_terrain_cannot():
    wall = plane(
        [[0, 0, 0], [1, 0, 0], [4, 0, 0], [4, 0, 3], [0, 0, 3]],
        kind="wall",
        normal=(0, -1, 0),
    )
    terrain = plane(
        [[0, 0, 0], [0, -2, 0], [4, -2, 0], [4, 0, 0], [1, 0, 0]],
        source=None,
        kind="terrain",
    )
    _, a = geometric_graph([wall, terrain])
    assert a["removed_redundant_vertices"] == 1
    fs = [
        plane(
            [[0, 0, 0], [4, 0, 0], [4, 0, 3], [0, 0, 3]], kind="wall", normal=(0, -1, 0)
        ),
        plane(
            [[0, 0, 0], [0, -2, 0], [1, -2, 0], [1, 0, 0]], source=None, kind="terrain"
        ),
        plane(
            [[1, 0, 0], [1, -2, 0], [4, -2, 1]],
            source=None,
            kind="terrain",
            normal=(-2, 0, 6),
        ),
        plane(
            [[1, 0, 0], [4, -2, 1], [4, 0, 0]],
            source=None,
            kind="terrain",
            normal=(0, 3, 6),
        ),
    ]
    _, a = geometric_graph(fs)
    assert len(fs[0].rings[0]) == 5
    assert any(
        r["reason"] == "necessary terrain contact" for r in a["retained_vertices"]
    )


def reference_base(patches):
    # Only the established block base is needed by domain construction.
    return Mesh(
        vertices=np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        faces=np.array([[0, 1, 2]] * len(patches)),
        markers=np.arange(len(patches)),
    )


def patch(p, coeff, source="a", i=0):
    return RoofPatch(Polygon(p), np.zeros(2), np.array(coeff, float), source, i, 0)


def test_different_height_sources_and_crossing_exposed_intervals():
    from dtcc_core.builder.meshing.roof_interfaces import node_interfaces

    patches = [
        patch([(0, 0), (2, 0), (2, 4), (0, 4)], [0, 1, 3]),
        patch([(2, 0), (4, 0), (4, 4), (2, 4)], [0, -1, 7], "b", 1),
    ]
    patches = node_interfaces(patches)
    fs = building_domains(patches, reference_base(patches))
    _, a = geometric_graph(fs)
    steps = [
        f
        for f in fs
        if f.kind == "wall" and abs(f.origin[0] - 2) < 1e-8 and abs(f.normal[0]) > 0.99
    ]
    assert {f.source_key for f in steps} == {"a", "b"}
    assert sum(f.polygon.area for f in steps) == pytest.approx(8)
    assert any(
        abs(r["xyz"][1] - 2) < 1e-8 and abs(r["xyz"][2] - 5) < 1e-8
        for r in a["retained_vertices"]
    )


def test_simple_hip_domains_keep_four_planes_and_minimal_walls():
    patches = [
        patch([(0, 0), (4, 0), (2, 2)], [0, 1, 3]),
        patch([(4, 0), (4, 4), (2, 2)], [-1, 0, 7]),
        patch([(4, 4), (0, 4), (2, 2)], [0, -1, 7]),
        patch([(0, 4), (0, 0), (2, 2)], [1, 0, 3]),
    ]
    fs = building_domains(patches, reference_base(patches))
    _, a = geometric_graph(fs)
    assert sum(f.kind == "roof" for f in fs) == 4
    assert sum(f.kind == "wall" for f in fs) == 4
    assert len(a["retained_vertices"]) == 9


def test_equivalent_same_owner_roof_patches_merge_without_losing_provenance():
    patches = [
        patch([(0, 0), (2, 0), (2, 2), (0, 2)], [0, 0, 3]),
        patch([(2, 0), (4, 0), (4, 2), (2, 2)], [0, 0, 3]),
    ]
    fs = building_domains(patches, reference_base(patches))
    roofs = [f for f in fs if f.kind == "roof"]
    assert len(roofs) == 1 and roofs[0].markers == [2, 3]
    geometric_graph(fs)
    assert len(roofs[0].rings[0]) == 4


def test_budget_limits_and_reject_optional_vertex_that_moves_sliver():
    for name, value in [
        ("max_face_vertices", 3),
        ("max_source_vertices", 9),
        ("max_group_vertices", 33),
    ]:
        with pytest.raises(ValueError):
            CompactSurfaceOptions(**{name: value})
    p = Polygon([(0, 0), (10, 0), (10, 0.001), (0, 0.001)])
    points, edges, t = triangulate_polygon(p)
    assert point_trial(points, t, p, edges) is None


def test_explicit_mode_validation():
    from dtcc_core.builder.geometry_builders.meshes import build_city_surface_mesh

    with pytest.raises(ValueError, match="requires lod='lod2'"):
        build_city_surface_mesh(None, compact_surface=CompactSurfaceOptions())


@pytest.mark.parametrize("slope", [0.0, 2.0])
@pytest.mark.parametrize("budgeted", [False, True])
def test_compact_shell_terrain_ownership_and_separate_round_trip(slope, budgeted):
    from dtcc_core.builder.geometry_builders.meshes import (
        _shell_closure_audit,
        _welded_mesh,
    )
    from dtcc_core.builder.meshing.compact_surface import compact_surface
    from dtcc_core.builder.meshing.roof_validation import validate_shell
    from dtcc_core.builder.meshing.wall_faces import separate_wall_mesh

    from .test_roof_catalogue import assemble

    patches = [
        patch([(0, 0), (8, 0), (8, 6), (0, 6)], [0, 0.1, 4]),
        patch([(8, 0), (14, 0), (14, 6), (8, 6)], [0, 0, 6], "b", 1),
    ]
    old, patches = assemble(patches, slope=slope)
    mesh, a = compact_surface(old, patches, CompactSurfaceOptions(budgeted=budgeted))
    assert a["topology"]["valid"] and a["shared_boundary_mismatches"] == 0
    assert all(
        sum(v.get("kind", "interior") == "interior" for v in r["optional_vertices"])
        <= 2
        for r in a["faces"]
    )
    assert all(n <= 8 for n in a["source_optional_vertices"].values())
    assert len(a["optional_vertices"]) <= 16
    assert all(r["geometry"]["plane_error_m"] <= 1e-7 for r in a["faces"])
    assert (
        a["unchanged_terrain_triangles"]
        >= a["terrain_locality"]["outside_star_triangles"]
        > 0
    )
    assert a["terrain_locality"]["outside_star_bitwise_preserved"]
    separate = separate_wall_mesh(mesh, len(patches))
    welded = _welded_mesh(separate)
    assert validate_shell(
        welded, num_patches=len(patches), closure_audit=_shell_closure_audit
    )["valid"]
    assert {
        r.attributes["source_key"]
        for r in mesh.regions
        if r.attributes["kind"] != "terrain"
    } == {"a", "b"}


def test_courtyard_is_preserved_by_direct_domain_construction():
    from dtcc_core.builder.meshing.compact_surface import compact_surface

    from .test_roof_catalogue import assemble

    p = RoofPatch(
        Polygon(
            [(0, 0), (12, 0), (12, 12), (0, 12)],
            holes=[[(4, 4), (4, 8), (8, 8), (8, 4)]],
        ),
        np.zeros(2),
        np.array([0.0, 0, 5]),
        "court",
    )
    old, patches = assemble([p])
    mesh, a = compact_surface(old, patches, CompactSurfaceOptions(budgeted=False))
    roofs = [r for r in a["faces"] if r["kind"] == "roof"]
    assert len(roofs) == 1 and roofs[0]["holes"] == 1
    assert roofs[0]["triangles"] == roofs[0]["boundary_only_expected"]
    assert a["topology"]["valid"]
