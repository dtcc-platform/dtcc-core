"""Numerical policies retain physical steps and exact DEM boundary extrema."""

import numpy as np
import pytest
import shapely
from affine import Affine
from shapely.geometry import Polygon, box
from shapely import affinity

from dtcc_core.model import Raster
from dtcc_core.builder.meshing.terrain_attachment import minimum_dem
from dtcc_core.builder.meshing.roof_interfaces import node_interfaces
from dtcc_core.builder.meshing.compact_surface import (
    CompactSurfaceOptions,
    compact_surface,
)
from dtcc_core.builder.roofs import RoofPatch


@pytest.mark.parametrize("offset", [(0.0, 0.0), (647700.0, 6638400.0)])
def test_foundation_minimum_includes_boundary_quadratic_extremum(offset):
    # z=-x*y on [0,1]^2. On x+y=1 its minimum is -.25 at the
    # edge midpoint, although every polygon vertex and contained DEM node is 0.
    x, y = offset
    raster = Raster(
        data=np.array([[0.0, 0.0], [0.0, -1.0]]),
        georef=Affine(1.0, 0.0, x - 0.5, 0.0, 1.0, y - 0.5),
    )
    polygon = Polygon([(x, y), (x + 1, y), (x, y + 1)])
    assert minimum_dem(polygon, raster) == pytest.approx(-0.25, abs=1e-12)


@pytest.mark.parametrize("slope", [0.0, 0.5])
@pytest.mark.parametrize("offset", [(0.0, 0.0), (647700.0, 6638400.0)])
def test_numerical_height_reconciliation_retains_a_real_small_roof_step(slope, offset):
    patches = [
        RoofPatch(
            affinity.translate(box(0, 0, 1, 1), *offset),
            np.asarray(offset),
            np.array([slope, 0.0, 1.0]),
            "a",
            0,
            0,
        ),
        RoofPatch(
            affinity.translate(box(1, 0, 2, 1), *offset),
            np.asarray(offset),
            np.array([slope, 0.0, 1.0 + 1e-7]),
            "b",
            1,
            0,
        ),
    ]
    mesh, audit = compact_surface(
        None,
        node_interfaces(patches),
        CompactSurfaceOptions(budgeted=False),
        base_z=0.0,
    )
    assert audit["topology"]["valid"]
    interface = mesh.vertices[
        np.isclose(mesh.vertices[:, 0], offset[0] + 1.0, atol=1e-12, rtol=0)
    ]
    assert np.any(np.isclose(interface[:, 2], 1.0 + slope, atol=1e-12, rtol=0))
    assert np.any(np.isclose(interface[:, 2], 1.0 + slope + 1e-7, atol=1e-12, rtol=0))


@pytest.mark.parametrize("offset", [(0.0, 0.0), (647700.0, 6638400.0)])
def test_projected_outline_samples_share_one_canonical_breakpoint(offset):
    from scipy.spatial import distance
    from shapely import affinity
    from shapely.geometry import LineString
    from dtcc_core.builder.cleaning.subdivision import canonical_boundary_graph

    domain = Polygon([(0, 0), (10, 1), (10, 10), (0, 10)])
    # Independently rounded samples on a sloping boundary project only 10 nm
    # apart. Both must enter the *same split outline*, rather than compete with
    # an unchanged full-length outline segment.
    lines = [
        domain.boundary,
        LineString([(5.0, 0.5000001), (5.0, 5.0)]),
        LineString([(5.0, 0.5), (4.0, 5.0)]),
    ]
    graph = shapely.union_all([affinity.translate(line, *offset) for line in lines])
    audit = {}
    result = canonical_boundary_graph(
        graph, {0: affinity.translate(domain, *offset)}, audit
    )
    coordinates = np.unique(
        np.concatenate([np.asarray(line.coords) for line in result.geoms]), axis=0
    )
    assert distance.pdist(coordinates).min() > 1e-7
    assert (
        audit["maximum_boundary_projection_m"] <= audit["boundary_projection_budget_m"]
    )


def test_shared_steep_roof_junction_minimizes_normal_plane_error():
    from dtcc_core.builder.meshing.compact_geometry import (
        building_domains,
        face_coordinates,
    )

    cut = 5e-8
    patches = [
        RoofPatch(
            box(-1, 0, cut, 1), np.zeros(2), np.array([10.0, 0.0, 20.0]), "a", 0, 0
        ),
        RoofPatch(
            box(cut, 0, 1, 1), np.zeros(2), np.array([-10.0, 0.0, 20.0]), "b", 1, 0
        ),
    ]
    faces = building_domains(patches, None, base_z=0.0)
    roofs = [face for face in faces if face.kind == "roof"]
    for face in roofs:
        assert face.original_plane_error <= 1e-7
        points = np.concatenate(face_coordinates(face))
        ridge = points[points[:, 0] == cut]
        np.testing.assert_allclose(ridge[:, 2], 20.0, rtol=0, atol=1e-12)


@pytest.mark.parametrize("offset", [(0.0, 0.0), (647700.0, 6638400.0)])
def test_vectorized_evidence_predicate_matches_scalar_selection(offset):
    from types import SimpleNamespace
    from shapely import affinity
    from shapely.geometry import Point
    from dtcc_core.builder.geometry_builders.lod2 import (
        _level_point_counts_in_region,
        EDGE_TOLERANCE,
    )

    polygon = affinity.translate(
        Polygon([(0, 0), (4, 0), (4, 4), (0, 4)], [[(1, 1), (2, 1), (2, 2), (1, 2)]]),
        *offset,
    )
    xy = (
        np.vstack(
            [
                np.random.default_rng(5).uniform(-1, 5, (1000, 2)),
                [[0, 0], [4, 4], [1, 1], [1.5, 1.5], [-EDGE_TOLERANCE, 0]],
            ]
        )
        + offset
    )
    points = np.c_[xy, np.arange(len(xy))]
    levels = [
        SimpleNamespace(point_indices=indices)
        for indices in (
            np.arange(len(xy)),
            np.arange(0, len(xy), 3),
            np.empty(0, dtype=int),
        )
    ]
    buffered = polygon.buffer(EDGE_TOLERANCE)
    scalar = [
        sum(buffered.covers(Point(*xy[index])) for index in level.point_indices)
        for level in levels
    ]
    assert _level_point_counts_in_region(points, polygon, levels) == scalar


@pytest.mark.parametrize("offset", [(0.0, 0.0), (647700.0, 6638400.0)])
def test_wall_winding_uses_oriented_incidence_at_a_thin_roof_tip(offset):
    from shapely import affinity

    polygon = affinity.translate(Polygon([(0, 0), (4, 0), (2, 5e-7)]), *offset)
    patch = RoofPatch(
        polygon, np.asarray(offset), np.array([0.0, 0.0, 1.0]), "thin-tip", 0, 0
    )
    mesh, audit = compact_surface(
        None,
        node_interfaces([patch]),
        CompactSurfaceOptions(budgeted=False),
        base_z=0.0,
    )
    # A fixed-distance inside probe can cross this tip. Ring incidence has no
    # probe distance and must preserve winding despite the constrained angle.
    assert audit["topology"]["valid"]
    assert audit["topology"]["closure"]["inconsistent_edges"] == 0
    assert len(mesh.faces) == 7


def test_repeated_collinear_removals_bound_total_original_boundary_error():
    from dtcc_core.builder.meshing.compact_geometry import planar_face, geometric_graph

    # Removing the 18 nm peak looks locally like a 9 nm change. Removing both
    # neighbouring 9 nm samples afterwards must not erase the entire wave.
    xyz = np.array(
        [
            [2.0, 0.0, 18e-9],
            [3.0, 0.0, 9e-9],
            [4.0, 0.0, 0.0],
            [4.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 9e-9],
        ]
    )
    domain = planar_face(
        [xyz], kind="wall", source="wavy", markers=[0], normal=[0.0, -1.0, 0.0]
    )
    registry, _ = geometric_graph([domain])
    actual = Polygon(
        (np.asarray(registry.vertices)[domain.rings[0]] - domain.origin)
        @ domain.basis.T
    )
    assert domain.polygon.boundary.hausdorff_distance(actual.boundary) <= 1e-8
    assert len(domain.rings[0]) > 4
