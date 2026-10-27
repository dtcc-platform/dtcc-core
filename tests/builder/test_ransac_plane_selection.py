"""Translation invariance and supported multi-plane RANSAC reconstruction."""

import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import box

from dtcc_core.builder.geometry_builders.lod2 import (
    RoofPlane,
    _candidate_lod2_from_parts,
    _fit_plane,
    _points_on_edge,
    is_watertight,
)
from dtcc_core.builder.geometry_builders.roof_plane_selection import envelope_candidates
from dtcc_core.builder.roofs.ransac import fit_ransac


def grid():
    x, y = np.meshgrid(np.linspace(-9.8, 9.8, 60), np.linspace(-5.8, 5.8, 40))
    return np.c_[x.ravel(), y.ravel()]


@pytest.mark.parametrize(
    "offset", [[0.0, 0.0], [91600.0, 6231250.0], [386500.0, 6175000.0]]
)
def test_dense_plane_fit_preserves_both_slopes_in_map_coordinates(offset):
    xy = grid()
    points = np.c_[xy + offset, 10.0 + 0.3 * xy[:, 0] - 0.6 * xy[:, 1]]
    before = points.copy()
    plane = _fit_plane(points)
    np.testing.assert_allclose([plane.a, plane.b], [0.3, -0.6], atol=1e-9)
    np.testing.assert_allclose(
        [plane.z_at(*p) for p in points[:, :2]], points[:, 2], atol=1e-8
    )
    np.testing.assert_array_equal(points, before)


@pytest.mark.parametrize("kind", ["hip", "valley"])
def test_continuous_multiplane_roof_is_watertight_and_accurate(kind):
    xy = grid()
    coefficients = np.array(
        [[0.6, 0, 12], [-0.6, 0, 12], [0, 0.6, 9.6], [0, -0.6, 9.6]]
    )
    z = (xy @ coefficients[:, :2].T + coefficients[:, 2]).min(axis=1)
    if kind == "valley":
        z = 18 - z
    points = np.c_[xy, z]
    polygon = box(-10, -6, 10, 6)
    shell = _candidate_lod2_from_parts(polygon, points, 0.0)
    assert shell is not None and is_watertight(shell)
    roof = fit_ransac(polygon, points, 0.0, source_key="test")
    assert roof.accepted, roof.diagnostics
    assert roof.diagnostics["plane_selection"]["selected_planes"] >= 3
    np.testing.assert_allclose(roof.model.predict(xy), z, atol=0.12)
    offset = np.array([91600.0, 6231250.0])
    shifted = fit_ransac(
        translate(polygon, *offset), points + [*offset, 0], 0.0, source_key="map"
    )
    assert shifted.accepted, shifted.diagnostics
    np.testing.assert_allclose(
        shifted.model.predict(xy + offset), roof.model.predict(xy), atol=1e-6
    )


def test_minor_plane_does_not_discard_supported_gable():
    xy = grid()
    z = 12 - 0.6 * abs(xy[:, 0])
    # Small elevated patch supplies a third detected plane without changing
    # the majority roof. It must not force a flat roof or an unsupported face.
    artifact_xy = xy[:80]
    points = np.vstack([np.c_[xy, z], np.c_[artifact_xy, np.full(80, 16.0)]])
    fit = fit_ransac(box(-10, -6, 10, 6), points, 0.0, source_key="minor")
    assert fit.accepted, fit.diagnostics
    selection = fit.diagnostics["plane_selection"]
    assert selection["detected_planes"] > selection["selected_planes"] == 2
    np.testing.assert_allclose(fit.model.predict(xy), z, atol=0.12)


def test_selection_does_not_bridge_balanced_height_steps():
    xy = grid()
    z = np.where(xy[:, 0] < 0, 8.0, 12.0)
    points = np.c_[xy, z]
    planes = [
        RoofPlane(0, 0, 8, np.flatnonzero(xy[:, 0] < 0)),
        RoofPlane(0, 0, 12, np.flatnonzero(xy[:, 0] >= 0)),
        RoofPlane(0.6, 0, 10, np.arange(20)),
    ]
    assert not list(
        envelope_candidates(
            box(-10, -6, 10, 6),
            points,
            planes,
            tolerance=0.4,
            min_points=20,
            min_area=2,
        )
    )


def test_selection_rejects_planes_without_spatial_support():
    xy = grid()
    z = 12 - 0.6 * abs(xy[:, 0])
    points = np.c_[xy, z]
    # Correct equations but swapped evidence: a low residual alone is not enough.
    planes = [
        RoofPlane(0.6, 0, 12, np.flatnonzero(xy[:, 0] > 0)),
        RoofPlane(-0.6, 0, 12, np.flatnonzero(xy[:, 0] < 0)),
    ]
    assert not list(
        envelope_candidates(
            box(-10, -6, 10, 6),
            points,
            planes,
            tolerance=0.4,
            min_points=20,
            min_area=2,
        )
    )


def test_wall_keeps_roof_corner_junctions_with_millimetre_tolerance():
    # On a 20 m edge a dimensionless .001 tolerance would erase this 1 cm edge.
    points = [
        np.array(p) for p in [[0.0, 0.0, 6.0], [0.01, 0.0, 6.006], [20.0, 0.0, 6.0]]
    ]
    selected = _points_on_edge(
        points, np.array([0.0, 0.0]), np.array([20.0, 0.0]), for_wall=True
    )
    np.testing.assert_array_equal(selected, points)
