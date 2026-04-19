import numpy as np

from dtcc_core.builder.geometry_builders.roof_config import RoofPlane
from dtcc_core.builder.evaluation.metrics import point_coverage_ratio


def _plane_with_inliers(n):
    boundary = np.array([[0.0, 0.0, 5.0], [1.0, 0.0, 5.0], [0.0, 1.0, 5.0]])
    return RoofPlane(
        normal=np.array([0.0, 0.0, 1.0]),
        offset=0.0,
        inliers=np.arange(n),
        boundary_3d=boundary,
    )


def test_coverage_ratio_is_sum_inliers_over_filtered():
    planes = [_plane_with_inliers(40), _plane_with_inliers(30)]
    rec = {"planes": planes, "classification": None}
    filtered_n = 100
    assert point_coverage_ratio(rec, filtered_n) == 0.70


def test_coverage_ratio_empty_planes():
    rec = {"planes": [], "classification": None}
    assert point_coverage_ratio(rec, 100) == 0.0


def test_coverage_ratio_zero_filtered_returns_none():
    rec = {"planes": [_plane_with_inliers(10)], "classification": None}
    assert point_coverage_ratio(rec, 0) is None


def test_coverage_ratio_clamps_at_one():
    planes = [_plane_with_inliers(80), _plane_with_inliers(80)]
    rec = {"planes": planes, "classification": None}
    assert point_coverage_ratio(rec, 100) == 1.0
