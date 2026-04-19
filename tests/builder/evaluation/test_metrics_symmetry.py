import math

import numpy as np
import pytest

from dtcc_core.builder.geometry_builders.roof_config import RoofPlane
from dtcc_core.builder.evaluation.metrics import slope_symmetry_error


def _plane_with_tilt(tilt_deg: float, side: float) -> RoofPlane:
    rad = math.radians(tilt_deg)
    n = np.array([math.sin(rad), 0.0, math.cos(rad)])
    boundary = np.array([
        [0.0, 0.0, 5.0],
        [side, 0.0, 5.0],
        [side, side, 5.0],
        [0.0, side, 5.0],
    ])
    return RoofPlane(
        normal=n,
        offset=0.0,
        inliers=np.arange(10),
        boundary_3d=boundary,
    )


def test_perfectly_symmetric_returns_zero():
    planes = [_plane_with_tilt(30.0, 3.0), _plane_with_tilt(30.0, 3.0)]
    assert slope_symmetry_error({"planes": planes}) == pytest.approx(0.0, abs=1e-9)


def test_asymmetric_returns_abs_difference_in_degrees():
    planes = [_plane_with_tilt(30.0, 3.0), _plane_with_tilt(20.0, 2.5)]
    err = slope_symmetry_error({"planes": planes})
    assert err == pytest.approx(10.0, abs=0.01)


def test_single_plane_returns_none():
    planes = [_plane_with_tilt(30.0, 3.0)]
    assert slope_symmetry_error({"planes": planes}) is None


def test_no_planes_returns_none():
    assert slope_symmetry_error({"planes": []}) is None
