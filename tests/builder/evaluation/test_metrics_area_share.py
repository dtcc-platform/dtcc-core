import math

import numpy as np
import pytest

from dtcc_core.builder.geometry_builders.roof_config import RoofPlane
from dtcc_core.builder.evaluation.metrics import plane_areas, top2_area_share


def _square_plane(side: float) -> RoofPlane:
    boundary = np.array([
        [0.0, 0.0, 5.0],
        [side, 0.0, 5.0],
        [side, side, 5.0],
        [0.0, side, 5.0],
    ])
    return RoofPlane(
        normal=np.array([0.0, 0.0, 1.0]),
        offset=0.0,
        inliers=np.arange(10),
        boundary_3d=boundary,
    )


def _plane(area: float) -> RoofPlane:
    return _square_plane(math.sqrt(area) if area > 0 else 0.0)


def test_plane_areas_extracts_areas():
    planes = [_plane(3.0), _plane(5.0), _plane(2.0)]
    areas = plane_areas(planes)
    assert len(areas) == 3
    assert areas[0] == pytest.approx(3.0, abs=1e-9)
    assert areas[1] == pytest.approx(5.0, abs=1e-9)
    assert areas[2] == pytest.approx(2.0, abs=1e-9)


def test_top2_area_share_with_three_planes():
    planes = [_plane(3.0), _plane(5.0), _plane(2.0)]  # total=10, top2=8
    rec = {"planes": planes}
    assert top2_area_share(rec) == pytest.approx(0.8, abs=1e-9)


def test_top2_area_share_single_plane_is_one():
    planes = [_plane(4.0)]
    rec = {"planes": planes}
    assert top2_area_share(rec) == pytest.approx(1.0, abs=1e-9)


def test_top2_area_share_no_planes_returns_none():
    rec = {"planes": []}
    assert top2_area_share(rec) is None


def test_top2_area_share_zero_total_area_returns_none():
    planes = [_plane(0.0), _plane(0.0)]
    rec = {"planes": planes}
    assert top2_area_share(rec) is None
