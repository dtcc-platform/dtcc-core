"""Tests for the shared bounding box helpers."""

from __future__ import annotations

import numpy as np
import pytest

from dtcc_core.reproject import bounds as reproject_bounds


def test_bounds_to_wgs84_returns_validated_bounds_for_wgs84_aliases(monkeypatch):
    def fail_reproject(*args, **kwargs):
        raise AssertionError("WGS84 aliases should not be reprojected")

    monkeypatch.setattr(reproject_bounds, "reproject_array", fail_reproject)

    assert reproject_bounds.bounds_to_wgs84((11, 57, 12, 58), "CRS84") == (
        11.0,
        57.0,
        12.0,
        58.0,
    )
    assert reproject_bounds.bounds_to_wgs84([11, 57, 12, 58], "EPSG:4326") == (
        11.0,
        57.0,
        12.0,
        58.0,
    )
    assert reproject_bounds.bounds_to_wgs84((11, 57, 12, 58), "WGS84") == (
        11.0,
        57.0,
        12.0,
        58.0,
    )


def test_bounds_to_wgs84_uses_all_four_bbox_corners(monkeypatch):
    captured = {}

    def fake_reproject(points, source_crs, target_crs):
        captured["points"] = np.asarray(points)
        captured["source_crs"] = source_crs
        captured["target_crs"] = target_crs
        x = points[:, 0]
        y = points[:, 1]
        return np.column_stack((x + y, y - x, np.zeros(points.shape[0])))

    monkeypatch.setattr(reproject_bounds, "reproject_array", fake_reproject)

    assert reproject_bounds.bounds_to_wgs84((0, 0, 2, 3), "EPSG:3006") == (
        0.0,
        -2.0,
        5.0,
        3.0,
    )
    assert captured["source_crs"] == "EPSG:3006"
    assert captured["target_crs"] == "EPSG:4326"
    assert captured["points"].tolist() == [
        [0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 3.0, 0.0],
    ]


@pytest.mark.parametrize(
    ("bounds", "message"),
    [
        ((0, 0, 1), "exactly four"),
        ((1, 0, 1, 2), "xmin < xmax"),
        ((0, 2, 1, 2), "ymin < ymax"),
        ((0, 0, float("nan"), 2), "finite"),
        (("x", 0, 1, 2), "numeric"),
    ],
)
def test_bounds_to_wgs84_rejects_invalid_bounds(bounds, message):
    with pytest.raises(ValueError, match=message):
        reproject_bounds.bounds_to_wgs84(bounds, "EPSG:3006")


def test_bounds_to_wgs84_rejects_missing_crs():
    with pytest.raises(ValueError, match="source_crs"):
        reproject_bounds.bounds_to_wgs84((0, 0, 1, 1), "")
