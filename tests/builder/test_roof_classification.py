import numpy as np
import pytest

from dtcc_core.model.enums import RoofType
from dtcc_core.builder.geometry_builders.roof_config import (
    RoofDetectionConfig,
    RoofPlane,
)
from dtcc_core.builder.geometry_builders.roof_classification import (
    merge_planes,
    classify_roof,
    ClassificationResult,
)


def _make_plane(normal, offset, n_inliers=50, boundary_z=10.0):
    normal = np.array(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    inliers = np.arange(n_inliers)
    boundary = np.array(
        [[0, 0, boundary_z], [5, 0, boundary_z],
         [5, 5, boundary_z], [0, 5, boundary_z]],
        dtype=float,
    )
    return RoofPlane(
        normal=normal, offset=offset, inliers=inliers, boundary_3d=boundary,
    )


def test_merge_similar_planes():
    p1 = _make_plane([0, 0, 1], 10.0, n_inliers=50)
    p2 = _make_plane([0.01, 0.01, 1], 10.0, n_inliers=30)
    config = RoofDetectionConfig(merge_normal_threshold=10.0, merge_distance=5.0)

    merged = merge_planes([p1, p2], config)
    assert len(merged) == 1


def test_no_merge_different_normals():
    p1 = _make_plane([0, 0, 1], 10.0)
    p2 = _make_plane([0, 1, 0.5], 5.0)
    config = RoofDetectionConfig(merge_normal_threshold=10.0, merge_distance=5.0)

    merged = merge_planes([p1, p2], config)
    assert len(merged) == 2


def test_classify_flat_roof():
    p = _make_plane([0, 0, 1], 10.0, n_inliers=100)
    config = RoofDetectionConfig()

    result = classify_roof([p], footprint_area=100.0, total_points=100, config=config)

    assert result.roof_type == RoofType.FLAT
    assert result.confidence >= 0.5


def test_classify_gabled_roof():
    n1 = [0.0, -np.sin(np.radians(30)), np.cos(np.radians(30))]
    n2 = [0.0, np.sin(np.radians(30)), np.cos(np.radians(30))]
    p1 = _make_plane(n1, 10.0, n_inliers=100)
    p2 = _make_plane(n2, 10.0, n_inliers=100)
    config = RoofDetectionConfig()

    result = classify_roof([p1, p2], footprint_area=50.0, total_points=200, config=config)

    assert result.roof_type == RoofType.GABLED
    assert result.confidence >= 0.5
    assert result.ridge_line is not None


def test_classify_unknown_low_confidence():
    p = _make_plane([0.3, 0.3, 0.7], 5.0, n_inliers=20)
    config = RoofDetectionConfig(min_roof_confidence=0.8)

    result = classify_roof([p], footprint_area=100.0, total_points=100, config=config)

    assert result.roof_type == RoofType.UNKNOWN
    assert result.confidence < 0.8


def test_classification_result_has_all_fields():
    p = _make_plane([0, 0, 1], 10.0, n_inliers=100)
    config = RoofDetectionConfig()

    result = classify_roof([p], footprint_area=100.0, total_points=100, config=config)

    assert isinstance(result, ClassificationResult)
    assert isinstance(result.roof_type, RoofType)
    assert isinstance(result.confidence, float)
    assert isinstance(result.roof_planes, list)
    assert isinstance(result.eave_height, float)
    assert isinstance(result.ridge_height, float)
