import json

import numpy as np
import pytest

from dtcc_core.model.enums import RoofType, SurfaceSemantic
from dtcc_core.builder.geometry_builders.roof_config import (
    RoofDetectionConfig,
    RoofPlane,
)


def test_roof_type_enum_values():
    assert RoofType.FLAT.name == "FLAT"
    assert RoofType.GABLED.name == "GABLED"
    assert RoofType.HIPPED.name == "HIPPED"
    assert RoofType.UNKNOWN.name == "UNKNOWN"


def test_surface_semantic_enum_values():
    assert SurfaceSemantic.GROUND.value == 0
    assert SurfaceSemantic.WALL.value == 1
    assert SurfaceSemantic.ROOF.value == 2


def test_roof_type_json_serializable():
    """Enum names must be JSON-serializable strings for building.attributes."""
    attrs = {"roof_type": RoofType.GABLED.name}
    serialized = json.dumps(attrs)
    deserialized = json.loads(serialized)
    assert deserialized["roof_type"] == "GABLED"
    assert RoofType[deserialized["roof_type"]] == RoofType.GABLED


def test_config_defaults():
    config = RoofDetectionConfig()
    assert config.ransac_distance_threshold == 0.15
    assert config.min_plane_points == 20
    assert config.max_roof_planes == 6
    assert config.flat_angle_threshold == 15.0
    assert config.min_roof_confidence == 0.5
    assert config.min_convexity_ratio == 0.85
    assert config.min_rectangularity_ratio == 0.80


def test_config_override():
    config = RoofDetectionConfig(ransac_distance_threshold=0.25, min_plane_points=30)
    assert config.ransac_distance_threshold == 0.25
    assert config.min_plane_points == 30
    assert config.flat_angle_threshold == 15.0  # unchanged default


def test_roof_plane_creation():
    normal = np.array([0.0, 0.0, 1.0])
    offset = 10.0
    inliers = np.array([0, 1, 2, 3, 4])
    boundary = np.array(
        [[0, 0, 10], [5, 0, 10], [5, 5, 10], [0, 5, 10]], dtype=float
    )
    plane = RoofPlane(
        normal=normal,
        offset=offset,
        inliers=inliers,
        boundary_3d=boundary,
    )
    assert plane.normal is normal
    assert plane.offset == 10.0
    assert len(plane.inliers) == 5
    assert plane.centroid is not None
    assert plane.area > 0
    assert plane.slope_deg == pytest.approx(0.0, abs=0.1)
    assert plane.role is None


def test_roof_plane_slope_and_azimuth():
    # 45-degree slope facing south (negative Y normal component)
    normal = np.array([0.0, -np.sin(np.radians(45)), np.cos(np.radians(45))])
    normal = normal / np.linalg.norm(normal)
    plane = RoofPlane(
        normal=normal,
        offset=5.0,
        inliers=np.array([0, 1, 2]),
        boundary_3d=np.array([[0, 0, 5], [5, 0, 5], [2.5, 5, 5]], dtype=float),
    )
    assert plane.slope_deg == pytest.approx(45.0, abs=1.0)
    # normal points toward -Y, azimuth should be ~270 deg
    assert 260 < plane.azimuth_deg < 280
