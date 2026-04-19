from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.spatial import ConvexHull


@dataclass
class RoofDetectionConfig:
    """Configuration for the LoD2 roof detection and construction pipeline."""

    # Roof point filtering (Spec Section 1)
    facade_angle_threshold: float = 75.0       # deg
    min_roof_points: int = 20

    # Normal estimation & RANSAC (Spec Section 2)
    normal_estimation_k: int = 20
    ransac_distance_threshold: float = 0.15    # m
    ransac_iterations: int = 1000
    max_roof_planes: int = 6
    min_plane_points: int = 20
    min_remaining_ratio: float = 0.05

    # Region-growing (Spec Section 2)
    region_normal_threshold: float = 15.0      # deg
    region_neighbor_k: int = 10

    # Dominant plane merging (Spec Section 3)
    merge_normal_threshold: float = 10.0       # deg
    merge_distance: float = 1.0                # m
    min_plane_area_ratio: float = 0.1

    # Classification (Spec Section 3)
    flat_angle_threshold: float = 15.0         # deg
    slope_symmetry_tolerance: float = 10.0     # deg
    azimuth_opposition_tolerance: float = 20.0 # deg
    hip_main_area_ratio: float = 0.30
    hip_perpendicular_tolerance: float = 30.0  # deg
    hip_inward_tolerance: float = 45.0         # deg
    min_roof_confidence: float = 0.5

    # Geometry (Spec Section 4)
    edge_snap_tolerance: float = 0.01          # m
    min_convexity_ratio: float = 0.85
    min_rectangularity_ratio: float = 0.80

    # Pipeline control
    rebuild: bool = True
    fallback_to_flat: bool = True


@dataclass
class RoofPlane:
    """A detected roof plane with geometric properties."""

    normal: np.ndarray          # (3,) unit normal of the fitted plane
    offset: float               # signed distance from origin to plane
    inliers: np.ndarray         # (K,) indices into the building point cloud
    boundary_3d: np.ndarray     # (M, 3) convex hull vertices in world coordinates
    role: Optional[str] = None  # assigned during classification

    @property
    def centroid(self) -> np.ndarray:
        """Centroid of the boundary polygon in 3D world coordinates."""
        return self.boundary_3d.mean(axis=0)

    @property
    def area(self) -> float:
        """Area of the convex hull boundary in world units (m^2)."""
        if len(self.boundary_3d) < 3:
            return 0.0
        try:
            hull = ConvexHull(self.boundary_3d[:, :2])
            return hull.volume  # ConvexHull.volume is area in 2D
        except Exception:
            return 0.0

    @property
    def slope_deg(self) -> float:
        """Angle between the plane normal and the vertical (Z-axis), in degrees.
        0 = horizontal (flat), 90 = vertical (wall)."""
        cos_angle = abs(self.normal[2]) / np.linalg.norm(self.normal)
        return float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))

    @property
    def azimuth_deg(self) -> float:
        """Horizontal direction the normal points toward, in degrees.
        0=+X, 90=+Y, 180=-X, 270=-Y. Only meaningful for non-horizontal planes."""
        return float(np.degrees(np.arctan2(self.normal[1], self.normal[0])) % 360)
