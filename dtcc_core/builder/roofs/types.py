"""Ephemeral roof evidence, planar patches and JSON-ready fit diagnostics."""
from dataclasses import dataclass, field

import numpy as np
from shapely.geometry import Polygon


@dataclass(frozen=True)
class RoofFitOptions:
    """Robust fitting thresholds in metres; these are not confidence probabilities.

    ``shapes`` selects the ordered candidates. Slopes are bounded to 60 degrees.
    Substantial faces (at least 10% of the area) require observed inliers.
    """
    shapes: tuple[str, ...] = ("flat", "shed", "gable", "hip")
    min_points: int = 12
    median_error: float = .35
    p90_error: float = 1.5
    face_inlier_distance: float = .6
    minimum_clearance: float = .25
    max_nfev: int = 100

    def __post_init__(self):
        from .catalogue import REGISTRY
        if not self.shapes or len(set(self.shapes)) != len(self.shapes):
            raise ValueError("shapes must contain distinct supported roof keys")
        unknown = set(self.shapes) - REGISTRY.keys()
        if unknown:
            raise ValueError(f"Unknown roof shapes: {sorted(unknown)}; supported: {tuple(REGISTRY)}")
        if self.min_points < 12 or self.max_nfev < 1:
            raise ValueError("min_points must be >= 12 and max_nfev must be positive")
        values = [self.median_error, self.p90_error, self.face_inlier_distance, self.minimum_clearance]
        if not all(np.isfinite(v) and v > 0 for v in values):
            raise ValueError("Roof fit thresholds must be positive and finite")


@dataclass
class RoofPatch:
    """A polygon (including holes) and z=a(x-ox)+b(y-oy)+c.

    Ownership identifies a source cell; patch_id identifies a plane piece, and
    block_id the cleaned output. These identities are independent of mesh markers.
    """
    polygon: Polygon
    origin: np.ndarray
    coefficients: np.ndarray
    source_key: str
    patch_id: int = 0
    block_id: int = 0

    def height(self, xy):
        return (np.asarray(xy)[..., :2] - self.origin) @ self.coefficients[:2] + self.coefficients[2]

    def native_plane(self):
        return [*map(float, self.origin), *map(float, self.coefficients)]


@dataclass
class RoofModel:
    """Source roof represented by a catalogue family and its lower-envelope planes."""
    source_key: str
    shape: str
    parameters: np.ndarray
    origin: np.ndarray
    planes: np.ndarray
    diagnostics: dict = field(default_factory=dict)

    def predict(self, xy):
        local = np.asarray(xy)[..., :2] - self.origin
        return np.min(local @ self.planes[:, :2].T + self.planes[:, 2], axis=-1)

    def patches(self, target, *, block_id=0):
        from .partition import planar_patches
        return planar_patches(self, target, block_id=block_id)


@dataclass
class RoofFitResult:
    """Accepted roof model or rejection, together with fitting diagnostics."""
    model: RoofModel | None
    diagnostics: dict

    @property
    def accepted(self):
        return bool(self.diagnostics.get("accepted", False))
