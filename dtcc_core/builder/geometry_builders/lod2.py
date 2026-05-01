from __future__ import annotations

from collections import Counter

import numpy as np

from ...model import Building, GeometryType, MultiSurface, Surface


MIN_ROOF_POINTS = 24
MIN_PLANE_INLIERS = 8
MAX_PLANES = 6
RANSAC_ITERATIONS = 200
RANSAC_SEED = 0
RANSAC_DISTANCE_THRESHOLD = 0.2
MIN_PATCH_AREA = 2.0
MAX_UNCOVERED_FOOTPRINT_FRACTION = 0.01
COPLANAR_NORMAL_TOLERANCE = 1e-3
COPLANAR_OFFSET_TOLERANCE = 0.2
NEAR_PARALLEL_NORMAL_TOLERANCE = 1e-2
EDGE_TOLERANCE = 1e-3


class RoofPlane:
    def __init__(self, a: float, b: float, c: float, inliers: np.ndarray):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.inliers = np.asarray(inliers, dtype=int)

    @property
    def normal(self) -> np.ndarray:
        normal = np.array([-self.a, -self.b, 1.0], dtype=float)
        return normal / np.linalg.norm(normal)

    def z_at(self, x: float, y: float) -> float:
        return self.a * float(x) + self.b * float(y) + self.c

    def distances(self, points: np.ndarray) -> np.ndarray:
        predicted = self.a * points[:, 0] + self.b * points[:, 1] + self.c
        return np.abs(points[:, 2] - predicted) / np.linalg.norm([-self.a, -self.b, 1.0])


def _fit_plane(points: np.ndarray, inliers: np.ndarray | None = None) -> RoofPlane:
    if inliers is None:
        sample = np.asarray(points, dtype=float)
        inlier_indices = np.arange(len(sample))
    else:
        inlier_indices = np.asarray(inliers, dtype=int)
        sample = np.asarray(points, dtype=float)[inlier_indices]
    matrix = np.column_stack([sample[:, 0], sample[:, 1], np.ones(len(sample))])
    a, b, c = np.linalg.lstsq(matrix, sample[:, 2], rcond=None)[0]
    return RoofPlane(a, b, c, inlier_indices)


def _fit_plane_from_three(points: np.ndarray, indices: np.ndarray) -> RoofPlane | None:
    sample = points[indices]
    v1 = sample[1] - sample[0]
    v2 = sample[2] - sample[0]
    normal = np.cross(v1, v2)
    if abs(normal[2]) < 1e-12:
        return None
    a = -normal[0] / normal[2]
    b = -normal[1] / normal[2]
    c = sample[0, 2] - a * sample[0, 0] - b * sample[0, 1]
    return RoofPlane(a, b, c, np.array([], dtype=int))


def _ransac_planes(
    points: np.ndarray,
    *,
    seed: int = RANSAC_SEED,
    max_planes: int = MAX_PLANES,
    iterations: int = RANSAC_ITERATIONS,
    distance_threshold: float = RANSAC_DISTANCE_THRESHOLD,
    min_inliers: int = MIN_PLANE_INLIERS,
) -> list[RoofPlane]:
    points = np.asarray(points, dtype=float)
    remaining = np.arange(len(points))
    planes: list[RoofPlane] = []
    rng = np.random.default_rng(seed)

    while len(remaining) >= min_inliers and len(planes) < max_planes:
        best_inliers = np.array([], dtype=int)
        for _ in range(iterations):
            sample_indices = rng.choice(remaining, size=3, replace=False)
            candidate = _fit_plane_from_three(points, sample_indices)
            if candidate is None:
                continue
            distances = candidate.distances(points[remaining])
            inliers = remaining[distances <= distance_threshold]
            if len(inliers) > len(best_inliers):
                best_inliers = inliers
        if len(best_inliers) < min_inliers:
            break
        refined = _fit_plane(points, best_inliers)
        planes.append(refined)
        remaining_set = set(best_inliers)
        remaining = np.array([idx for idx in remaining if idx not in remaining_set], dtype=int)

    return planes


def _vertex_key(vertex: np.ndarray, tolerance: float = EDGE_TOLERANCE) -> tuple[int, int, int]:
    return tuple(np.round(np.asarray(vertex, dtype=float) / tolerance).astype(int))


def _surface_edges(surface: Surface, tolerance: float = EDGE_TOLERANCE):
    vertices = surface.vertices
    if len(vertices) < 3:
        return
    keys = [_vertex_key(vertex, tolerance) for vertex in vertices]
    for index, start in enumerate(keys):
        end = keys[(index + 1) % len(keys)]
        yield tuple(sorted((start, end)))


def _edge_counts(multisurface: MultiSurface, tolerance: float = EDGE_TOLERANCE) -> Counter:
    counts = Counter()
    for surface in multisurface.surfaces:
        counts.update(_surface_edges(surface, tolerance))
    return counts


def is_watertight(multisurface: MultiSurface, tolerance: float = EDGE_TOLERANCE) -> bool:
    if not isinstance(multisurface, MultiSurface) or len(multisurface.surfaces) == 0:
        return False
    counts = _edge_counts(multisurface, tolerance)
    return len(counts) > 0 and all(count == 2 for count in counts.values())
