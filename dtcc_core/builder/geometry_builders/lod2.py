from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass

import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import split, unary_union

from dtcc_core.logging import info

from ...model import Building, GeometryType, MultiSurface, Surface
from .buildings import build_lod1_buildings


MIN_ROOF_POINTS = 24
MIN_PLANE_INLIERS = 20
# RANSAC may inspect richer roofs; shell assembly currently accepts one or two planes.
MAX_PLANES = 6
RANSAC_ITERATIONS = 200
RANSAC_SEED = 0
RANSAC_DISTANCE_THRESHOLD = 0.2
INLIER_SHARE_RECOVERY_THRESHOLDS = (0.05, 0.10, 0.15, 0.20)
MIN_PATCH_AREA = 2.0
MAX_UNCOVERED_FOOTPRINT_FRACTION = 0.01
COPLANAR_NORMAL_TOLERANCE = 1e-3
COPLANAR_OFFSET_TOLERANCE = 0.2
NEAR_PARALLEL_NORMAL_TOLERANCE = 1e-2
EDGE_TOLERANCE = 1e-3
SHARED_VERTEX_KEY_DECIMALS = 6
RECTANGULAR_FOOTPRINT_MIN_RATIO = 0.85
HIP_PLANE_INLIER_SHARE = 0.10
NEAR_SQUARE_SIDE_RATIO = 0.75
MAX_PYRAMID_FOOTPRINT_AREA = 200.0
FOOTPRINT_DECOMPOSITION_SIMPLIFY_TOLERANCE = 1.0
MAX_DECOMPOSITION_CONCAVE_VERTICES = 2
MIN_DECOMPOSITION_REGION_AREA = MIN_PATCH_AREA
DECOMPOSITION_RECTANGULARITY_MIN_RATIO = 0.80
DECOMPOSITION_AXIS_ANGLE_TOLERANCE = 15.0
INVALID_FOOTPRINT = "invalid_footprint"
INSUFFICIENT_ROOF_POINTS = "insufficient_roof_points"
NO_PLANES_FOUND = "no_planes_found"
UNSUPPORTED_PLANE_COUNT = "unsupported_plane_count"
SHELL_ASSEMBLY_FAILED = "shell_assembly_failed"
SPLIT_FAILED = "split_failed"
PATCH_ASSIGNMENT_FAILED = "patch_assignment_failed"
PATCH_COVERAGE_FAILED = "patch_coverage_failed"
NEAR_PARALLEL_PLANES = "near_parallel_planes"
MISSING_SHARED_RIDGE = "missing_shared_ridge"
RIDGE_HEIGHT_MISMATCH = "ridge_height_mismatch"
ROOF_SURFACE_COVERAGE_FAILED = "roof_surface_coverage_failed"
WATERTIGHT_SHELL_FAILED = "watertight_shell_failed"
UNSUPPORTED_RECT_4PLANE = "unsupported_rect_4plane"
UNSUPPORTED_NEAR_SQUARE_4PLANE = "unsupported_near_square_4plane"
UNSUPPORTED_IRREGULAR_OR_OTHER = "unsupported_irregular_or_other"
UNSUPPORTED_IRREGULAR_FOOTPRINT = "unsupported_irregular_footprint"
UNSUPPORTED_RECT_FEW_DOMINANT = "unsupported_rect_few_dominant"
UNSUPPORTED_OTHER = "unsupported_other"
DECOMPOSITION_CANDIDATE = "decomposition_candidate"
DECOMPOSITION_SUCCESS = "decomposition_success"
DECOMPOSITION_UNSUPPORTED_SHAPE = "decomposition_unsupported_shape"
DECOMPOSITION_NO_VALID_SLICES = "decomposition_no_valid_slices"
DECOMPOSITION_SPARSE_REGION_POINTS = "decomposition_sparse_region_points"
DECOMPOSITION_REGION_ROOF_FAILED = "decomposition_region_roof_failed"
DECOMPOSITION_JUNCTION_FAILED = "decomposition_junction_failed"
DECOMPOSITION_COVERAGE_FAILED = "decomposition_coverage_failed"
DECOMPOSITION_WATERTIGHT_FAILED = "decomposition_watertight_failed"
NO_VALID_SLICES_AXIS_MISALIGNED = "no_valid_slices_axis_misaligned"
NO_VALID_SLICES_PIECE_TOO_SMALL = "no_valid_slices_piece_too_small"
NO_VALID_SLICES_PIECE_NOT_RECTANGULAR = "no_valid_slices_piece_not_rectangular"
NO_VALID_SLICES_COVERAGE_FAILED = "no_valid_slices_coverage_failed"
NO_VALID_SLICES_OTHER = "no_valid_slices_other"
AXIS_MISALIGNED_ANGLE_15_20 = "axis_misaligned_angle_15_20"
AXIS_MISALIGNED_ANGLE_20_25 = "axis_misaligned_angle_20_25"
AXIS_MISALIGNED_ANGLE_25_30 = "axis_misaligned_angle_25_30"
AXIS_MISALIGNED_ANGLE_30_35 = "axis_misaligned_angle_30_35"
AXIS_MISALIGNED_ANGLE_35_40 = "axis_misaligned_angle_35_40"
AXIS_MISALIGNED_ANGLE_40_45 = "axis_misaligned_angle_40_45"
AXIS_MISALIGNED_ANGLE_OVER_45 = "axis_misaligned_angle_over_45"
REGION_ROOF_NO_PLANES = "region_roof_no_planes"
REGION_ROOF_UNSUPPORTED_COUNT = "region_roof_unsupported_count"
REGION_ROOF_SPLIT_FAILED = "region_roof_split_failed"
REGION_ROOF_ASSIGNMENT_FAILED = "region_roof_assignment_failed"
REGION_ROOF_OTHER = "region_roof_other"
DECOMPOSITION_L_LIKE = "decomposition_l_like"
DECOMPOSITION_T_OR_U_LIKE = "decomposition_t_or_u_like"
DECOMPOSITION_OTHER_SHAPE = "decomposition_other_shape"
TEMPLATE_GATE_REASONS = (
    UNSUPPORTED_RECT_4PLANE,
    UNSUPPORTED_NEAR_SQUARE_4PLANE,
    UNSUPPORTED_IRREGULAR_OR_OTHER,
    UNSUPPORTED_IRREGULAR_FOOTPRINT,
    UNSUPPORTED_RECT_FEW_DOMINANT,
    UNSUPPORTED_OTHER,
)
NO_VALID_SLICE_REASONS = (
    NO_VALID_SLICES_AXIS_MISALIGNED,
    NO_VALID_SLICES_PIECE_TOO_SMALL,
    NO_VALID_SLICES_PIECE_NOT_RECTANGULAR,
    NO_VALID_SLICES_COVERAGE_FAILED,
    NO_VALID_SLICES_OTHER,
)
AXIS_MISALIGNED_ANGLE_REASONS = (
    AXIS_MISALIGNED_ANGLE_15_20,
    AXIS_MISALIGNED_ANGLE_20_25,
    AXIS_MISALIGNED_ANGLE_25_30,
    AXIS_MISALIGNED_ANGLE_30_35,
    AXIS_MISALIGNED_ANGLE_35_40,
    AXIS_MISALIGNED_ANGLE_40_45,
    AXIS_MISALIGNED_ANGLE_OVER_45,
)
REGION_ROOF_REASONS = (
    REGION_ROOF_NO_PLANES,
    REGION_ROOF_UNSUPPORTED_COUNT,
    REGION_ROOF_SPLIT_FAILED,
    REGION_ROOF_ASSIGNMENT_FAILED,
    REGION_ROOF_OTHER,
)
DECOMPOSITION_REASONS = (
    DECOMPOSITION_CANDIDATE,
    DECOMPOSITION_SUCCESS,
    DECOMPOSITION_UNSUPPORTED_SHAPE,
    DECOMPOSITION_NO_VALID_SLICES,
    DECOMPOSITION_SPARSE_REGION_POINTS,
    DECOMPOSITION_REGION_ROOF_FAILED,
    DECOMPOSITION_JUNCTION_FAILED,
    DECOMPOSITION_COVERAGE_FAILED,
    DECOMPOSITION_WATERTIGHT_FAILED,
    *NO_VALID_SLICE_REASONS,
    *AXIS_MISALIGNED_ANGLE_REASONS,
    *REGION_ROOF_REASONS,
    DECOMPOSITION_L_LIKE,
    DECOMPOSITION_T_OR_U_LIKE,
    DECOMPOSITION_OTHER_SHAPE,
)
DECOMPOSITION_FAILURE_REASONS = (
    DECOMPOSITION_UNSUPPORTED_SHAPE,
    DECOMPOSITION_NO_VALID_SLICES,
    DECOMPOSITION_SPARSE_REGION_POINTS,
    DECOMPOSITION_REGION_ROOF_FAILED,
    DECOMPOSITION_JUNCTION_FAILED,
    DECOMPOSITION_COVERAGE_FAILED,
    DECOMPOSITION_WATERTIGHT_FAILED,
)
REJECTION_REASONS = (
    INVALID_FOOTPRINT,
    INSUFFICIENT_ROOF_POINTS,
    NO_PLANES_FOUND,
    UNSUPPORTED_PLANE_COUNT,
    SPLIT_FAILED,
    PATCH_ASSIGNMENT_FAILED,
    PATCH_COVERAGE_FAILED,
    NEAR_PARALLEL_PLANES,
    MISSING_SHARED_RIDGE,
    RIDGE_HEIGHT_MISMATCH,
    ROOF_SURFACE_COVERAGE_FAILED,
    WATERTIGHT_SHELL_FAILED,
    SHELL_ASSEMBLY_FAILED,
)


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


@dataclass(frozen=True)
class FootprintDecomposition:
    pieces: list[Polygon]
    slice_lines: list[LineString]
    family_reason: str
    concave_vertex_count: int


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


def _footprint_polygon(building: Building) -> Polygon | None:
    footprint = building.lod0
    if footprint is None or len(footprint.vertices) < 3 or len(footprint.holes) > 0:
        return None
    polygon = footprint.to_polygon()
    if polygon.is_empty or not polygon.is_valid or polygon.area <= 0:
        return None
    if len(polygon.interiors) > 0:
        return None
    return polygon


def _roof_points(building: Building) -> np.ndarray | None:
    point_cloud = building.point_cloud
    if point_cloud is None or len(point_cloud.points) < MIN_ROOF_POINTS:
        return None
    return np.asarray(point_cloud.points, dtype=float)


def _surface_from_xy(poly: Polygon, plane: RoofPlane) -> Surface:
    coords = np.array(poly.exterior.coords[:-1], dtype=float)
    vertices = np.column_stack([coords[:, 0], coords[:, 1], [plane.z_at(x, y) for x, y in coords]])
    return Surface(vertices=vertices)


def _ground_surface(footprint: Polygon, ground_height: float) -> Surface:
    coords = np.array(footprint.exterior.coords[:-1], dtype=float)[::-1]
    vertices = np.column_stack([coords[:, 0], coords[:, 1], np.full(len(coords), ground_height)])
    return Surface(vertices=vertices)


def _edge_parameter(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    direction = end - start
    denom = float(np.dot(direction, direction))
    if denom == 0.0:
        return 0.0
    return float(np.dot(point - start, direction) / denom)


def _points_on_edge(points: list[np.ndarray], start: np.ndarray, end: np.ndarray, *, for_wall: bool = False) -> list[np.ndarray]:
    edge = end - start
    length = np.linalg.norm(edge)
    if length == 0:
        return []
    selected = []
    for point in points:
        point_xy = point[:2]
        offset = point_xy - start
        distance = abs(edge[0] * offset[1] - edge[1] * offset[0]) / length
        parameter = _edge_parameter(point_xy, start, end)
        if distance <= EDGE_TOLERANCE and -EDGE_TOLERANCE <= parameter <= 1.0 + EDGE_TOLERANCE:
            selected.append(point)
    selected.sort(key=lambda point: _edge_parameter(point[:2], start, end))
    unique = []
    for point in selected:
        if not unique or np.linalg.norm(point - unique[-1]) > EDGE_TOLERANCE:
            unique.append(point)
    if not for_wall:
        return unique
    grouped: list[list[np.ndarray]] = []
    for point in unique:
        if not grouped:
            grouped.append([point])
            continue
        previous_parameter = _edge_parameter(grouped[-1][0][:2], start, end)
        parameter = _edge_parameter(point[:2], start, end)
        if abs(parameter - previous_parameter) <= EDGE_TOLERANCE:
            grouped[-1].append(point)
        else:
            grouped.append([point])
    ordered: list[np.ndarray] = []
    for index, group in enumerate(grouped):
        if len(group) == 1:
            ordered.extend(group)
            continue
        group = sorted(group, key=lambda point: point[2])
        if index == 0 and len(grouped) > 1:
            target_z = float(np.mean([point[2] for point in grouped[1]]))
            ordered.append(min(group, key=lambda point: abs(point[2] - target_z)))
        elif index == len(grouped) - 1 and ordered:
            target_z = ordered[-1][2]
            closest = min(group, key=lambda point: abs(point[2] - target_z))
            remaining = [point for point in group if np.linalg.norm(point - closest) > EDGE_TOLERANCE]
            remaining.sort(key=lambda point: abs(point[2] - target_z))
            ordered.extend([closest, *remaining])
        else:
            if ordered:
                target_z = ordered[-1][2]
                closest = min(group, key=lambda point: abs(point[2] - target_z))
                remaining = [point for point in group if np.linalg.norm(point - closest) > EDGE_TOLERANCE]
                remaining.sort(key=lambda point: point[2])
                ordered.extend([closest, *remaining])
            else:
                ordered.extend(group)
    return ordered


def _wall_surfaces(footprint: Polygon, roof_surfaces: list[Surface], ground_height: float) -> list[Surface]:
    roof_vertices = [vertex for surface in roof_surfaces for vertex in surface.vertices]
    coords = np.array(footprint.exterior.coords[:-1], dtype=float)
    walls = []
    for index, start in enumerate(coords):
        end = coords[(index + 1) % len(coords)]
        top = _points_on_edge(roof_vertices, start, end, for_wall=True)
        if len(top) < 2:
            return []
        bottom_end = np.array([end[0], end[1], ground_height], dtype=float)
        bottom_start = np.array([start[0], start[1], ground_height], dtype=float)
        walls.append(Surface(vertices=np.vstack([top, bottom_end, bottom_start])))
    return walls


def _build_shell(footprint: Polygon, roof_surfaces: list[Surface], ground_height: float) -> MultiSurface | None:
    walls = _wall_surfaces(footprint, roof_surfaces, ground_height)
    if len(walls) == 0:
        return None
    shell = MultiSurface(surfaces=[*roof_surfaces, *walls, _ground_surface(footprint, ground_height)])
    if not is_watertight(shell):
        return None
    return shell


def _internal_junction_surfaces(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> list[Surface] | None:
    junctions: list[Surface] = []
    roof_vertices = [vertex for surface in roof_surfaces for vertex in surface.vertices]
    handled_lines: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for line in slice_lines:
        coords = np.asarray(line.coords, dtype=float)
        start = coords[0]
        end = coords[-1]
        line_key = tuple(sorted((
            tuple(np.round(start / EDGE_TOLERANCE).astype(int)),
            tuple(np.round(end / EDGE_TOLERANCE).astype(int)),
        )))
        if line_key in handled_lines:
            continue
        handled_lines.add(line_key)
        top = _points_on_edge(roof_vertices, start, end)
        if len(top) < 2:
            continue
        by_xy: dict[tuple[int, int], list[np.ndarray]] = {}
        for vertex in top:
            key = tuple(np.round(vertex[:2] / EDGE_TOLERANCE).astype(int))
            by_xy.setdefault(key, []).append(vertex)
        paired = [vertices for vertices in by_xy.values() if len(vertices) >= 2]
        if not paired:
            continue
        low_high = []
        for vertices in paired:
            vertices = sorted(vertices, key=lambda vertex: vertex[2])
            if abs(vertices[-1][2] - vertices[0][2]) > EDGE_TOLERANCE:
                low_high.append((vertices[0], vertices[-1]))
        if len(low_high) < 2:
            continue
        low_high.sort(key=lambda pair: _edge_parameter(pair[0][:2], start, end))
        for first, second in zip(low_high, low_high[1:]):
            first_xy = first[0][:2]
            second_xy = second[0][:2]
            if np.linalg.norm(second_xy - first_xy) <= EDGE_TOLERANCE:
                continue
            midpoint = 0.5 * (first_xy + second_xy)
            if footprint.boundary.distance(Point(midpoint)) <= EDGE_TOLERANCE:
                continue
            if not footprint.buffer(EDGE_TOLERANCE).covers(LineString([first_xy, second_xy])):
                continue
            junctions.append(Surface(vertices=np.asarray([first[0], second[0], second[1], first[1]], dtype=float)))
    return junctions


def _planes_are_coplanar(first: RoofPlane, second: RoofPlane) -> bool:
    normals_close = abs(1.0 - abs(float(np.dot(first.normal, second.normal)))) <= COPLANAR_NORMAL_TOLERANCE
    offset_close = abs(first.c - second.c) <= COPLANAR_OFFSET_TOLERANCE
    return normals_close and offset_close


def _planes_are_near_parallel(first: RoofPlane, second: RoofPlane) -> bool:
    return abs(1.0 - abs(float(np.dot(first.normal, second.normal)))) <= NEAR_PARALLEL_NORMAL_TOLERANCE


def _coplanar_pair_count(planes: list[RoofPlane]) -> int:
    count = 0
    for i in range(len(planes)):
        for j in range(i + 1, len(planes)):
            if _planes_are_coplanar(planes[i], planes[j]):
                count += 1
    return count


def _record_inlier_share_diagnostics(
    planes: list[RoofPlane],
    recovery_counts: Counter | None,
    trailing_plane_counts: Counter | None,
) -> None:
    if len(planes) <= 2:
        return
    inlier_counts = [len(plane.inliers) for plane in planes]
    if not inlier_counts:
        return
    dominant_count = max(inlier_counts)
    if dominant_count == 0:
        return
    if trailing_plane_counts is not None:
        trailing_plane_counts["trailing_planes_at_min_inliers"] += sum(
            MIN_PLANE_INLIERS <= count <= MIN_PLANE_INLIERS + 4
            for count in inlier_counts
            if count != dominant_count
        )
    if recovery_counts is None:
        return
    for threshold in INLIER_SHARE_RECOVERY_THRESHOLDS:
        kept_count = sum(count >= threshold * dominant_count for count in inlier_counts)
        if 1 <= kept_count <= 2:
            recovery_counts[threshold] += 1


def _minimum_rotated_rectangle_metrics(footprint: Polygon) -> tuple[float, float] | None:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, module="shapely.constructive")
        rectangle = footprint.minimum_rotated_rectangle
    if rectangle.is_empty or rectangle.area <= 0:
        return None
    coords = np.asarray(rectangle.exterior.coords[:-1], dtype=float)
    if coords.shape != (4, 2):
        return None
    edge_lengths = [
        float(np.linalg.norm(coords[(index + 1) % 4] - coords[index]))
        for index in range(4)
    ]
    longest = max(edge_lengths)
    shortest = min(edge_lengths)
    if longest <= 0:
        return None
    rectangularity = float(footprint.area / rectangle.area)
    side_ratio = float(shortest / longest)
    return rectangularity, side_ratio


def _signed_area(coords: np.ndarray) -> float:
    shifted = np.roll(coords, -1, axis=0)
    return 0.5 * float(np.sum(coords[:, 0] * shifted[:, 1] - shifted[:, 0] * coords[:, 1]))


def _concave_vertex_indices(footprint: Polygon) -> list[int]:
    coords = np.asarray(footprint.exterior.coords[:-1], dtype=float)
    if len(coords) < 4:
        return []
    orientation = 1.0 if _signed_area(coords) >= 0.0 else -1.0
    concave: list[int] = []
    for index, current in enumerate(coords):
        previous = coords[index - 1]
        following = coords[(index + 1) % len(coords)]
        incoming = current - previous
        outgoing = following - current
        denom = float(np.linalg.norm(incoming) * np.linalg.norm(outgoing))
        if denom == 0.0:
            continue
        cross = float(incoming[0] * outgoing[1] - incoming[1] * outgoing[0])
        turn_sine = cross / denom
        if abs(turn_sine) <= 1e-3:
            continue
        if turn_sine * orientation < 0.0:
            concave.append(index)
    return concave


def _simplified_footprint_for_decomposition(footprint: Polygon) -> Polygon | None:
    simplified = footprint.simplify(
        FOOTPRINT_DECOMPOSITION_SIMPLIFY_TOLERANCE,
        preserve_topology=True,
    )
    if simplified.geom_type != "Polygon" or simplified.is_empty or not simplified.is_valid or simplified.area <= 0:
        return None
    return simplified


def _decomposition_shape_reason(footprint: Polygon) -> str:
    simplified = _simplified_footprint_for_decomposition(footprint)
    if simplified is None:
        return DECOMPOSITION_OTHER_SHAPE
    concave_count = len(_concave_vertex_indices(simplified))
    if concave_count == 1:
        return DECOMPOSITION_L_LIKE
    if concave_count == 2:
        return DECOMPOSITION_T_OR_U_LIKE
    return DECOMPOSITION_OTHER_SHAPE


def _minimum_rotated_rectangle_axes(footprint: Polygon) -> tuple[np.ndarray, np.ndarray] | None:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, module="shapely.constructive")
        rectangle = footprint.minimum_rotated_rectangle
    coords = np.asarray(rectangle.exterior.coords[:-1], dtype=float)
    if coords.shape != (4, 2):
        return None
    edges = [coords[(index + 1) % 4] - coords[index] for index in range(4)]
    lengths = [float(np.linalg.norm(edge)) for edge in edges]
    if max(lengths) <= 0:
        return None
    u = edges[int(np.argmax(lengths))]
    u = u / np.linalg.norm(u)
    v = np.array([-u[1], u[0]], dtype=float)
    return u, v


def _angle_to_axis_degrees(vector: np.ndarray, axes: tuple[np.ndarray, np.ndarray]) -> float:
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        return 90.0
    direction = vector / norm
    best_alignment = max(abs(float(np.dot(direction, axes[0]))), abs(float(np.dot(direction, axes[1]))))
    best_alignment = min(1.0, max(-1.0, best_alignment))
    return float(np.degrees(np.arccos(best_alignment)))


def _axis_misalignment_angle_reason(angle: float) -> str:
    if angle < 20.0:
        return AXIS_MISALIGNED_ANGLE_15_20
    if angle < 25.0:
        return AXIS_MISALIGNED_ANGLE_20_25
    if angle < 30.0:
        return AXIS_MISALIGNED_ANGLE_25_30
    if angle < 35.0:
        return AXIS_MISALIGNED_ANGLE_30_35
    if angle < 40.0:
        return AXIS_MISALIGNED_ANGLE_35_40
    if angle <= 45.0:
        return AXIS_MISALIGNED_ANGLE_40_45
    return AXIS_MISALIGNED_ANGLE_OVER_45


def _max_concave_edge_angle_to_axes(footprint: Polygon, concave_indices: list[int], axes: tuple[np.ndarray, np.ndarray]) -> float:
    coords = np.asarray(footprint.exterior.coords[:-1], dtype=float)
    max_angle = 0.0
    for index in concave_indices:
        previous_edge = coords[index] - coords[index - 1]
        next_edge = coords[(index + 1) % len(coords)] - coords[index]
        max_angle = max(
            max_angle,
            _angle_to_axis_degrees(previous_edge, axes),
            _angle_to_axis_degrees(next_edge, axes),
        )
    return max_angle


def _concave_edges_align_to_axes(footprint: Polygon, concave_indices: list[int], axes: tuple[np.ndarray, np.ndarray]) -> bool:
    return _max_concave_edge_angle_to_axes(
        footprint,
        concave_indices,
        axes,
    ) <= DECOMPOSITION_AXIS_ANGLE_TOLERANCE


def _slice_line_through(point: np.ndarray, direction: np.ndarray, footprint: Polygon) -> LineString:
    xmin, ymin, xmax, ymax = footprint.bounds
    span = max(xmax - xmin, ymax - ymin) * 4.0
    start = point - direction * span
    end = point + direction * span
    return LineString([(float(start[0]), float(start[1])), (float(end[0]), float(end[1]))])


def _flatten_polygons(geometry) -> list[Polygon]:
    if geometry.geom_type == "Polygon":
        return [geometry]
    if hasattr(geometry, "geoms"):
        polygons: list[Polygon] = []
        for geom in geometry.geoms:
            polygons.extend(_flatten_polygons(geom))
        return polygons
    return []


def _split_by_slice_lines(footprint: Polygon, lines: list[LineString], *, filter_small: bool = True) -> list[Polygon]:
    pieces = [footprint]
    for line in lines:
        next_pieces: list[Polygon] = []
        for piece in pieces:
            split_result = split(piece, line)
            split_polygons = _flatten_polygons(split_result)
            if filter_small:
                split_polygons = [
                    polygon for polygon in split_polygons
                    if polygon.area >= MIN_DECOMPOSITION_REGION_AREA
                ]
            if len(split_polygons) == 0:
                next_pieces.append(piece)
            else:
                next_pieces.extend(split_polygons)
        pieces = next_pieces
    pieces.sort(key=lambda polygon: (round(polygon.bounds[0], 6), round(polygon.bounds[1], 6), round(polygon.area, 6)))
    return pieces


def _piece_rectangularity(piece: Polygon) -> float:
    metrics = _minimum_rotated_rectangle_metrics(piece)
    if metrics is None:
        return 0.0
    rectangularity, _ = metrics
    return rectangularity


def _decomposition_score(footprint: Polygon, pieces: list[Polygon]) -> tuple[int, float, float, float]:
    min_rectangularity = min(_piece_rectangularity(piece) for piece in pieces)
    uncovered_fraction = footprint.difference(unary_union(pieces)).area / footprint.area
    areas = [piece.area for piece in pieces]
    area_imbalance = max(areas) / min(areas)
    return (len(pieces), -min_rectangularity, uncovered_fraction, area_imbalance)


def _decomposition_piece_failure_reason(footprint: Polygon, pieces: list[Polygon], expected_count: int) -> str | None:
    if len(pieces) != expected_count:
        return NO_VALID_SLICES_OTHER
    if any(piece.area < MIN_DECOMPOSITION_REGION_AREA for piece in pieces):
        return NO_VALID_SLICES_PIECE_TOO_SMALL
    if not _patches_cover_footprint(footprint, pieces):
        return NO_VALID_SLICES_COVERAGE_FAILED
    if any(_piece_rectangularity(piece) < DECOMPOSITION_RECTANGULARITY_MIN_RATIO for piece in pieces):
        return NO_VALID_SLICES_PIECE_NOT_RECTANGULAR
    return None


def _most_common_no_valid_slice_reason(reasons: Counter) -> str:
    for reason in NO_VALID_SLICE_REASONS:
        if reasons[reason] == max(reasons.values()):
            return reason
    return NO_VALID_SLICES_OTHER


def _decompose_footprint(
    footprint: Polygon,
    decomposition_counts: Counter | None = None,
) -> FootprintDecomposition | None:
    simplified = _simplified_footprint_for_decomposition(footprint)
    if simplified is None:
        _record_decomposition_subreason(decomposition_counts, NO_VALID_SLICES_OTHER)
        return None
    concave_indices = _concave_vertex_indices(simplified)
    concave_count = len(concave_indices)
    if concave_count < 1 or concave_count > MAX_DECOMPOSITION_CONCAVE_VERTICES:
        _record_decomposition_subreason(decomposition_counts, NO_VALID_SLICES_OTHER)
        return None
    axes = _minimum_rotated_rectangle_axes(simplified)
    if axes is None:
        _record_decomposition_subreason(decomposition_counts, NO_VALID_SLICES_OTHER)
        return None
    if not _concave_edges_align_to_axes(simplified, concave_indices, axes):
        _record_decomposition_subreason(decomposition_counts, NO_VALID_SLICES_AXIS_MISALIGNED)
        _record_decomposition_subreason(
            decomposition_counts,
            _axis_misalignment_angle_reason(
                _max_concave_edge_angle_to_axes(simplified, concave_indices, axes)
            ),
        )
        return None
    family_reason = DECOMPOSITION_L_LIKE if concave_count == 1 else DECOMPOSITION_T_OR_U_LIKE
    expected_count = concave_count + 1
    coords = np.asarray(simplified.exterior.coords[:-1], dtype=float)
    candidate_lines: list[LineString] = []
    for index in concave_indices:
        point = coords[index]
        for direction in axes:
            candidate_lines.append(_slice_line_through(point, direction, footprint))
            candidate_lines.append(_slice_line_through(point, -direction, footprint))
    candidates: list[tuple[tuple[int, float, float, float], list[Polygon], list[LineString]]] = []
    if concave_count == 1:
        line_sets = [[line] for line in candidate_lines]
    else:
        line_sets = [
            [first, second]
            for first_index, first in enumerate(candidate_lines)
            for second in candidate_lines[first_index + 1:]
        ]
    failure_reasons = Counter()
    for lines in line_sets:
        pieces = _split_by_slice_lines(footprint, lines)
        failure_reason = _decomposition_piece_failure_reason(footprint, pieces, expected_count)
        if failure_reason == NO_VALID_SLICES_OTHER:
            raw_pieces = _split_by_slice_lines(footprint, lines, filter_small=False)
            raw_failure_reason = _decomposition_piece_failure_reason(footprint, raw_pieces, expected_count)
            if raw_failure_reason == NO_VALID_SLICES_PIECE_TOO_SMALL:
                failure_reason = raw_failure_reason
        if failure_reason is not None:
            failure_reasons[failure_reason] += 1
            continue
        candidates.append((_decomposition_score(footprint, pieces), pieces, lines))
    if not candidates:
        _record_decomposition_subreason(
            decomposition_counts,
            _most_common_no_valid_slice_reason(failure_reasons) if failure_reasons else NO_VALID_SLICES_OTHER,
        )
        return None
    candidates.sort(key=lambda item: item[0])
    _, pieces, lines = candidates[0]
    return FootprintDecomposition(
        pieces=pieces,
        slice_lines=lines,
        family_reason=family_reason,
        concave_vertex_count=concave_count,
    )


def _dominant_plane_count(planes: list[RoofPlane], inlier_share: float = HIP_PLANE_INLIER_SHARE) -> int:
    if not planes:
        return 0
    inlier_counts = np.asarray([len(plane.inliers) for plane in planes], dtype=float)
    dominant_count = float(np.max(inlier_counts))
    if dominant_count <= 0:
        return 0
    return int(np.sum(inlier_counts >= inlier_share * dominant_count))


def _record_template_gate_diagnostics(
    footprint: Polygon,
    planes: list[RoofPlane],
    template_gate_counts: Counter | None,
) -> None:
    if template_gate_counts is None or len(planes) <= 2:
        return
    metrics = _minimum_rotated_rectangle_metrics(footprint)
    if metrics is None:
        template_gate_counts[UNSUPPORTED_IRREGULAR_OR_OTHER] += 1
        template_gate_counts[UNSUPPORTED_OTHER] += 1
        return
    rectangularity, side_ratio = metrics
    if rectangularity < RECTANGULAR_FOOTPRINT_MIN_RATIO:
        template_gate_counts[UNSUPPORTED_IRREGULAR_OR_OTHER] += 1
        template_gate_counts[UNSUPPORTED_IRREGULAR_FOOTPRINT] += 1
        return
    dominant_plane_count = _dominant_plane_count(planes)
    if dominant_plane_count < 4:
        template_gate_counts[UNSUPPORTED_IRREGULAR_OR_OTHER] += 1
        template_gate_counts[UNSUPPORTED_RECT_FEW_DOMINANT] += 1
        return
    template_gate_counts[UNSUPPORTED_RECT_4PLANE] += 1
    if side_ratio >= NEAR_SQUARE_SIDE_RATIO and footprint.area <= MAX_PYRAMID_FOOTPRINT_AREA:
        template_gate_counts[UNSUPPORTED_NEAR_SQUARE_4PLANE] += 1


def _is_irregular_footprint(footprint: Polygon) -> bool:
    metrics = _minimum_rotated_rectangle_metrics(footprint)
    if metrics is None:
        return False
    rectangularity, _ = metrics
    return rectangularity < RECTANGULAR_FOOTPRINT_MIN_RATIO


def _record_decomposition_candidate(decomposition_counts: Counter | None, family_reason: str) -> None:
    if decomposition_counts is None:
        return
    decomposition_counts[DECOMPOSITION_CANDIDATE] += 1
    decomposition_counts[family_reason] += 1


def _record_decomposition_failure(decomposition_counts: Counter | None, reason: str) -> None:
    if decomposition_counts is not None:
        if reason not in DECOMPOSITION_FAILURE_REASONS:
            reason = DECOMPOSITION_REGION_ROOF_FAILED
        decomposition_counts[reason] += 1


def _record_decomposition_subreason(decomposition_counts: Counter | None, reason: str) -> None:
    if decomposition_counts is not None:
        decomposition_counts[reason] += 1


def _patches_cover_footprint(footprint: Polygon, patches: list[Polygon]) -> bool:
    covered = unary_union(patches)
    missing = footprint.difference(covered)
    return missing.area / footprint.area <= MAX_UNCOVERED_FOOTPRINT_FRACTION


def _roof_surfaces_cover_footprint(footprint: Polygon, roof_surfaces: list[Surface]) -> bool:
    projected_patches = []
    for surface in roof_surfaces:
        if len(surface.vertices) < 3:
            continue
        patch = Polygon(surface.vertices[:, :2])
        if not patch.is_empty and patch.is_valid and patch.area > 0:
            projected_patches.append(patch)
    if not projected_patches:
        return False
    return _patches_cover_footprint(footprint, projected_patches)


def _plane_equality_line(first: RoofPlane, second: RoofPlane, footprint: Polygon) -> LineString | None:
    a = first.a - second.a
    b = first.b - second.b
    c = first.c - second.c
    if abs(a) + abs(b) < 1e-12:
        return None
    xmin, ymin, xmax, ymax = footprint.bounds
    span = max(xmax - xmin, ymax - ymin) * 4.0
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    if abs(b) >= abs(a):
        x0 = cx - span
        x1 = cx + span
        y0 = -(a * x0 + c) / b
        y1 = -(a * x1 + c) / b
    else:
        y0 = cy - span
        y1 = cy + span
        x0 = -(b * y0 + c) / a
        x1 = -(b * y1 + c) / a
    return LineString([(x0, y0), (x1, y1)])


def _split_footprint_for_two_planes(footprint: Polygon, first: RoofPlane, second: RoofPlane) -> list[Polygon] | None:
    line = _plane_equality_line(first, second, footprint)
    if line is None:
        return None
    pieces = split(footprint, line)
    polygons = [geom for geom in pieces.geoms if geom.area >= MIN_PATCH_AREA]
    if len(polygons) != 2:
        return None
    return polygons


def _assign_patches_to_planes(points: np.ndarray, planes: list[RoofPlane], patches: list[Polygon]) -> list[Polygon] | None:
    assigned: list[Polygon | None] = [None] * len(planes)
    used_patch_indices: set[int] = set()
    for plane_index, plane in enumerate(planes):
        best_patch_index = -1
        best_count = -1
        for patch_index, patch in enumerate(patches):
            if patch_index in used_patch_indices:
                continue
            buffered = patch.buffer(EDGE_TOLERANCE)
            count = sum(
                buffered.contains(Point(x, y)) or buffered.touches(Point(x, y))
                for x, y in points[plane.inliers][:, :2]
            )
            if count > best_count:
                best_count = count
                best_patch_index = patch_index
        if best_patch_index < 0 or best_count == 0:
            return None
        assigned[plane_index] = patches[best_patch_index]
        used_patch_indices.add(best_patch_index)
    return [patch for patch in assigned if patch is not None]


def _ridge_xy(first: Polygon, second: Polygon) -> list[tuple[float, float]]:
    shared = first.boundary.intersection(second.boundary)
    if shared.is_empty:
        return []
    if shared.geom_type == "LineString":
        return list(shared.coords)
    if shared.geom_type == "MultiLineString":
        longest = max(shared.geoms, key=lambda geom: geom.length)
        return list(longest.coords)
    return []


def _surface_from_patch_with_shared_edges(
    patch: Polygon,
    plane: RoofPlane,
    shared_edges: dict[tuple[float, float], np.ndarray],
) -> Surface:
    vertices = []
    for x, y in patch.exterior.coords[:-1]:
        key = (
            round(float(x), SHARED_VERTEX_KEY_DECIMALS),
            round(float(y), SHARED_VERTEX_KEY_DECIMALS),
        )
        if key in shared_edges:
            vertices.append(shared_edges[key])
        else:
            vertices.append(np.array([x, y, plane.z_at(x, y)], dtype=float))
    return Surface(vertices=np.asarray(vertices, dtype=float))


def _multi_plane_roof_surfaces(
    points: np.ndarray,
    footprint: Polygon,
    planes: list[RoofPlane],
) -> tuple[list[Surface] | None, str | None]:
    if len(planes) != 2:
        return None, SHELL_ASSEMBLY_FAILED
    split_patches = _split_footprint_for_two_planes(footprint, planes[0], planes[1])
    if split_patches is None:
        return None, SPLIT_FAILED
    patches = _assign_patches_to_planes(points, planes, split_patches)
    if patches is None:
        return None, PATCH_ASSIGNMENT_FAILED
    kept_planes = planes
    if len(patches) < 2 or not _patches_cover_footprint(footprint, patches):
        return None, PATCH_COVERAGE_FAILED

    shared_vertices_by_patch: list[dict[tuple[float, float], np.ndarray]] = [dict() for _ in patches]
    for i in range(len(patches)):
        for j in range(i + 1, len(patches)):
            if _planes_are_coplanar(kept_planes[i], kept_planes[j]):
                continue
            if _planes_are_near_parallel(kept_planes[i], kept_planes[j]):
                return None, NEAR_PARALLEL_PLANES
            coords = _ridge_xy(patches[i], patches[j])
            if len(coords) < 2:
                return None, MISSING_SHARED_RIDGE
            for x, y in coords:
                z_i = kept_planes[i].z_at(x, y)
                z_j = kept_planes[j].z_at(x, y)
                if abs(z_i - z_j) > RANSAC_DISTANCE_THRESHOLD:
                    return None, RIDGE_HEIGHT_MISMATCH
                vertex = np.array([x, y, 0.5 * (z_i + z_j)], dtype=float)
                key = (
                    round(float(x), SHARED_VERTEX_KEY_DECIMALS),
                    round(float(y), SHARED_VERTEX_KEY_DECIMALS),
                )
                shared_vertices_by_patch[i][key] = vertex
                shared_vertices_by_patch[j][key] = vertex

    roof_surfaces = [
        _surface_from_patch_with_shared_edges(patch, plane, shared)
        for patch, plane, shared in zip(patches, kept_planes, shared_vertices_by_patch)
    ]
    if not _roof_surfaces_cover_footprint(footprint, roof_surfaces):
        return None, ROOF_SURFACE_COVERAGE_FAILED
    return roof_surfaces, None


def _record_rejection(rejections: Counter | None, reason: str) -> None:
    if rejections is not None:
        rejections[reason] += 1


def _points_in_polygon(points: np.ndarray, polygon: Polygon) -> np.ndarray:
    buffered = polygon.buffer(EDGE_TOLERANCE)
    selected = [
        point for point in points
        if buffered.contains(Point(point[0], point[1])) or buffered.touches(Point(point[0], point[1]))
    ]
    return np.asarray(selected, dtype=float)


def _region_roof_failure_reason(reason: str | None) -> str:
    if reason == NO_PLANES_FOUND:
        return REGION_ROOF_NO_PLANES
    if reason == UNSUPPORTED_PLANE_COUNT:
        return REGION_ROOF_UNSUPPORTED_COUNT
    if reason == SPLIT_FAILED:
        return REGION_ROOF_SPLIT_FAILED
    if reason == PATCH_ASSIGNMENT_FAILED:
        return REGION_ROOF_ASSIGNMENT_FAILED
    return REGION_ROOF_OTHER


def _region_roof_surfaces(
    points: np.ndarray,
    region: Polygon,
    decomposition_counts: Counter | None = None,
) -> tuple[list[Surface] | None, str | None]:
    region_points = _points_in_polygon(points, region)
    if len(region_points) < MIN_ROOF_POINTS:
        return None, DECOMPOSITION_SPARSE_REGION_POINTS
    planes = _ransac_planes(region_points)
    roof_surfaces, reason = _roof_surfaces_for_planes(region_points, region, planes)
    if roof_surfaces is None:
        _record_decomposition_subreason(
            decomposition_counts,
            _region_roof_failure_reason(reason),
        )
        return None, reason or DECOMPOSITION_REGION_ROOF_FAILED
    return roof_surfaces, None


def _roof_surfaces_for_planes(
    points: np.ndarray,
    footprint: Polygon,
    planes: list[RoofPlane],
) -> tuple[list[Surface] | None, str | None]:
    if len(planes) == 0:
        return None, NO_PLANES_FOUND
    if len(planes) > 2:
        return None, UNSUPPORTED_PLANE_COUNT
    if len(planes) == 1:
        return [_surface_from_xy(footprint, planes[0])], None
    return _multi_plane_roof_surfaces(points, footprint, planes)


def _build_decomposed_shell(
    footprint: Polygon,
    roof_points: np.ndarray,
    ground_height: float,
    decomposition_counts: Counter | None = None,
) -> tuple[MultiSurface | None, str]:
    decomposition = _decompose_footprint(footprint, decomposition_counts)
    if decomposition is None:
        return None, DECOMPOSITION_NO_VALID_SLICES
    roof_surfaces: list[Surface] = []
    for region in decomposition.pieces:
        region_surfaces, reason = _region_roof_surfaces(roof_points, region, decomposition_counts)
        if region_surfaces is None:
            return None, reason or DECOMPOSITION_REGION_ROOF_FAILED
        roof_surfaces.extend(region_surfaces)
    if not _roof_surfaces_cover_footprint(footprint, roof_surfaces):
        return None, DECOMPOSITION_COVERAGE_FAILED
    junctions = _internal_junction_surfaces(footprint, roof_surfaces, decomposition.slice_lines)
    if junctions is None:
        return None, DECOMPOSITION_JUNCTION_FAILED
    walls = _wall_surfaces(footprint, roof_surfaces, ground_height)
    if len(walls) == 0:
        return None, DECOMPOSITION_WATERTIGHT_FAILED
    shell = MultiSurface(surfaces=[*roof_surfaces, *junctions, *walls, _ground_surface(footprint, ground_height)])
    if not is_watertight(shell):
        return None, DECOMPOSITION_WATERTIGHT_FAILED
    if decomposition_counts is not None:
        decomposition_counts[DECOMPOSITION_SUCCESS] += 1
    return shell, DECOMPOSITION_SUCCESS


def _candidate_lod2_from_parts(
    footprint: Polygon,
    roof_points: np.ndarray,
    ground_height: float,
    rejections: Counter | None = None,
    plane_counts: Counter | None = None,
    max_plane_coplanar_pair_counts: Counter | None = None,
    inlier_share_recovery_counts: Counter | None = None,
    trailing_plane_counts: Counter | None = None,
    template_gate_counts: Counter | None = None,
    decomposition_counts: Counter | None = None,
) -> MultiSurface | None:
    planes = _ransac_planes(roof_points)
    if plane_counts is not None:
        plane_counts[len(planes)] += 1
    if len(planes) == MAX_PLANES and max_plane_coplanar_pair_counts is not None:
        max_plane_coplanar_pair_counts[_coplanar_pair_count(planes)] += 1
    _record_inlier_share_diagnostics(
        planes,
        inlier_share_recovery_counts,
        trailing_plane_counts,
    )
    if len(planes) > 2:
        _record_template_gate_diagnostics(footprint, planes, template_gate_counts)
        if _is_irregular_footprint(footprint):
            family_reason = _decomposition_shape_reason(footprint)
            _record_decomposition_candidate(decomposition_counts, family_reason)
            if family_reason == DECOMPOSITION_OTHER_SHAPE:
                _record_decomposition_failure(decomposition_counts, DECOMPOSITION_UNSUPPORTED_SHAPE)
            else:
                shell, decomposition_reason = _build_decomposed_shell(
                    footprint,
                    roof_points,
                    ground_height,
                    decomposition_counts,
                )
                if shell is not None:
                    return shell
                _record_decomposition_failure(decomposition_counts, decomposition_reason)
        _record_rejection(rejections, UNSUPPORTED_PLANE_COUNT)
        return None

    roof_surfaces, rejection_reason = _roof_surfaces_for_planes(roof_points, footprint, planes)
    if roof_surfaces is None:
        if rejection_reason == SPLIT_FAILED:
            family_reason = _decomposition_shape_reason(footprint)
            _record_decomposition_candidate(decomposition_counts, family_reason)
            if family_reason == DECOMPOSITION_OTHER_SHAPE:
                _record_decomposition_failure(decomposition_counts, DECOMPOSITION_UNSUPPORTED_SHAPE)
            else:
                shell, decomposition_reason = _build_decomposed_shell(
                    footprint,
                    roof_points,
                    ground_height,
                    decomposition_counts,
                )
                if shell is not None:
                    return shell
                _record_decomposition_failure(decomposition_counts, decomposition_reason)
        _record_rejection(rejections, rejection_reason or SHELL_ASSEMBLY_FAILED)
        return None
    shell = _build_shell(footprint, roof_surfaces, float(ground_height))
    if shell is None:
        _record_rejection(rejections, WATERTIGHT_SHELL_FAILED)
    return shell


def _log_lod2_summary(
    *,
    total: int,
    lod2_count: int,
    fallback_count: int,
    skipped_existing: int,
    rejections: Counter,
    plane_counts: Counter,
    max_plane_coplanar_pair_counts: Counter,
    inlier_share_recovery_counts: Counter,
    trailing_plane_counts: Counter,
    template_gate_counts: Counter,
    decomposition_counts: Counter,
) -> None:
    info(
        "LOD2 build summary: "
        f"total={total} lod2={lod2_count} fallback={fallback_count} "
        f"skipped_existing={skipped_existing}"
    )
    plane_count_parts = [
        f"planes_{count}={plane_counts[count]}"
        for count in sorted(plane_counts)
        if plane_counts[count] > 0
    ]
    if plane_count_parts:
        info("LOD2 plane count summary: " + " ".join(plane_count_parts))
    coplanar_pair_parts = [
        f"coplanar_pairs_{count}={max_plane_coplanar_pair_counts[count]}"
        for count in sorted(max_plane_coplanar_pair_counts)
        if max_plane_coplanar_pair_counts[count] > 0
    ]
    if coplanar_pair_parts:
        info("LOD2 max-plane coplanar pair summary: " + " ".join(coplanar_pair_parts))
    recovery_parts = [
        f"would_recover_at_{int(threshold * 100)}pct={inlier_share_recovery_counts[threshold]}"
        for threshold in INLIER_SHARE_RECOVERY_THRESHOLDS
        if inlier_share_recovery_counts[threshold] > 0
    ]
    if recovery_parts:
        info("LOD2 inlier-share recovery simulation: " + " ".join(recovery_parts))
    trailing_parts = [
        f"{reason}={trailing_plane_counts[reason]}"
        for reason in sorted(trailing_plane_counts)
        if trailing_plane_counts[reason] > 0
    ]
    if trailing_parts:
        info("LOD2 trailing plane summary: " + " ".join(trailing_parts))
    template_gate_parts = [
        f"{reason}={template_gate_counts[reason]}"
        for reason in TEMPLATE_GATE_REASONS
        if template_gate_counts[reason] > 0
    ]
    if template_gate_parts:
        info("LOD2 template gate summary: " + " ".join(template_gate_parts))
    decomposition_parts = [
        f"{reason}={decomposition_counts[reason]}"
        for reason in DECOMPOSITION_REASONS
        if decomposition_counts[reason] > 0
    ]
    if decomposition_parts:
        info("LOD2 decomposition summary: " + " ".join(decomposition_parts))
    rejection_parts = [
        f"{reason}={rejections[reason]}"
        for reason in REJECTION_REASONS
        if rejections[reason] > 0
    ]
    if rejection_parts:
        info("LOD2 rejection summary: " + " ".join(rejection_parts))


def _candidate_lod2(
    building: Building,
    default_ground_height: float,
    always_use_default_ground: bool,
    rejections: Counter | None = None,
    plane_counts: Counter | None = None,
    max_plane_coplanar_pair_counts: Counter | None = None,
    inlier_share_recovery_counts: Counter | None = None,
    trailing_plane_counts: Counter | None = None,
    template_gate_counts: Counter | None = None,
    decomposition_counts: Counter | None = None,
) -> MultiSurface | None:
    footprint = _footprint_polygon(building)
    if footprint is None:
        _record_rejection(rejections, INVALID_FOOTPRINT)
        return None
    roof_points = _roof_points(building)
    if roof_points is None:
        _record_rejection(rejections, INSUFFICIENT_ROOF_POINTS)
        return None
    ground_height = default_ground_height if always_use_default_ground else building.attributes.get("ground_height", default_ground_height)
    return _candidate_lod2_from_parts(
        footprint,
        roof_points,
        float(ground_height),
        rejections,
        plane_counts,
        max_plane_coplanar_pair_counts,
        inlier_share_recovery_counts,
        trailing_plane_counts,
        template_gate_counts,
        decomposition_counts,
    )


def build_lod2_buildings(
    buildings: list[Building],
    *,
    default_ground_height: float = 0.0,
    always_use_default_ground: bool = False,
    rebuild: bool = True,
    build_lod1_fallback: bool = True,
    log_rejections: bool = False,
) -> list[Building]:
    rejections = Counter()
    plane_counts = Counter()
    max_plane_coplanar_pair_counts = Counter()
    inlier_share_recovery_counts = Counter()
    trailing_plane_counts = Counter()
    template_gate_counts = Counter()
    decomposition_counts = Counter()
    lod2_count = 0
    fallback_count = 0
    skipped_existing = 0

    for building in buildings:
        if building.lod2 is not None and not rebuild:
            skipped_existing += 1
            continue
        if rebuild:
            building.remove_geometry(GeometryType.LOD2)
        candidate = _candidate_lod2(
            building,
            default_ground_height,
            always_use_default_ground,
            rejections if log_rejections else None,
            plane_counts if log_rejections else None,
            max_plane_coplanar_pair_counts if log_rejections else None,
            inlier_share_recovery_counts if log_rejections else None,
            trailing_plane_counts if log_rejections else None,
            template_gate_counts if log_rejections else None,
            decomposition_counts if log_rejections else None,
        )
        if candidate is not None:
            building.add_geometry(candidate, GeometryType.LOD2)
            lod2_count += 1
        else:
            fallback_count += 1
            if build_lod1_fallback and building.lod1 is None:
                build_lod1_buildings(
                    [building],
                    default_ground_height=default_ground_height,
                    always_use_default_ground=always_use_default_ground,
                    rebuild=False,
                )
    if log_rejections:
        _log_lod2_summary(
            total=len(buildings),
            lod2_count=lod2_count,
            fallback_count=fallback_count,
            skipped_existing=skipped_existing,
            rejections=rejections,
            plane_counts=plane_counts,
            max_plane_coplanar_pair_counts=max_plane_coplanar_pair_counts,
            inlier_share_recovery_counts=inlier_share_recovery_counts,
            trailing_plane_counts=trailing_plane_counts,
            template_gate_counts=template_gate_counts,
            decomposition_counts=decomposition_counts,
        )
    return buildings


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
