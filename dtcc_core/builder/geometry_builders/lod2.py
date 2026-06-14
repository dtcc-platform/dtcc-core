from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass

import numpy as np
from shapely.geometry import LineString, MultiPoint, Point, Polygon, box
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
BOUNDARY_ENDPOINT_Z_WELD_TOLERANCE = 0.02
SHARED_VERTEX_KEY_DECIMALS = 6
RECTANGULAR_FOOTPRINT_MIN_RATIO = 0.85
HIP_PLANE_INLIER_SHARE = 0.10
MAX_TEMPLATE_DROPPED_INLIER_SHARE = 0.30
NEAR_SQUARE_SIDE_RATIO = 0.75
MAX_PYRAMID_FOOTPRINT_AREA = 200.0
FOOTPRINT_DECOMPOSITION_SIMPLIFY_TOLERANCE = 1.0
MAX_DECOMPOSITION_CONCAVE_VERTICES = 2
MIN_DECOMPOSITION_REGION_AREA = MIN_PATCH_AREA
DECOMPOSITION_RECTANGULARITY_MIN_RATIO = 0.80
DECOMPOSITION_AXIS_ANGLE_TOLERANCE = 15.0
STEPPED_FLAT_MAX_SLOPE = 0.12
STEPPED_FLAT_MIN_LEVELS = 2
STEPPED_FLAT_MIN_LEVEL_INLIERS = MIN_PLANE_INLIERS
STEPPED_FLAT_MIN_PATCH_SUPPORT_SHARE = 0.70
STEPPED_FLAT_MAX_PATCH_CONTAMINATION_SHARE = 0.10
STEPPED_FLAT_REGION_MAX_MULTIPART = 1
STEPPED_FLAT_REGION_GRID_TARGET_CELLS = 24
STEPPED_FLAT_REGION_GRID_MIN_SIZE = 1.0
STEPPED_FLAT_REGION_MAX_PATCH_SUPPORT_AREA_RATIO = 1.5
FLAT_COLLAPSE_MIN_DOMINANT_SHARE = 0.90
FLAT_COLLAPSE_BAND_TOLERANCE = 2.0 * RANSAC_DISTANCE_THRESHOLD
FLAT_COLLAPSE_MAX_DOMINANT_RMSE = RANSAC_DISTANCE_THRESHOLD
STEPPED_FLAT_CANDIDATE = "stepped_flat_candidate"
STEPPED_FLAT_SUCCESS = "stepped_flat_success"
STEPPED_FLAT_FLAT_COLLAPSE_SUCCESS = "stepped_flat_flat_collapse_success"
STEPPED_FLAT_SLOPED_MAJOR_EVIDENCE = "stepped_flat_sloped_major_evidence"
STEPPED_FLAT_INSUFFICIENT_LEVELS = "stepped_flat_insufficient_levels"
STEPPED_FLAT_DROPPED_TOO_MANY = "stepped_flat_dropped_too_many"
STEPPED_FLAT_PATCH_FAILED = "stepped_flat_patch_failed"
STEPPED_FLAT_STEP_WALL_FAILED = "stepped_flat_step_wall_failed"
STEPPED_FLAT_WATERTIGHT_FAILED = "stepped_flat_watertight_failed"
STEPPED_FLAT_REASONS = (
    STEPPED_FLAT_CANDIDATE,
    STEPPED_FLAT_SUCCESS,
    STEPPED_FLAT_FLAT_COLLAPSE_SUCCESS,
    STEPPED_FLAT_SLOPED_MAJOR_EVIDENCE,
    STEPPED_FLAT_INSUFFICIENT_LEVELS,
    STEPPED_FLAT_DROPPED_TOO_MANY,
    STEPPED_FLAT_PATCH_FAILED,
    STEPPED_FLAT_STEP_WALL_FAILED,
    STEPPED_FLAT_WATERTIGHT_FAILED,
)
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
WATERTIGHT_TOO_FEW_PAIRED_VERTICES = "watertight_too_few_paired_vertices"
WATERTIGHT_UNPAIRED_RIDGE_ON_SLICE = "watertight_unpaired_ridge_on_slice"
WATERTIGHT_EDGE_COUNT_MISMATCH = "watertight_edge_count_mismatch"
WATERTIGHT_OTHER = "watertight_other"
EDGE_COUNT_UNMATCHED_ON_SLICE = "edge_count_unmatched_on_slice"
EDGE_COUNT_UNMATCHED_ON_FOOTPRINT_EXTERIOR = "edge_count_unmatched_on_footprint_exterior"
EDGE_COUNT_UNMATCHED_INTERIOR = "edge_count_unmatched_interior"
EDGE_COUNT_EXCESS_COUNT = "edge_count_excess_count"
EDGE_COUNT_OTHER = "edge_count_other"
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
REGION_ROOF_DROPPED_TOO_MANY = "region_roof_dropped_too_many"
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
    REGION_ROOF_DROPPED_TOO_MANY,
    REGION_ROOF_OTHER,
)
WATERTIGHT_FAILURE_REASONS = (
    WATERTIGHT_TOO_FEW_PAIRED_VERTICES,
    WATERTIGHT_UNPAIRED_RIDGE_ON_SLICE,
    WATERTIGHT_EDGE_COUNT_MISMATCH,
    WATERTIGHT_OTHER,
)
EDGE_COUNT_MISMATCH_REASONS = (
    EDGE_COUNT_UNMATCHED_ON_SLICE,
    EDGE_COUNT_UNMATCHED_ON_FOOTPRINT_EXTERIOR,
    EDGE_COUNT_UNMATCHED_INTERIOR,
    EDGE_COUNT_EXCESS_COUNT,
    EDGE_COUNT_OTHER,
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
    *WATERTIGHT_FAILURE_REASONS,
    *EDGE_COUNT_MISMATCH_REASONS,
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


@dataclass(frozen=True)
class SteppedFlatLevel:
    plane: RoofPlane
    point_indices: np.ndarray


@dataclass(frozen=True)
class AssignedSteppedFlatRegion:
    polygon: Polygon
    level_index: int


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


def _surface_is_simple(surface: Surface) -> bool:
    vertices = np.asarray(surface.vertices, dtype=float)
    if len(vertices) < 3:
        return False
    centered = vertices - vertices.mean(axis=0)
    _, _, rotation = np.linalg.svd(centered)
    projected = centered @ rotation[:2].T
    polygon = Polygon(projected)
    return polygon.is_valid and polygon.area > EDGE_TOLERANCE * EDGE_TOLERANCE


def _surfaces_are_simple(surfaces: list[Surface]) -> bool:
    return all(_surface_is_simple(surface) for surface in surfaces)


def _internal_junction_surfaces(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> list[Surface] | None:
    junctions: list[Surface] = []
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
        group_xys = _slice_line_group_xys(roof_surfaces, line)
        if len(group_xys) < 2:
            continue
        for first_xy, second_xy in zip(group_xys, group_xys[1:]):
            if np.linalg.norm(second_xy - first_xy) <= EDGE_TOLERANCE:
                continue
            midpoint = 0.5 * (first_xy + second_xy)
            if footprint.boundary.distance(Point(midpoint)) <= EDGE_TOLERANCE:
                continue
            if not footprint.buffer(EDGE_TOLERANCE).covers(LineString([first_xy, second_xy])):
                continue
            edges = _slice_line_surface_edges(roof_surfaces, first_xy, second_xy)
            if len(edges) < 2:
                continue
            first_edge, second_edge = edges[:2]
            start_gap = first_edge[0][2] - second_edge[0][2]
            end_gap = first_edge[1][2] - second_edge[1][2]
            if abs(start_gap) <= EDGE_TOLERANCE and abs(end_gap) <= EDGE_TOLERANCE:
                # Roof edges coincide and pair with each other directly.
                continue
            if _crossing_parameter(first_edge, second_edge) is not None:
                # Roof planes still cross inside this segment; a quad here
                # would self-intersect. Crossings are split beforehand by
                # _split_crossing_junction_edges, so refuse the shell rather
                # than emit invalid geometry.
                return None
            if abs(start_gap) <= EDGE_TOLERANCE:
                corners = [first_edge[0], first_edge[1], second_edge[1]]
            elif abs(end_gap) <= EDGE_TOLERANCE:
                corners = [first_edge[0], first_edge[1], second_edge[0]]
            else:
                corners = [first_edge[0], first_edge[1], second_edge[1], second_edge[0]]
            junctions.append(Surface(vertices=np.asarray(corners, dtype=float)))
    return junctions


def _slice_line_group_xys(
    roof_surfaces: list[Surface],
    line: LineString,
) -> list[np.ndarray]:
    coords = np.asarray(line.coords, dtype=float)
    start = coords[0]
    end = coords[-1]
    top, _, _ = _slice_line_vertex_groups(roof_surfaces, line)
    by_xy: dict[tuple[int, int], np.ndarray] = {}
    for vertex in top:
        key = tuple(np.round(vertex[:2] / EDGE_TOLERANCE).astype(int))
        by_xy.setdefault(key, vertex[:2])
    group_xys = list(by_xy.values())
    group_xys.sort(key=lambda xy: _edge_parameter(xy, start, end))
    return group_xys


def _crossing_parameter(
    first_edge: tuple[np.ndarray, np.ndarray],
    second_edge: tuple[np.ndarray, np.ndarray],
) -> float | None:
    start_gap = first_edge[0][2] - second_edge[0][2]
    end_gap = first_edge[1][2] - second_edge[1][2]
    if abs(start_gap) <= EDGE_TOLERANCE or abs(end_gap) <= EDGE_TOLERANCE:
        return None
    if start_gap * end_gap >= 0.0:
        return None
    return start_gap / (start_gap - end_gap)


def _snap_crossing_vertex_heights(
    first_surface: Surface,
    second_surface: Surface,
    xy: np.ndarray,
) -> None:
    indices = []
    for surface in (first_surface, second_surface):
        index = next(
            index
            for index, vertex in enumerate(surface.vertices)
            if np.linalg.norm(np.asarray(vertex[:2], dtype=float) - xy) <= EDGE_TOLERANCE
        )
        indices.append(index)
    mean_z = 0.5 * (
        first_surface.vertices[indices[0]][2] + second_surface.vertices[indices[1]][2]
    )
    first_surface.vertices[indices[0]][2] = mean_z
    second_surface.vertices[indices[1]][2] = mean_z


def _split_crossing_junction_edges(
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> list[Surface]:
    """Insert a shared vertex where adjacent region roof planes cross.

    Two nearly coincident roof planes can swap height order along a slice
    line. A single junction quad between them would self-intersect, so the
    crossing point is inserted into both roof edges; the junction builder
    then emits two simple triangles instead. When the crossing falls too
    close to an existing vertex to insert into both edges, the geometry is
    left untouched and the junction builder rejects the shell.
    """
    surfaces = list(roof_surfaces)
    for line in slice_lines:
        group_xys = _slice_line_group_xys(surfaces, line)
        for first_xy, second_xy in zip(group_xys, group_xys[1:]):
            if np.linalg.norm(second_xy - first_xy) <= EDGE_TOLERANCE:
                continue
            owners = [
                index
                for index, surface in enumerate(surfaces)
                if _surface_edge_between_xys(surface, first_xy, second_xy) is not None
            ]
            if len(owners) < 2:
                continue
            first_owner, second_owner = owners[:2]
            parameter = _crossing_parameter(
                _surface_edge_between_xys(surfaces[first_owner], first_xy, second_xy),
                _surface_edge_between_xys(surfaces[second_owner], first_xy, second_xy),
            )
            if parameter is None:
                continue
            crossing_xy = first_xy + parameter * (second_xy - first_xy)
            first_split = _insert_xy_on_surface_edge(surfaces[first_owner], crossing_xy)
            second_split = _insert_xy_on_surface_edge(surfaces[second_owner], crossing_xy)
            if first_split is surfaces[first_owner] or second_split is surfaces[second_owner]:
                continue
            _snap_crossing_vertex_heights(first_split, second_split, crossing_xy)
            surfaces[first_owner] = first_split
            surfaces[second_owner] = second_split
    return surfaces


def _slice_line_vertex_groups(
    roof_surfaces: list[Surface],
    line: LineString,
) -> tuple[list[np.ndarray], list[list[np.ndarray]], list[tuple[np.ndarray, np.ndarray]]]:
    coords = np.asarray(line.coords, dtype=float)
    start = coords[0]
    end = coords[-1]
    roof_vertices = [vertex for surface in roof_surfaces for vertex in surface.vertices]
    top = _points_on_edge(roof_vertices, start, end)
    by_xy: dict[tuple[int, int], list[np.ndarray]] = {}
    for vertex in top:
        key = tuple(np.round(vertex[:2] / EDGE_TOLERANCE).astype(int))
        by_xy.setdefault(key, []).append(vertex)
    unpaired = [vertices for vertices in by_xy.values() if len(vertices) == 1]
    low_high = []
    for vertices in by_xy.values():
        if len(vertices) < 2:
            continue
        vertices = sorted(vertices, key=lambda vertex: vertex[2])
        if abs(vertices[-1][2] - vertices[0][2]) > EDGE_TOLERANCE:
            low_high.append((vertices[0], vertices[-1]))
    return top, unpaired, low_high


def _decomposed_watertight_failure_reason(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
    shell: MultiSurface | None,
) -> str:
    for line in slice_lines:
        top, unpaired, low_high = _slice_line_vertex_groups(roof_surfaces, line)
        for vertices in unpaired:
            point = Point(vertices[0][:2])
            if footprint.boundary.distance(point) > EDGE_TOLERANCE:
                return WATERTIGHT_UNPAIRED_RIDGE_ON_SLICE
        if len(top) >= 2 and len(low_high) < 2:
            return WATERTIGHT_TOO_FEW_PAIRED_VERTICES
    if shell is not None:
        edge_counts = _edge_counts(shell)
        if edge_counts and any(count != 2 for count in edge_counts.values()):
            return WATERTIGHT_EDGE_COUNT_MISMATCH
    return WATERTIGHT_OTHER


def _edge_xy_midpoint(edge: tuple[tuple[int, int, int], tuple[int, int, int]]) -> np.ndarray:
    start, end = edge
    return np.array(
        [
            0.5 * (start[0] + end[0]) * EDGE_TOLERANCE,
            0.5 * (start[1] + end[1]) * EDGE_TOLERANCE,
        ],
        dtype=float,
    )


def _edge_count_mismatch_reason(
    footprint: Polygon,
    slice_lines: list[LineString],
    shell: MultiSurface,
) -> str:
    bad_edges = [
        (edge, count)
        for edge, count in _edge_counts(shell).items()
        if count != 2
    ]
    if any(count > 2 for _, count in bad_edges):
        return EDGE_COUNT_EXCESS_COUNT
    for edge, _ in bad_edges:
        midpoint = Point(_edge_xy_midpoint(edge))
        if any(line.distance(midpoint) <= EDGE_TOLERANCE for line in slice_lines):
            return EDGE_COUNT_UNMATCHED_ON_SLICE
    for edge, _ in bad_edges:
        midpoint = Point(_edge_xy_midpoint(edge))
        if footprint.boundary.distance(midpoint) <= EDGE_TOLERANCE:
            return EDGE_COUNT_UNMATCHED_ON_FOOTPRINT_EXTERIOR
    for edge, _ in bad_edges:
        midpoint = Point(_edge_xy_midpoint(edge))
        if footprint.contains(midpoint):
            return EDGE_COUNT_UNMATCHED_INTERIOR
    return EDGE_COUNT_OTHER


def _xy_close(first: np.ndarray, second: np.ndarray) -> bool:
    return np.linalg.norm(first - second) <= EDGE_TOLERANCE


def _surface_edge_between_xys(
    surface: Surface,
    first_xy: np.ndarray,
    second_xy: np.ndarray,
) -> tuple[np.ndarray, np.ndarray] | None:
    vertices = np.asarray(surface.vertices, dtype=float)
    for index, start in enumerate(vertices):
        end = vertices[(index + 1) % len(vertices)]
        if _xy_close(start[:2], first_xy) and _xy_close(end[:2], second_xy):
            return start, end
        if _xy_close(start[:2], second_xy) and _xy_close(end[:2], first_xy):
            return end, start
    return None


def _slice_line_surface_edges(
    roof_surfaces: list[Surface],
    first_xy: np.ndarray,
    second_xy: np.ndarray,
) -> list[tuple[np.ndarray, np.ndarray]]:
    edges = []
    for surface in roof_surfaces:
        edge = _surface_edge_between_xys(surface, first_xy, second_xy)
        if edge is not None:
            edges.append(edge)
    return edges


def _surface_has_xy(surface: Surface, xy: np.ndarray) -> bool:
    return any(np.linalg.norm(vertex[:2] - xy) <= EDGE_TOLERANCE for vertex in surface.vertices)


def _edge_contains_xy(xy: np.ndarray, start: np.ndarray, end: np.ndarray) -> bool:
    edge = end[:2] - start[:2]
    length = float(np.linalg.norm(edge))
    if length == 0.0:
        return False
    offset = xy - start[:2]
    distance = abs(edge[0] * offset[1] - edge[1] * offset[0]) / length
    parameter = _edge_parameter(xy, start[:2], end[:2])
    return distance <= EDGE_TOLERANCE and EDGE_TOLERANCE < parameter < 1.0 - EDGE_TOLERANCE


def _insert_xy_on_surface_edge(surface: Surface, xy: np.ndarray) -> Surface:
    vertices = np.asarray(surface.vertices, dtype=float)
    if _surface_has_xy(surface, xy):
        return surface
    for index, start in enumerate(vertices):
        end = vertices[(index + 1) % len(vertices)]
        if not _edge_contains_xy(xy, start, end):
            continue
        parameter = _edge_parameter(xy, start[:2], end[:2])
        z = start[2] + parameter * (end[2] - start[2])
        vertex = np.array([xy[0], xy[1], z], dtype=float)
        new_vertices = np.vstack([vertices[: index + 1], vertex, vertices[index + 1 :]])
        return Surface(vertices=new_vertices)
    return surface


def _surface_has_vertex(surface: Surface, vertex: np.ndarray) -> bool:
    return any(np.linalg.norm(existing - vertex) <= EDGE_TOLERANCE for existing in surface.vertices)


def _insert_vertex_on_surface_edge(surface: Surface, vertex: np.ndarray) -> Surface:
    vertices = np.asarray(surface.vertices, dtype=float)
    if _surface_has_vertex(surface, vertex):
        return surface
    xy = vertex[:2]
    for index, start in enumerate(vertices):
        end = vertices[(index + 1) % len(vertices)]
        if not _edge_contains_xy(xy, start, end):
            continue
        parameter = _edge_parameter(xy, start[:2], end[:2])
        z = start[2] + parameter * (end[2] - start[2])
        if abs(z - vertex[2]) > EDGE_TOLERANCE:
            continue
        inserted = np.array([xy[0], xy[1], z], dtype=float)
        new_vertices = np.vstack([vertices[: index + 1], inserted, vertices[index + 1 :]])
        return Surface(vertices=new_vertices)
    return surface


def _reconcile_slice_line_vertices(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> list[Surface]:
    reconciled = [
        Surface(vertices=np.asarray(surface.vertices, dtype=float).copy())
        for surface in roof_surfaces
    ]
    for line in slice_lines:
        _, unpaired, _ = _slice_line_vertex_groups(reconciled, line)
        for vertices in unpaired:
            xy = vertices[0][:2]
            if footprint.boundary.distance(Point(xy)) <= EDGE_TOLERANCE:
                continue
            reconciled = [_insert_xy_on_surface_edge(surface, xy) for surface in reconciled]
    return reconciled


def _weld_near_equal_boundary_endpoint_heights(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> list[Surface]:
    welded = [
        Surface(vertices=np.asarray(surface.vertices, dtype=float).copy())
        for surface in roof_surfaces
    ]
    endpoint_xys: dict[tuple[int, int], np.ndarray] = {}
    for line in slice_lines:
        top, _, _ = _slice_line_vertex_groups(welded, line)
        for vertex in top:
            xy = vertex[:2]
            if footprint.boundary.distance(Point(xy)) > EDGE_TOLERANCE:
                continue
            key = tuple(np.round(xy / EDGE_TOLERANCE).astype(int))
            endpoint_xys[key] = xy
    for xy in endpoint_xys.values():
        vertices: list[tuple[Surface, int, float]] = []
        for surface in welded:
            for index, vertex in enumerate(surface.vertices):
                if np.linalg.norm(vertex[:2] - xy) <= EDGE_TOLERANCE:
                    vertices.append((surface, index, float(vertex[2])))
        vertices.sort(key=lambda item: item[2])
        cluster: list[tuple[Surface, int, float]] = []
        for item in vertices:
            if cluster and item[2] - cluster[0][2] > BOUNDARY_ENDPOINT_Z_WELD_TOLERANCE:
                if len(cluster) > 1:
                    mean_z = float(np.mean([entry[2] for entry in cluster]))
                    for surface, index, _ in cluster:
                        surface.vertices[index][2] = mean_z
                cluster = []
            cluster.append(item)
        if len(cluster) > 1:
            mean_z = float(np.mean([entry[2] for entry in cluster]))
            for surface, index, _ in cluster:
                surface.vertices[index][2] = mean_z
    return welded


def _reconcile_shell_edge_vertices(surfaces: list[Surface]) -> list[Surface]:
    reconciled = [
        Surface(vertices=np.asarray(surface.vertices, dtype=float).copy())
        for surface in surfaces
    ]
    vertices: list[np.ndarray] = []
    seen: set[tuple[int, int, int]] = set()
    for surface in reconciled:
        for vertex in surface.vertices:
            vertex = np.asarray(vertex, dtype=float)
            key = tuple(np.round(vertex / EDGE_TOLERANCE).astype(int))
            if key in seen:
                continue
            seen.add(key)
            vertices.append(vertex)
    for vertex in vertices:
        reconciled = [_insert_vertex_on_surface_edge(surface, vertex) for surface in reconciled]
    return reconciled


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


def _record_stepped_flat_reason(stepped_flat_counts: Counter | None, reason: str) -> None:
    if stepped_flat_counts is not None:
        stepped_flat_counts[reason] += 1


def _stepped_flat_evidence_indices(planes: list[RoofPlane]) -> np.ndarray | None:
    if not planes:
        return None
    all_inliers: set[int] = set()
    near_horizontal_inliers: set[int] = set()
    for plane in planes:
        plane_inliers = {int(index) for index in plane.inliers}
        all_inliers.update(plane_inliers)
        if float(np.hypot(plane.a, plane.b)) <= STEPPED_FLAT_MAX_SLOPE:
            near_horizontal_inliers.update(plane_inliers)
    if not all_inliers or not near_horizontal_inliers:
        return None
    sloped_inliers = all_inliers - near_horizontal_inliers
    if len(sloped_inliers) / len(all_inliers) > MAX_TEMPLATE_DROPPED_INLIER_SHARE:
        return None
    return np.asarray(sorted(near_horizontal_inliers), dtype=int)


def _major_planes_are_near_horizontal(planes: list[RoofPlane]) -> bool:
    return _stepped_flat_evidence_indices(planes) is not None


def _height_cluster_indices(points: np.ndarray) -> list[np.ndarray]:
    if len(points) == 0:
        return []
    order = np.argsort(points[:, 2])
    clusters: list[list[int]] = [[int(order[0])]]
    previous_z = float(points[order[0], 2])
    for point_index in order[1:]:
        z = float(points[point_index, 2])
        if abs(z - previous_z) > RANSAC_DISTANCE_THRESHOLD:
            clusters.append([int(point_index)])
        else:
            clusters[-1].append(int(point_index))
        previous_z = z
    return [np.asarray(cluster, dtype=int) for cluster in clusters]


def _stepped_flat_height_levels(
    points: np.ndarray,
    planes: list[RoofPlane],
) -> tuple[list[SteppedFlatLevel] | None, str | None]:
    points = np.asarray(points, dtype=float)
    evidence_indices = _stepped_flat_evidence_indices(planes)
    if evidence_indices is None:
        return None, STEPPED_FLAT_SLOPED_MAJOR_EVIDENCE
    evidence_indices = evidence_indices[
        (evidence_indices >= 0) & (evidence_indices < len(points))
    ]
    evidence_index_set = {int(index) for index in evidence_indices}
    clusters = _height_cluster_indices(points)
    if not clusters:
        return None, STEPPED_FLAT_INSUFFICIENT_LEVELS
    strongest = max(len(cluster) for cluster in clusters)
    minimum_level_count = max(
        STEPPED_FLAT_MIN_LEVEL_INLIERS,
        int(np.ceil(HIP_PLANE_INLIER_SHARE * strongest)),
    )
    kept: list[np.ndarray] = [
        cluster for cluster in clusters if len(cluster) >= minimum_level_count
    ]
    if len(kept) < STEPPED_FLAT_MIN_LEVELS:
        return None, STEPPED_FLAT_INSUFFICIENT_LEVELS
    kept_evidence_indices = {
        int(index)
        for cluster in kept
        for index in cluster
        if int(index) in evidence_index_set
    }
    dropped_count = len(evidence_index_set - kept_evidence_indices)
    if (
        len(evidence_index_set) > 0
        and dropped_count / len(evidence_index_set) > MAX_TEMPLATE_DROPPED_INLIER_SHARE
    ):
        return None, STEPPED_FLAT_DROPPED_TOO_MANY
    levels = []
    for cluster in kept:
        level_z = float(np.median(points[cluster, 2]))
        levels.append(
            SteppedFlatLevel(
                plane=RoofPlane(0.0, 0.0, level_z, cluster),
                point_indices=cluster,
            )
        )
    levels.sort(key=lambda level: level.plane.c)
    return levels, None


def _dominant_height_band(heights: np.ndarray, tolerance: float) -> tuple[float, np.ndarray]:
    order = np.argsort(heights)
    sorted_heights = heights[order]
    best_start = 0
    best_end = 0
    start = 0
    for end, height in enumerate(sorted_heights):
        while height - sorted_heights[start] > 2.0 * tolerance:
            start += 1
        if end - start > best_end - best_start:
            best_start = start
            best_end = end
    band_indices = order[best_start : best_end + 1]
    level_z = float(np.median(heights[band_indices]))
    band_mask = np.abs(heights - level_z) <= tolerance
    if np.any(band_mask):
        level_z = float(np.median(heights[band_mask]))
        band_mask = np.abs(heights - level_z) <= tolerance
    return level_z, band_mask


def _flat_collapse_plane(points: np.ndarray, planes: list[RoofPlane]) -> RoofPlane | None:
    points = np.asarray(points, dtype=float)
    evidence_indices = _stepped_flat_evidence_indices(planes)
    if evidence_indices is None:
        return None
    evidence_indices = evidence_indices[
        (evidence_indices >= 0) & (evidence_indices < len(points))
    ]
    if len(evidence_indices) == 0:
        return None
    heights = points[evidence_indices, 2]
    level_z, band_mask = _dominant_height_band(
        heights,
        FLAT_COLLAPSE_BAND_TOLERANCE,
    )
    if np.count_nonzero(band_mask) / len(heights) < FLAT_COLLAPSE_MIN_DOMINANT_SHARE:
        return None
    residuals = heights[band_mask] - level_z
    if len(residuals) == 0:
        return None
    rmse = float(np.sqrt(np.mean(residuals * residuals)))
    if rmse > FLAT_COLLAPSE_MAX_DOMINANT_RMSE:
        return None
    return RoofPlane(0.0, 0.0, level_z, evidence_indices[band_mask])


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


def _cap_region_roof_planes(planes: list[RoofPlane]) -> tuple[list[RoofPlane], str | None]:
    if len(planes) <= 2:
        return planes, None
    planes_by_inliers = sorted(planes, key=lambda plane: len(plane.inliers), reverse=True)
    total_inliers = sum(len(plane.inliers) for plane in planes_by_inliers)
    dropped_inliers = sum(len(plane.inliers) for plane in planes_by_inliers[2:])
    if total_inliers > 0 and dropped_inliers / total_inliers > MAX_TEMPLATE_DROPPED_INLIER_SHARE:
        return [], REGION_ROOF_DROPPED_TOO_MANY
    return planes_by_inliers[:2], None


def _region_roof_surfaces(
    points: np.ndarray,
    region: Polygon,
    decomposition_counts: Counter | None = None,
) -> tuple[list[Surface] | None, str | None]:
    region_points = _points_in_polygon(points, region)
    if len(region_points) < MIN_ROOF_POINTS:
        return None, DECOMPOSITION_SPARSE_REGION_POINTS
    planes = _ransac_planes(region_points)
    planes, cap_reason = _cap_region_roof_planes(planes)
    if cap_reason is not None:
        _record_decomposition_subreason(decomposition_counts, cap_reason)
        return None, cap_reason
    roof_surfaces, reason = _roof_surfaces_for_planes(region_points, region, planes)
    if roof_surfaces is None:
        _record_decomposition_subreason(
            decomposition_counts,
            _region_roof_failure_reason(reason),
        )
        return None, reason or DECOMPOSITION_REGION_ROOF_FAILED
    return roof_surfaces, None


def _project_xy_to_axes(
    xy: np.ndarray,
    axes: tuple[np.ndarray, np.ndarray],
) -> np.ndarray:
    xy = np.asarray(xy, dtype=float)
    return np.column_stack((xy @ axes[0], xy @ axes[1]))


def _local_rect_to_world_polygon(
    strip_axis: int,
    strip_min: float,
    strip_max: float,
    cross_min: float,
    cross_max: float,
    axes: tuple[np.ndarray, np.ndarray],
) -> Polygon:
    if strip_axis == 0:
        local_coords = [
            (strip_min, cross_min),
            (strip_max, cross_min),
            (strip_max, cross_max),
            (strip_min, cross_max),
        ]
    else:
        local_coords = [
            (cross_min, strip_min),
            (cross_max, strip_min),
            (cross_max, strip_max),
            (cross_min, strip_max),
        ]
    coords = [
        tuple(float(value) for value in local[0] * axes[0] + local[1] * axes[1])
        for local in local_coords
    ]
    return Polygon(coords)


def _single_polygon_or_none(geometry) -> Polygon | None:
    polygons = [
        polygon
        for polygon in _flatten_polygons(geometry)
        if not polygon.is_empty and polygon.is_valid and polygon.area >= MIN_PATCH_AREA
    ]
    if len(polygons) != 1:
        return None
    return polygons[0]


def _single_grid_cell_polygon_or_none(geometry) -> Polygon | None:
    polygons = [
        polygon
        for polygon in _flatten_polygons(geometry)
        if not polygon.is_empty and polygon.is_valid and polygon.area > EDGE_TOLERANCE
    ]
    if len(polygons) != 1:
        return None
    return polygons[0]


def _level_patch_support_is_clean(
    points: np.ndarray,
    patch: Polygon,
    level: SteppedFlatLevel,
    all_level_indices: set[int],
) -> bool:
    level_indices = set(int(index) for index in level.point_indices)
    if not level_indices:
        return False
    other_level_indices = all_level_indices - level_indices
    buffered = patch.buffer(EDGE_TOLERANCE)
    support_count = 0
    contamination_count = 0
    for point_index in all_level_indices:
        point = points[point_index]
        point_xy = Point(point[0], point[1])
        if not (buffered.contains(point_xy) or buffered.touches(point_xy)):
            continue
        if point_index in level_indices:
            support_count += 1
        elif point_index in other_level_indices:
            contamination_count += 1
    if support_count / len(level_indices) < STEPPED_FLAT_MIN_PATCH_SUPPORT_SHARE:
        return False
    accepted_patch_count = support_count + contamination_count
    if accepted_patch_count == 0:
        return False
    if not _region_area_has_level_support_extent(
        points,
        patch,
        level,
        _region_grid_cell_size(patch),
    ):
        return False
    return contamination_count / accepted_patch_count <= STEPPED_FLAT_MAX_PATCH_CONTAMINATION_SHARE


def _level_point_counts_in_region(
    points: np.ndarray,
    region: Polygon,
    levels: list[SteppedFlatLevel],
) -> list[int]:
    buffered = region.buffer(EDGE_TOLERANCE)
    counts: list[int] = []
    for level in levels:
        count = 0
        for point_index in level.point_indices:
            point = points[int(point_index)]
            point_xy = Point(point[0], point[1])
            if buffered.contains(point_xy) or buffered.touches(point_xy):
                count += 1
        counts.append(count)
    return counts


def _assign_candidate_regions_to_levels(
    points: np.ndarray,
    footprint: Polygon,
    candidate_regions: list[Polygon],
    levels: list[SteppedFlatLevel],
) -> list[Polygon] | None:
    if len(candidate_regions) < len(levels):
        return None
    points = np.asarray(points, dtype=float)
    assignments: list[AssignedSteppedFlatRegion] = []
    for region in candidate_regions:
        if region.is_empty or not region.is_valid or region.area < MIN_PATCH_AREA:
            return None
        counts = _level_point_counts_in_region(points, region, levels)
        total = sum(counts)
        if total == 0:
            return None
        level_index = int(np.argmax(counts))
        contamination = (total - counts[level_index]) / total
        if contamination > STEPPED_FLAT_MAX_PATCH_CONTAMINATION_SHARE:
            return None
        assignments.append(AssignedSteppedFlatRegion(region, level_index))

    all_level_indices = {
        int(index)
        for level in levels
        for index in level.point_indices
    }
    patches_by_level: list[Polygon] = []
    for level_index, level in enumerate(levels):
        regions = [
            assignment.polygon
            for assignment in assignments
            if assignment.level_index == level_index
        ]
        if not regions:
            return None
        merged = unary_union(regions)
        polygons = [
            polygon
            for polygon in _flatten_polygons(merged)
            if not polygon.is_empty and polygon.is_valid and polygon.area >= MIN_PATCH_AREA
        ]
        if len(polygons) != STEPPED_FLAT_REGION_MAX_MULTIPART:
            return None
        patch = _single_polygon_or_none(polygons[0].intersection(footprint))
        if patch is None:
            return None
        if len(patch.interiors) > 0:
            return None
        if not _level_patch_support_is_clean(points, patch, level, all_level_indices):
            return None
        patches_by_level.append(patch)

    covered = unary_union(patches_by_level)
    overlap_area = sum(patch.area for patch in patches_by_level) - covered.area
    if overlap_area / footprint.area > MAX_UNCOVERED_FOOTPRINT_FRACTION:
        return None
    if not _patches_cover_footprint(footprint, patches_by_level):
        return None
    return patches_by_level


def _decomposition_candidate_regions(footprint: Polygon) -> list[Polygon] | None:
    decomposition = _decompose_footprint(footprint)
    if decomposition is None:
        return None
    regions = [
        region
        for region in decomposition.pieces
        if not region.is_empty and region.is_valid and region.area >= MIN_PATCH_AREA
    ]
    if len(regions) != len(decomposition.pieces):
        return None
    if not _patches_cover_footprint(footprint, regions):
        return None
    return regions


def _region_grid_cell_size(footprint: Polygon) -> float:
    minx, miny, maxx, maxy = footprint.bounds
    span = max(maxx - minx, maxy - miny)
    if span <= 0.0:
        return STEPPED_FLAT_REGION_GRID_MIN_SIZE
    return max(
        STEPPED_FLAT_REGION_GRID_MIN_SIZE,
        span / STEPPED_FLAT_REGION_GRID_TARGET_CELLS,
    )


def _nearest_level_index_to_region(
    points: np.ndarray,
    region: Polygon,
    levels: list[SteppedFlatLevel],
) -> int | None:
    representative = region.representative_point()
    xy = np.array([representative.x, representative.y], dtype=float)
    distances: list[float] = []
    for level in levels:
        if len(level.point_indices) == 0:
            return None
        level_xy = points[level.point_indices, :2]
        distances.append(float(np.min(np.linalg.norm(level_xy - xy, axis=1))))
    if not distances:
        return None
    return int(np.argmin(distances))


def _level_support_extent_area(
    points: np.ndarray,
    level: SteppedFlatLevel,
    buffer_distance: float,
) -> float | None:
    if len(level.point_indices) == 0:
        return None
    level_xy = points[level.point_indices, :2]
    support = MultiPoint([tuple(xy) for xy in level_xy]).convex_hull
    extent = support.buffer(buffer_distance)
    if extent.is_empty or not extent.is_valid or extent.area <= EDGE_TOLERANCE:
        return None
    return float(extent.area)


def _region_area_has_level_support_extent(
    points: np.ndarray,
    region: Polygon,
    level: SteppedFlatLevel,
    buffer_distance: float,
) -> bool:
    extent_area = _level_support_extent_area(points, level, buffer_distance)
    if extent_area is None:
        return False
    return (
        region.area / extent_area
        <= STEPPED_FLAT_REGION_MAX_PATCH_SUPPORT_AREA_RATIO
    )


def _grid_candidate_regions(
    points: np.ndarray,
    footprint: Polygon,
    levels: list[SteppedFlatLevel],
) -> list[Polygon] | None:
    cell_size = _region_grid_cell_size(footprint)
    minx, miny, maxx, maxy = footprint.bounds
    x_values = np.arange(minx, maxx, cell_size)
    y_values = np.arange(miny, maxy, cell_size)
    if len(x_values) == 0 or len(y_values) == 0:
        return None

    cells_by_level: dict[int, list[Polygon]] = {
        index: []
        for index in range(len(levels))
    }
    for x in x_values:
        for y in y_values:
            cell = box(
                float(x),
                float(y),
                float(min(x + cell_size, maxx)),
                float(min(y + cell_size, maxy)),
            )
            clipped = _single_grid_cell_polygon_or_none(cell.intersection(footprint))
            if clipped is None:
                continue
            counts = _level_point_counts_in_region(
                points,
                clipped.buffer(0.5 * cell_size),
                levels,
            )
            total = sum(counts)
            level_index = None
            if total > 0:
                candidate_level_index = int(np.argmax(counts))
                contamination = (total - counts[candidate_level_index]) / total
                if contamination <= STEPPED_FLAT_MAX_PATCH_CONTAMINATION_SHARE:
                    level_index = candidate_level_index
            if level_index is None:
                level_index = _nearest_level_index_to_region(points, clipped, levels)
            if level_index is None:
                continue
            cells_by_level[level_index].append(clipped)

    candidate_regions: list[Polygon] = []
    for level_index, regions in cells_by_level.items():
        if not regions:
            return None
        merged = unary_union(regions).buffer(0.0)
        polygons = [
            polygon
            for polygon in _flatten_polygons(merged)
            if not polygon.is_empty and polygon.is_valid and polygon.area >= MIN_PATCH_AREA
        ]
        if len(polygons) != STEPPED_FLAT_REGION_MAX_MULTIPART:
            return None
        if any(len(polygon.interiors) > 0 for polygon in polygons):
            return None
        candidate_regions.extend(polygons)
    covered = unary_union(candidate_regions)
    overlap_area = sum(region.area for region in candidate_regions) - covered.area
    if overlap_area / footprint.area > MAX_UNCOVERED_FOOTPRINT_FRACTION:
        return None
    if not _patches_cover_footprint(footprint, candidate_regions):
        return None
    return candidate_regions


def _stepped_flat_strip_level_patches(
    points: np.ndarray,
    footprint: Polygon,
    levels: list[SteppedFlatLevel] | None,
) -> tuple[list[Polygon] | None, str | None]:
    if levels is None or len(levels) < STEPPED_FLAT_MIN_LEVELS:
        return None, STEPPED_FLAT_INSUFFICIENT_LEVELS
    footprint_metrics = _minimum_rotated_rectangle_metrics(footprint)
    if footprint_metrics is None or footprint_metrics[0] < RECTANGULAR_FOOTPRINT_MIN_RATIO:
        return None, STEPPED_FLAT_PATCH_FAILED
    axes = _minimum_rotated_rectangle_axes(footprint)
    if axes is None:
        return None, STEPPED_FLAT_PATCH_FAILED

    points = np.asarray(points, dtype=float)
    footprint_local = _project_xy_to_axes(
        np.asarray(footprint.exterior.coords[:-1], dtype=float),
        axes,
    )
    points_local = _project_xy_to_axes(points[:, :2], axes)
    level_centers = np.asarray(
        [
            np.mean(points_local[level.point_indices], axis=0)
            for level in levels
        ],
        dtype=float,
    )
    strip_axis = int(np.argmax(np.ptp(level_centers, axis=0)))
    cross_axis = 1 - strip_axis
    ordered = sorted(
        enumerate(levels),
        key=lambda item: float(level_centers[item[0], strip_axis]),
    )
    ordered_ranges = [
        (
            float(np.min(points_local[level.point_indices, strip_axis])),
            float(np.max(points_local[level.point_indices, strip_axis])),
        )
        for _, level in ordered
    ]
    strip_edges = [float(np.min(footprint_local[:, strip_axis]))]
    for index in range(len(ordered_ranges) - 1):
        current_max = ordered_ranges[index][1]
        next_min = ordered_ranges[index + 1][0]
        if next_min <= current_max + EDGE_TOLERANCE:
            return None, STEPPED_FLAT_PATCH_FAILED
        strip_edges.append(0.5 * (current_max + next_min))
    strip_edges.append(float(np.max(footprint_local[:, strip_axis])))
    cross_min = float(np.min(footprint_local[:, cross_axis])) - EDGE_TOLERANCE
    cross_max = float(np.max(footprint_local[:, cross_axis])) + EDGE_TOLERANCE

    all_level_indices = {
        int(index)
        for level in levels
        for index in level.point_indices
    }
    patches_by_level: list[Polygon | None] = [None] * len(levels)
    for ordered_index, (level_index, level) in enumerate(ordered):
        strip_min = strip_edges[ordered_index]
        strip_max = strip_edges[ordered_index + 1]
        if strip_max - strip_min <= EDGE_TOLERANCE:
            return None, STEPPED_FLAT_PATCH_FAILED
        strip = _local_rect_to_world_polygon(
            strip_axis,
            strip_min,
            strip_max,
            cross_min,
            cross_max,
            axes,
        )
        patch = _single_polygon_or_none(strip.intersection(footprint))
        if patch is None:
            return None, STEPPED_FLAT_PATCH_FAILED
        patch_metrics = _minimum_rotated_rectangle_metrics(patch)
        if patch_metrics is None or patch_metrics[0] < DECOMPOSITION_RECTANGULARITY_MIN_RATIO:
            return None, STEPPED_FLAT_PATCH_FAILED
        if not _level_patch_support_is_clean(points, patch, level, all_level_indices):
            return None, STEPPED_FLAT_PATCH_FAILED
        patches_by_level[level_index] = patch

    patches = [patch for patch in patches_by_level if patch is not None]
    if len(patches) != len(levels) or not _patches_cover_footprint(footprint, patches):
        return None, STEPPED_FLAT_PATCH_FAILED
    return patches, None


def _stepped_flat_region_level_patches(
    points: np.ndarray,
    footprint: Polygon,
    levels: list[SteppedFlatLevel] | None,
) -> tuple[list[Polygon] | None, str | None]:
    if levels is None or len(levels) < STEPPED_FLAT_MIN_LEVELS:
        return None, STEPPED_FLAT_INSUFFICIENT_LEVELS
    decomposition_regions = _decomposition_candidate_regions(footprint)
    if decomposition_regions is not None:
        patches = _assign_candidate_regions_to_levels(
            points,
            footprint,
            decomposition_regions,
            levels,
        )
        if patches is not None:
            return patches, None
    grid_regions = _grid_candidate_regions(points, footprint, levels)
    if grid_regions is not None:
        patches = _assign_candidate_regions_to_levels(
            points,
            footprint,
            grid_regions,
            levels,
        )
        if patches is not None:
            return patches, None
    return None, STEPPED_FLAT_PATCH_FAILED


def _stepped_flat_level_patches(
    points: np.ndarray,
    footprint: Polygon,
    levels: list[SteppedFlatLevel] | None,
) -> tuple[list[Polygon] | None, str | None]:
    patches, reason = _stepped_flat_strip_level_patches(points, footprint, levels)
    if patches is not None:
        return patches, reason
    if reason != STEPPED_FLAT_PATCH_FAILED:
        return None, reason
    return _stepped_flat_region_level_patches(points, footprint, levels)


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
    roof_surfaces = _reconcile_slice_line_vertices(
        footprint,
        roof_surfaces,
        decomposition.slice_lines,
    )
    roof_surfaces = _weld_near_equal_boundary_endpoint_heights(
        footprint,
        roof_surfaces,
        decomposition.slice_lines,
    )
    roof_surfaces = _split_crossing_junction_edges(roof_surfaces, decomposition.slice_lines)
    junctions = _internal_junction_surfaces(footprint, roof_surfaces, decomposition.slice_lines)
    if junctions is None:
        return None, DECOMPOSITION_JUNCTION_FAILED
    walls = _wall_surfaces(footprint, roof_surfaces, ground_height)
    if len(walls) == 0:
        _record_decomposition_subreason(decomposition_counts, WATERTIGHT_OTHER)
        return None, DECOMPOSITION_WATERTIGHT_FAILED
    shell = MultiSurface(surfaces=[*roof_surfaces, *junctions, *walls, _ground_surface(footprint, ground_height)])
    if not is_watertight(shell):
        watertight_reason = _decomposed_watertight_failure_reason(
            footprint,
            roof_surfaces,
            decomposition.slice_lines,
            shell,
        )
        _record_decomposition_subreason(decomposition_counts, watertight_reason)
        if watertight_reason == WATERTIGHT_EDGE_COUNT_MISMATCH:
            _record_decomposition_subreason(
                decomposition_counts,
                _edge_count_mismatch_reason(footprint, decomposition.slice_lines, shell),
            )
        return None, DECOMPOSITION_WATERTIGHT_FAILED
    if decomposition_counts is not None:
        decomposition_counts[DECOMPOSITION_SUCCESS] += 1
    return shell, DECOMPOSITION_SUCCESS


def _line_strings_from_boundary_intersection(geometry) -> list[LineString]:
    if geometry.is_empty:
        return []
    geom_type = getattr(geometry, "geom_type", "")
    if geom_type in ("LineString", "LinearRing"):
        coords = list(geometry.coords)
        return [
            LineString([start, end])
            for start, end in zip(coords, coords[1:])
            if LineString([start, end]).length > EDGE_TOLERANCE
        ]
    if hasattr(geometry, "geoms"):
        lines: list[LineString] = []
        for part in geometry.geoms:
            lines.extend(_line_strings_from_boundary_intersection(part))
        return lines
    return []


def _stepped_flat_slice_lines(patches: list[Polygon]) -> list[LineString]:
    slice_lines: list[LineString] = []
    handled_lines: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for first_index, first_patch in enumerate(patches):
        for second_patch in patches[first_index + 1 :]:
            intersection = first_patch.boundary.intersection(second_patch.boundary)
            for line in _line_strings_from_boundary_intersection(intersection):
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
                slice_lines.append(line)
    return slice_lines


def _node_stepped_flat_roof_surfaces(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> list[Surface]:
    noded = [
        Surface(vertices=np.asarray(surface.vertices, dtype=float).copy())
        for surface in roof_surfaces
    ]
    for line in slice_lines:
        coords = np.asarray(line.coords, dtype=float)
        for xy in (coords[0], coords[-1]):
            if footprint.boundary.distance(Point(xy)) <= EDGE_TOLERANCE:
                continue
            noded = [_insert_xy_on_surface_edge(surface, xy) for surface in noded]
    return noded


def _has_multi_height_slice_endpoint(
    footprint: Polygon,
    roof_surfaces: list[Surface],
    slice_lines: list[LineString],
) -> bool:
    endpoints: dict[tuple[int, int], np.ndarray] = {}
    for line in slice_lines:
        coords = np.asarray(line.coords, dtype=float)
        for xy in (coords[0], coords[-1]):
            if footprint.boundary.distance(Point(xy)) <= EDGE_TOLERANCE:
                continue
            key = tuple(np.round(xy / EDGE_TOLERANCE).astype(int))
            endpoints[key] = xy
    for xy in endpoints.values():
        heights = sorted(
            float(vertex[2])
            for surface in roof_surfaces
            for vertex in surface.vertices
            if np.linalg.norm(vertex[:2] - xy) <= EDGE_TOLERANCE
        )
        unique_heights: list[float] = []
        for height in heights:
            if not unique_heights or abs(height - unique_heights[-1]) > EDGE_TOLERANCE:
                unique_heights.append(height)
        if len(unique_heights) > 2:
            return True
    return False


def _build_stepped_flat_shell(
    footprint: Polygon,
    patches: list[Polygon] | None,
    levels: list[SteppedFlatLevel] | None,
    ground_height: float,
) -> tuple[MultiSurface | None, str]:
    if patches is None or levels is None or len(patches) == 0 or len(patches) != len(levels):
        return None, STEPPED_FLAT_PATCH_FAILED

    roof_surfaces = [
        _surface_from_xy(patch, level.plane)
        for patch, level in zip(patches, levels)
    ]
    slice_lines = _stepped_flat_slice_lines(patches)
    roof_surfaces = _node_stepped_flat_roof_surfaces(footprint, roof_surfaces, slice_lines)
    if _has_multi_height_slice_endpoint(footprint, roof_surfaces, slice_lines):
        return None, STEPPED_FLAT_STEP_WALL_FAILED
    junctions = _internal_junction_surfaces(footprint, roof_surfaces, slice_lines)
    if junctions is None:
        return None, STEPPED_FLAT_STEP_WALL_FAILED

    walls = _wall_surfaces(footprint, roof_surfaces, ground_height)
    if len(walls) == 0:
        return None, STEPPED_FLAT_WATERTIGHT_FAILED
    surfaces = _reconcile_shell_edge_vertices(
        [*roof_surfaces, *junctions, *walls, _ground_surface(footprint, ground_height)]
    )
    shell = MultiSurface(surfaces=surfaces)
    if not _surfaces_are_simple(shell.surfaces):
        return None, STEPPED_FLAT_WATERTIGHT_FAILED
    if not is_watertight(shell):
        return None, STEPPED_FLAT_WATERTIGHT_FAILED
    return shell, STEPPED_FLAT_SUCCESS


def _candidate_stepped_flat_lod2(
    footprint: Polygon,
    roof_points: np.ndarray,
    ground_height: float,
    planes: list[RoofPlane],
    stepped_flat_counts: Counter | None = None,
) -> tuple[MultiSurface | None, str | None]:
    _record_stepped_flat_reason(stepped_flat_counts, STEPPED_FLAT_CANDIDATE)
    levels, reason = _stepped_flat_height_levels(roof_points, planes)
    if levels is None:
        if reason == STEPPED_FLAT_INSUFFICIENT_LEVELS:
            flat_plane = _flat_collapse_plane(roof_points, planes)
            if flat_plane is not None:
                shell = _build_shell(
                    footprint,
                    [_surface_from_xy(footprint, flat_plane)],
                    ground_height,
                )
                if shell is not None:
                    _record_stepped_flat_reason(
                        stepped_flat_counts,
                        STEPPED_FLAT_FLAT_COLLAPSE_SUCCESS,
                    )
                    return shell, None
                reason = STEPPED_FLAT_WATERTIGHT_FAILED
        if reason is not None:
            _record_stepped_flat_reason(stepped_flat_counts, reason)
        return None, reason
    patches, reason = _stepped_flat_level_patches(roof_points, footprint, levels)
    if patches is None:
        if reason is not None:
            _record_stepped_flat_reason(stepped_flat_counts, reason)
        return None, reason
    shell, reason = _build_stepped_flat_shell(footprint, patches, levels, ground_height)
    if shell is None:
        if reason is not None:
            _record_stepped_flat_reason(stepped_flat_counts, reason)
        return None, reason
    _record_stepped_flat_reason(stepped_flat_counts, STEPPED_FLAT_SUCCESS)
    return shell, None


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
    stepped_flat_counts: Counter | None = None,
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
        shell, _ = _candidate_stepped_flat_lod2(
            footprint,
            roof_points,
            ground_height,
            planes,
            stepped_flat_counts,
        )
        if shell is not None:
            return shell
        _record_rejection(rejections, UNSUPPORTED_PLANE_COUNT)
        return None

    roof_surfaces, rejection_reason = _roof_surfaces_for_planes(roof_points, footprint, planes)
    if roof_surfaces is None:
        if rejection_reason == SPLIT_FAILED:
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
            shell, _ = _candidate_stepped_flat_lod2(
                footprint,
                roof_points,
                ground_height,
                planes,
                stepped_flat_counts,
            )
            if shell is not None:
                return shell
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
    stepped_flat_counts: Counter,
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
    stepped_flat_parts = [
        f"{reason}={stepped_flat_counts[reason]}"
        for reason in STEPPED_FLAT_REASONS
        if stepped_flat_counts[reason] > 0
    ]
    if stepped_flat_parts:
        info("LOD2 stepped-flat summary: " + " ".join(stepped_flat_parts))
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
    stepped_flat_counts: Counter | None = None,
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
        stepped_flat_counts,
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
    stepped_flat_counts = Counter()
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
            stepped_flat_counts if log_rejections else None,
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
            stepped_flat_counts=stepped_flat_counts,
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
