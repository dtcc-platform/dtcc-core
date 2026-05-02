from __future__ import annotations

from collections import Counter

import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import split, unary_union

from dtcc_core.logging import info

from ...model import Building, GeometryType, MultiSurface, Surface
from .buildings import build_lod1_buildings


MIN_ROOF_POINTS = 24
MIN_PLANE_INLIERS = 8
# RANSAC may inspect richer roofs; shell assembly currently accepts one or two planes.
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
SHARED_VERTEX_KEY_DECIMALS = 6
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


def _points_on_edge(points: list[np.ndarray], start: np.ndarray, end: np.ndarray) -> list[np.ndarray]:
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
    return unique


def _wall_surfaces(footprint: Polygon, roof_surfaces: list[Surface], ground_height: float) -> list[Surface]:
    roof_vertices = [vertex for surface in roof_surfaces for vertex in surface.vertices]
    coords = np.array(footprint.exterior.coords[:-1], dtype=float)
    walls = []
    for index, start in enumerate(coords):
        end = coords[(index + 1) % len(coords)]
        top = _points_on_edge(roof_vertices, start, end)
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


def _log_lod2_summary(
    *,
    total: int,
    lod2_count: int,
    fallback_count: int,
    skipped_existing: int,
    rejections: Counter,
    plane_counts: Counter,
    max_plane_coplanar_pair_counts: Counter,
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
    planes = _ransac_planes(roof_points)
    if plane_counts is not None:
        plane_counts[len(planes)] += 1
    if len(planes) == MAX_PLANES and max_plane_coplanar_pair_counts is not None:
        max_plane_coplanar_pair_counts[_coplanar_pair_count(planes)] += 1
    if len(planes) == 0:
        _record_rejection(rejections, NO_PLANES_FOUND)
        return None
    if len(planes) > 2:
        _record_rejection(rejections, UNSUPPORTED_PLANE_COUNT)
        return None
    if len(planes) == 1:
        roof_surfaces = [_surface_from_xy(footprint, planes[0])]
        shell = _build_shell(footprint, roof_surfaces, float(ground_height))
        if shell is None:
            _record_rejection(rejections, WATERTIGHT_SHELL_FAILED)
        return shell
    roof_surfaces, rejection_reason = _multi_plane_roof_surfaces(roof_points, footprint, planes)
    if roof_surfaces is None:
        _record_rejection(rejections, rejection_reason or SHELL_ASSEMBLY_FAILED)
        return None
    shell = _build_shell(footprint, roof_surfaces, float(ground_height))
    if shell is None:
        _record_rejection(rejections, WATERTIGHT_SHELL_FAILED)
    return shell


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
