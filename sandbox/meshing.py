# Copyright (C) 2025 Vasilis Naserentin
# Licensed under the MIT License
# Codex cli, gpt-5-codex high
# Specs requirements: 1) Manifold geometry 2) No overlapping and intersecting triangles 
# 3) No disconnections 4) Not (very) large variations in triangle size
# 4) Being able to enforce min/max quality constraints (e.g. angle, AR)
# Usage: python3 meshing.py --min-angle X --threads X --obj mesh.obj
                                                             
import math
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple, TypeVar

import gmsh

Point2D = Tuple[float, float]
T = TypeVar("T")


@dataclass(frozen=True)
class LengthStats:
    minimum: float
    median: float
    maximum: float
    mean: float
    p90: float


def read_polygons(path: Path) -> Tuple[List[Point2D], List[List[Point2D]]]:
    """Parse polygonal loops from the testcase file."""
    polygons: List[List[Point2D]] = []
    with path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            tokens = line.split()
            if len(tokens) % 2 != 0:
                raise ValueError(f"Line has an odd number of coordinates: '{line}'")
            coords = [(float(tokens[i]), float(tokens[i + 1])) for i in range(0, len(tokens), 2)]
            polygons.append(_normalise_polygon(coords))
    if not polygons:
        raise ValueError("No polygons found in testcase")
    exterior, *holes = polygons
    exterior = _ensure_orientation(exterior, clockwise=False)
    holes = [_ensure_orientation(hole, clockwise=True) for hole in holes]
    return exterior, holes


def _ensure_orientation(coords: Sequence[Point2D], *, clockwise: bool) -> List[Point2D]:
    """Return coords with requested orientation."""
    area = _signed_area(coords)
    if clockwise and area > 0:
        return list(reversed(coords))
    if not clockwise and area < 0:
        return list(reversed(coords))
    return list(coords)


def _signed_area(coords: Sequence[Point2D]) -> float:
    area = 0.0
    for (x1, y1), (x2, y2) in _pairwise_closed(coords):
        area += x1 * y2 - x2 * y1
    return 0.5 * area


def _normalise_polygon(coords: Sequence[Point2D]) -> List[Point2D]:
    cleaned: List[Point2D] = []
    for x, y in coords:
        if not cleaned or _dist2((x, y), cleaned[-1]) > 1e-12:
            cleaned.append((x, y))
    if len(cleaned) >= 3 and _dist2(cleaned[0], cleaned[-1]) <= 1e-12:
        cleaned.pop()
    if len(cleaned) < 3:
        raise ValueError("Polygon must have at least three distinct vertices")
    return cleaned


def _dist2(p: Point2D, q: Point2D) -> float:
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2


def _pairwise_closed(items: Sequence[T]) -> Iterable[Tuple[T, T]]:
    for i in range(len(items)):
        yield items[i], items[(i + 1) % len(items)]


def _edge_lengths(coords: Sequence[Point2D]) -> List[float]:
    return [math.hypot(x2 - x1, y2 - y1) for (x1, y1), (x2, y2) in _pairwise_closed(coords)]


def _length_statistics(polygons: Sequence[Sequence[Point2D]]) -> LengthStats:
    lengths: List[float] = []
    for poly in polygons:
        lengths.extend(_edge_lengths(poly))
    filtered = [l for l in lengths if l > 1e-8]
    if not filtered:
        raise ValueError("Could not determine characteristic lengths")
    filtered.sort()
    min_len = filtered[0]
    max_len = filtered[-1]
    median_len = statistics.median(filtered)
    mean_len = statistics.fmean(filtered) if hasattr(statistics, "fmean") else statistics.mean(filtered)
    if len(filtered) == 1:
        p90 = filtered[0]
    else:
        idx = int(round(0.9 * (len(filtered) - 1)))
        p90 = filtered[idx]
    return LengthStats(min_len, median_len, max_len, mean_len, p90)


def _vertex_angle_deg(coords: Sequence[Point2D], index: int) -> float:
    n = len(coords)
    prev_pt = coords[index - 1]
    curr_pt = coords[index]
    next_pt = coords[(index + 1) % n]
    v1x, v1y = prev_pt[0] - curr_pt[0], prev_pt[1] - curr_pt[1]
    v2x, v2y = next_pt[0] - curr_pt[0], next_pt[1] - curr_pt[1]
    n1 = math.hypot(v1x, v1y)
    n2 = math.hypot(v2x, v2y)
    if n1 < 1e-12 or n2 < 1e-12:
        return 180.0
    cos_theta = (v1x * v2x + v1y * v2y) / (n1 * n2)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))


def _point_mesh_size(
    coords: Sequence[Point2D],
    index: int,
    *,
    min_size: float,
    max_size: float,
    corner_scale: float,
) -> float:
    prev_pt = coords[index - 1]
    curr_pt = coords[index]
    next_pt = coords[(index + 1) % len(coords)]
    prev_len = math.hypot(curr_pt[0] - prev_pt[0], curr_pt[1] - prev_pt[1])
    next_len = math.hypot(next_pt[0] - curr_pt[0], next_pt[1] - curr_pt[1])
    local_len = max(min(prev_len, next_len), 1e-6)
    angle = _vertex_angle_deg(coords, index)
    angle_factor = 1.0
    if angle < 75.0:
        angle_factor = max(angle / 75.0, 0.25)
    size = local_len * corner_scale * angle_factor
    size = min(size, local_len)
    size = max(min_size, min(max_size, size))
    return size


def _create_curve_loop(
    coords: Sequence[Point2D],
    *,
    min_size: float,
    max_size: float,
    corner_scale: float,
) -> int:
    point_tags: List[int] = []
    for idx, (x, y) in enumerate(coords):
        mesh_size = _point_mesh_size(
            coords,
            idx,
            min_size=min_size,
            max_size=max_size,
            corner_scale=corner_scale,
        )
        point_tags.append(gmsh.model.geo.addPoint(x, y, 0.0, mesh_size))
    line_tags: List[int] = []
    for start, end in _pairwise_closed(point_tags):
        line_tags.append(gmsh.model.geo.addLine(start, end))
    return gmsh.model.geo.addCurveLoop(line_tags)


def _maybe_set_option(name: str, value: float) -> bool:
    try:
        gmsh.option.setNumber(name, float(value))
        return True
    except Exception as exc:  # pragma: no cover - gmsh handles logging
        sys.stderr.write(f"[gmsh] Warning: could not set option '{name}': {exc}\n")
        return False


def _safe_optimize(method: str) -> None:
    try:
        gmsh.model.mesh.optimize(method)
    except Exception as exc:  # pragma: no cover - gmsh handles logging
        sys.stderr.write(f"[gmsh] Warning: optimize '{method}' failed: {exc}\n")


def _triangle_angles_deg(points: Sequence[Point2D]) -> List[float]:
    (ax, ay), (bx, by), (cx, cy) = points

    def angle(p: Point2D, q: Point2D, r: Point2D) -> float:
        v1x, v1y = p[0] - q[0], p[1] - q[1]
        v2x, v2y = r[0] - q[0], r[1] - q[1]
        n1 = math.hypot(v1x, v1y)
        n2 = math.hypot(v2x, v2y)
        if n1 < 1e-12 or n2 < 1e-12:
            return 180.0
        cos_theta = (v1x * v2x + v1y * v2y) / (n1 * n2)
        cos_theta = max(-1.0, min(1.0, cos_theta))
        return math.degrees(math.acos(cos_theta))

    return [
        angle((ax, ay), (bx, by), (cx, cy)),
        angle((bx, by), (cx, cy), (ax, ay)),
        angle((cx, cy), (ax, ay), (bx, by)),
    ]


def _triangle_aspect_ratio(points: Sequence[Point2D]) -> float:
    (ax, ay), (bx, by), (cx, cy) = points
    ab = math.hypot(bx - ax, by - ay)
    bc = math.hypot(cx - bx, cy - by)
    ca = math.hypot(ax - cx, ay - cy)
    longest = max(ab, bc, ca)
    if longest < 1e-12:
        return float("inf")
    area2 = abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
    if area2 < 1e-16:
        return float("inf")
    return (longest ** 2) / area2


def _triangle_metrics() -> List[Tuple[int, float, float]]:
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    if not len(node_tags):
        return []
    coords = {
        int(tag): (node_coords[i * 3], node_coords[i * 3 + 1])
        for i, tag in enumerate(node_tags)
    }

    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    result: List[Tuple[int, float, float]] = []
    for elem_type, tags, nodes in zip(elem_types, elem_tags, elem_node_tags):
        name, dim, *_ = gmsh.model.mesh.getElementProperties(elem_type)
        if dim != 2 or "triangle" not in name.lower():
            continue
        tag_list = [int(t) for t in tags]
        if not tag_list:
            continue
        num_nodes = len(nodes) // len(tag_list)
        for idx, tag in enumerate(tag_list):
            start = idx * num_nodes
            tri_nodes = nodes[start : start + 3]
            try:
                tri_points = [coords[int(tn)] for tn in tri_nodes]
            except KeyError:
                continue
            angles = _triangle_angles_deg(tri_points)
            aspect = _triangle_aspect_ratio(tri_points)
            result.append((tag, min(angles), aspect))
    return result


def _min_triangle_angle_deg() -> float:
    metrics = _triangle_metrics()
    return min((ang for _, ang, _ in metrics), default=180.0)


def _collect_bad_triangles(threshold: float, limit: int = 2000) -> List[int]:
    candidates = [(tag, angle) for tag, angle, _ in _triangle_metrics() if angle < threshold]
    if not candidates:
        return []
    candidates.sort(key=lambda item: item[1])
    if limit and len(candidates) > limit:
        candidates = candidates[:limit]
    return [tag for tag, _ in candidates]


def _refine_bad_triangles(tags: Sequence[int]) -> bool:
    if not tags:
        return False
    splitter = getattr(gmsh.model.mesh, "splitElements", None)
    if callable(splitter):
        try:
            splitter(2, list(tags))
            return True
        except Exception as exc:  # pragma: no cover - gmsh handles logging
            sys.stderr.write(f"[gmsh] Warning: splitElements failed: {exc}\n")
    try:
        gmsh.model.mesh.refine()
        return True
    except Exception as exc:  # pragma: no cover - gmsh handles logging
        sys.stderr.write(f"[gmsh] Warning: global refine failed: {exc}\n")
        return False


def _enforce_minimum_angle(
    target_angle_deg: float,
    *,
    max_optimize_passes: int = 8,
    max_refinement_rounds: int = 5,
) -> float:
    min_angle = _min_triangle_angle_deg()
    opt_passes = 0
    while min_angle < target_angle_deg and opt_passes < max_optimize_passes:
        _safe_optimize("Netgen")
        _safe_optimize("Relocate2D")
        _safe_optimize("Laplace2D")
        min_angle = _min_triangle_angle_deg()
        opt_passes += 1
    if opt_passes:
        sys.stderr.write(
            f"[gmsh] Applied {opt_passes} smoothing pass(es); min angle is {min_angle:.2f}°\n"
        )
    if min_angle >= target_angle_deg:
        return min_angle

    refinement_round = 0
    while min_angle < target_angle_deg and refinement_round < max_refinement_rounds:
        bad_tags = _collect_bad_triangles(max(target_angle_deg - 0.5, target_angle_deg * 0.8))
        if not _refine_bad_triangles(bad_tags):
            break
        gmsh.model.mesh.removeDuplicateNodes()
        _safe_optimize("Netgen")
        _safe_optimize("Relocate2D")
        _safe_optimize("Laplace2D")
        refinement_round += 1
        min_angle = _min_triangle_angle_deg()
        sys.stderr.write(
            f"[gmsh] Refined round {refinement_round}; min angle now {min_angle:.2f}°\n"
        )
    return min_angle


def _min_polygon_angle(polygons: Sequence[Sequence[Point2D]]) -> float:
    min_angle = 180.0
    for poly in polygons:
        n = len(poly)
        for i in range(n):
            p_prev = poly[i - 1]
            p_curr = poly[i]
            p_next = poly[(i + 1) % n]
            v1x, v1y = p_prev[0] - p_curr[0], p_prev[1] - p_curr[1]
            v2x, v2y = p_next[0] - p_curr[0], p_next[1] - p_curr[1]
            n1 = math.hypot(v1x, v1y)
            n2 = math.hypot(v2x, v2y)
            if n1 < 1e-12 or n2 < 1e-12:
                continue
            cos_theta = (v1x * v2x + v1y * v2y) / (n1 * n2)
            cos_theta = max(-1.0, min(1.0, cos_theta))
            angle = math.degrees(math.acos(cos_theta))
            min_angle = min(min_angle, angle)
    return min_angle


def build_mesh(
    case_path: Path,
    *,
    msh_path: Path,
    geo_path: Path | None = None,
    obj_path: Path | None = None,
    min_angle_deg: float = 30.0,
    threads: int | None = None,
) -> None:
    min_angle_deg = max(1.0, min(min_angle_deg, 89.0))
    if threads is not None:
        threads = max(1, int(threads))
    exterior, holes = read_polygons(case_path)
    stats = _length_statistics([exterior, *holes])
    geom_min_angle = _min_polygon_angle([exterior, *holes])
    target_angle = min_angle_deg
    geom_slack = 1.0
    if geom_min_angle + 1e-9 < target_angle:
        adjusted_target = max(1.0, geom_min_angle - geom_slack)
        if adjusted_target < target_angle:
            sys.stderr.write(
                f"[gmsh] Geometry corner of {geom_min_angle:.2f}° detected; lowering target to {adjusted_target:.2f}°.\n"
            )
            target_angle = adjusted_target

    outer_min_size = max(0.02, stats.minimum * 0.3, stats.median * 0.02)
    hole_min_size = max(0.02, outer_min_size * 0.6, stats.minimum * 0.2)
    outer_max_size = min(stats.maximum, max(stats.p90, stats.median * 2.0))
    hole_max_size = min(stats.maximum, max(hole_min_size * 3.0, stats.median * 1.2))
    outer_max_size = max(outer_max_size, outer_min_size * 1.2)
    hole_max_size = max(hole_max_size, hole_min_size * 1.2)
    global_min_size = min(outer_min_size, hole_min_size)
    global_max_size = max(outer_max_size, hole_max_size)
    outer_corner_scale = 0.85
    inner_corner_scale = 0.65

    gmsh_args: List[str] = []
    if threads:
        gmsh_args.extend(["-nt", str(threads)])

    gmsh.initialize(gmsh_args)
    gmsh.option.setNumber("General.Terminal", 1)
    if threads:
        gmsh.option.setNumber("General.NumThreads", float(threads))
        _maybe_set_option("Mesh.MaxNumThreads", float(threads))
        _maybe_set_option("Mesh.OptimizeMaxNumThreads", float(threads))
    try:
        gmsh.model.add("site_mesh")

        outer_loop = _create_curve_loop(
            exterior,
            min_size=outer_min_size,
            max_size=outer_max_size,
            corner_scale=outer_corner_scale,
        )
        hole_loops = [
            _create_curve_loop(
                hole,
                min_size=hole_min_size,
                max_size=hole_max_size,
                corner_scale=inner_corner_scale,
            )
            for hole in holes
        ]

        gmsh.model.geo.addPlaneSurface([outer_loop, *hole_loops])
        gmsh.model.geo.synchronize()

        _maybe_set_option("Mesh.Algorithm", 6)  # Frontal-Delaunay for 2D
        if not _maybe_set_option("Mesh.MinimumElementQuality", 0.4):
            _maybe_set_option("Mesh.OptimizeThreshold", 0.4)
        _maybe_set_option("Mesh.SmoothSteps", 15)
        _maybe_set_option("Mesh.Optimize", 1)
        _maybe_set_option("Mesh.OptimizeNetgen", 1)
        _maybe_set_option("Mesh.CharacteristicLengthMin", global_min_size)
        _maybe_set_option("Mesh.CharacteristicLengthMax", global_max_size)
        _maybe_set_option("Mesh.MeshSizeFromPoints", 1)
        _maybe_set_option("Mesh.MeshSizeFromCurvature", 1)
        _maybe_set_option("Mesh.MeshSizeExtendFromBoundary", 1)

        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.removeDuplicateNodes()
        _safe_optimize("Netgen")
        _safe_optimize("Relocate2D")
        _safe_optimize("Laplace2D")
        gmsh.model.mesh.setOrder(1)

        min_angle = _enforce_minimum_angle(target_angle)
        if min_angle < target_angle - 1.0:
            raise RuntimeError(
                f"Unable to reach requested minimum angle: {min_angle:.2f}° < {target_angle:.2f}°."
            )

        metrics = _triangle_metrics()
        aspects = [aspect for _, _, aspect in metrics if aspect > 0 and math.isfinite(aspect)]
        if aspects:
            mean_func = getattr(statistics, "fmean", None)
            mean_ar = mean_func(aspects) if callable(mean_func) else sum(aspects) / len(aspects)
            print(
                f"[mesh] Triangle aspect ratio -> min: {min(aspects):.3f}, max: {max(aspects):.3f}, mean: {mean_ar:.3f}"
            )
        else:
            print("[mesh] Triangle aspect ratio -> unavailable (no valid triangles)")

        gmsh.write(str(msh_path))
        if geo_path is not None:
            gmsh.write(str(geo_path))
        if obj_path is not None:
            try:
                gmsh.write(str(obj_path))
            except Exception as exc:  # pragma: no cover - gmsh handles logging
                sys.stderr.write(f"[gmsh] Warning: could not write OBJ '{obj_path}': {exc}\n")
    finally:
        gmsh.finalize()


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Generate a high-quality 2D mesh from testcase polygons.")
    parser.add_argument("case", nargs="?", default="testcase.txt", type=Path, help="Input testcase file")
    parser.add_argument("--msh", default=Path("mesh.msh"), type=Path, help="Output .msh path")
    parser.add_argument("--geo", default=None, type=Path, help="Optional .geo_unrolled output path")
    parser.add_argument("--obj", default=None, type=Path, help="Optional .obj output path")
    parser.add_argument(
        "--min-angle",
        default=30.0,
        type=float,
        help="Requested minimum triangle angle in degrees (clamped by geometry corners)",
    )
    parser.add_argument(
        "--threads",
        default=None,
        type=int,
        help="Number of parallel threads for gmsh (leave unset for gmsh default)",
    )
    args = parser.parse_args()

    build_mesh(
        args.case,
        msh_path=args.msh,
        geo_path=args.geo,
        obj_path=args.obj,
        min_angle_deg=args.min_angle,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
