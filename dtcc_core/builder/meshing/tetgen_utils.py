from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, Mapping, Sequence

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import triangulate
from shapely.validation import explain_validity

from ...model import Mesh
from ...model.mixins.mesh.quality import tri_aspect_ratio, tri_element_quality


_MIN_EDGE_RATIO_WARNING = 1.0e-3
_MIN_AREA_RATIO_WARNING = 1.0e-6
_MIN_TRI_QUALITY_WARNING = 2.0e-2
_MAX_TRI_ASPECT_RATIO_WARNING = 1.0e2
_BOUNDARY_PINCH_RATIO_WARNING = 1.0e-3
_MAX_VERTICAL_FACE_ASPECT_RATIO = 15.0
_POLYGON_SIDEWALL_MIN_TRI_QUALITY = 2.5e-1
_POLYGON_SIDEWALL_MAX_TRI_ASPECT_RATIO = 4.0


def _vertical_quad_strip_count(
    horizontal_length: float,
    vertical_height: float,
    *,
    target_aspect_ratio: float = _MAX_VERTICAL_FACE_ASPECT_RATIO,
) -> int:
    if horizontal_length <= 0.0 or vertical_height <= 0.0:
        return 1

    target_height = horizontal_length * float(target_aspect_ratio)
    if target_height <= 0.0:
        return 1

    return max(1, int(np.ceil(vertical_height / target_height)))


def _sidewall_strip_count(
    vertices: np.ndarray,
    bottom_loop: np.ndarray,
    top_loop: np.ndarray,
) -> int:
    if len(bottom_loop) < 2 or len(top_loop) != len(bottom_loop):
        return 1

    bottom_points = vertices[np.asarray(bottom_loop, dtype=np.int64)]
    top_points = vertices[np.asarray(top_loop, dtype=np.int64)]
    edge_vectors = np.diff(bottom_points[:, :2], axis=0)
    edge_lengths = np.linalg.norm(edge_vectors, axis=1)
    positive_lengths = edge_lengths[edge_lengths > 0.0]
    if positive_lengths.size == 0:
        return 1

    min_horizontal_length = float(np.min(positive_lengths))
    max_vertical_height = float(np.max(np.abs(top_points[:, 2] - bottom_points[:, 2])))
    return _vertical_quad_strip_count(min_horizontal_length, max_vertical_height)


def _boundary_sidewall_strip_count(
    vertices: np.ndarray,
    bottom_loops: Mapping[str, np.ndarray],
    top_loops: Mapping[str, np.ndarray],
) -> int:
    return max(
        _sidewall_strip_count(
            vertices,
            np.asarray(bottom_loops[name], dtype=np.int64),
            np.asarray(top_loops[name], dtype=np.int64),
        )
        for name in ("south", "east", "north", "west")
    )


def _prefer_polygon_sidewalls(vertices: np.ndarray, faces: np.ndarray) -> bool:
    F = np.asarray(faces, dtype=np.int64)
    if F.size == 0:
        return True

    V = np.asarray(vertices, dtype=float)
    triangle_quality = tri_element_quality(V, F)
    triangle_aspect_ratio = tri_aspect_ratio(V, F)
    if triangle_quality.size == 0 or triangle_aspect_ratio.size == 0:
        return True

    return (
        float(triangle_quality.min()) >= _POLYGON_SIDEWALL_MIN_TRI_QUALITY
        and float(triangle_aspect_ratio.max()) <= _POLYGON_SIDEWALL_MAX_TRI_ASPECT_RATIO
    )


def _append_vertical_sidewall_facet(
    vertices_out: list[list[float]],
    boundary_facets: list[list[int]],
    *,
    bottom_v0: int,
    bottom_v1: int,
    top_v0: int,
    top_v1: int,
    strip_count: int,
    vertical_vertex_cache: dict[tuple[int, int, int, int], int],
) -> None:
    def _get_strip_vertex(bottom_v: int, top_v: int, step: int) -> int:
        if step <= 0:
            return int(bottom_v)
        if step >= strip_count:
            return int(top_v)

        cache_key = (int(bottom_v), int(top_v), int(step), int(strip_count))
        cached = vertical_vertex_cache.get(cache_key)
        if cached is not None:
            return cached

        bottom_point = np.asarray(vertices_out[int(bottom_v)], dtype=float)
        top_point = np.asarray(vertices_out[int(top_v)], dtype=float)
        t = float(step) / float(strip_count)
        point = bottom_point + t * (top_point - bottom_point)
        vertex_index = len(vertices_out)
        vertices_out.append(point.tolist())
        vertical_vertex_cache[cache_key] = vertex_index
        return vertex_index

    lower_v0 = int(bottom_v0)
    lower_v1 = int(bottom_v1)
    for step in range(1, max(1, int(strip_count)) + 1):
        upper_v0 = _get_strip_vertex(bottom_v0, top_v0, step)
        upper_v1 = _get_strip_vertex(bottom_v1, top_v1, step)
        boundary_facets.append([lower_v0, lower_v1, upper_v1])
        boundary_facets.append([lower_v0, upper_v1, upper_v0])
        lower_v0 = upper_v0
        lower_v1 = upper_v1


def _triangle_edge_orientation_stats(
    shell_faces: np.ndarray,
    boundary_facets: Sequence[Sequence[int]],
) -> dict[str, int]:
    shell_faces = np.asarray(shell_faces, dtype=np.int64)
    boundary_triangles = [
        np.asarray(facet, dtype=np.int64).reshape(-1)
        for facet in boundary_facets
        if len(np.asarray(facet, dtype=np.int64).reshape(-1)) == 3
    ]
    if shell_faces.ndim != 2 or shell_faces.shape[1] != 3 or not boundary_triangles:
        return {
            "shared_edge_count": 0,
            "same_direction_shared_edge_count": 0,
            "opposite_direction_shared_edge_count": 0,
            "shell_boundary_same_direction_edge_count": 0,
            "shell_boundary_opposite_direction_edge_count": 0,
            "boundary_boundary_same_direction_edge_count": 0,
            "boundary_boundary_opposite_direction_edge_count": 0,
        }

    edge_entries: defaultdict[tuple[int, int], list[tuple[str, int, int, int]]] = defaultdict(list)
    for source, faces in (
        ("shell", shell_faces),
        ("boundary", np.vstack(boundary_triangles)),
    ):
        for face_index, face in enumerate(np.asarray(faces, dtype=np.int64)):
            for u, v in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
                edge_entries[tuple(sorted((int(u), int(v))))].append(
                    (source, face_index, int(u), int(v))
                )

    stats = {
        "shared_edge_count": 0,
        "same_direction_shared_edge_count": 0,
        "opposite_direction_shared_edge_count": 0,
        "shell_boundary_same_direction_edge_count": 0,
        "shell_boundary_opposite_direction_edge_count": 0,
        "boundary_boundary_same_direction_edge_count": 0,
        "boundary_boundary_opposite_direction_edge_count": 0,
    }
    for entries in edge_entries.values():
        if len(entries) != 2:
            continue
        stats["shared_edge_count"] += 1
        (source_a, _, u_a, v_a), (source_b, _, u_b, v_b) = entries
        same_direction = u_a == u_b and v_a == v_b
        if same_direction:
            stats["same_direction_shared_edge_count"] += 1
        else:
            stats["opposite_direction_shared_edge_count"] += 1

        if source_a != source_b:
            if same_direction:
                stats["shell_boundary_same_direction_edge_count"] += 1
            else:
                stats["shell_boundary_opposite_direction_edge_count"] += 1
        elif source_a == "boundary":
            if same_direction:
                stats["boundary_boundary_same_direction_edge_count"] += 1
            else:
                stats["boundary_boundary_opposite_direction_edge_count"] += 1

    return stats


def _triangle_signed_volume_sum(
    vertices: np.ndarray,
    faces: np.ndarray,
) -> float:
    if len(faces) == 0:
        return 0.0
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces, dtype=np.int64)
    p0 = points[triangles[:, 0], :3]
    p1 = points[triangles[:, 1], :3]
    p2 = points[triangles[:, 2], :3]
    signed_volume = np.einsum("ij,ij->i", p0, np.cross(p1, p2)) / 6.0
    return float(np.sum(signed_volume))


def _orient_boundary_triangle_facets_to_shell(
    shell_faces: np.ndarray,
    boundary_facets: Sequence[Sequence[int]],
) -> list[list[int]]:
    shell_faces = np.asarray(shell_faces, dtype=np.int64)
    boundary_triangles = [
        np.asarray(facet, dtype=np.int64).reshape(-1)
        for facet in boundary_facets
        if len(np.asarray(facet, dtype=np.int64).reshape(-1)) == 3
    ]
    if (
        shell_faces.ndim != 2
        or shell_faces.shape[1] != 3
        or len(boundary_triangles) != len(boundary_facets)
        or len(boundary_triangles) == 0
    ):
        return [
            np.asarray(facet, dtype=np.int64).reshape(-1).tolist()
            for facet in boundary_facets
        ]

    shell_edge_entries: defaultdict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for face in shell_faces:
        for u, v in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            shell_edge_entries[tuple(sorted((int(u), int(v))))].append((int(u), int(v)))

    shell_boundary_orientations = {
        edge: entries[0]
        for edge, entries in shell_edge_entries.items()
        if len(entries) == 1
    }

    boundary_edges: defaultdict[tuple[int, int], list[tuple[int, tuple[int, int]]]] = defaultdict(list)
    for facet_index, facet in enumerate(boundary_triangles):
        for u, v in ((facet[0], facet[1]), (facet[1], facet[2]), (facet[2], facet[0])):
            boundary_edges[tuple(sorted((int(u), int(v))))].append(
                (facet_index, (int(u), int(v)))
            )

    adjacency: list[list[tuple[int, bool]]] = [[] for _ in boundary_triangles]
    facet_constraints: list[bool | None] = [None] * len(boundary_triangles)

    for edge, entries in boundary_edges.items():
        if len(entries) == 2:
            (left_index, left_edge), (right_index, right_edge) = entries
            stored_same_direction = left_edge == right_edge
            adjacency[left_index].append((right_index, stored_same_direction))
            adjacency[right_index].append((left_index, stored_same_direction))

        shell_edge = shell_boundary_orientations.get(edge)
        if shell_edge is None:
            continue
        for facet_index, facet_edge in entries:
            required_flip = facet_edge == shell_edge
            existing = facet_constraints[facet_index]
            if existing is None:
                facet_constraints[facet_index] = required_flip
            elif existing != required_flip:
                raise ValueError(
                    "Boundary triangle closure has inconsistent shell-edge orientation constraints."
                )

    visited = np.zeros(len(boundary_triangles), dtype=bool)
    flip = np.zeros(len(boundary_triangles), dtype=bool)

    for start in range(len(boundary_triangles)):
        if visited[start]:
            continue
        root_flip = facet_constraints[start]
        if root_flip is None:
            root_flip = False
        flip[start] = bool(root_flip)
        visited[start] = True
        queue: deque[int] = deque([start])

        while queue:
            current = queue.popleft()
            current_flip = bool(flip[current])
            for neighbor, stored_same_direction in adjacency[current]:
                neighbor_flip = current_flip ^ bool(stored_same_direction)
                if not visited[neighbor]:
                    constraint = facet_constraints[neighbor]
                    if constraint is not None and bool(constraint) != bool(neighbor_flip):
                        raise ValueError(
                            "Boundary triangle closure cannot be oriented consistently with the shell."
                        )
                    flip[neighbor] = neighbor_flip
                    visited[neighbor] = True
                    queue.append(neighbor)
                elif bool(flip[neighbor]) != bool(neighbor_flip):
                    raise ValueError(
                        "Boundary triangle closure contains inconsistent shared-edge winding."
                    )

    oriented_boundary_facets: list[list[int]] = []
    for facet, should_flip in zip(boundary_triangles, flip):
        triangle = np.asarray(facet, dtype=np.int64)
        if bool(should_flip):
            triangle = triangle[[0, 2, 1]]
        oriented_boundary_facets.append(triangle.tolist())

    return oriented_boundary_facets


def orient_closed_triangle_plc(
    vertices: np.ndarray,
    shell_faces: np.ndarray,
    boundary_facets: Sequence[Sequence[int]],
) -> tuple[np.ndarray, list[list[int]]]:
    shell_faces = np.asarray(shell_faces, dtype=np.int64)
    oriented_boundary_facets = _orient_boundary_triangle_facets_to_shell(
        shell_faces,
        boundary_facets,
    )
    combined_faces = np.vstack(
        [shell_faces, np.asarray(oriented_boundary_facets, dtype=np.int64)]
    )
    signed_volume = _triangle_signed_volume_sum(vertices, combined_faces)
    if signed_volume >= 0.0:
        return shell_faces, oriented_boundary_facets

    return (
        shell_faces[:, [0, 2, 1]],
        [
            np.asarray(facet, dtype=np.int64)[[0, 2, 1]].tolist()
            for facet in oriented_boundary_facets
        ],
    )


@dataclass
class TetgenPLCDiagnostics:
    num_vertices: int
    num_faces: int
    num_boundary_facets: int
    min_edge_length: float
    median_edge_length: float
    min_face_area: float
    median_face_area: float
    min_triangle_quality: float
    max_triangle_aspect_ratio: float
    degenerate_face_count: int
    duplicate_face_count: int
    nonmanifold_edge_count: int
    open_edge_count: int
    boundary_facets: Dict[str, Dict[str, float | int | bool | str]]
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _mesh_scale(vertices: np.ndarray) -> float:
    if vertices.size == 0:
        return 1.0
    extent = np.ptp(vertices, axis=0)
    scale = float(np.linalg.norm(extent))
    return scale if scale > 0.0 else 1.0


def _geometry_tolerance(vertices: np.ndarray) -> float:
    return max(_mesh_scale(vertices) * 1.0e-12, 1.0e-12)


def _edge_multiplicity(faces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    return np.unique(edges, axis=0, return_counts=True)


def _triangle_areas(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)


def _minimum_nonadjacent_vertex_distance(coords: np.ndarray) -> float:
    if len(coords) < 4:
        return float("inf")

    best = float("inf")
    for i in range(len(coords)):
        for j in range(i + 2, len(coords)):
            if i == 0 and j == len(coords) - 1:
                continue
            best = min(best, float(np.linalg.norm(coords[i] - coords[j])))
    return best


def _project_boundary_facet(name: str, points: np.ndarray) -> tuple[np.ndarray, float]:
    if name in {"south", "north"}:
        return points[:, [0, 2]], float(np.ptp(points[:, 1]))
    if name in {"east", "west"}:
        return points[:, [1, 2]], float(np.ptp(points[:, 0]))
    if name == "top":
        return points[:, [0, 1]], float(np.ptp(points[:, 2]))

    if len(points) >= 3:
        reference = points[0]
        normal = np.zeros(3, dtype=float)
        for i in range(1, len(points) - 1):
            normal += np.cross(points[i] - reference, points[i + 1] - reference)
        if np.linalg.norm(normal) > 0.0:
            drop_axis = int(np.argmax(np.abs(normal)))
            keep_axes = [axis for axis in range(3) if axis != drop_axis]
            return points[:, keep_axes], float(np.ptp(points[:, drop_axis]))

    extent = np.ptp(points, axis=0)
    drop_axis = int(np.argmin(extent))
    keep_axes = [axis for axis in range(3) if axis != drop_axis]
    return points[:, keep_axes], float(extent[drop_axis])


def _normalize_boundary_facets(
    boundary_facets: Mapping[str, Sequence[int]] | Sequence[Sequence[int]],
) -> Dict[str, np.ndarray]:
    if isinstance(boundary_facets, Mapping):
        return {
            str(name): np.asarray(indices, dtype=np.int64)
            for name, indices in boundary_facets.items()
        }
    return {
        f"facet_{i}": np.asarray(indices, dtype=np.int64)
        for i, indices in enumerate(boundary_facets)
    }


def inspect_tetgen_plc(
    vertices: np.ndarray,
    faces: np.ndarray,
    boundary_facets: Mapping[str, Sequence[int]] | Sequence[Sequence[int]],
) -> TetgenPLCDiagnostics:
    """Inspect the PLC handed to TetGen and flag hard-invalid and risky shells."""

    V = np.asarray(vertices, dtype=float)
    F = np.asarray(faces, dtype=np.int64)
    named_facets = _normalize_boundary_facets(boundary_facets)

    errors: list[str] = []
    warnings: list[str] = []
    facet_reports: Dict[str, Dict[str, float | int | bool | str]] = {}

    if V.ndim != 2 or V.shape[1] != 3:
        raise ValueError("vertices must have shape (N, 3)")
    if F.ndim != 2 or F.shape[1] != 3:
        raise ValueError("faces must have shape (M, 3)")

    tol = _geometry_tolerance(V)

    if not np.isfinite(V).all():
        errors.append("Surface shell contains non-finite vertex coordinates.")
    if F.size and ((F < 0).any() or (F >= len(V)).any()):
        errors.append("Surface shell contains face indices outside the vertex range.")

    face_vertex_repeats = np.count_nonzero(
        (F[:, 0] == F[:, 1]) | (F[:, 1] == F[:, 2]) | (F[:, 2] == F[:, 0])
    )
    if face_vertex_repeats:
        errors.append(
            f"Surface shell contains {face_vertex_repeats} faces with repeated vertex indices."
        )

    areas = _triangle_areas(V, F) if len(F) else np.empty(0, dtype=float)
    degenerate_face_count = int(np.count_nonzero(areas <= tol * tol))
    if degenerate_face_count:
        errors.append(
            f"Surface shell contains {degenerate_face_count} degenerate triangles."
        )

    duplicate_face_count = 0
    if len(F):
        sorted_faces = np.sort(F, axis=1)
        _, counts = np.unique(sorted_faces, axis=0, return_counts=True)
        duplicate_face_count = int(np.count_nonzero(counts > 1))
        if duplicate_face_count:
            errors.append(
                f"Surface shell contains {duplicate_face_count} duplicate triangles."
            )

    unique_edges = np.empty((0, 2), dtype=np.int64)
    edge_counts = np.empty(0, dtype=np.int64)
    if len(F):
        unique_edges, edge_counts = _edge_multiplicity(F)

    nonmanifold_edge_count = int(np.count_nonzero(edge_counts > 2))
    open_edge_count = int(np.count_nonzero(edge_counts == 1))
    if nonmanifold_edge_count:
        errors.append(
            f"Surface shell contains {nonmanifold_edge_count} non-manifold edges shared by more than two triangles."
        )

    referenced_vertices = set()
    if F.size:
        referenced_vertices.update(np.unique(F.reshape(-1)).tolist())

    edge_lengths = np.empty(0, dtype=float)
    if len(unique_edges):
        edge_lengths = np.linalg.norm(
            V[unique_edges[:, 0]] - V[unique_edges[:, 1]], axis=1
        )

    min_edge_length = float(edge_lengths.min()) if edge_lengths.size else 0.0
    median_edge_length = float(np.median(edge_lengths)) if edge_lengths.size else 0.0
    min_face_area = float(areas.min()) if areas.size else 0.0
    median_face_area = float(np.median(areas)) if areas.size else 0.0

    triangle_quality = tri_element_quality(V, F) if len(F) else np.empty(0, dtype=float)
    triangle_aspect_ratio = (
        tri_aspect_ratio(V, F) if len(F) else np.empty(0, dtype=float)
    )
    min_triangle_quality = (
        float(triangle_quality.min()) if triangle_quality.size else 1.0
    )
    max_triangle_aspect_ratio = (
        float(triangle_aspect_ratio.max()) if triangle_aspect_ratio.size else 1.0
    )

    reference_edge = median_edge_length if median_edge_length > 0.0 else _mesh_scale(V)
    reference_area = median_face_area if median_face_area > 0.0 else reference_edge**2

    if edge_lengths.size and min_edge_length < reference_edge * _MIN_EDGE_RATIO_WARNING:
        warnings.append(
            "Surface shell contains edges that are orders of magnitude smaller than the median edge length."
        )
    if areas.size and min_face_area < reference_area * _MIN_AREA_RATIO_WARNING:
        warnings.append(
            "Surface shell contains triangles that are orders of magnitude smaller than the median triangle area."
        )
    if triangle_quality.size and min_triangle_quality < _MIN_TRI_QUALITY_WARNING:
        warnings.append(
            f"Surface shell minimum triangle quality is very low ({min_triangle_quality:.3g})."
        )
    if (
        triangle_aspect_ratio.size
        and max_triangle_aspect_ratio > _MAX_TRI_ASPECT_RATIO_WARNING
    ):
        warnings.append(
            f"Surface shell maximum triangle aspect ratio is very high ({max_triangle_aspect_ratio:.3g})."
        )

    for name, indices in named_facets.items():
        indices = np.asarray(indices, dtype=np.int64).reshape(-1)
        if len(indices) < 3:
            errors.append(
                f"Boundary facet '{name}' has fewer than three vertices."
            )
            continue
        if ((indices < 0) | (indices >= len(V))).any():
            errors.append(
                f"Boundary facet '{name}' references vertices outside the valid range."
            )
            continue

        referenced_vertices.update(np.unique(indices).tolist())

        points = V[indices]
        projected, plane_spread = _project_boundary_facet(name, points)
        if len(np.unique(indices)) < 3:
            errors.append(
                f"Boundary facet '{name}' does not contain three distinct vertices."
            )
            continue

        closed_projected = np.vstack([projected, projected[0]])
        projected_edge_lengths = np.linalg.norm(
            np.diff(closed_projected, axis=0), axis=1
        )
        consecutive_duplicates = int(np.count_nonzero(projected_edge_lengths <= tol))
        if consecutive_duplicates:
            errors.append(
                f"Boundary facet '{name}' contains repeated projected vertices."
            )
            continue

        polygon = Polygon(projected)
        projected_area = float(abs(polygon.area))
        valid = bool(polygon.is_valid and projected_area > tol * tol)
        if not valid:
            errors.append(
                f"Boundary facet '{name}' is not a simple projected polygon ({explain_validity(polygon)})."
            )

        min_nonadjacent_distance = _minimum_nonadjacent_vertex_distance(projected)
        if (
            np.isfinite(min_nonadjacent_distance)
            and min_nonadjacent_distance < reference_edge * _BOUNDARY_PINCH_RATIO_WARNING
        ):
            warnings.append(
                f"Boundary facet '{name}' contains very close nonadjacent projected vertices."
            )

        facet_reports[name] = {
            "vertex_count": int(len(indices)),
            "projected_area": projected_area,
            "projected_min_edge_length": float(projected_edge_lengths.min()),
            "projected_min_nonadjacent_distance": float(min_nonadjacent_distance),
            "support_plane_spread": plane_spread,
            "valid": valid,
        }

    unreferenced_vertex_count = len(V) - len(referenced_vertices)
    if unreferenced_vertex_count:
        errors.append(
            f"PLC contains {unreferenced_vertex_count} unreferenced vertices that are not used by any shell face or boundary facet."
        )

    return TetgenPLCDiagnostics(
        num_vertices=int(len(V)),
        num_faces=int(len(F)),
        num_boundary_facets=int(len(named_facets)),
        min_edge_length=min_edge_length,
        median_edge_length=median_edge_length,
        min_face_area=min_face_area,
        median_face_area=median_face_area,
        min_triangle_quality=min_triangle_quality,
        max_triangle_aspect_ratio=max_triangle_aspect_ratio,
        degenerate_face_count=degenerate_face_count,
        duplicate_face_count=duplicate_face_count,
        nonmanifold_edge_count=nonmanifold_edge_count,
        open_edge_count=open_edge_count,
        boundary_facets=facet_reports,
        errors=errors,
        warnings=warnings,
    )


def format_tetgen_plc_diagnostics(diagnostics: TetgenPLCDiagnostics) -> str:
    return (
        "TetGen PLC precheck: "
        f"{diagnostics.num_vertices} vertices, "
        f"{diagnostics.num_faces} faces, "
        f"min_edge={diagnostics.min_edge_length:.3g}, "
        f"median_edge={diagnostics.median_edge_length:.3g}, "
        f"min_area={diagnostics.min_face_area:.3g}, "
        f"median_area={diagnostics.median_face_area:.3g}, "
        f"tri_quality_min={diagnostics.min_triangle_quality:.3g}, "
        f"tri_aspect_max={diagnostics.max_triangle_aspect_ratio:.3g}, "
        f"nonmanifold_edges={diagnostics.nonmanifold_edge_count}, "
        f"open_edges={diagnostics.open_edge_count}"
    )


def get_east_boundary_vertices(vertices, xmax=None, tol=1e-3):
    """
    Return indices on the east boundary (x near ``xmax``), sorted south to north.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    xmax : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 0].max()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted south to north.

    Notes
    -----
    Sorting:
    - 2D: by ``y`` ascending
    - 3D: by ``(y, then z)`` ascending

    Building the east wall polygon with outward normal ``+x``:
    let ``g`` be indices on ground (``z = zmin``), ``r`` the corresponding roof
    indices (``z = zmax``). Use CCW order as seen from ``+x``:
    ``[ g (south to north), r (north to south) ]``.
    """
    V = np.asarray(vertices)
    if xmax is None:
        xmax = V[:, 0].max()
    mask = V[:, 0] >= (xmax - tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            order = np.lexsort((V[idx, 2], V[idx, 1]))  # by y, then z
        else:
            order = np.argsort(V[idx, 1])               # by y
        idx = idx[order]
    return idx


def get_west_boundary_vertices(vertices, ymin=None, ymax=None, xmin=None, tol=1e-3):
    """
    Return indices on the west boundary (x near ``xmin``), sorted north to south.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    xmin : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 0].min()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted north to south.

    Notes
    -----
    Sorting:
    - 2D: by ``y`` descending
    - 3D: by ``(y`` descending, then ``z`` ascending)

    For the west wall (outward normal ``-x``), a CCW ordering seen from ``-x`` is
    ``[ g (north to south), r (south to north) ]``, which yields outward normal
    ``-x``.
    """
    V = np.asarray(vertices)
    if xmin is None:
        xmin = V[:, 0].min()
    mask = V[:, 0] <= (xmin + tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            # sort by y DESC, then z ASC
            order = np.lexsort((V[idx, 2], -V[idx, 1]))
        else:
            order = np.argsort(-V[idx, 1])  # y descending
        idx = idx[order]
    return idx


def get_south_boundary_vertices(vertices, ymin=None, tol=1e-3):
    """
    Return indices on the south boundary (y near ``ymin``), sorted west to east.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    ymin : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 1].min()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted west to east.

    Notes
    -----
    Sorting:
    - 2D: by ``x`` ascending
    - 3D: by ``(x, then z)`` ascending

    For the south wall (outward normal ``-y``), viewed from ``-y``, a CCW ordering is
    ``[ g (west to east), r (east to west) ]``; reversing the roof indices closes
    the quad and keeps the normal pointing outward.
    """
    V = np.asarray(vertices)
    if ymin is None:
        ymin = V[:, 1].min()
    mask = V[:, 1] <= (ymin + tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            order = np.lexsort((V[idx, 2], V[idx, 0]))  # by x, then z
        else:
            order = np.argsort(V[idx, 0])               # by x
        idx = idx[order]
    return idx


def get_north_boundary_vertices(vertices, ymax=None, tol=1e-3):
    """
    Return indices on the north boundary (y near ``ymax``), sorted east to west.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    ymax : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 1].max()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted east to west.

    Notes
    -----
    Sorting:
    - 2D: by ``x`` descending
    - 3D: by ``(x`` descending, then ``z`` ascending)

    For the north wall (outward normal ``+y``), viewed from ``+y``, a CCW ordering is
    ``[ g (east to west), r (west to east) ]``; reversing the roof segment closes
    the quad and yields outward normal ``+y``.
    """
    V = np.asarray(vertices)
    if ymax is None:
        ymax = V[:, 1].max()
    mask = V[:, 1] >= (ymax - tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            # sort by x DESC, then z ASC
            order = np.lexsort((V[idx, 2], -V[idx, 0]))
        else:
            order = np.argsort(-V[idx, 0])  # x descending
        idx = idx[order]
    return idx
def _remove_duplicate_consecutive_indices(
    indices: Sequence[int],
    vertices: np.ndarray,
    tol: float,
) -> list[int]:
    if not indices:
        return []

    deduped: list[int] = [int(indices[0])]
    for index in indices[1:]:
        current = int(index)
        if np.linalg.norm(vertices[current] - vertices[deduped[-1]]) <= tol:
            continue
        deduped.append(current)

    if len(deduped) > 1 and np.linalg.norm(vertices[deduped[0]] - vertices[deduped[-1]]) <= tol:
        deduped.pop()
    return deduped


def _project_boundary_loop(name: str, points: np.ndarray) -> np.ndarray:
    if name in {"south", "north"}:
        return points[:, [0, 2]]
    if name in {"east", "west"}:
        return points[:, [1, 2]]
    return points[:, [0, 1]]


def _compute_top_plane(vertices: np.ndarray, top_height: float) -> tuple[float, float]:
    xmin, ymin, zmin = np.min(vertices, axis=0)
    xmax, ymax, zmax = np.max(vertices, axis=0)
    del xmin, ymin, xmax, ymax

    domain_h = float(zmax - zmin)
    height = float(top_height)
    if height <= domain_h:
        height = 1.5 * domain_h if domain_h > 0 else max(1.0, top_height)

    return float(zmin), float(zmin + height)


def _simplify_boundary_loop(
    name: str,
    indices: Sequence[int],
    vertices: np.ndarray,
    tol: float,
) -> np.ndarray:
    del name
    return np.asarray(
        _remove_duplicate_consecutive_indices(indices, vertices, tol),
        dtype=np.int64,
    )


def _boundary_loops(vertices: np.ndarray, tol: float) -> dict[str, np.ndarray]:
    xmin, ymin, _ = np.min(vertices, axis=0)
    xmax, ymax, _ = np.max(vertices, axis=0)
    return {
        "south": get_south_boundary_vertices(vertices, ymin=ymin, tol=tol),
        "east": get_east_boundary_vertices(vertices, xmax=xmax, tol=tol),
        "north": get_north_boundary_vertices(vertices, ymax=ymax, tol=tol),
        "west": get_west_boundary_vertices(vertices, xmin=xmin, tol=tol),
    }


def _validate_boundary_loop_alignment(
    bottom_vertices: np.ndarray,
    top_source_vertices: np.ndarray,
    bottom_loops: Mapping[str, np.ndarray],
    top_loops: Mapping[str, np.ndarray],
    tol: float,
) -> None:
    xy_tol = max(tol, 1.0e-6)
    for name in ("south", "east", "north", "west"):
        bottom_loop = np.asarray(bottom_loops[name], dtype=np.int64)
        top_loop = np.asarray(top_loops[name], dtype=np.int64)
        if len(bottom_loop) != len(top_loop):
            raise ValueError(
                f"Boundary loop '{name}' length mismatch between surface mesh "
                f"({len(bottom_loop)}) and closure mesh ({len(top_loop)})."
            )
        if len(bottom_loop) == 0:
            continue
        if not np.allclose(
            bottom_vertices[bottom_loop, :2],
            top_source_vertices[top_loop, :2],
            atol=xy_tol,
            rtol=0.0,
        ):
            raise ValueError(
                f"Boundary loop '{name}' in closure mesh does not align with "
                "the surface mesh boundary."
            )


def _realign_closure_boundary_loops(
    bottom_vertices: np.ndarray,
    bottom_loops: Mapping[str, np.ndarray],
    top_source_vertices: np.ndarray,
    top_loops: Mapping[str, np.ndarray],
    tol: float,
) -> dict[str, np.ndarray] | None:
    xy_tol = max(tol, 1.0e-6)
    aligned: dict[str, np.ndarray] = {}
    for name in ("south", "east", "north", "west"):
        bottom_loop = np.asarray(bottom_loops[name], dtype=np.int64)
        top_loop = np.asarray(top_loops[name], dtype=np.int64)
        if len(top_loop) < len(bottom_loop):
            return None

        matched: list[int] = []
        cursor = 0
        for bottom_index in bottom_loop:
            bottom_xy = bottom_vertices[int(bottom_index), :2]
            found = None
            for offset in range(cursor, len(top_loop)):
                candidate_index = int(top_loop[offset])
                if np.allclose(
                    top_source_vertices[candidate_index, :2],
                    bottom_xy,
                    atol=xy_tol,
                    rtol=0.0,
                ):
                    found = offset
                    matched.append(candidate_index)
                    cursor = offset + 1
                    break
            if found is None:
                return None
        aligned[name] = np.asarray(matched, dtype=np.int64)
    return aligned


def _outer_boundary_ring_indices(boundary_loops: Mapping[str, np.ndarray]) -> np.ndarray:
    south = np.asarray(boundary_loops["south"], dtype=np.int64)
    east = np.asarray(boundary_loops["east"], dtype=np.int64)
    north = np.asarray(boundary_loops["north"], dtype=np.int64)
    west = np.asarray(boundary_loops["west"], dtype=np.int64)

    ring = np.concatenate(
        [
            south,
            east[1:],
            north[1:],
            west[1:-1],
        ]
    )
    if ring.size == 0:
        return np.empty((0,), dtype=np.int64)
    deduped = [int(ring[0])]
    for index in ring[1:]:
        current = int(index)
        if current == deduped[-1]:
            continue
        deduped.append(current)
    if len(deduped) > 1 and deduped[0] == deduped[-1]:
        deduped.pop()
    return np.asarray(deduped, dtype=np.int64)


def _build_lifted_top_boundary_vertices(
    boundary_vertices: np.ndarray,
    boundary_loops: Mapping[str, np.ndarray],
    z_top: float,
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    ring_indices = _outer_boundary_ring_indices(boundary_loops)
    if len(ring_indices) < 4:
        raise ValueError("Top-cap boundary must contain at least four vertices.")

    top_vertices = np.asarray(boundary_vertices[ring_indices], dtype=float).copy()
    top_vertices[:, 2] = float(z_top)

    source_to_local = {int(source_index): local_index for local_index, source_index in enumerate(ring_indices)}
    top_loops: dict[str, np.ndarray] = {}
    for name in ("south", "east", "north", "west"):
        source_loop = np.asarray(boundary_loops[name], dtype=np.int64)
        try:
            top_loops[name] = np.asarray(
                [source_to_local[int(source_index)] for source_index in source_loop],
                dtype=np.int64,
            )
        except KeyError as exc:
            raise ValueError(
                f"Boundary loop '{name}' references a vertex outside the top-cap boundary ring."
            ) from exc

    top_ring = np.arange(len(ring_indices), dtype=np.int64)
    return top_vertices, top_loops, top_ring


def _triangulate_boundary_ring_polygon(
    boundary_vertices: np.ndarray,
    ring_indices: np.ndarray,
    tol: float,
) -> np.ndarray:
    ring_indices = np.asarray(ring_indices, dtype=np.int64)
    if len(ring_indices) < 4:
        raise ValueError("Top-cap boundary must contain at least four vertices.")

    ring_points = np.asarray(boundary_vertices[ring_indices], dtype=float)
    polygon = Polygon(ring_points[:, :2])
    if polygon.is_empty or not polygon.is_valid or polygon.area <= tol * tol:
        raise ValueError(
            f"Top-cap boundary polygon is invalid ({explain_validity(polygon)})."
        )

    coord_lookup: dict[tuple[float, float], int] = {}
    for local_index, point in enumerate(ring_points[:, :2]):
        coord_lookup[(round(float(point[0]), 12), round(float(point[1]), 12))] = local_index

    triangles: list[list[int]] = []
    for triangle in triangulate(polygon):
        if not polygon.covers(triangle.representative_point()):
            continue
        coords = np.asarray(triangle.exterior.coords[:-1], dtype=float)
        if len(coords) != 3:
            continue

        tri_indices: list[int] = []
        for coord in coords:
            key = (round(float(coord[0]), 12), round(float(coord[1]), 12))
            local_index = coord_lookup.get(key)
            if local_index is None:
                raise ValueError(
                    "Top-cap triangulation introduced a vertex outside the prescribed boundary ring."
                )
            tri_indices.append(int(ring_indices[local_index]))

        if len(set(tri_indices)) < 3:
            continue

        points = np.asarray(boundary_vertices[np.asarray(tri_indices, dtype=np.int64)], dtype=float)
        if np.cross(points[1] - points[0], points[2] - points[0])[2] < 0.0:
            tri_indices = [tri_indices[0], tri_indices[2], tri_indices[1]]
        triangles.append(tri_indices)

    if not triangles:
        raise ValueError("Top-cap polygon triangulation produced no valid triangles.")

    return np.asarray(triangles, dtype=np.int64)


def _build_top_cap_mesh(
    *,
    boundary_vertices: np.ndarray,
    boundary_loops: Mapping[str, np.ndarray],
    backend: str,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    tol: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    from .backends import available_2d_meshers
    from .flat_mesh_backends import build_city_flat_mesh_from_coverage

    ring_indices = _outer_boundary_ring_indices(boundary_loops)
    if len(ring_indices) < 4:
        raise ValueError("Top-cap boundary must contain at least four vertices.")

    ring_xy = boundary_vertices[ring_indices, :2]
    polygon = Polygon(ring_xy)
    if polygon.is_empty or not polygon.is_valid or polygon.area <= tol * tol:
        raise ValueError(
            f"Top-cap boundary polygon is invalid ({explain_validity(polygon)})."
        )

    minx, miny, maxx, maxy = polygon.bounds

    def _build_with_backend(selected_backend: str) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
        top_cap_mesh = build_city_flat_mesh_from_coverage(
            region_polygons=[polygon],
            region_markers=[0],
            bounds=(float(minx), float(miny), float(maxx), float(maxy)),
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            backend=selected_backend,
            sort_triangles=False,
        )
        top_vertices = np.asarray(top_cap_mesh.vertices, dtype=float)
        top_faces = np.asarray(top_cap_mesh.faces, dtype=np.int64)
        if top_vertices.ndim != 2 or top_vertices.shape[1] != 3:
            raise ValueError("Top-cap mesh vertices must have shape (N, 3).")
        if top_faces.ndim != 2 or top_faces.shape[1] != 3 or top_faces.shape[0] == 0:
            raise ValueError("Top-cap mesh must provide triangle faces.")
        top_loops = _boundary_loops(top_vertices, tol)
        _validate_boundary_loop_alignment(
            boundary_vertices,
            top_vertices,
            boundary_loops,
            top_loops,
            tol,
        )
        return top_vertices, top_faces, top_loops

    def _is_retryable_top_cap_error(exc: Exception) -> bool:
        message = str(exc)
        return (
            "Boundary loop" in message
            or "Triangle support not built" in message
            or "Requested 2D mesher" in message
            or "No supported 2D mesher backend is available" in message
        )

    fallback_backends = [backend]
    for candidate in available_2d_meshers():
        if candidate not in fallback_backends:
            fallback_backends.append(candidate)

    first_retryable_error: Exception | None = None
    for candidate in fallback_backends:
        try:
            return _build_with_backend(candidate)
        except (RuntimeError, ValueError) as exc:
            if not _is_retryable_top_cap_error(exc):
                raise
            if first_retryable_error is None:
                first_retryable_error = exc

    if first_retryable_error is not None:
        raise first_retryable_error
    raise RuntimeError("Failed to build a valid top-cap mesh with any available 2D backend.")


def compute_boundary_triangle_facets(
    mesh: Mesh,
    closure_mesh: Mesh,
    top_height: float = 100.0,
    tol: float = 1e-3,
    top_cap_backend: str = "auto",
    top_cap_max_mesh_size: float | None = None,
    top_cap_min_mesh_angle: float = 25.0,
) -> tuple[np.ndarray, list[list[int]]]:
    """
    Build a fully triangulated PLC closure for TetGen using the 2D ground mesh boundary.

    The surface mesh contributes the terrain/building shell. The closure mesh
    contributes the outer-domain boundary loop ordering. The top cap is then
    re-triangulated from that rectangular boundary, and the four side walls are
    built as stacked triangle strips between matching boundary loops.
    """

    shell_vertices = np.asarray(mesh.vertices, dtype=float)
    cap_source_vertices = np.asarray(closure_mesh.vertices, dtype=float)
    if shell_vertices.ndim != 2 or shell_vertices.shape[1] != 3:
        raise ValueError("Surface mesh vertices must have shape (N, 3).")
    if cap_source_vertices.ndim != 2 or cap_source_vertices.shape[1] != 3:
        raise ValueError("Closure mesh vertices must have shape (N, 3).")

    _, z_top = _compute_top_plane(shell_vertices, top_height)

    bottom_loops = _boundary_loops(shell_vertices, tol)
    cap_source_loops = _boundary_loops(cap_source_vertices, tol)
    try:
        _validate_boundary_loop_alignment(
            shell_vertices,
            cap_source_vertices,
            bottom_loops,
            cap_source_loops,
            tol,
        )
    except ValueError:
        realigned_cap_source_loops = _realign_closure_boundary_loops(
            shell_vertices,
            bottom_loops,
            cap_source_vertices,
            cap_source_loops,
            tol,
        )
        if realigned_cap_source_loops is None:
            raise
        cap_source_loops = realigned_cap_source_loops
        _validate_boundary_loop_alignment(
            shell_vertices,
            cap_source_vertices,
            bottom_loops,
            cap_source_loops,
            tol,
        )

    top_vertices, top_faces, top_loops = _build_top_cap_mesh(
        boundary_vertices=cap_source_vertices,
        boundary_loops=cap_source_loops,
        backend=top_cap_backend,
        max_mesh_size=top_cap_max_mesh_size,
        min_mesh_angle=top_cap_min_mesh_angle,
        tol=tol,
    )
    top_vertices = top_vertices.copy()
    top_vertices[:, 2] = z_top

    offset = shell_vertices.shape[0]
    vertices_out = np.vstack([shell_vertices, top_vertices]).tolist()
    boundary_facets: list[list[int]] = []
    top_boundary_loops = {
        name: np.asarray(top_loops[name], dtype=np.int64) + offset
        for name in ("south", "east", "north", "west")
    }
    sidewall_strip_count = _boundary_sidewall_strip_count(
        np.asarray(vertices_out, dtype=float),
        bottom_loops,
        top_boundary_loops,
    )
    vertical_vertex_cache: dict[tuple[int, int, int, int], int] = {}

    for name in ("south", "east", "north", "west"):
        bottom_loop = np.asarray(bottom_loops[name], dtype=np.int64)
        top_loop = top_boundary_loops[name]
        for b0, b1, t0, t1 in zip(
            bottom_loop[:-1],
            bottom_loop[1:],
            top_loop[:-1],
            top_loop[1:],
        ):
            _append_vertical_sidewall_facet(
                vertices_out,
                boundary_facets,
                bottom_v0=int(b0),
                bottom_v1=int(b1),
                top_v0=int(t0),
                top_v1=int(t1),
                strip_count=sidewall_strip_count,
                vertical_vertex_cache=vertical_vertex_cache,
            )

    for face in top_faces:
        tri = np.asarray(face, dtype=np.int64) + offset
        points = np.asarray(vertices_out, dtype=float)[tri]
        if np.cross(points[1] - points[0], points[2] - points[0])[2] < 0.0:
            tri = np.array([tri[0], tri[2], tri[1]], dtype=np.int64)
        boundary_facets.append([int(tri[0]), int(tri[1]), int(tri[2])])

    return np.asarray(vertices_out, dtype=float), boundary_facets


def compute_oriented_boundary_triangle_plc(
    mesh: Mesh,
    closure_mesh: Mesh,
    top_height: float = 100.0,
    tol: float = 1e-3,
    top_cap_backend: str = "auto",
    top_cap_max_mesh_size: float | None = None,
    top_cap_min_mesh_angle: float = 25.0,
) -> tuple[np.ndarray, np.ndarray, list[list[int]]]:
    vertices_out, boundary_facets = compute_boundary_triangle_facets(
        mesh,
        closure_mesh,
        top_height=top_height,
        tol=tol,
        top_cap_backend=top_cap_backend,
        top_cap_max_mesh_size=top_cap_max_mesh_size,
        top_cap_min_mesh_angle=top_cap_min_mesh_angle,
    )
    oriented_shell_faces, oriented_boundary_facets = orient_closed_triangle_plc(
        vertices_out,
        np.asarray(mesh.faces, dtype=np.int64),
        boundary_facets,
    )
    return vertices_out, oriented_shell_faces, oriented_boundary_facets


def compute_oriented_boundary_plc(
    mesh: Mesh,
    closure_mesh: Mesh,
    top_height: float = 100.0,
    tol: float = 1e-3,
    top_cap_backend: str = "auto",
    top_cap_max_mesh_size: float | None = None,
    top_cap_min_mesh_angle: float = 25.0,
) -> tuple[np.ndarray, np.ndarray, list[list[int]], np.ndarray]:
    """
    Build the PLC used for TetGen with a single polygon top cap.

    The top cap is represented as one polygon facet over the lifted outer boundary
    ring. A triangle-only copy of that cap is returned separately for auditing and
    orientation checks.
    """

    shell_vertices = np.asarray(mesh.vertices, dtype=float)
    shell_faces = np.asarray(mesh.faces, dtype=np.int64)
    cap_source_vertices = np.asarray(closure_mesh.vertices, dtype=float)
    if shell_vertices.ndim != 2 or shell_vertices.shape[1] != 3:
        raise ValueError("Surface mesh vertices must have shape (N, 3).")
    if cap_source_vertices.ndim != 2 or cap_source_vertices.shape[1] != 3:
        raise ValueError("Closure mesh vertices must have shape (N, 3).")

    _, z_top = _compute_top_plane(shell_vertices, top_height)

    bottom_loops = _boundary_loops(shell_vertices, tol)
    cap_source_loops = _boundary_loops(cap_source_vertices, tol)
    try:
        _validate_boundary_loop_alignment(
            shell_vertices,
            cap_source_vertices,
            bottom_loops,
            cap_source_loops,
            tol,
        )
    except ValueError:
        realigned_cap_source_loops = _realign_closure_boundary_loops(
            shell_vertices,
            bottom_loops,
            cap_source_vertices,
            cap_source_loops,
            tol,
        )
        if realigned_cap_source_loops is None:
            raise
        cap_source_loops = realigned_cap_source_loops
        _validate_boundary_loop_alignment(
            shell_vertices,
            cap_source_vertices,
            bottom_loops,
            cap_source_loops,
            tol,
        )

    top_boundary_vertices, top_boundary_loops_local, top_boundary_ring_local = (
        _build_lifted_top_boundary_vertices(
            cap_source_vertices,
            cap_source_loops,
            z_top,
        )
    )
    top_cap_triangles_local = _triangulate_boundary_ring_polygon(
        top_boundary_vertices,
        top_boundary_ring_local,
        tol,
    )

    offset = shell_vertices.shape[0]
    vertices_out = np.vstack([shell_vertices, top_boundary_vertices]).tolist()
    top_boundary_loops = {
        name: np.asarray(top_boundary_loops_local[name], dtype=np.int64) + offset
        for name in ("south", "east", "north", "west")
    }
    top_boundary_ring = np.asarray(top_boundary_ring_local, dtype=np.int64) + offset
    top_cap_triangles = np.asarray(top_cap_triangles_local, dtype=np.int64) + offset

    polygon_sidewall_facets: list[list[int]] = []
    simple_sidewall_triangles: list[list[int]] = []
    for name in ("south", "east", "north", "west"):
        bottom_loop = np.asarray(bottom_loops[name], dtype=np.int64)
        top_loop = top_boundary_loops[name]
        polygon_sidewall_facets.append(
            np.concatenate([bottom_loop, top_loop[::-1]]).astype(np.int64).tolist()
        )
        for b0, b1, t0, t1 in zip(
            bottom_loop[:-1],
            bottom_loop[1:],
            top_loop[:-1],
            top_loop[1:],
        ):
            simple_sidewall_triangles.append([int(b0), int(b1), int(t1)])
            simple_sidewall_triangles.append([int(b0), int(t1), int(t0)])

    use_polygon_sidewalls = _prefer_polygon_sidewalls(
        np.asarray(vertices_out, dtype=float),
        shell_faces,
    )
    if use_polygon_sidewalls:
        sidewall_triangles = simple_sidewall_triangles
    else:
        sidewall_triangles = []
        sidewall_strip_count = _boundary_sidewall_strip_count(
            np.asarray(vertices_out, dtype=float),
            bottom_loops,
            top_boundary_loops,
        )
        vertical_vertex_cache: dict[tuple[int, int, int, int], int] = {}
        for name in ("south", "east", "north", "west"):
            bottom_loop = np.asarray(bottom_loops[name], dtype=np.int64)
            top_loop = top_boundary_loops[name]
            for b0, b1, t0, t1 in zip(
                bottom_loop[:-1],
                bottom_loop[1:],
                top_loop[:-1],
                top_loop[1:],
            ):
                _append_vertical_sidewall_facet(
                    vertices_out,
                    sidewall_triangles,
                    bottom_v0=int(b0),
                    bottom_v1=int(b1),
                    top_v0=int(t0),
                    top_v1=int(t1),
                    strip_count=sidewall_strip_count,
                    vertical_vertex_cache=vertical_vertex_cache,
                )

    oriented_sidewall_triangles = _orient_boundary_triangle_facets_to_shell(
        shell_faces,
        sidewall_triangles,
    )
    combined_faces = np.vstack(
        [
            shell_faces,
            np.asarray(oriented_sidewall_triangles, dtype=np.int64),
            np.asarray(top_cap_triangles, dtype=np.int64),
        ]
    )
    oriented_shell_faces = shell_faces
    oriented_top_cap_ring = top_boundary_ring.tolist()
    oriented_top_cap_triangles = np.asarray(top_cap_triangles, dtype=np.int64)
    if _triangle_signed_volume_sum(np.asarray(vertices_out, dtype=float), combined_faces) < 0.0:
        oriented_shell_faces = shell_faces[:, [0, 2, 1]]
        oriented_sidewall_triangles = [
            np.asarray(facet, dtype=np.int64)[[0, 2, 1]].tolist()
            for facet in oriented_sidewall_triangles
        ]
        polygon_sidewall_facets = [list(reversed(facet)) for facet in polygon_sidewall_facets]
        oriented_top_cap_ring = list(reversed(oriented_top_cap_ring))
        oriented_top_cap_triangles = np.asarray(top_cap_triangles, dtype=np.int64)[:, [0, 2, 1]]

    if use_polygon_sidewalls:
        boundary_facets = list(polygon_sidewall_facets)
    else:
        boundary_facets = [
            list(np.asarray(facet, dtype=np.int64))
            for facet in oriented_sidewall_triangles
        ]
    boundary_facets.append(oriented_top_cap_ring)
    audit_boundary_triangles = np.vstack(
        [
            np.asarray(oriented_sidewall_triangles, dtype=np.int64),
            np.asarray(oriented_top_cap_triangles, dtype=np.int64),
        ]
    )
    return (
        np.asarray(vertices_out, dtype=float),
        np.asarray(oriented_shell_faces, dtype=np.int64),
        boundary_facets,
        audit_boundary_triangles,
    )


def compute_boundary_facets(mesh: Mesh, top_height=100.0, tol=1e-3):
    """
    Build PLC facet polygons for a rectangular box around the mesh domain.

    Polygons are wound CCW as seen from outside so normals point outward
    (right-hand rule):
    - South -> outward ``-y``
    - East  -> outward ``+x``
    - North -> outward ``+y``
    - West  -> outward ``-x``
    - Top   -> outward ``+z``

    Parameters
    ----------
    mesh : Mesh
        Input mesh used to derive the domain bounds.
    top_height : float, optional
        Height of the top cap above ``zmin``; defaults to 100.0.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertices with the four added top points appended at the end.
    dict[str, numpy.ndarray]
        Indices of the five polygons with outward normals:
        ``{"south": [...], "east": [...], "north": [...], "west": [...], "top": [...]}``.
    """
    V = np.asarray(mesh.vertices, dtype=float)
    xmin, ymin, _ = np.min(V, axis=0)
    xmax, ymax, zmax = np.max(V, axis=0)

    # 1) Height check / adjust
    zmin, z_top = _compute_top_plane(V, top_height)

    # 2) Grab boundary indices (sorted for correct ground-edge order)
    east_idx  = get_east_boundary_vertices (V, xmax=xmax, tol=tol)   # south->north
    west_idx  = get_west_boundary_vertices (V, xmin=xmin, tol=tol)   # north->south
    south_idx = get_south_boundary_vertices(V, ymin=ymin, tol=tol)   # west->east
    north_idx = get_north_boundary_vertices(V, ymax=ymax, tol=tol)   # east->west

    # 3) Create 4 top-corner points (appended to the vertex array)
    top_points = np.array([
        [xmin, ymin, z_top],  # t_sw: south-west  (index = N + 0)
        [xmin, ymax, z_top],  # t_nw: north-west  (index = N + 1)
        [xmax, ymin, z_top],  # t_se: south-east  (index = N + 2)
        [xmax, ymax, z_top],  # t_ne: north-east  (index = N + 3)
    ], dtype=float)

    N0 = V.shape[0]
    t_sw, t_nw, t_se, t_ne = N0 + 0, N0 + 1, N0 + 2, N0 + 3
    V_out = np.vstack([V, top_points])

    # 4) Build five polygons with CCW winding as seen from outside

    # SOUTH wall (y = ymin, outward -y): CCW seen from -y
    #   ground: west->east  then top: east->west (reverse)
    south_poly = list(south_idx) + [t_se, t_sw]

    # EAST wall (x = xmax, outward +x): CCW seen from +x
    #   ground: south->north then top: north->south (reverse)
    east_poly  = list(east_idx)  + [t_ne, t_se]

    # NORTH wall (y = ymax, outward +y): CCW seen from +y
    #   ground: east->west then top: west->east
    north_poly = list(north_idx) + [t_nw, t_ne]

    # WEST wall (x = xmin, outward -x): CCW seen from -x
    #   ground: north->south then top: south->north
    west_poly  = list(west_idx)  + [t_sw, t_nw]

    # TOP cap (z = z_top, outward +z): CCW in XY as seen from above (+z)
    #   CCW rectangle: (xmin,ymin)->(xmax,ymin)->(xmax,ymax)->(xmin,ymax)
    top_poly   = [t_sw, t_se, t_ne, t_nw]

    facets = {
        "south": _simplify_boundary_loop("south", south_poly, V_out, tol),
        "east": _simplify_boundary_loop("east", east_poly, V_out, tol),
        "north": _simplify_boundary_loop("north", north_poly, V_out, tol),
        "west": _simplify_boundary_loop("west", west_poly, V_out, tol),
        "top": _simplify_boundary_loop("top", top_poly, V_out, tol),
    }

    return V_out, facets
