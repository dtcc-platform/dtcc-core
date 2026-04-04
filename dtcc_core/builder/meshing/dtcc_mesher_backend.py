from __future__ import annotations

import importlib
import math

import numpy as np
from shapely.geometry import GeometryCollection, MultiPolygon, Point, Polygon, box
from shapely.ops import unary_union
from shapely.prepared import prep

from .. import _dtcc_builder
from ..logging import warning
from ..model_conversion import builder_mesh_to_mesh, create_builder_surface
from ...model import Mesh, Surface


_DTCC_MESHER_MODULE = None


def _load_dtcc_mesher():
    global _DTCC_MESHER_MODULE
    if _DTCC_MESHER_MODULE is None:
        _DTCC_MESHER_MODULE = importlib.import_module("dtcc_mesher")
    return _DTCC_MESHER_MODULE


def _sanitize_loop(points: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    if len(points) == 0:
        return np.empty((0, 2), dtype=np.float64)

    clean: list[np.ndarray] = []
    for point in np.asarray(points, dtype=np.float64):
        if not clean or np.linalg.norm(point - clean[-1]) > tol:
            clean.append(point)

    if len(clean) > 1 and np.linalg.norm(clean[0] - clean[-1]) <= tol:
        clean.pop()

    return np.asarray(clean, dtype=np.float64)


class _PlanarDomainBuilder:
    def __init__(self, tolerance: float = 1e-9):
        self._scale = 1.0 / tolerance
        self._points: list[np.ndarray] = []
        self._segments: list[tuple[int, int]] = []
        self._point_index: dict[tuple[int, int], int] = {}
        self._segment_index: set[tuple[int, int]] = set()

    def _key(self, point: np.ndarray) -> tuple[int, int]:
        return (
            int(round(float(point[0]) * self._scale)),
            int(round(float(point[1]) * self._scale)),
        )

    def add_point(self, point: np.ndarray, *, reuse_existing: bool = True) -> int:
        key = self._key(point)
        if reuse_existing and key in self._point_index:
            return self._point_index[key]

        index = len(self._points)
        self._points.append(np.asarray(point, dtype=np.float64))
        if reuse_existing:
            self._point_index[key] = index
        return index

    def add_segment(self, start: int, end: int) -> None:
        if start == end:
            return

        canonical = (start, end) if start < end else (end, start)
        if canonical in self._segment_index:
            return

        self._segment_index.add(canonical)
        self._segments.append((start, end))

    def add_loop(self, loop: np.ndarray) -> None:
        clean_loop = _sanitize_loop(loop)
        if len(clean_loop) < 3:
            return

        # Keep coincident vertices from different loops distinct. dtcc_mesher's
        # PSLG validator requires loop vertices to have degree 2, so collapsing
        # shared corner coordinates across separate rings can create invalid
        # degree-4 nodes even when the geometry is otherwise valid.
        indices = [self.add_point(point, reuse_existing=False) for point in clean_loop]
        for index, next_index in zip(indices, indices[1:] + indices[:1]):
            self.add_segment(index, next_index)

    def to_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        points = np.asarray(self._points, dtype=np.float64)
        segments = np.asarray(self._segments, dtype=np.uint32)

        if points.size == 0:
            points = np.empty((0, 2), dtype=np.float64)
        if segments.size == 0:
            segments = np.empty((0, 2), dtype=np.uint32)

        return points, segments


def _polygon_interior_seeds(polygon: Polygon, spacing: float | None) -> list[np.ndarray]:
    if spacing is None or spacing <= 0:
        return []

    minx, miny, maxx, maxy = polygon.bounds
    if maxx - minx <= spacing or maxy - miny <= spacing:
        return []

    prepared = prep(polygon)
    vertical_step = spacing * math.sqrt(3.0) / 2.0
    seeds: list[np.ndarray] = []

    row = 0
    y = miny + 0.5 * vertical_step
    while y < maxy:
        x_offset = 0.5 * spacing if row % 2 else 0.0
        x = minx + 0.5 * spacing + x_offset
        while x < maxx:
            point = Point(x, y)
            if prepared.contains(point):
                seeds.append(np.array([x, y], dtype=np.float64))
            x += spacing
        y += vertical_step
        row += 1

    return seeds


def _subdivide_loop(loop: np.ndarray, max_edge_length: float | None) -> np.ndarray:
    clean_loop = _sanitize_loop(loop)
    if len(clean_loop) < 3 or max_edge_length is None or max_edge_length <= 0:
        return clean_loop

    subdivided: list[np.ndarray] = []
    for start, end in zip(clean_loop, np.vstack([clean_loop[1:], clean_loop[:1]])):
        edge = end - start
        length = float(np.linalg.norm(edge))
        segment_count = max(1, int(math.ceil(length / max_edge_length)))

        if not subdivided:
            subdivided.append(start)

        for step in range(1, segment_count):
            t = step / segment_count
            subdivided.append(start + t * edge)
        subdivided.append(end)

    if len(subdivided) > 1 and np.allclose(subdivided[0], subdivided[-1]):
        subdivided.pop()

    return np.asarray(subdivided, dtype=np.float64)


def _surface_basis(surface: Surface) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if surface.normal.shape != (3,):
        surface.calculate_normal()

    normal = np.asarray(surface.normal, dtype=np.float64)
    normal /= np.linalg.norm(normal)

    origin = np.asarray(surface.centroid, dtype=np.float64)

    tangent = None
    for index in range(len(surface.vertices)):
        start = np.asarray(surface.vertices[index], dtype=np.float64)
        end = np.asarray(surface.vertices[(index + 1) % len(surface.vertices)], dtype=np.float64)
        candidate = end - start
        candidate = candidate - np.dot(candidate, normal) * normal
        if np.linalg.norm(candidate) > 1e-10:
            tangent = candidate / np.linalg.norm(candidate)
            break

    if tangent is None:
        reference = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if abs(np.dot(reference, normal)) > 0.9:
            reference = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        tangent = np.cross(normal, reference)
        tangent /= np.linalg.norm(tangent)

    bitangent = np.cross(normal, tangent)
    bitangent /= np.linalg.norm(bitangent)

    return origin, tangent, bitangent, normal


def _project_surface(surface: Surface) -> tuple[np.ndarray, list[np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    origin, tangent, bitangent, normal = _surface_basis(surface)

    def project(points: np.ndarray) -> np.ndarray:
        vectors = np.asarray(points, dtype=np.float64) - origin
        return np.column_stack(
            [vectors @ tangent, vectors @ bitangent]
        ).astype(np.float64, copy=False)

    outer = _sanitize_loop(project(surface.vertices))
    holes = [_sanitize_loop(project(hole)) for hole in surface.holes if len(hole) >= 3]
    holes = [hole for hole in holes if len(hole) >= 3]

    return outer, holes, (origin, tangent, bitangent, normal)


def _build_surface_domain(
    outer_loop: np.ndarray,
    hole_loops: list[np.ndarray],
    triangle_size: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    builder = _PlanarDomainBuilder()
    outer_loop = _subdivide_loop(outer_loop, triangle_size)
    hole_loops = [_subdivide_loop(loop, triangle_size) for loop in hole_loops]

    builder.add_loop(outer_loop)
    for hole_loop in hole_loops:
        builder.add_loop(hole_loop)

    polygon = Polygon(outer_loop, [loop for loop in hole_loops])
    for seed in _polygon_interior_seeds(polygon, triangle_size):
        builder.add_point(seed)

    hole_points = None
    if hole_loops:
        hole_points = np.asarray(
            [np.array(Polygon(loop).representative_point().coords[0]) for loop in hole_loops],
            dtype=np.float64,
        )

    return (*builder.to_arrays(), hole_points)


def _lift_vertices(
    vertices_2d: np.ndarray,
    transform: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    origin, tangent, bitangent, _normal = transform
    return origin + np.outer(vertices_2d[:, 0], tangent) + np.outer(vertices_2d[:, 1], bitangent)


def _orient_faces(vertices: np.ndarray, faces: np.ndarray, normal: np.ndarray) -> np.ndarray:
    if len(faces) == 0:
        return faces

    face_normals = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]],
        vertices[faces[:, 2]] - vertices[faces[:, 0]],
    )
    flip = np.einsum("ij,j->i", face_normals, normal) < 0
    if np.any(flip):
        oriented = faces.copy()
        oriented[flip, 1] = faces[flip, 2]
        oriented[flip, 2] = faces[flip, 1]
        return oriented
    return faces


def mesh_surface_with_dtcc_mesher(
    surface: Surface,
    *,
    triangle_size: float | None,
    min_mesh_angle: float,
) -> Mesh:
    if len(surface.vertices) < 3:
        return Mesh()

    dtcc_mesher = _load_dtcc_mesher()
    outer_loop, hole_loops, transform = _project_surface(surface)
    points, segments, hole_points = _build_surface_domain(
        outer_loop,
        hole_loops,
        triangle_size if triangle_size is not None and triangle_size > 0 else None,
    )

    raw_mesh = dtcc_mesher.generate(
        points,
        segments=segments if len(segments) > 0 else None,
        holes=hole_points,
        min_angle=min_mesh_angle,
        refine=triangle_size is not None and triangle_size > 0,
    )

    vertices = _lift_vertices(np.asarray(raw_mesh.points, dtype=np.float64), transform)
    _, _, _, normal = transform
    faces = _orient_faces(
        vertices,
        np.asarray(raw_mesh.triangles, dtype=np.int64),
        normal,
    )

    return Mesh(vertices=vertices, faces=faces)


def _polygon_loops(polygon: Polygon) -> tuple[np.ndarray, list[np.ndarray]]:
    outer = np.asarray(polygon.exterior.coords[:-1], dtype=np.float64)
    holes = [np.asarray(ring.coords[:-1], dtype=np.float64) for ring in polygon.interiors]
    return outer, holes


def _build_flat_domain(
    *,
    bounds: tuple[float, float, float, float],
    hole_polygons: list[Polygon],
    max_mesh_size: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    xmin, ymin, xmax, ymax = bounds
    outer_loop = np.array(
        [
            [xmin, ymin],
            [xmax, ymin],
            [xmax, ymax],
            [xmin, ymax],
        ],
        dtype=np.float64,
    )

    builder = _PlanarDomainBuilder()
    builder.add_loop(_subdivide_loop(outer_loop, max_mesh_size))

    hole_points: list[np.ndarray] = []
    excluded_loops: list[np.ndarray] = []
    for polygon in hole_polygons:
        if polygon.is_empty:
            continue
        outer, holes = _polygon_loops(polygon)
        excluded_loops.append(_subdivide_loop(outer, max_mesh_size))
        excluded_loops.extend(_subdivide_loop(hole, max_mesh_size) for hole in holes)

    for loop in excluded_loops:
        builder.add_loop(loop)
        hole_points.append(np.array(Polygon(loop).representative_point().coords[0]))

    domain_polygon = Polygon(
        outer_loop,
        [loop for loop in excluded_loops],
    )
    for seed in _polygon_interior_seeds(domain_polygon, max_mesh_size):
        builder.add_point(seed)

    points, segments = builder.to_arrays()
    if not hole_points:
        return points, segments, None

    return points, segments, np.asarray(hole_points, dtype=np.float64)


def _mesh_planar_polygon(
    dtcc_mesher,
    polygon: Polygon,
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
) -> Mesh:
    if polygon.is_empty:
        return Mesh()

    normalized_geometry = polygon
    if polygon.interiors:
        shell_polygon = Polygon(np.asarray(polygon.exterior.coords, dtype=np.float64))
        hole_geometries = [
            Polygon(np.asarray(ring.coords, dtype=np.float64))
            for ring in polygon.interiors
        ]
        normalized_geometry = shell_polygon.difference(unary_union(hole_geometries))

    meshes: list[Mesh] = []
    for component in _iter_polygons(normalized_geometry):
        meshes.append(
            _mesh_planar_polygon_component(
                dtcc_mesher,
                component,
                max_mesh_size=max_mesh_size,
                min_mesh_angle=min_mesh_angle,
            )
        )

    return _merge_planar_meshes(meshes)


def _mesh_planar_polygon_component(
    dtcc_mesher,
    polygon: Polygon,
    *,
    max_mesh_size: float | None,
    min_mesh_angle: float,
) -> Mesh:
    outer_loop, hole_loops = _polygon_loops(polygon)
    points, segments, hole_points = _build_surface_domain(
        outer_loop,
        hole_loops,
        max_mesh_size,
    )

    try:
        raw_mesh = dtcc_mesher.generate(
            points,
            segments=segments if len(segments) > 0 else None,
            holes=hole_points,
            min_angle=min_mesh_angle,
            refine=False,
        )
    except RuntimeError as exc:
        warning(
            "dtcc_mesher could not mesh one planar polygon component; "
            "falling back to builder earcut without refinement."
        )
        return _mesh_planar_polygon_with_builder_earcut(
            polygon,
            min_mesh_angle=min_mesh_angle,
            original_error=exc,
        )

    vertices_2d = np.asarray(raw_mesh.points, dtype=np.float64)
    vertices = np.column_stack(
        [vertices_2d, np.zeros(len(vertices_2d), dtype=np.float64)]
    )
    faces = np.asarray(raw_mesh.triangles, dtype=np.int64)
    return Mesh(vertices=vertices, faces=faces)


def _mesh_planar_polygon_with_builder_earcut(
    polygon: Polygon,
    *,
    min_mesh_angle: float,
    original_error: RuntimeError,
) -> Mesh:
    surface = Surface()
    surface.from_polygon(polygon, 0.0)

    try:
        builder_surface = create_builder_surface(surface)
        builder_mesh = _dtcc_builder.mesh_surface(builder_surface, -1, min_mesh_angle)
    except RuntimeError as fallback_error:
        raise RuntimeError(
            "dtcc_mesher failed to triangulate a planar polygon component, and "
            "the builder earcut fallback was also unavailable."
        ) from fallback_error

    mesh = builder_mesh_to_mesh(builder_mesh)
    mesh.vertices = np.asarray(mesh.vertices, dtype=np.float64)
    mesh.faces = np.asarray(mesh.faces, dtype=np.int64)
    return mesh


def _merge_planar_meshes(meshes: list[Mesh], tolerance: float = 1e-9) -> Mesh:
    if not meshes:
        return Mesh()

    scale = 1.0 / tolerance
    vertex_index: dict[tuple[int, int, int], int] = {}
    vertices: list[np.ndarray] = []
    faces: list[np.ndarray] = []
    markers: list[int] = []

    for mesh in meshes:
        local_to_global: dict[int, int] = {}
        for local_index, vertex in enumerate(np.asarray(mesh.vertices, dtype=np.float64)):
            key = (
                int(round(float(vertex[0]) * scale)),
                int(round(float(vertex[1]) * scale)),
                int(round(float(vertex[2]) * scale)),
            )
            if key not in vertex_index:
                vertex_index[key] = len(vertices)
                vertices.append(vertex)
            local_to_global[local_index] = vertex_index[key]

        mesh_markers = np.asarray(mesh.markers, dtype=np.int64)
        if mesh_markers.size == 0:
            mesh_markers = np.full(len(mesh.faces), -2, dtype=np.int64)

        for face_index, face in enumerate(np.asarray(mesh.faces, dtype=np.int64)):
            faces.append(
                np.array(
                    [local_to_global[int(face[0])], local_to_global[int(face[1])], local_to_global[int(face[2])]],
                    dtype=np.int64,
                )
            )
            markers.append(int(mesh_markers[face_index]))

    return Mesh(
        vertices=np.asarray(vertices, dtype=np.float64),
        faces=np.asarray(faces, dtype=np.int64),
        markers=np.asarray(markers, dtype=np.int64),
    )


def _add_halo_markers(mesh: Mesh) -> Mesh:
    if len(mesh.faces) == 0:
        return mesh

    markers = np.asarray(mesh.markers, dtype=np.int64).copy()
    is_building_vertex = np.zeros(len(mesh.vertices), dtype=bool)

    for face, marker in zip(mesh.faces, markers):
        if marker >= 0:
            is_building_vertex[face] = True

    for face_index, face in enumerate(mesh.faces):
        if markers[face_index] == -2 and np.any(is_building_vertex[face]):
            markers[face_index] = -1

    mesh.markers = markers
    return mesh


def _iter_polygons(geometry) -> list[Polygon]:
    if geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return [polygon for polygon in geometry.geoms if not polygon.is_empty]
    if isinstance(geometry, GeometryCollection):
        polygons: list[Polygon] = []
        for item in geometry.geoms:
            polygons.extend(_iter_polygons(item))
        return polygons
    return []


def build_city_flat_mesh_with_dtcc_mesher(
    *,
    building_polygons: list[Polygon],
    building_markers: list[int] | None = None,
    hole_polygons: list[Polygon],
    bounds: tuple[float, float, float, float],
    max_mesh_size: float,
    min_mesh_angle: float,
) -> Mesh:
    dtcc_mesher = _load_dtcc_mesher()
    ground_domain = box(*bounds)
    excluded_polygons = [polygon for polygon in [*building_polygons, *hole_polygons] if not polygon.is_empty]
    if excluded_polygons:
        ground_domain = ground_domain.difference(unary_union(excluded_polygons))

    meshes: list[Mesh] = []
    for ground_polygon in _iter_polygons(ground_domain):
        ground_mesh = _mesh_planar_polygon(
            dtcc_mesher,
            ground_polygon,
            max_mesh_size=max_mesh_size if max_mesh_size > 0 else None,
            min_mesh_angle=min_mesh_angle,
        )
        if len(ground_mesh.faces) == 0:
            continue
        ground_mesh.markers = np.full(
            len(ground_mesh.faces),
            -2,
            dtype=np.int64,
        )
        meshes.append(ground_mesh)

    if building_markers is not None and len(building_markers) != len(building_polygons):
        raise ValueError("building_markers length must match building_polygons length")

    for polygon_index, polygon in enumerate(building_polygons):
        building_mesh = _mesh_planar_polygon(
            dtcc_mesher,
            polygon,
            max_mesh_size=max_mesh_size if max_mesh_size > 0 else None,
            min_mesh_angle=min_mesh_angle,
        )
        if len(building_mesh.faces) == 0:
            continue
        marker = (
            int(building_markers[polygon_index])
            if building_markers is not None
            else polygon_index
        )
        building_mesh.markers = np.full(
            len(building_mesh.faces),
            marker,
            dtype=np.int64,
        )
        meshes.append(building_mesh)

    merged = _merge_planar_meshes(meshes)
    return _add_halo_markers(merged)
