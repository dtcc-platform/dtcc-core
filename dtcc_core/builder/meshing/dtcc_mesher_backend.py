from __future__ import annotations

import importlib

import numpy as np

from shapely.geometry import Polygon

from ...model import Mesh, Surface


_DTCC_MESHER_MODULE = None


def _load_dtcc_mesher():
    global _DTCC_MESHER_MODULE
    if _DTCC_MESHER_MODULE is None:
        _DTCC_MESHER_MODULE = importlib.import_module("dtcc_mesher")
    return _DTCC_MESHER_MODULE


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

    outer = project(surface.vertices)
    holes = [project(hole) for hole in surface.holes if len(hole) >= 3]

    return outer, holes, (origin, tangent, bitangent, normal)


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


def _meshing_options(
    dtcc_mesher,
    *,
    min_mesh_angle: float,
    max_edge_length: float | None,
):
    return dtcc_mesher.MeshingOptions(
        min_angle=min_mesh_angle,
        max_edge_length=max_edge_length,
        refine=False,
    )


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
    max_edge_length = triangle_size if triangle_size is not None and triangle_size > 0 else None

    raw_mesh = dtcc_mesher.mesh(
        dtcc_mesher.Domain.from_loops(outer_loop, hole_loops),
        options=_meshing_options(
            dtcc_mesher,
            min_mesh_angle=min_mesh_angle,
            max_edge_length=max_edge_length,
        ),
    )

    vertices = _lift_vertices(np.asarray(raw_mesh.points, dtype=np.float64), transform)
    _, _, _, normal = transform
    faces = _orient_faces(
        vertices,
        np.asarray(raw_mesh.triangles, dtype=np.int64),
        normal,
    )

    return Mesh(vertices=vertices, faces=faces)


def build_city_flat_mesh_with_dtcc_mesher(
    *,
    region_polygons: list[Polygon],
    region_markers: list[int],
    max_mesh_size: float,
    min_mesh_angle: float,
) -> Mesh:
    dtcc_mesher = _load_dtcc_mesher()
    max_edge_length = max_mesh_size if max_mesh_size > 0 else None

    if len(region_polygons) != len(region_markers):
        raise ValueError("region_markers length must match region_polygons length")
    if not region_polygons:
        return Mesh()

    raw_mesh = dtcc_mesher.mesh(
        dtcc_mesher.Coverage(region_polygons, region_markers),
        options=_meshing_options(
            dtcc_mesher,
            min_mesh_angle=min_mesh_angle,
            max_edge_length=max_edge_length,
        ),
    )

    vertices_2d = np.asarray(raw_mesh.points, dtype=np.float64)
    vertices = np.column_stack(
        [vertices_2d, np.zeros(len(vertices_2d), dtype=np.float64)]
    )
    faces = np.asarray(raw_mesh.triangles, dtype=np.int64)
    markers = (
        np.asarray(raw_mesh.markers, dtype=np.int64)
        if raw_mesh.markers is not None
        else np.empty((0,), dtype=np.int64)
    )
    return Mesh(vertices=vertices, faces=faces, markers=markers)
