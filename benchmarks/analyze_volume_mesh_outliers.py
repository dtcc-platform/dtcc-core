#!/usr/bin/env python3
"""
Inspect the worst tetrahedra in a DTCC volume mesh and export local neighborhoods.

This script is meant for debugging the remaining 3D quality tail after the
upstream footprint conditioning / 2D meshing stages are already trusted.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Iterable

import meshio
import numpy as np

from dtcc_core.io.meshes import load_volume_mesh
from dtcc_core.model import VolumeMesh
from dtcc_core.model.mixins.mesh.quality import (
    tet_aspect_ratio,
    tet_element_quality,
    tet_radius_ratio,
)


def _sorted_face(face: np.ndarray) -> tuple[int, int, int]:
    a, b, c = sorted(map(int, face.tolist()))
    return (a, b, c)


def _tet_faces(cell: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return (
        cell[[1, 2, 3]],
        cell[[0, 2, 3]],
        cell[[0, 1, 3]],
        cell[[0, 1, 2]],
    )


def _build_face_maps(
    cells: np.ndarray,
    boundary_faces: np.ndarray,
    boundary_markers: np.ndarray,
) -> tuple[dict[tuple[int, int, int], list[int]], dict[tuple[int, int, int], int]]:
    face_to_cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for cell_index, cell in enumerate(cells):
        for face in _tet_faces(cell):
            face_to_cells[_sorted_face(face)].append(int(cell_index))

    boundary_marker_map = {
        _sorted_face(face): int(marker)
        for face, marker in zip(boundary_faces, boundary_markers)
    }
    return face_to_cells, boundary_marker_map


def _build_cell_neighbors(cells: np.ndarray) -> list[list[int]]:
    face_to_cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for cell_index, cell in enumerate(cells):
        for face in _tet_faces(cell):
            face_to_cells[_sorted_face(face)].append(int(cell_index))

    neighbors: list[set[int]] = [set() for _ in range(len(cells))]
    for owners in face_to_cells.values():
        if len(owners) < 2:
            continue
        for owner in owners:
            neighbors[owner].update(other for other in owners if other != owner)
    return [sorted(cell_neighbors) for cell_neighbors in neighbors]


def _threshold_components(
    selected_cells: np.ndarray,
    neighbors: list[list[int]],
) -> list[list[int]]:
    selected = set(map(int, selected_cells.tolist()))
    components: list[list[int]] = []
    visited: set[int] = set()

    for start in selected:
        if start in visited:
            continue
        component: list[int] = []
        queue: deque[int] = deque([start])
        visited.add(start)
        while queue:
            cell_index = queue.popleft()
            component.append(cell_index)
            for neighbor in neighbors[cell_index]:
                if neighbor not in selected or neighbor in visited:
                    continue
                visited.add(neighbor)
                queue.append(neighbor)
        components.append(sorted(component))

    return sorted(components, key=len, reverse=True)


def _n_hop_neighborhood(
    seeds: Iterable[int],
    neighbors: list[list[int]],
    hops: int,
) -> list[int]:
    visited = set(map(int, seeds))
    frontier = set(visited)
    for _ in range(max(int(hops), 0)):
        next_frontier: set[int] = set()
        for cell_index in frontier:
            next_frontier.update(neighbors[cell_index])
        next_frontier -= visited
        visited.update(next_frontier)
        frontier = next_frontier
        if not frontier:
            break
    return sorted(visited)


def _point_to_triangles_distance(point: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    # Vectorized point-triangle distance based on closest-point queries.
    a = triangles[:, 0]
    b = triangles[:, 1]
    c = triangles[:, 2]

    ab = b - a
    ac = c - a
    ap = point - a

    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)

    distances = np.full(len(triangles), np.inf, dtype=float)

    mask = (d1 <= 0.0) & (d2 <= 0.0)
    distances[mask] = np.linalg.norm(ap[mask], axis=1)

    bp = point - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    mask = (d3 >= 0.0) & (d4 <= d3)
    distances[mask] = np.linalg.norm(bp[mask], axis=1)

    vc = d1 * d4 - d3 * d2
    mask = (vc <= 0.0) & (d1 >= 0.0) & (d3 <= 0.0)
    if np.any(mask):
        v = d1[mask] / (d1[mask] - d3[mask])
        projection = a[mask] + v[:, None] * ab[mask]
        distances[mask] = np.linalg.norm(point - projection, axis=1)

    cp = point - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    mask = (d6 >= 0.0) & (d5 <= d6)
    distances[mask] = np.linalg.norm(cp[mask], axis=1)

    vb = d5 * d2 - d1 * d6
    mask = (vb <= 0.0) & (d2 >= 0.0) & (d6 <= 0.0)
    if np.any(mask):
        w = d2[mask] / (d2[mask] - d6[mask])
        projection = a[mask] + w[:, None] * ac[mask]
        distances[mask] = np.linalg.norm(point - projection, axis=1)

    va = d3 * d6 - d5 * d4
    mask = (va <= 0.0) & ((d4 - d3) >= 0.0) & ((d5 - d6) >= 0.0)
    if np.any(mask):
        denom = (d4[mask] - d3[mask]) + (d5[mask] - d6[mask])
        w = (d4[mask] - d3[mask]) / denom
        projection = b[mask] + w[:, None] * (c[mask] - b[mask])
        distances[mask] = np.linalg.norm(point - projection, axis=1)

    mask = ~np.isfinite(distances)
    if np.any(mask):
        denom = va[mask] + vb[mask] + vc[mask]
        denom = np.where(np.abs(denom) > 0.0, denom, 1.0)
        v = vb[mask] / denom
        w = vc[mask] / denom
        projection = a[mask] + ab[mask] * v[:, None] + ac[mask] * w[:, None]
        distances[mask] = np.linalg.norm(point - projection, axis=1)

    return distances


def _extract_cell_submesh(
    vertices: np.ndarray,
    cells: np.ndarray,
    selected_cells: list[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = np.asarray(sorted(set(selected_cells)), dtype=np.int64)
    used_vertices, inverse = np.unique(cells[selected].ravel(), return_inverse=True)
    sub_vertices = vertices[used_vertices]
    sub_cells = inverse.reshape(len(selected), 4)
    return sub_vertices, sub_cells, selected


def _extract_boundary_surface(
    vertices: np.ndarray,
    boundary_faces: np.ndarray,
    boundary_markers: np.ndarray,
    face_to_cells: dict[tuple[int, int, int], list[int]],
    selected_cells: set[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    keep_indices: list[int] = []
    for face_index, face in enumerate(boundary_faces):
        owners = face_to_cells.get(_sorted_face(face), [])
        if any(owner in selected_cells for owner in owners):
            keep_indices.append(face_index)

    if not keep_indices:
        return (
            np.empty((0, 3), dtype=float),
            np.empty((0, 3), dtype=np.int64),
            np.empty((0,), dtype=np.int64),
        )

    kept_faces = boundary_faces[np.asarray(keep_indices, dtype=np.int64)]
    kept_markers = boundary_markers[np.asarray(keep_indices, dtype=np.int64)]
    used_vertices, inverse = np.unique(kept_faces.ravel(), return_inverse=True)
    sub_vertices = vertices[used_vertices]
    sub_faces = inverse.reshape(len(keep_indices), 3)
    return sub_vertices, sub_faces, kept_markers


def _boundary_orientation_summary(
    vertices: np.ndarray,
    faces: np.ndarray,
    markers: np.ndarray,
) -> dict[str, dict[str, int]]:
    if len(faces) == 0:
        return {}

    face_points = vertices[faces]
    normals = np.cross(
        face_points[:, 1] - face_points[:, 0],
        face_points[:, 2] - face_points[:, 0],
    )
    normal_norms = np.linalg.norm(normals, axis=1)
    normal_norms = np.where(normal_norms > 0.0, normal_norms, 1.0)
    z_alignment = np.abs(normals[:, 2]) / normal_norms

    summary: dict[str, dict[str, int]] = {}
    for marker in np.unique(markers):
        mask = markers == marker
        summary[str(int(marker))] = {
            "vertical": int(np.count_nonzero(z_alignment[mask] < 0.2)),
            "horizontal": int(np.count_nonzero(z_alignment[mask] > 0.8)),
            "mixed": int(np.count_nonzero((z_alignment[mask] >= 0.2) & (z_alignment[mask] <= 0.8))),
        }
    return summary


def _component_boundary_markers(
    component_cells: list[int],
    cells: np.ndarray,
    boundary_marker_map: dict[tuple[int, int, int], int],
) -> Counter[int]:
    markers: Counter[int] = Counter()
    for cell_index in component_cells:
        for face in _tet_faces(cells[cell_index]):
            marker = boundary_marker_map.get(_sorted_face(face))
            if marker is not None:
                markers[int(marker)] += 1
    return markers


def analyze_volume_mesh(
    mesh: VolumeMesh,
    *,
    top_k: int,
    aspect_threshold: float,
    neighborhood_hops: int,
    output_dir: Path,
) -> dict:
    vertices = np.asarray(mesh.vertices[:, :3], dtype=np.float64)
    cells = np.asarray(mesh.cells, dtype=np.int64)
    boundary_faces = np.asarray(getattr(mesh, "boundary_faces", []), dtype=np.int64)
    boundary_markers = np.asarray(getattr(mesh, "boundary_markers", []), dtype=np.int64)

    aspect_ratio = tet_aspect_ratio(vertices, cells)
    radius_ratio = tet_radius_ratio(vertices, cells)
    element_quality = tet_element_quality(vertices, cells)
    centroids = vertices[cells].mean(axis=1)

    neighbors = _build_cell_neighbors(cells)
    face_to_cells, boundary_marker_map = _build_face_maps(
        cells,
        boundary_faces,
        boundary_markers,
    )

    top_indices = np.argsort(aspect_ratio)[-top_k:][::-1]
    threshold_cells = np.flatnonzero(aspect_ratio > float(aspect_threshold))
    components = _threshold_components(threshold_cells, neighbors)
    top_index_set = set(map(int, top_indices.tolist()))

    boundary_triangles = (
        vertices[boundary_faces]
        if len(boundary_faces) > 0
        else np.empty((0, 3, 3), dtype=np.float64)
    )

    component_reports: list[dict] = []
    top_cell_to_component: dict[int, int] = {}
    for component_id, component_cells in enumerate(components, start=1):
        for cell_index in component_cells:
            top_cell_to_component[int(cell_index)] = component_id
        export_component = bool(top_index_set.intersection(component_cells))

        component_aspect = aspect_ratio[np.asarray(component_cells, dtype=np.int64)]
        component_radius = radius_ratio[np.asarray(component_cells, dtype=np.int64)]
        component_quality = element_quality[np.asarray(component_cells, dtype=np.int64)]
        component_centroids = centroids[np.asarray(component_cells, dtype=np.int64)]
        boundary_counter = _component_boundary_markers(
            component_cells,
            cells,
            boundary_marker_map,
        )

        cells_path = None
        boundary_path = None
        neighborhood_boundary_counter: dict[str, int] = {}
        neighborhood_boundary_orientation: dict[str, dict[str, int]] = {}
        if export_component:
            neighborhood = _n_hop_neighborhood(component_cells, neighbors, neighborhood_hops)
            neighborhood_vertices, neighborhood_cells, original_cell_indices = _extract_cell_submesh(
                vertices,
                cells,
                neighborhood,
            )
            neighborhood_faces_vertices, neighborhood_faces, neighborhood_face_markers = (
                _extract_boundary_surface(
                    vertices,
                    boundary_faces,
                    boundary_markers,
                    face_to_cells,
                    set(neighborhood),
                )
            )

            component_name = f"component_{component_id:02d}"
            cells_path = output_dir / f"{component_name}_cells.vtu"
            meshio.write_points_cells(
                str(cells_path),
                neighborhood_vertices,
                [("tetra", neighborhood_cells)],
                cell_data={
                    "original_cell_index": [original_cell_indices.astype(np.int64)],
                    "aspect_ratio": [aspect_ratio[original_cell_indices]],
                    "radius_ratio": [radius_ratio[original_cell_indices]],
                    "element_quality": [element_quality[original_cell_indices]],
                    "is_component_core": [
                        np.isin(original_cell_indices, np.asarray(component_cells, dtype=np.int64)).astype(np.int32)
                    ],
                },
            )
            if len(neighborhood_faces) > 0:
                neighborhood_boundary_counter = {
                    str(marker): int(count)
                    for marker, count in sorted(
                        Counter(map(int, neighborhood_face_markers.tolist())).items()
                    )
                }
                neighborhood_boundary_orientation = _boundary_orientation_summary(
                    neighborhood_faces_vertices,
                    neighborhood_faces,
                    neighborhood_face_markers,
                )
                boundary_path = output_dir / f"{component_name}_boundary.vtu"
                meshio.write_points_cells(
                    str(boundary_path),
                    neighborhood_faces_vertices,
                    [("triangle", neighborhood_faces)],
                    cell_data={
                        "marker": [neighborhood_face_markers.astype(np.int64)],
                    },
                )

        component_reports.append(
            {
                "component_id": component_id,
                "num_core_cells": int(len(component_cells)),
                "num_neighborhood_cells": (
                    0 if not export_component else int(len(neighborhood))
                ),
                "aspect_ratio_max": float(component_aspect.max()),
                "aspect_ratio_mean": float(component_aspect.mean()),
                "radius_ratio_max": float(component_radius.max()),
                "element_quality_min": float(component_quality.min()),
                "centroid_mean": component_centroids.mean(axis=0).tolist(),
                "boundary_markers": {str(marker): int(count) for marker, count in sorted(boundary_counter.items())},
                "neighborhood_boundary_markers": neighborhood_boundary_counter,
                "neighborhood_boundary_orientation": neighborhood_boundary_orientation,
                "cells_vtu": None if cells_path is None else cells_path.name,
                "boundary_vtu": None if boundary_path is None else boundary_path.name,
                "top_cell_indices": [
                    int(idx)
                    for idx in np.asarray(component_cells, dtype=np.int64)[np.argsort(component_aspect)[-5:][::-1]]
                ],
            }
        )

    exported_component_ids = {
        int(item["component_id"]) for item in component_reports if item["cells_vtu"] is not None
    }

    worst_cells: list[dict] = []
    for rank, cell_index in enumerate(top_indices, start=1):
        cell = cells[int(cell_index)]
        direct_markers = []
        for face in _tet_faces(cell):
            marker = boundary_marker_map.get(_sorted_face(face))
            if marker is not None:
                direct_markers.append(int(marker))

        nearest_marker = None
        nearest_distance = None
        if len(boundary_triangles) > 0:
            distances = _point_to_triangles_distance(centroids[int(cell_index)], boundary_triangles)
            nearest_face_index = int(np.argmin(distances))
            nearest_marker = int(boundary_markers[nearest_face_index])
            nearest_distance = float(distances[nearest_face_index])

        worst_cells.append(
            {
                "rank": rank,
                "cell_index": int(cell_index),
                "component_id": top_cell_to_component.get(int(cell_index)),
                "aspect_ratio": float(aspect_ratio[int(cell_index)]),
                "radius_ratio": float(radius_ratio[int(cell_index)]),
                "element_quality": float(element_quality[int(cell_index)]),
                "centroid": centroids[int(cell_index)].tolist(),
                "direct_boundary_markers": direct_markers,
                "nearest_boundary_marker": nearest_marker,
                "nearest_boundary_distance": nearest_distance,
            }
        )

    return {
        "num_cells": int(len(cells)),
        "boundary_marker_histogram": {
            str(marker): int(count)
            for marker, count in sorted(Counter(map(int, boundary_markers.tolist())).items())
        },
        "aspect_threshold": float(aspect_threshold),
        "num_cells_above_threshold": int(len(threshold_cells)),
        "num_threshold_components": int(len(components)),
        "num_exported_components": int(len(exported_component_ids)),
        "top_cells": worst_cells,
        "components": [
            item for item in component_reports if int(item["component_id"]) in exported_component_ids
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect worst tetrahedra and export local neighborhoods."
    )
    parser.add_argument("mesh", type=Path, help="Input DTCC volume mesh (.xdmf/.vtu).")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for the JSON report and exported neighborhoods.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=12,
        help="Number of worst cells to report.",
    )
    parser.add_argument(
        "--aspect-threshold",
        type=float,
        default=20.0,
        help="Aspect-ratio threshold used to form bad-cell components.",
    )
    parser.add_argument(
        "--neighborhood-hops",
        type=int,
        default=2,
        help="Number of face-adjacency hops exported around each bad-cell component.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mesh_path = args.mesh.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else mesh_path.with_suffix("")
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh = load_volume_mesh(mesh_path)
    report = analyze_volume_mesh(
        mesh,
        top_k=args.top_k,
        aspect_threshold=args.aspect_threshold,
        neighborhood_hops=args.neighborhood_hops,
        output_dir=output_dir,
    )
    report["mesh"] = str(mesh_path)

    report_path = output_dir / "outlier_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
