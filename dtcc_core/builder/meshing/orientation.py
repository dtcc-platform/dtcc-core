"""Consistent face orientation for triangle meshes.

A mesh assembled from independently oriented surfaces, such as an extruded
building, can end up with neighbouring triangles wound in opposite directions.
Exchange formats derive facet normals from the vertex order (STL stores no
usable normals of its own), so mixed winding reaches other tools as inverted
normals and is reported as a defective mesh.
"""

from __future__ import annotations

import numpy as np

from ...model import Mesh


def _corner_edges(faces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the (3F,) start and end vertex of every face corner edge."""
    starts = faces.reshape(-1)
    ends = np.roll(faces, -1, axis=1).reshape(-1)
    return starts, ends


def _edge_ids(faces: np.ndarray) -> np.ndarray:
    """Give every undirected edge an id, as a (3F,) array of corner edges."""
    starts, ends = _corner_edges(faces)
    keys = np.stack([np.minimum(starts, ends), np.maximum(starts, ends)], axis=1)
    _, inverse = np.unique(keys, axis=0, return_inverse=True)
    return inverse.reshape(-1)


def _faces_per_edge(edge_ids: np.ndarray, face_count: int) -> dict[int, np.ndarray]:
    face_of_corner = np.repeat(np.arange(face_count), 3)
    order = np.argsort(edge_ids, kind="stable")
    sorted_ids = edge_ids[order]
    sorted_faces = face_of_corner[order]
    splits = np.flatnonzero(np.diff(sorted_ids)) + 1
    groups = np.split(sorted_faces, splits)
    starts = np.concatenate(([0], splits))
    return {int(sorted_ids[start]): group for start, group in zip(starts, groups)}


def _directed_edges(face: np.ndarray) -> set[tuple[int, int]]:
    a, b, c = (int(index) for index in face)
    return {(a, b), (b, c), (c, a)}


def _signed_volume(vertices: np.ndarray, faces: np.ndarray) -> float:
    v0 = vertices[faces[:, 0]]
    cross = np.cross(vertices[faces[:, 1]] - v0, vertices[faces[:, 2]] - v0)
    return float(np.einsum("ij,ij->i", v0, cross).sum() / 6.0)


def _is_closed(edge_ids: np.ndarray, faces_in_component: np.ndarray) -> bool:
    corner_ids = edge_ids.reshape(-1, 3)[faces_in_component].reshape(-1)
    _, counts = np.unique(corner_ids, return_counts=True)
    return bool(np.all(counts == 2))


def orient_faces_consistently(mesh: Mesh) -> Mesh:
    """Wind every face of ``mesh`` consistently with its neighbours.

    Faces sharing an edge are consistent when each traverses that edge in the
    opposite direction. Each connected component is made consistent, then
    turned the right way out: a closed component is oriented to a positive
    signed volume (normals pointing outwards), while an open one keeps the
    orientation most of its faces already had.

    The mesh is modified in place and returned.
    """
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if faces.size == 0:
        return mesh

    faces = faces.copy()
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    edge_ids = _edge_ids(faces)
    edge_faces = _faces_per_edge(edge_ids, len(faces))
    corner_ids = edge_ids.reshape(-1, 3)

    visited = np.zeros(len(faces), dtype=bool)
    flipped = np.zeros(len(faces), dtype=bool)

    for seed in range(len(faces)):
        if visited[seed]:
            continue
        visited[seed] = True
        component = [seed]
        stack = [seed]
        while stack:
            current = stack.pop()
            current_edges = _directed_edges(faces[current])
            for edge_id in corner_ids[current]:
                for neighbour in edge_faces[int(edge_id)]:
                    neighbour = int(neighbour)
                    if visited[neighbour]:
                        continue
                    # Consistent neighbours walk the shared edge the other way.
                    if current_edges & _directed_edges(faces[neighbour]):
                        faces[neighbour] = faces[neighbour][::-1]
                        flipped[neighbour] = True
                    visited[neighbour] = True
                    component.append(neighbour)
                    stack.append(neighbour)

        members = np.asarray(component, dtype=np.int64)
        if _is_closed(edge_ids, members):
            turn_over = _signed_volume(vertices, faces[members]) < 0.0
        else:
            turn_over = bool(flipped[members].sum() * 2 > len(members))
        if turn_over:
            faces[members] = faces[members][:, ::-1]
            flipped[members] = ~flipped[members]

    mesh.faces = faces
    normals = np.asarray(mesh.normals) if mesh.normals is not None else None
    if normals is not None and normals.shape == (len(faces), 3):
        normals[flipped] = -normals[flipped]
        mesh.normals = normals
    return mesh
