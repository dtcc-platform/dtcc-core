"""Planar domains extracted from an existing, source-labelled city shell."""

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import shapely
from shapely.geometry import Polygon

from .wall_faces import _face_geometry, _wall_groups


@dataclass
class SurfaceFace:
    """Source-owned planar domain with a 3D frame and shared boundary references."""

    id: int
    kind: str
    source_key: str | None
    markers: list[int]
    original_faces: np.ndarray
    origin: np.ndarray
    basis: np.ndarray
    normal: np.ndarray
    polygon: Polygon
    boundary: list[tuple[int, int]]
    rings: list[list[int]]
    original_plane_error: float


def boundary_rings(edges):
    """Follow oriented boundary edges without dropping collinear vertices."""
    following = {}
    incoming = set()
    for a, b in edges:
        if a in following or b in incoming:
            raise ValueError("Planar face has a branching boundary")
        following[a] = b
        incoming.add(b)
    if set(following) != incoming:
        raise ValueError("Planar face has an open boundary")
    rings = []
    while following:
        first = min(following)
        ring, current = [], first
        while True:
            ring.append(current)
            current = following.pop(current)
            if current == first:
                break
        rings.append(ring)
    return rings


def _make_face(mesh, indices, origin, normal, *, face_id, kind, source_key):
    normal = normal / np.linalg.norm(normal)
    # Choose a long projected edge to avoid amplifying short-edge roundoff.
    tri = mesh.vertices[mesh.faces[indices]]
    differences = np.concatenate(
        [tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0], tri[:, 2] - tri[:, 1]]
    )
    tangent = differences[np.argmax(np.linalg.norm(differences, axis=1))]
    tangent = tangent - normal * np.dot(tangent, normal)
    tangent /= np.linalg.norm(tangent)
    basis = np.array([tangent, np.cross(normal, tangent)])
    if kind == "wall":
        polygon, boundary, basis = _face_geometry(mesh, indices, origin, normal)
    else:
        vertex_ids = np.unique(mesh.faces[indices])
        xy = (mesh.vertices[vertex_ids] - origin) @ basis.T
        inverse = {int(v): i for i, v in enumerate(vertex_ids)}
        local = np.array([[inverse[int(v)] for v in t] for t in mesh.faces[indices]])
        polygon = shapely.union_all(shapely.polygons(xy[local]))
        if polygon.geom_type != "Polygon" or not polygon.is_valid:
            raise ValueError(f"Invalid {kind} planar face union {face_id}")
        edges = defaultdict(list)
        for t in mesh.faces[indices]:
            for a, b in zip(t, np.roll(t, -1)):
                edges[tuple(sorted((int(a), int(b))))].append((int(a), int(b)))
        boundary = [items[0] for items in edges.values() if len(items) == 1]
    rings = boundary_rings(boundary)
    # The outside loop encloses the largest area. Preserve its original IDs.
    rings.sort(key=lambda r: -Polygon((mesh.vertices[r] - origin) @ basis.T).area)
    return SurfaceFace(
        face_id,
        kind,
        source_key,
        sorted(set(map(int, np.asarray(mesh.markers)[indices]))),
        np.asarray(indices, dtype=np.int64),
        origin,
        basis,
        normal,
        polygon,
        boundary,
        rings,
        float(np.max(abs((tri - origin) @ normal))),
    )


def extract_surface_faces(mesh, owners):
    """Remove artificial wall seams; keep roof patches and terrain planes.

    Terrain unions require equal oriented normals and a 1 nm root-plane
    residual. All other terrain edges remain constraints of the original
    piecewise planar surface. Roof groups never cross patch markers.
    """
    v, f, markers = mesh.vertices, mesh.faces, np.asarray(mesh.markers)
    normals = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    wall_groups, _ = _wall_groups(mesh, owners)
    grouped = [(g, o, n, "wall", owners[int(markers[g[0]])]) for g, o, n in wall_groups]
    nonwalls = np.flatnonzero((markers < 0) | (markers >= len(owners)))
    edges = defaultdict(list)
    for i in nonwalls:
        for a, b in zip(f[i], np.roll(f[i], -1)):
            edges[tuple(sorted((int(a), int(b))))].append(int(i))
    adjacent = defaultdict(list)
    for items in edges.values():
        if len(items) != 2:
            continue
        i, j = items
        if markers[i] == markers[j] and np.dot(normals[i], normals[j]) >= 1 - 1e-14:
            adjacent[i].append(j)
            adjacent[j].append(i)
    seen = set()
    for i in nonwalls:
        if i in seen:
            continue
        origin, normal = v[f[i, 0]], normals[i]
        todo, indices = [int(i)], []
        seen.add(int(i))
        while todo:
            j = todo.pop()
            indices.append(j)
            for k in adjacent[j]:
                if k not in seen and np.max(abs((v[f[k]] - origin) @ normal)) <= 1e-9:
                    seen.add(k)
                    todo.append(k)
        kind = "terrain" if markers[i] < 0 else "roof"
        owner = None if kind == "terrain" else owners[int(markers[i]) - len(owners)]
        grouped.append((indices, origin, normal, kind, owner))
    return [
        _make_face(mesh, g, o, n, face_id=i, kind=kind, source_key=owner)
        for i, (g, o, n, kind, owner) in enumerate(grouped)
    ]
