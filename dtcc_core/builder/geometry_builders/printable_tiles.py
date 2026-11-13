"""Cut an oriented city surface into closed printable solids."""

from itertools import pairwise

import numpy as np
import shapely
from shapely.geometry import Polygon

from ...model import Mesh
from ..geometry_policy import XY_PRECISION
from ..meshing.compact_triangulation import orient_triangles, triangulate
from ..meshing.meshing import disjoint_meshes
from ..meshing.roof_validation import vertex_link_audit
from ..meshing.surface_faces import boundary_rings

BASE_SIZE_MM = 200.0
BASE_HEIGHT_MM = 20.0
MAGNET_INSET_MM = 8.0
MAGNET_RADIUS_MM = 5.0
MAGNET_DEPTH_MM = 2.5
MINIMUM_MATERIAL_MM = 1.0
MAGNET_SEGMENTS = 32
MINIMUM_TILE_SIZE_MM = 140.0
MAXIMUM_TILE_SIZE_MM = 250.0


def _clip_polygon(polygon, axis, limit, direction):
    distances = direction * (polygon[:, axis] - limit)
    result = []
    for i, a in enumerate(polygon):
        j = (i + 1) % len(polygon)
        if distances[i] >= 0:
            result.append(a)
        if (distances[i] < 0 < distances[j]) or (distances[j] < 0 < distances[i]):
            b = polygon[j]
            if a[axis] <= b[axis]:
                low, high = a, b
            else:
                low, high = b, a
            point = low + (limit - low[axis]) / (high[axis] - low[axis]) * (high - low)
            point[axis] = limit
            result.append(point)
    return np.asarray(result)


def _clip_triangles(triangles, markers, planes):
    """Clip indexed source facets in XYZ, preserving their orientation."""
    clipped, owners = [], []
    for triangle, marker in zip(triangles, markers, strict=True):
        polygon = triangle
        for axis, limit, direction in planes:
            d = direction * (polygon[:, axis] - limit)
            if np.any(d < 0):
                polygon = _clip_polygon(polygon, axis, limit, direction)
                if len(polygon) < 3:
                    break
        for i in range(1, len(polygon) - 1):
            clipped.append([polygon[0], polygon[i], polygon[i + 1]])
            owners.append(marker)
    return np.asarray(clipped).reshape(-1, 3, 3), np.asarray(owners)


def _indexed(triangles, markers=None):
    vertices, inverse = np.unique(
        np.asarray(triangles).reshape(-1, 3), axis=0, return_inverse=True
    )
    faces = inverse.reshape(-1, 3)
    keep = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2])
    keep &= faces[:, 2] != faces[:, 0]
    used, indices = np.unique(faces[keep], return_inverse=True)
    return Mesh(
        vertices=vertices[used],
        faces=indices.reshape(-1, 3),
        markers=np.empty(0, dtype=int) if markers is None else markers[keep],
    )


def _edges(faces):
    directed = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    _, first, inverse, counts = np.unique(
        np.sort(directed, axis=1),
        axis=0,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    directions = np.bincount(
        inverse.reshape(-1), weights=directed[:, 0] < directed[:, 1]
    )
    if np.any((counts > 2) | ((counts == 2) & (directions != 1))):
        raise ValueError("inconsistent winding or nonmanifold edges")
    return directed[first[counts == 1]]


def _clip_surface(triangles, markers, bounds):
    low, high = np.asarray(bounds).reshape(2, 2)
    xyz = triangles.copy()
    xy = xyz[:, :, :2]
    for edge in (low, high):
        close = np.abs(xy - edge) <= XY_PRECISION
        np.copyto(xy, edge, where=close)
    selected = np.all((xy.max(axis=1) > low) & (xy.min(axis=1) < high), axis=1)
    planes = [(0, low[0], 1), (0, high[0], -1), (1, low[1], 1), (1, high[1], -1)]
    clipped, owners = _clip_triangles(xyz[selected], markers[selected], planes)
    clipped[:, :, :2] -= low
    upper_edge = clipped[:, :, :2] == high - low
    clipped = np.round(clipped, 7)
    np.copyto(clipped[:, :, :2], high - low, where=upper_edge)
    return clipped, owners


def _ring_edges(ring):
    return np.column_stack((ring, np.roll(ring, -1)))


def _side_walls(mesh, boundary):
    vertices = mesh.vertices.tolist()
    profiles = {}
    for index in np.unique(boundary):
        xy = tuple(mesh.vertices[index, :2])
        profiles.setdefault(xy, []).append(int(index))
    columns = {}
    for xy, indices in profiles.items():
        base = len(vertices)
        vertices.append([*xy, 0.0])
        ordered = sorted(indices, key=lambda i: mesh.vertices[i, 2], reverse=True)
        for i, vertex in enumerate(ordered):
            columns[vertex] = [*ordered[i:], base]

    faces = []
    for a, b in boundary:
        left, right = columns[a], columns[b]
        if left[-1] != right[-1]:
            faces.extend((b, p, q) for p, q in pairwise(left))
            faces.extend((left[-1], q, p) for p, q in pairwise(right))
    return vertices, np.asarray(faces).reshape(-1, 3)


def _add_underside(vertices, walls, size):
    walls = np.asarray(walls)
    bottom = np.asarray(vertices)[walls, 2] == 0
    sides = walls[bottom.sum(axis=1) == 2]
    floor_edges = [
        (a, b) for a, b in _edges(sides) if vertices[a][2] == vertices[b][2] == 0
    ]
    [outline] = boundary_rings(floor_edges)
    points = np.asarray(vertices)[outline, :2]
    footprint = Polygon(points)
    if not footprint.is_valid or abs(footprint.area - size**2) > size * 1e-4:
        raise ValueError("surface does not cover the tile footprint")

    angles = np.arange(MAGNET_SEGMENTS) * (2 * np.pi / MAGNET_SEGMENTS)
    circle = MAGNET_RADIUS_MM * np.column_stack((np.cos(angles), np.sin(angles)))
    centers = [
        (x, y)
        for y in (MAGNET_INSET_MM, size - MAGNET_INSET_MM)
        for x in (MAGNET_INSET_MM, size - MAGNET_INSET_MM)
    ]
    openings = [circle + center for center in centers]
    xy = np.vstack([points, *openings])
    rings = np.split(
        np.arange(len(xy)), [len(points) + n * MAGNET_SEGMENTS for n in range(4)]
    )
    polygon = Polygon(points, openings)
    segments = np.vstack([_ring_edges(ring) for ring in rings])
    cap, _ = triangulate(xy, segments, polygon)
    first = len(vertices)
    mapping = np.r_[outline, np.arange(first, first + 4 * MAGNET_SEGMENTS)]
    vertices.extend([[x, y, 0.0] for x, y in xy[len(points) :]])
    faces = [mapping[orient_triangles(xy, cap)[:, ::-1]]]
    for ring, center in zip(rings[1:], centers, strict=True):
        ceiling = len(vertices) + np.arange(MAGNET_SEGMENTS)
        hub = len(vertices) + MAGNET_SEGMENTS
        vertices.extend([[*point, MAGNET_DEPTH_MM] for point in xy[ring]])
        vertices.append([*center, MAGNET_DEPTH_MM])
        a, b = _ring_edges(mapping[ring]).T
        c, d = _ring_edges(ceiling).T
        faces.extend(
            (
                np.column_stack((b, a, d)),
                np.column_stack((d, a, c)),
                np.column_stack((d, c, np.full(len(c), hub))),
            )
        )
    return np.concatenate(faces)


def _validate_solid(mesh):
    triangles = mesh.vertices[mesh.faces]
    relative = triangles - mesh.vertices[0]
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    if not np.all(np.isfinite(triangles)) or np.any(np.all(normals == 0, axis=1)):
        raise ValueError("nonfinite coordinates or zero-area triangles")
    if len(_edges(mesh.faces)):
        raise ValueError("open edges")
    if (
        vertex_link_audit(mesh)["nonmanifold_vertices"]
        or len(disjoint_meshes(mesh)) != 1
    ):
        raise ValueError("expected one manifold solid")
    if (
        np.einsum("ij,ij->", relative[:, 0], np.cross(relative[:, 1], relative[:, 2]))
        <= 0
    ):
        raise ValueError("nonpositive volume")


def _resize_base(mesh, size):
    """Resize the central span while preserving the corner and connector dimensions."""
    if size == BASE_SIZE_MM:
        return mesh
    triangles = np.asarray(mesh.vertices, dtype=float)[mesh.faces]
    triangles[np.abs(triangles) <= XY_PRECISION] = 0
    for axis in (0, 1):
        triangles = _stretch_span(triangles, axis, size - BASE_SIZE_MM)
    result = _indexed(np.round(triangles, 8).astype(np.float32).astype(float))
    _validate_solid(result)
    return result


def _stretch_span(triangles, axis, extra):
    start, stop = 60.0, 140.0
    intervals = (
        (triangles[:, :, axis].min(axis=1) < start, [(axis, start, -1)]),
        (np.ones(len(triangles), dtype=bool), [(axis, start, 1), (axis, stop, -1)]),
        (triangles[:, :, axis].max(axis=1) > stop, [(axis, stop, 1)]),
    )
    sections = []
    for selected, planes in intervals:
        chosen = triangles[selected]
        clipped, _ = _clip_triangles(chosen, np.zeros(len(chosen), dtype=int), planes)
        clipped[:, :, axis] += extra * np.clip(
            (clipped[:, :, axis] - start) / (stop - start), 0, 1
        )
        sections.append(clipped)
    return np.concatenate(sections)


def _building_polygons(surface):
    triangles = surface.vertices[surface.faces, :2]
    a, b = triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    roofs = (surface.markers >= 0) & (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0] > 0)
    polygons = shapely.union_all(shapely.polygons(triangles[roofs]))
    return [
        [np.asarray(ring.coords).tolist() for ring in (p.exterior, *p.interiors)]
        for p in shapely.get_parts(polygons)
        if isinstance(p, Polygon) and not p.is_empty
    ]


def _build_printable_tile(triangles, markers, bounds, size, scale, anchor):
    clipped, owners = _clip_surface(triangles, markers, bounds)
    if not len(clipped):
        raise ValueError("city surface does not cover the tile")
    mm_per_m = 1000.0 / scale
    bottom = np.floor(
        (
            (clipped[:, :, 2].min() - anchor) * mm_per_m
            - MAGNET_DEPTH_MM
            - MINIMUM_MATERIAL_MM
        )
        / BASE_HEIGHT_MM
    )
    base_z = anchor + bottom * BASE_HEIGHT_MM / mm_per_m
    clipped[:, :, 2] -= base_z
    clipped *= mm_per_m
    for axis in (0, 1):
        edge = (bounds[axis + 2] - bounds[axis]) * mm_per_m
        clipped[:, :, axis][np.abs(clipped[:, :, axis] - edge) < 1e-5] = size
    surface = _indexed(clipped.astype(np.float32).astype(float), owners)
    boundary = _edges(surface.faces)
    xy = surface.vertices[boundary, :2]
    if np.any(
        ~np.any(
            np.all(xy == 0, axis=1) | np.all(xy == np.float32(size), axis=1), axis=1
        )
    ):
        raise ValueError("surface has an opening away from the tile boundary")
    vertices, walls = _side_walls(surface, boundary)
    underside = _add_underside(vertices, walls, float(np.float32(size)))
    faces = np.concatenate((surface.faces, walls, underside))
    solid = _indexed(np.asarray(vertices, dtype=np.float32)[faces].astype(float))
    _validate_solid(solid)
    return solid, int(bottom), _building_polygons(surface)
