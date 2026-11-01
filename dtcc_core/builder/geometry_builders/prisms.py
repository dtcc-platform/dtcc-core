"""Private decomposition of closed, ground-supported column solids.

All coordinates remain in one local metre frame. Geometry is copied; the
kernel neither repairs open shells nor infers occupancy from roof area alone.
Internal caps and touching independent shells are unsupported.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from math import pi, radians, sin

import numpy as np
from shapely import normalize, union_all
from shapely.affinity import translate
from shapely.geometry import MultiPolygon, Polygon, box
from shapely.ops import polygonize
from shapely.validation import explain_validity

__all__ = []


class _PrismError(ValueError):
    """Invalid or unsupported solid, with its source location."""

    def __init__(
        self, reason, message, *, component=None, face_index=None, location=None
    ):
        self.reason = reason
        self.component = component
        self.face_index = face_index
        self.location = location
        super().__init__(f"{reason}: {message}")


@dataclass
class _PrismFace:
    rings: tuple[np.ndarray, ...]
    source_index: int = 0


@dataclass
class _PrismComponent:
    ground_z: float
    base: Polygon | MultiPolygon
    levels: tuple[tuple[float, Polygon | MultiPolygon], ...]
    diagnostics: dict = field(default_factory=dict)
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0)
    srs: str = ""


def _check_parameters(tolerance, wall_angle_tolerance):
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if not np.isfinite(wall_angle_tolerance) or not 0 <= wall_angle_tolerance < 90:
        raise ValueError("wall_angle_tolerance must be finite and in [0, 90) degrees")


def _frame(geometry):
    """Read a translation-only frame without applying it to local coordinates."""
    affine = np.asarray(geometry.transform.affine)
    if (
        affine.shape != (4, 4)
        or not np.issubdtype(affine.dtype, np.number)
        or np.iscomplexobj(affine)
        or not np.array_equal(affine[:3, :3], np.eye(3))
        or not np.array_equal(affine[3], [0, 0, 0, 1])
        or not np.isfinite(affine).all()
    ):
        raise _PrismError(
            "unsupported_frame", "only translation-only metre frames are supported"
        )
    return affine[:3, 3].copy(), geometry.transform.srs


def _check_metadata(geometry):
    if (
        geometry.fields
        or geometry.regions
        or geometry.schema_id is not None
        or geometry.schema_version is not None
        or geometry.dataset_context is not None
    ):
        raise _PrismError(
            "unsupported_metadata",
            "geometry fields, regions, context and schema declarations require explicit transfer",
        )


def _extract_faces(geometry):
    if isinstance(geometry, (list, tuple)):
        if not all(isinstance(face, _PrismFace) for face in geometry):
            raise TypeError("Ring records must be _PrismFace instances")
        return list(geometry), (np.zeros(3), "")

    from dtcc_core.model import Mesh, MultiSurface, Solid

    if not isinstance(geometry, (Mesh, MultiSurface, Solid)):
        raise _PrismError(
            "missing_geometry", "expected Mesh, MultiSurface or one-shell Solid"
        )
    _check_metadata(geometry)
    offset, srs = _frame(geometry)
    if isinstance(geometry, Mesh):
        return _mesh_records(geometry), (offset, srs)
    if isinstance(geometry, Solid):
        if len(geometry.shells) != 1 or sorted(geometry.shells[0]) != list(
            range(len(geometry.surfaces))
        ):
            raise _PrismError(
                "not_2_5d",
                "one complete exterior shell is required; cavities are unsupported",
            )
    frames = [_frame(surface) for surface in geometry.surfaces]
    translations = {tuple(value) for value, _ in frames}
    if len(translations) > 1:
        raise _PrismError("unsupported_frame", "surfaces must share one translation")
    child_offset = frames[0][0] if frames else np.zeros(3)
    for surface, (_, child_srs) in zip(geometry.surfaces, frames, strict=True):
        _check_metadata(surface)
        if srs and child_srs and srs != child_srs:
            raise _PrismError("unsupported_frame", "surfaces must share one SRS")
    nonempty_srs = {value for _, value in frames if value}
    if len(nonempty_srs) > 1:
        raise _PrismError("unsupported_frame", "surfaces must share one SRS")
    records = [
        _PrismFace((surface.vertices, *surface.holes), i)
        for i, surface in enumerate(geometry.surfaces)
    ]
    return records, (offset + child_offset, srs or next(iter(nonempty_srs), ""))


def _mesh_records(geometry):
    vertices, faces = np.asarray(geometry.vertices), np.asarray(geometry.faces)
    if (
        vertices.ndim != 2
        or vertices.shape[1] != 3
        or faces.ndim != 2
        or faces.shape[1] != 3
    ):
        raise _PrismError(
            "degenerate", "mesh requires (N, 3) vertices and triangular faces"
        )
    if not np.issubdtype(vertices.dtype, np.number) or not np.isfinite(vertices).all():
        raise _PrismError("degenerate", "all mesh vertices must be finite numbers")
    if (
        not np.issubdtype(faces.dtype, np.integer)
        or np.any(faces < 0)
        or np.any(faces >= len(vertices))
    ):
        raise _PrismError("degenerate", "mesh face indices must reference vertices")
    if geometry.markers.size:
        raise _PrismError(
            "unsupported_metadata", "mesh markers require explicit transfer"
        )
    return [_PrismFace((vertices[face],), i) for i, face in enumerate(faces)]


def _normal(ring):
    local = ring - ring[0]
    return np.cross(local, np.roll(local, -1, axis=0)).sum(axis=0)


def _clean_ring(ring, face_index, tolerance):
    values = np.asarray(ring)
    if (
        values.ndim != 2
        or values.shape[1] != 3
        or np.iscomplexobj(values)
        or not np.issubdtype(values.dtype, np.number)
    ):
        raise _PrismError(
            "degenerate", "rings require real (N, 3) coordinates", face_index=face_index
        )
    if not np.isfinite(values).all():
        raise _PrismError(
            "degenerate", "nonfinite ring coordinates", face_index=face_index
        )
    source_values = (
        values
        if np.issubdtype(values.dtype, np.floating)
        else values.astype(np.float64)
    )
    if values.size and np.max(abs(np.spacing(source_values))) >= tolerance / 4:
        raise _PrismError(
            "precision_loss",
            "stored coordinate precision cannot resolve the requested tolerance",
            face_index=face_index,
        )
    values = np.asarray(values, dtype=np.float64).copy()
    if len(values):
        values = values[np.r_[True, np.any(values[1:] != values[:-1], axis=1)]]
        if len(values) > 1 and np.array_equal(values[0], values[-1]):
            values = values[:-1]
    if len(np.unique(values, axis=0)) < 3 or np.linalg.norm(_normal(values)) == 0:
        raise _PrismError(
            "degenerate", "ring has no finite positive area", face_index=face_index
        )
    return values


def _clean_faces(records, tolerance):
    result = []
    for face in records:
        if not face.rings:
            raise _PrismError(
                "degenerate", "face has no exterior ring", face_index=face.source_index
            )
        rings = [_clean_ring(ring, face.source_index, tolerance) for ring in face.rings]
        normal = _normal(rings[0])
        rings[1:] = [
            ring[::-1].copy() if np.dot(normal, _normal(ring)) > 0 else ring
            for ring in rings[1:]
        ]
        result.append(_PrismFace(tuple(rings), face.source_index))
    if not result:
        raise _PrismError("missing_geometry", "no faces")
    return result


def _clusters(points, tolerance):
    """Distance neighbourhoods with a whole-cluster diameter guard."""
    from scipy.spatial import cKDTree

    points, inverse = np.unique(points, axis=0, return_inverse=True)
    parents = np.arange(len(points))

    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    for a, b in sorted(cKDTree(points).query_pairs(tolerance)):
        a, b = root(a), root(b)
        parents[max(a, b)] = min(a, b)
    groups = defaultdict(list)
    for i in range(len(points)):
        groups[root(i)].append(i)
    representatives = points.copy()
    for members in groups.values():
        values = points[members]
        # The diameter bound prevents unbounded single-link chains.
        if len(values) > 1 and np.any(
            np.linalg.norm(values[:, None] - values[None, :], axis=2) > tolerance
        ):
            raise _PrismError(
                "precision_loss",
                "vertex neighbourhood exceeds the whole-cluster diameter",
                location=values[0].tolist(),
            )
        representatives[members] = (values.min(axis=0) + values.max(axis=0)) / 2
    return representatives[inverse]


def _replace_rings(faces, values):
    result, start = [], 0
    for face in faces:
        rings = []
        for ring in face.rings:
            rings.append(values[start : start + len(ring)].copy())
            start += len(ring)
        result.append(_PrismFace(tuple(rings), face.source_index))
    return result


def _vertex_groups(faces):
    """Provisional components, including point contacts that topology must reject."""
    owners = defaultdict(set)
    for i, face in enumerate(faces):
        for ring in face.rings:
            for point in ring:
                owners[tuple(point)].add(i)
    neighbours = [set() for _ in faces]
    for members in owners.values():
        seed = min(members)
        for i in members:
            neighbours[seed].add(i)
            neighbours[i].add(seed)
    groups, visited = [], set()
    for seed in range(len(faces)):
        if seed in visited:
            continue
        stack, members = [seed], []
        visited.add(seed)
        while stack:
            i = stack.pop()
            members.append(i)
            for j in sorted(neighbours[i] - visited):
                visited.add(j)
                stack.append(j)
        groups.append(sorted(members))
    return groups


def _classify(face, tolerance, angle):
    points = np.vstack(face.rings)
    normal = _normal(face.rings[0])
    length = np.linalg.norm(normal)
    if length == 0 or not np.isfinite(length):
        raise _PrismError(
            "precision_loss",
            "normalization collapsed a face",
            face_index=face.source_index,
        )
    residual = float(np.max(np.abs((points - points[0]) @ (normal / length))))
    if residual > tolerance:
        raise _PrismError(
            "sloped_face",
            "face is not planar within tolerance",
            face_index=face.source_index,
        )
    if np.ptp(points[:, 2]) <= tolerance:
        return "horizontal", residual
    if abs(normal[2]) / length > sin(radians(angle)):
        raise _PrismError(
            "sloped_face",
            "face is neither horizontal nor a vertical wall",
            face_index=face.source_index,
        )
    return "wall", residual


def _height_groups(faces, kinds, tolerance):
    intervals = sorted(
        (min(r[:, 2].min() for r in face.rings), max(r[:, 2].max() for r in face.rings))
        for face, kind in zip(faces, kinds, strict=True)
        if kind == "horizontal"
    )
    groups = []
    for low, high in intervals:
        if groups and low - groups[-1][1] <= tolerance:
            if high - groups[-1][0] > tolerance:
                raise _PrismError(
                    "ambiguous_levels",
                    "horizontal height chain exceeds the whole-cluster span",
                )
            groups[-1][1] = max(groups[-1][1], high)
        else:
            groups.append([low, high])
    if len(groups) < 2:
        raise _PrismError(
            "ambiguous_levels", "base and roof levels cannot be separated"
        )
    return groups


def _normalize_faces(originals, welded, tolerance, angle, epsilon, origin):
    classified = [_classify(face, tolerance, angle) for face in originals]
    kinds = [kind for kind, _ in classified]
    heights = _height_groups(originals, kinds, tolerance)
    levels = np.array([(low + high) / 2 for low, high in heights])
    initial = np.vstack([ring for face in originals for ring in face.rings])
    points = np.vstack([ring for face in welded for ring in face.rings])
    path = np.linalg.norm(initial - points, axis=1)
    vertical = points.copy()
    vertical[:, :2] = _clusters(points[:, :2], tolerance)
    nearest = np.argmin(abs(points[:, 2, None] - levels), axis=1)
    close = abs(points[:, 2] - levels[nearest]) <= tolerance
    vertical[close, 2] = levels[nearest[close]]
    path += np.linalg.norm(vertical - points, axis=1)
    path += np.linalg.norm((vertical + origin) - origin - vertical, axis=1)
    movement = float(path.max())
    if movement > tolerance + epsilon:
        raise _PrismError(
            "precision_loss", "combined welding and verticalization exceeds tolerance"
        )
    faces = _replace_rings(welded, vertical)
    for face in faces:
        _classify(face, epsilon, angle)
        for ring in face.rings:
            if len(np.unique(ring, axis=0)) != len(ring):
                raise _PrismError(
                    "precision_loss",
                    "normalization collapsed a nonzero ring edge",
                    face_index=face.source_index,
                )
    return faces, kinds, heights, movement, max(residual for _, residual in classified)


def _canonical_ring(ring):
    values = [tuple(point) for point in ring]
    candidates = []
    for sequence in (values, values[::-1]):
        start = sequence.index(min(sequence))
        candidates.append(tuple(sequence[start:] + sequence[:start]))
    return min(candidates)


def _face_key(face):
    return _canonical_ring(face.rings[0]), tuple(
        sorted(_canonical_ring(ring) for ring in face.rings[1:])
    )


def _node_ring(ring, vertices, epsilon):
    values = []
    for a, b in zip(ring, np.roll(ring, -1, axis=0), strict=True):
        edge = b - a
        length2 = np.dot(edge, edge)
        if length2 == 0:
            raise _PrismError("precision_loss", "collapsed edge")
        t = (vertices - a) @ edge / length2
        distance = np.linalg.norm(vertices - (a + t[:, None] * edge), axis=1)
        indices = np.flatnonzero((t >= 0) & (t < 1) & (distance <= epsilon))
        values.extend(vertices[indices[np.argsort(t[indices], kind="stable")]])
    return np.asarray(values)


def _node_and_deduplicate(faces, originals, kinds, epsilon):
    vertices = np.unique(
        np.vstack([ring for face in faces for ring in face.rings]), axis=0
    )
    records = []
    for face, original, kind in zip(faces, originals, kinds, strict=True):
        noded = _PrismFace(
            tuple(_node_ring(r, vertices, epsilon) for r in face.rings),
            face.source_index,
        )
        records.append((_face_key(noded), _face_key(original), noded, original, kind))
    records.sort(key=lambda record: record[:2])
    unique, duplicates = [], 0
    for record in records:
        if unique and record[0] == unique[-1][0]:
            duplicates += 1
        else:
            unique.append(record)
    return (
        [r[2] for r in unique],
        [r[3] for r in unique],
        [r[4] for r in unique],
        duplicates,
    )


def _shell_orientation(faces):
    incidence = defaultdict(list)
    for i, face in enumerate(faces):
        for ring in face.rings:
            for a, b in zip(ring, np.roll(ring, -1, axis=0), strict=True):
                a, b = tuple(a), tuple(b)
                edge = tuple(sorted((a, b)))
                incidence[edge].append((i, 1 if a < b else -1))
    for edge, adjacent in incidence.items():
        if len(adjacent) != 2:
            reason = "open_shell" if len(adjacent) == 1 else "nonmanifold"
            raise _PrismError(
                reason,
                "each atomic shell edge requires exactly two incident faces; internal/shared faces are unsupported",
                face_index=faces[adjacent[0][0]].source_index,
                location=list(edge[0]),
            )
    neighbours = defaultdict(list)
    for (i, a), (j, b) in incidence.values():
        neighbours[i].append((j, -a * b))
        neighbours[j].append((i, -a * b))
    orientation, stack = {0: 1}, [0]
    while stack:
        i = stack.pop()
        for j, relative in neighbours[i]:
            sign = orientation[i] * relative
            if j in orientation:
                if orientation[j] != sign:
                    raise _PrismError("nonmanifold", "shell is not orientable")
            else:
                orientation[j] = sign
                stack.append(j)
    if len(orientation) != len(faces):
        raise _PrismError("nonmanifold", "shells touch only at vertices")
    return [orientation[i] for i in range(len(faces))], len(incidence)


def _polygon(rings, face_index=None):
    polygon = Polygon(rings[0], rings[1:])
    if (
        not polygon.is_valid
        or polygon.is_empty
        or not np.isfinite(polygon.area)
        or polygon.area <= 0
    ):
        raise _PrismError(
            "degenerate",
            f"invalid face polygon: {explain_validity(polygon)}",
            face_index=face_index,
        )
    return polygon


def _parts(polygon):
    if isinstance(polygon, Polygon):
        return [polygon]
    if isinstance(polygon, MultiPolygon):
        return list(polygon.geoms)
    raise _PrismError(
        "precision_loss", "positive polygon coverage did not survive normalization"
    )


def _geometry_key(polygon):
    return (*polygon.bounds, normalize(polygon).wkb_hex)


def _coverage(actual, expected, epsilon, reason):
    area = float(actual.symmetric_difference(expected).area)
    distance = float(actual.boundary.hausdorff_distance(expected.boundary))
    if (
        area > epsilon * (actual.length + expected.length)
        or not np.isfinite(distance)
        or distance > epsilon
    ):
        raise _PrismError(
            reason,
            f"surface coverage differs: area={area:g}, boundary distance={distance:g}",
        )
    return area, distance


def _roof_cells(faces, kinds, epsilon):
    by_height = defaultdict(list)
    for face, kind in zip(faces, kinds, strict=True):
        if kind == "horizontal":
            by_height[float(face.rings[0][0, 2])].append(
                _polygon([ring[:, :2] for ring in face.rings], face.source_index)
            )
    unions = {}
    for z, polygons in sorted(by_height.items()):
        merged = union_all(polygons)
        if sum(p.area for p in polygons) - merged.area > epsilon * sum(
            p.length for p in polygons
        ):
            raise _PrismError("not_2_5d", "overlapping coplanar caps are unsupported")
        unions[z] = merged
    ground = min(unions)
    base = unions.pop(ground)
    roofs = union_all(list(unions.values()))
    area, distance = _coverage(base, roofs, epsilon, "overhang")
    linework = union_all([base.boundary, *(p.boundary for p in unions.values())])
    cells = defaultdict(list)
    for cell in polygonize(linework):
        point = cell.representative_point()
        if not base.covers(point):
            continue
        labels = [z for z, polygon in unions.items() if polygon.covers(point)]
        if len(labels) != 1:
            raise _PrismError(
                "not_2_5d",
                "roof overlap or missing roof; internal caps are unsupported",
                location=list(point.coords[0]),
            )
        _coverage(cell.intersection(base), cell, epsilon, "not_2_5d")
        _coverage(cell.intersection(unions[labels[0]]), cell, epsilon, "not_2_5d")
        cells[labels[0]].append(cell)
    levels = [(z, union_all(polygons)) for z, polygons in sorted(cells.items())]
    _coverage(union_all([p for _, p in levels]), base, epsilon, "precision_loss")
    # Dissolving equal levels can remove a neighbouring level's subdivision.
    nodes = np.unique(
        np.vstack(
            [
                np.asarray(r.coords)[:-1]
                for _, polygon in levels
                for piece in _parts(polygon)
                for r in (piece.exterior, *piece.interiors)
            ]
        ),
        axis=0,
    )
    result = []
    for z, polygon in levels:
        pieces = []
        for piece in _parts(polygon):
            rings = [
                _node_ring(np.asarray(r.coords)[:-1], nodes, epsilon)
                for r in (piece.exterior, *piece.interiors)
            ]
            pieces.append(normalize(_polygon(rings)))
        pieces.sort(key=_geometry_key)
        result.append((z, pieces[0] if len(pieces) == 1 else MultiPolygon(pieces)))
    return ground, base, result, area, distance


def _wall_edges(levels, ground):
    edges = defaultdict(list)
    for z, polygon in levels:
        for piece in _parts(polygon):
            for ring in (piece.exterior, *piece.interiors):
                points = np.asarray(ring.coords)[:-1]
                for a, b in zip(points, np.roll(points, -1, axis=0), strict=True):
                    edges[tuple(sorted((tuple(a), tuple(b))))].append(z)
    result = []
    for (a, b), values in sorted(edges.items()):
        if len(values) > 2:
            raise _PrismError("nonmanifold", "roof cells meet on an ambiguous boundary")
        low = min(values) if len(values) == 2 else ground
        high = max(values)
        if high > low:
            result.append((np.asarray(a), np.asarray(b), low, high))
    return result


def _wall_planes(levels, ground, epsilon):
    planes = []
    for a, b, low, high in _wall_edges(levels, ground):
        direction = b - a
        direction /= np.linalg.norm(direction)
        for origin, axis, _expected, _ in planes:
            normal = np.array([-axis[1], axis[0]])
            if max(abs((a - origin) @ normal), abs((b - origin) @ normal)) <= epsilon:
                break
        else:
            origin, axis, _expected = a, direction, []
            planes.append((origin, axis, _expected, []))
        values = [(a - origin) @ axis, (b - origin) @ axis]
        _expected.append(box(min(values), low, max(values), high))
    return planes


def _validate_walls(faces, kinds, levels, ground, epsilon):
    planes = _wall_planes(levels, ground, epsilon)
    for face, kind in zip(faces, kinds, strict=True):
        if kind != "wall":
            continue
        points = np.vstack(face.rings)
        for origin, axis, _, _provided in planes:
            normal = np.array([-axis[1], axis[0]])
            if np.max(abs((points[:, :2] - origin) @ normal)) <= epsilon:
                break
        else:
            raise _PrismError(
                "not_2_5d",
                "wall does not belong to the expected exterior",
                face_index=face.source_index,
            )
        rings = [
            np.column_stack(((ring[:, :2] - origin) @ axis, ring[:, 2]))
            for ring in face.rings
        ]
        _provided.append(_polygon(rings, face.source_index))
    max_area, max_distance = 0.0, 0.0
    for _, _, expected, provided in planes:
        if not provided:
            raise _PrismError("open_shell", "missing exterior or step wall")
        merged = union_all(provided)
        if sum(p.area for p in provided) - merged.area > epsilon * sum(
            p.length for p in provided
        ):
            raise _PrismError("nonmanifold", "overlapping coplanar walls")
        area, distance = _coverage(merged, union_all(expected), epsilon, "not_2_5d")
        max_area, max_distance = max(max_area, area), max(max_distance, distance)
    return max_area, max_distance


def _shell_volume(faces, orientation, reference):
    # The oriented area-vector integral includes signed holes and concavities;
    # it does not approximate a cap by an unsigned triangle fan.
    return float(
        sum(
            sign * np.dot(ring[0] - reference, _normal(ring)) / 6
            for face, sign in zip(faces, orientation, strict=True)
            for ring in face.rings
        )
    )


def _volume_budget(levels, ground, displacement, epsilon):
    """Boundary-sweep area times height, plus two moving caps times area."""
    return float(
        sum(
            (
                p.length * displacement
                + pi
                * displacement**2
                * sum(1 + len(piece.interiors) for piece in _parts(p))
            )
            * (z - ground)
            + 2 * displacement * p.area
            for z, p in levels
        )
        + epsilon * sum(p.area + p.length * (z - ground) for z, p in levels)
    )


def _base_residual(originals, faces, kinds, ground, base, movement, epsilon):
    """Check source caps individually so soup seams are not exterior boundaries."""
    raw_bases, distance = [], 0.0
    for original, kind, face in zip(originals, kinds, faces, strict=True):
        if kind != "horizontal" or face.rings[0][0, 2] != ground:
            continue
        raw = _polygon([r[:, :2] for r in original.rings])
        normalized = _polygon([r[:, :2] for r in face.rings])
        distance = max(distance, raw.boundary.hausdorff_distance(normalized.boundary))
        raw_bases.append(raw)
    raw_base = union_all(raw_bases)
    area = float(raw_base.symmetric_difference(base).area)
    if distance > movement + epsilon or area > (raw_base.length + base.length) * (
        movement + epsilon
    ):
        raise _PrismError(
            "precision_loss", "normalized base exceeds measured displacement budget"
        )
    return area, float(distance)


def _component(originals, welded, tolerance, angle, epsilon, origin):
    faces, kinds, heights, movement, planarity = _normalize_faces(
        originals, welded, tolerance, angle, epsilon, origin
    )
    faces, originals, kinds, duplicates = _node_and_deduplicate(
        faces, originals, kinds, epsilon
    )
    orientation, edge_count = _shell_orientation(faces)
    ground, base, levels, coverage_area, boundary_distance = _roof_cells(
        faces, kinds, epsilon
    )
    wall_area, wall_distance = _validate_walls(faces, kinds, levels, ground, epsilon)
    volume = sum(p.area * (z - ground) for z, p in levels)
    reference = faces[0].rings[0][0]
    shell_volume = abs(_shell_volume(faces, orientation, reference))
    source_volume = abs(_shell_volume(originals, orientation, reference))
    budget = _volume_budget(levels, ground, movement, epsilon)
    if (
        abs(shell_volume - volume) > _volume_budget(levels, ground, 0, epsilon)
        or abs(source_volume - volume) > budget
    ):
        raise _PrismError(
            "not_2_5d", "independent shell volume differs from prism volume"
        )
    source_area_residual, source_boundary_residual = _base_residual(
        originals, faces, kinds, ground, base, movement, epsilon
    )
    diagnostics = dict(
        duplicate_faces=duplicates,
        max_displacement=movement,
        max_planarity_residual=planarity,
        height_cluster_spans=[float(high - low) for low, high in heights],
        coverage_area_residual=max(coverage_area, wall_area),
        boundary_residual=max(boundary_distance, wall_distance),
        source_base_area_residual=float(source_area_residual),
        source_base_face_boundary_residual=source_boundary_residual,
        shell_volume=shell_volume,
        source_shell_volume=source_volume,
        prism_volume=float(volume),
        volume_residual=float(abs(shell_volume - volume)),
        source_volume_residual=float(abs(source_volume - volume)),
        volume_budget=budget,
        atomic_edges=edge_count,
        redundant_cap_area=0.0,
    )
    return _PrismComponent(ground, base, tuple(levels), diagnostics)


def _reject_overlapping_components(components):
    for i, first in enumerate(components):
        for j, second in enumerate(components[i + 1 :], start=i + 1):
            for z1, p1 in first.levels:
                for z2, p2 in second.levels:
                    if min(z1, z2) >= max(
                        first.ground_z, second.ground_z
                    ) and p1.intersects(p2):
                        raise _PrismError(
                            "not_2_5d",
                            f"shells {i} and {j} touch or overlap; internal/shared faces and cavities are unsupported",
                            component=i,
                        )


def _decompose_prisms(geometry, *, tolerance=0.002, wall_angle_tolerance=0.1):
    """Return deterministic column components after independent shell validation.

    The caller assumes ground-supported buildings: an isolated box at z=3..4
    has inferred ground 3. Missing walls, overhangs, cavities, slopes, internal
    caps and shared-wall solids fail explicitly. Duplicate complete faces are
    deduplicated. Heights and polygons use the original local coordinates;
    translation/SRS are returned separately and are never applied twice.
    """
    _check_parameters(tolerance, wall_angle_tolerance)
    records, (translation, srs) = _extract_faces(geometry)
    originals = _clean_faces(records, tolerance)
    all_points = np.vstack([ring for face in originals for ring in face.rings])
    origin = all_points.min(axis=0)
    local = all_points - origin
    epsilon = np.finfo(float).eps * max(float(np.ptp(local, axis=0).max()), 1.0) * 64
    if epsilon >= tolerance / 4:
        raise _PrismError(
            "precision_loss",
            "local extent cannot resolve the requested absolute tolerance",
        )
    originals = _replace_rings(originals, local)
    try:
        welded = _replace_rings(originals, _clusters(local, tolerance))
    except _PrismError as error:
        if error.location is not None:
            error.location = (np.asarray(error.location) + origin).tolist()
        raise
    results = []
    for component_index, members in enumerate(_vertex_groups(welded)):
        try:
            result = _component(
                [originals[i] for i in members],
                [welded[i] for i in members],
                tolerance,
                wall_angle_tolerance,
                epsilon,
                origin,
            )
        except _PrismError as error:
            error.component = component_index
            if error.location is not None:
                error.location = (
                    np.asarray(error.location) + origin[: len(error.location)]
                ).tolist()
            raise
        result.translation = tuple(float(v) for v in translation)
        result.srs = srs
        results.append(result)
    _reject_overlapping_components(results)
    for result in results:
        result.ground_z += float(origin[2])
        result.base = normalize(translate(result.base, origin[0], origin[1]))
        result.levels = tuple(
            (float(z + origin[2]), normalize(translate(p, origin[0], origin[1])))
            for z, p in result.levels
        )
    return tuple(
        sorted(
            results,
            key=lambda result: (
                result.ground_z,
                tuple(z for z, _ in result.levels),
                _geometry_key(result.base),
            ),
        )
    )
