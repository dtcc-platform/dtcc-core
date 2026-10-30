"""Minimal source-owned planar domains and their jointly simplified edge graph."""

import heapq
from collections import defaultdict
from dataclasses import replace

import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

from ...common.progress import report_progress
from ..geometry_policy import COINCIDENCE, canonical_profile_heights
from .roof_interfaces import RoofAssemblyError
from .shared_edges import SharedEdges
from .surface_faces import SurfaceFace, boundary_rings

TOLERANCE = COINCIDENCE


def planar_face(rings, *, kind, source, markers, normal, face_id=0, plane_origin=None):
    """Construct a planar domain and retain its canonical shared 3D ring vertices."""
    normal = np.asarray(normal, dtype=float)
    normal /= np.linalg.norm(normal)
    origin = np.asarray(
        rings[0][0] if plane_origin is None else plane_origin, dtype=float
    )
    if kind == "wall":
        tangent = np.array([-normal[1], normal[0], 0.0])
        tangent /= np.linalg.norm(tangent)
        basis = np.array([tangent, [0.0, 0.0, 1.0]])
    else:
        directions = np.diff(np.asarray(rings[0]), axis=0)
        tangent = directions[np.argmax(np.linalg.norm(directions, axis=1))]
        tangent -= normal * (normal @ tangent)
        tangent /= np.linalg.norm(tangent)
        basis = np.array([tangent, np.cross(normal, tangent)])
    local = [(np.asarray(r) - origin) @ basis.T for r in rings]
    polygon = orient(Polygon(local[0], local[1:]), sign=1.0)
    if not polygon.is_valid or polygon.area <= 1e-12:
        raise ValueError("Invalid compact planar domain")
    face = SurfaceFace(
        face_id,
        kind,
        source,
        sorted(markers),
        np.array([], dtype=int),
        origin,
        basis,
        normal,
        polygon,
        [],
        [],
        float(max(np.max(abs((np.asarray(r) - origin) @ normal)) for r in rings)),
    )
    # Keep the canonical shared 3D endpoints; projecting a numerically reconciled
    # junction onto each incident plane independently would separate its copies.
    face._geometric_rings = [np.asarray(r).copy() for r in rings]
    return face


def face_coordinates(face):
    """Return canonical 3D rings, reconstructing them from the face plane if needed."""
    if hasattr(face, "_geometric_rings"):
        return face._geometric_rings
    return [
        face.origin + np.asarray(r.coords)[:-1] @ face.basis
        for r in [face.polygon.exterior, *face.polygon.interiors]
    ]


def merge_planar(faces):
    """Union edge-connected domains under one fixed supporting plane and owner.

    A point contact does not join planar faces. Overlaying such domains in an
    oblique chart can manufacture a tiny segment at the common vertex.
    """
    groups = []
    by_owner = defaultdict(list)
    for face in faces:
        for group in by_owner[(face.kind, face.source_key)]:
            root = group[0]
            if (
                root.normal @ face.normal > 1 - 1e-14
                and abs((face.origin - root.origin) @ root.normal) <= TOLERANCE
                and all(
                    np.max(abs((ring - root.origin) @ root.normal)) <= TOLERANCE
                    for ring in face_coordinates(face)
                )
            ):
                group.append(face)
                break
        else:
            group = [face]
            groups.append(group)
            by_owner[(face.kind, face.source_key)].append(group)
    connected = []
    for group in groups:
        root = group[0]
        polygons = []
        for f in group:
            rings = [(r - root.origin) @ root.basis.T for r in face_coordinates(f)]
            polygons.append(Polygon(rings[0], rings[1:]))
        tree = shapely.STRtree(polygons)
        adjacent = defaultdict(set)
        for i, polygon in enumerate(polygons):
            for j in tree.query(polygon, predicate="intersects"):
                j = int(j)
                if (
                    j < i
                    and polygon.boundary.intersection(polygons[j].boundary).length
                    > TOLERANCE
                ):
                    adjacent[i].add(j)
                    adjacent[j].add(i)
        remaining = set(range(len(group)))
        while remaining:
            first = min(remaining)
            component, todo = [], [first]
            remaining.remove(first)
            while todo:
                i = todo.pop()
                component.append(group[i])
                for j in sorted(adjacent[i]):
                    if j in remaining:
                        remaining.remove(j)
                        todo.append(j)
            connected.append(component)
    result = []
    for group in connected:
        root = group[0]
        if len(group) == 1:
            result.append(root)
            continue
        polygons = []
        for f in group:
            rings = [(r - root.origin) @ root.basis.T for r in face_coordinates(f)]
            polygons.append(Polygon(rings[0], rings[1:]))
        union = shapely.union_all(polygons)
        for p in getattr(union, "geoms", [union]):
            contributors = [
                f for f, q in zip(group, polygons) if q.intersection(p).area > 1e-12
            ]
            merged = replace(
                root,
                polygon=orient(p, sign=1.0),
                markers=sorted(set(m for f in contributors for m in f.markers)),
            )
            originals = np.concatenate(
                [r for f in contributors for r in face_coordinates(f)]
            )
            original_projection = (originals - merged.origin) @ merged.basis.T
            rings = []
            for ring in [merged.polygon.exterior, *merged.polygon.interiors]:
                uv = np.asarray(ring.coords)[:-1]
                xyz = merged.origin + uv @ merged.basis
                for k, x in enumerate(uv):
                    distances = np.linalg.norm(original_projection - x, axis=1)
                    nearest = int(np.argmin(distances))
                    if distances[nearest] <= TOLERANCE:
                        xyz[k] = originals[nearest]
                rings.append(xyz)
            merged._geometric_rings = rings
            result.append(merged)
    return result


def building_domains(
    patches, reference, *, allow_source_steps=False, base_z=None, block_bases=None
):
    """Use fitted domains and native or explicit bases, with exposed affine intervals.

    Canonical oriented roof edges supply exact incidence. No roof tessellation
    or vertical strip levels enter these domains.
    """
    n = len(patches)
    if block_bases is not None:
        bases = block_bases
    elif base_z is not None:
        if not np.isfinite(base_z):
            raise ValueError("Building base must be finite")
        bases = {p.block_id: float(base_z) for p in patches}
        below_base = {
            p.block_id
            for p in patches
            for ring in [p.polygon.exterior, *p.polygon.interiors]
            if np.any(p.height(np.asarray(ring.coords)) <= base_z)
        }
        if below_base:
            raise RoofAssemblyError(
                "Roof must lie above the building base", block_ids=below_base
            )
    else:
        bases = {}
        for i, patch in enumerate(patches):
            ids = reference.faces[np.asarray(reference.markers) == i]
            if len(ids):
                level = float(reference.vertices[ids, 2].min())
                bases[patch.block_id] = min(bases.get(patch.block_id, level), level)
    # Heights are canonical only at common XY junctions, as in the assembler.
    at_xy = defaultdict(list)
    for i, p in enumerate(patches):
        for ring in [p.polygon.exterior, *p.polygon.interiors]:
            for xy in np.asarray(ring.coords)[:-1]:
                at_xy[tuple(xy)].append((float(p.height(xy)), i))
    heights = {}
    slopes = {i: p.coefficients[:2] for i, p in enumerate(patches)}
    for xy, values in at_xy.items():
        try:
            canonical = canonical_profile_heights(values, slopes)
        except ValueError as exc:
            raise RoofAssemblyError(
                str(exc), block_ids=[patches[i].block_id for z, i in values]
            ) from exc
        heights.update({(xy, i): z for i, z in canonical.items()})

    def height(i, xy):
        return heights.get((tuple(xy), i), float(patches[i].height(xy)))

    roofs = []
    for i, p in enumerate(patches):
        rings = []
        for ring in [p.polygon.exterior, *p.polygon.interiors]:
            xy = np.asarray(ring.coords)[:-1]
            rings.append(np.c_[xy, [height(i, x) for x in xy]])
        roofs.append(
            planar_face(
                rings,
                kind="roof",
                source=p.source_key,
                markers=[n + i],
                normal=[-p.coefficients[0], -p.coefficients[1], 1.0],
                plane_origin=[*p.origin, p.coefficients[2]],
            )
        )
    # The subdivision stage supplies fully noded oriented rings. Incidence is
    # combinatorial: distance to a short edge midpoint is not ownership.
    incidence = defaultdict(list)
    for i, p in enumerate(patches):
        for ring in [p.polygon.exterior, *p.polygon.interiors]:
            xy = np.asarray(ring.coords)
            for a, b in zip(xy[:-1], xy[1:]):
                key = tuple(sorted((tuple(a), tuple(b))))
                if key[0] == key[1]:
                    raise RoofAssemblyError(
                        "Zero-length canonical roof edge", block_ids=[p.block_id]
                    )
                incidence[key].append((i, a, b))
    walls = []
    for key, incident in sorted(incidence.items()):
        if not 1 <= len(incident) <= 2:
            raise RoofAssemblyError(
                "Nonmanifold canonical roof edge",
                block_ids=[patches[i].block_id for i, a, b in incident],
            )
        if len(incident) == 2 and np.array_equal(incident[0][1], incident[1][1]):
            raise RoofAssemblyError(
                "Canonical roof edge has inconsistent oriented incidence",
                block_ids=[patches[i].block_id for i, a, b in incident],
            )
        high, a, b = max(
            incident,
            key=lambda item: (
                height(item[0], item[1]) + height(item[0], item[2]),
                -item[0],
            ),
        )
        low = next((i for i, _, _ in incident if i != high), None)
        top = np.array([height(high, a), height(high, b)])
        bottom = (
            np.full(2, bases[patches[high].block_id])
            if low is None
            else np.array([height(low, a), height(low, b)])
        )
        if low is not None:
            if np.max(abs(top - bottom)) <= COINCIDENCE:
                continue
            if np.min(top - bottom) < -TOLERANCE:
                raise RoofAssemblyError(
                    "Unnoded compact height crossing",
                    block_ids=[patches[high].block_id],
                )
            if (
                patches[high].source_key == patches[low].source_key
                and not allow_source_steps
            ):
                raise RoofAssemblyError(
                    f"Discontinuous source roof {patches[high].source_key}: "
                    f"patches {patches[high].patch_id}/{patches[low].patch_id}, "
                    f"endpoint height differences {(top - bottom).tolist()}",
                    block_ids=[patches[high].block_id],
                )
        direction = b - a
        # Every exterior and hole ring has occupied material on its left.
        # The exposed wall points to its right, without an offset-point probe.
        normal = np.array([direction[1], -direction[0], 0.0])
        ring = np.array([[*a, bottom[0]], [*b, bottom[1]], [*b, top[1]], [*a, top[0]]])
        ring = np.array(
            [
                x
                for j, x in enumerate(ring)
                if np.linalg.norm(x - ring[j - 1]) > TOLERANCE
            ]
        )
        walls.append(
            planar_face(
                [ring],
                kind="wall",
                source=patches[high].source_key,
                markers=[high],
                normal=normal,
            )
        )
    return merge_planar(roofs + walls)


def geometric_graph(faces):
    """Node all face loops, then remove a vertex only from every incident loop.

    Tests use numerical 3D collinearity and incidence, with no detail tolerance.
    Equal XY at different elevations remains distinct.
    """
    vertices, lookup, boundaries = [], {}, {}
    original_edges = {}

    def vertex(p):
        key = tuple(np.rint(p / TOLERANCE).astype(np.int64))
        # Inspect neighbouring hash cells to avoid grid-boundary mismatches.
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for i in lookup.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        if np.linalg.norm(vertices[i] - p) <= TOLERANCE:
                            return i
        i = len(vertices)
        vertices.append(p)
        lookup.setdefault(key, []).append(i)
        return i

    for i, face in enumerate(faces):
        face.id = i
        rings = []
        for coordinates in face_coordinates(face):
            ids = [vertex(x) for x in coordinates]
            rings.append(ids)
            for k, (a, b) in enumerate(zip(ids, ids[1:] + ids[:1])):
                original_edges[(i, a, b)] = np.asarray(
                    [coordinates[k], coordinates[(k + 1) % len(ids)]]
                )
        boundaries[i] = [(a, b) for r in rings for a, b in zip(r, r[1:] + r[:1])]
    v = np.asarray(vertices, dtype=float).reshape(-1, 3)
    # Node T junctions even when the endpoint belongs to a shorter vertical edge.
    from scipy.spatial import cKDTree

    tree = cKDTree(v)
    chains = {}
    support = {}
    for i, boundary in boundaries.items():
        edges = []
        for a, b in boundary:
            key = tuple(sorted((a, b)))
            if a == b:
                raise ValueError(
                    f"Shared geometry face {i} has a zero-length boundary edge"
                )
            if key not in chains:
                start, end = v[list(key)]
                d = end - start
                candidates = np.asarray(
                    tree.query_ball_point(
                        (start + end) / 2, np.linalg.norm(d) / 2 + TOLERANCE
                    ),
                    dtype=int,
                )
                t = (v[candidates] - start) @ d / (d @ d)
                residual = np.linalg.norm(
                    v[candidates] - (start + t[:, None] * d), axis=1
                )
                selected = (t >= -1e-12) & (t <= 1 + 1e-12) & (residual <= TOLERANCE)
                chains[key] = [
                    int(j) for _, j in sorted(zip(t[selected], candidates[selected]))
                ]
            ids = chains[key] if a == key[0] else chains[key][::-1]
            original_a, original_b = original_edges[(i, a, b)]
            direction = v[b] - v[a]
            parameters = np.clip(
                (v[ids] - v[a]) @ direction / (direction @ direction), 0.0, 1.0
            )
            samples = original_a + parameters[:, None] * (original_b - original_a)
            for k, (left, right) in enumerate(zip(ids[:-1], ids[1:])):
                support[(i, left, right)] = samples[k : k + 2].copy()
            edges.extend(zip(ids[:-1], ids[1:]))
        faces[i].rings = boundary_rings(edges)
        faces[i].rings.sort(
            key=lambda r: -Polygon((v[r] - faces[i].origin) @ faces[i].basis.T).area
        )
        if (i + 1) % 100 == 0 or i + 1 == len(faces):
            report_progress(
                percent=70 * (i + 1) / len(faces),
                message=f"Reconciled shared edges for {i + 1}/{len(faces)} faces",
            )
    report_progress(percent=70, message="Simplifying shared LOD2 edge graph...")
    raw_count = len(vertices)
    raw_per_face = {str(f.id): sum(map(len, f.rings)) for f in faces}
    removed = []
    support_limited = set()
    occurrences = defaultdict(list)
    for f in faces:
        for r in f.rings:
            for x in r:
                occurrences[x].append((f.id, r))
    pending = sorted(occurrences)
    heapq.heapify(pending)
    queued = set(pending)
    while pending:
        x = heapq.heappop(pending)
        queued.discard(x)
        uses = occurrences.get(x, [])
        if not uses or any(
            getattr(faces[i], "_fixed_boundary", False) for i, _ in uses
        ):
            continue
        neighbours = set()
        merged_support = []
        for fi, r in uses:
            k = r.index(x)
            a, b = r[k - 1], r[(k + 1) % len(r)]
            neighbours.update((a, b))
            d = v[b] - v[a]
            length2 = float(d @ d)
            t = float((v[x] - v[a]) @ d / length2) if length2 else -1.0
            if (
                len(r) <= 3
                or not 0 < t < 1
                or np.linalg.norm(v[x] - (v[a] + t * d)) > TOLERANCE
            ):
                break
            # Bound total displacement from the original admitted boundary,
            # including samples removed earlier. Local epsilon tests alone can
            # accumulate into a larger change as neighbouring vertices vanish.
            samples = np.vstack([support[(fi, a, x)], support[(fi, x, b)]])
            parameters = np.clip((samples - v[a]) @ d / length2, 0.0, 1.0)
            if (
                np.max(
                    np.linalg.norm(samples - (v[a] + parameters[:, None] * d), axis=1)
                )
                > TOLERANCE
            ):
                support_limited.add(x)
                break
            merged_support.append((fi, a, b, samples))
        else:
            if len(neighbours) != 2:
                continue
            for fi, a, b, samples in merged_support:
                del support[(fi, a, x)], support[(fi, x, b)]
                support[(fi, a, b)] = samples
            for _, r in uses:
                r.remove(x)
            del occurrences[x]
            removed.append(x)
            # Only immediate neighbours change their collinearity predicate.
            for neighbour in sorted(neighbours):
                if neighbour not in queued:
                    heapq.heappush(pending, neighbour)
                    queued.add(neighbour)
    report_progress(percent=90, message="Registering shared LOD2 edges...")
    reasons = []
    occurrences = defaultdict(list)
    for f in faces:
        for r in f.rings:
            for k, x in enumerate(r):
                occurrences[x].append((f, r[k - 1], r[(k + 1) % len(r)]))
    for x, uses in sorted(occurrences.items()):
        kinds = {f.kind for f, _, _ in uses}
        sources = {f.source_key for f, _, _ in uses if f.source_key is not None}
        corners, neighbours = [], set()
        for f, a, b in uses:
            d = v[b] - v[a]
            residual = np.linalg.norm(np.cross(v[x] - v[a], d)) / np.linalg.norm(d)
            if residual > TOLERANCE:
                corners.append(f.id)
            neighbours.update((a, b))
        if "terrain" in kinds and len(kinds) > 1:
            reason = "necessary terrain contact"
        elif len(sources) > 1:
            reason = "source interface"
        elif len(neighbours) > 2 or len(uses) > 2:
            reason = "height crossing/junction"
        else:
            reason = "geometric corner/crease"
        fixed = [f.id for f, _, _ in uses if getattr(f, "_fixed_boundary", False)]
        reasons.append(
            dict(
                vertex=x,
                xyz=v[x].tolist(),
                reason=reason,
                incident_faces=[f.id for f, _, _ in uses],
                corner_faces=corners,
                fixed_terrain_faces=fixed,
                neighbours=sorted(neighbours),
                collinear_retained_because=(
                    "fixed outside-contact terrain triangulation"
                    if fixed
                    else "total original boundary displacement bound"
                    if x in support_limited
                    else "junction changes edge incidence"
                    if not corners
                    else "corner in incident faces " + str(corners)
                ),
            )
        )
    for f in faces:
        f.boundary = [(a, b) for r in f.rings for a, b in zip(r, r[1:] + r[:1])]
    registry = SharedEdges(v, {f.id: f.boundary for f in faces}, tolerance=TOLERANCE)
    return registry, dict(
        raw_vertices=raw_count,
        raw_face_boundary_vertices=raw_per_face,
        removed_redundant_vertices=len(removed),
        removed=[
            dict(
                vertex=i, xyz=v[i].tolist(), reason="jointly redundant collinear sample"
            )
            for i in removed
        ],
        retained_vertices=reasons,
    )


def terrain_domains(reference):
    """Retain the sampled terrain triangles as immutable geometric faces.

    Coplanar triangles need not form a simple polygon when a contact star
    pinches at a vertex. Keeping their original faces avoids a competing terrain
    union and preserves every original terrain triangle outside that star.
    """
    from .surface_faces import _make_face

    markers = np.asarray(reference.markers)
    contact = getattr(reference, "terrain_contact_vertices", None)
    if contact is None:
        contact = np.intersect1d(
            np.unique(reference.faces[markers < 0]),
            np.unique(reference.faces[markers >= 0]),
        )
    local = np.any(np.isin(reference.faces, contact), axis=1)
    result = []
    for i in np.flatnonzero(markers < 0):
        triangle = reference.vertices[reference.faces[i]]
        normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
        normal /= np.linalg.norm(normal)
        face = _make_face(
            reference,
            [i],
            triangle[0],
            normal,
            face_id=0,
            kind="terrain",
            source_key=None,
        )
        face._fixed_boundary = not local[i]
        face._geometric_rings = [reference.vertices[r].copy() for r in face.rings]
        result.append(face)
    return result, dict(
        contact_vertices=contact.tolist(),
        outside_star_face_ids=np.flatnonzero(~local & (markers < 0)).tolist(),
        attachment_star_triangles=int(np.sum(local & (markers < 0))),
        outside_star_triangles=int(np.sum(~local & (markers < 0))),
    )
