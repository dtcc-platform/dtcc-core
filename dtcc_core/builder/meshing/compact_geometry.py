"""Minimal source-owned planar domains and their jointly simplified edge graph."""
from collections import defaultdict
from dataclasses import replace

import numpy as np
import shapely
from shapely.geometry import Polygon, Point
from shapely.geometry.polygon import orient

from .surface_faces import SurfaceFace, extract_surface_faces, boundary_rings
from .shared_edges import SharedEdges

TOLERANCE = 1e-8


def planar_face(rings, *, kind, source, markers, normal, face_id=0):
    """Construct a planar domain and retain its canonical shared 3D ring vertices."""
    normal = np.asarray(normal, dtype=float)
    normal /= np.linalg.norm(normal)
    origin = np.asarray(rings[0][0], dtype=float)
    directions = np.diff(np.asarray(rings[0]), axis=0)
    tangent = directions[np.argmax(np.linalg.norm(directions, axis=1))]
    tangent -= normal * (normal @ tangent)
    tangent /= np.linalg.norm(tangent)
    basis = np.array([tangent, np.cross(normal, tangent)])
    local = [(np.asarray(r)-origin) @ basis.T for r in rings]
    polygon = orient(Polygon(local[0], local[1:]), sign=1.)
    if not polygon.is_valid or polygon.area <= 1e-12:
        raise ValueError('Invalid compact planar domain')
    face = SurfaceFace(face_id, kind, source, sorted(markers), np.array([], dtype=int),
                       origin, basis, normal, polygon, [], [], 0.)
    # Keep the canonical shared 3D endpoints; projecting a numerically reconciled
    # junction onto each incident plane independently would separate its copies.
    face._geometric_rings = [np.asarray(r).copy() for r in rings]
    return face


def face_coordinates(face):
    """Return canonical 3D rings, reconstructing them from the face plane if needed."""
    if hasattr(face, "_geometric_rings"):
        return face._geometric_rings
    return [face.origin + np.asarray(r.coords)[:-1] @ face.basis
            for r in [face.polygon.exterior, *face.polygon.interiors]]


def merge_planar(faces):
    """Union compatible domains under one fixed supporting plane and owner."""
    groups = []
    for face in faces:
        for group in groups:
            root = group[0]
            if (root.kind == face.kind and root.source_key == face.source_key and
                root.normal @ face.normal > 1-1e-14 and
                max(np.max(abs((r-root.origin) @ root.normal)) for r in face_coordinates(face)) <= TOLERANCE):
                group.append(face)
                break
        else:
            groups.append([face])
    result = []
    for group in groups:
        root = group[0]
        if len(group) == 1:
            result.append(root)
            continue
        polygons = []
        for f in group:
            rings = [(r-root.origin) @ root.basis.T for r in face_coordinates(f)]
            polygons.append(Polygon(rings[0], rings[1:]))
        union = shapely.union_all(polygons)
        for p in getattr(union, 'geoms', [union]):
            contributors = [f for f, q in zip(group, polygons) if q.intersection(p).area > 1e-12]
            merged = replace(root, polygon=orient(p, sign=1.), markers=sorted(set(
                m for f in contributors for m in f.markers)))
            originals = np.concatenate([r for f in contributors for r in face_coordinates(f)])
            rings = []
            for ring in [merged.polygon.exterior, *merged.polygon.interiors]:
                xyz = merged.origin + np.asarray(ring.coords)[:-1] @ merged.basis
                for k, x in enumerate(xyz):
                    distances = np.linalg.norm(originals-x, axis=1)
                    nearest = int(np.argmin(distances))
                    if distances[nearest] <= 1e-7:
                        xyz[k] = originals[nearest]
                rings.append(xyz)
            merged._geometric_rings = rings
            result.append(merged)
    return result


def building_domains(patches, reference, *, allow_source_steps=False, base_z=None):
    """Use fitted domains and native or explicit bases, with exposed affine intervals.

    Shared edges use the same 1 micrometre height coincidence rule as the native
    assembler. No roof tessellation or vertical strip levels enter these domains.
    """
    n = len(patches)
    if base_z is not None:
        if not np.isfinite(base_z):
            raise ValueError('Building base must be finite')
        bases = {p.block_id: float(base_z) for p in patches}
        if any(p.height(xy) <= base_z for p in patches
               for ring in [p.polygon.exterior, *p.polygon.interiors] for xy in ring.coords):
            raise ValueError('Roof must lie above the building base')
    else:
        bases = {}
        for i, patch in enumerate(patches):
            ids = reference.faces[np.asarray(reference.markers) == i]
            if len(ids):
                level = float(reference.vertices[ids, 2].min())
                bases[patch.block_id] = min(bases.get(patch.block_id, level), level)
    polygons = [p.polygon for p in patches]
    boundaries = [p.boundary for p in polygons]
    tree = shapely.STRtree(boundaries)
    # Heights are canonical only at common XY junctions, as in the assembler.
    at_xy = defaultdict(list)
    for i, p in enumerate(patches):
        for ring in [p.polygon.exterior, *p.polygon.interiors]:
            for xy in np.asarray(ring.coords)[:-1]:
                at_xy[tuple(xy)].append((float(p.height(xy)), i))
    heights = {}
    for xy, values in at_xy.items():
        values.sort()
        canonical = values[0][0]
        for z, i in values:
            if z-canonical > 1e-6:
                canonical = z
            heights[(xy, i)] = canonical
    def height(i, xy):
        return heights.get((tuple(xy), i), float(patches[i].height(xy)))
    roofs = []
    for i, p in enumerate(patches):
        rings = []
        for ring in [p.polygon.exterior, *p.polygon.interiors]:
            xy = np.asarray(ring.coords)[:-1]
            rings.append(np.c_[xy, [height(i, x) for x in xy]])
        roofs.append(planar_face(rings, kind='roof', source=p.source_key,
                                 markers=[n+i], normal=[-p.coefficients[0], -p.coefficients[1], 1.]))
    graph = shapely.union_all(boundaries)
    lines = list(graph.geoms) if hasattr(graph, 'geoms') else [graph]
    walls = []
    for line in lines:
        xy = np.asarray(line.coords)
        for a, b in zip(xy[:-1], xy[1:]):
            midpoint = Point((a+b)/2)
            adjacent = [int(i) for i in tree.query(midpoint.buffer(TOLERANCE))
                        if boundaries[int(i)].distance(midpoint) <= TOLERANCE]
            if not 1 <= len(adjacent) <= 2:
                raise ValueError('Ambiguous compact source interface')
            high = max(adjacent, key=lambda i: (height(i, a)+height(i, b), -i))
            low = next((i for i in adjacent if i != high), None)
            top = np.array([height(high, a), height(high, b)])
            bottom = (np.full(2, bases[patches[high].block_id]) if low is None else
                      np.array([height(low, a), height(low, b)]))
            if low is not None:
                if np.max(abs(top-bottom)) <= 1e-6:
                    continue
                if np.min(top-bottom) < -TOLERANCE:
                    raise ValueError('Unnoded compact height crossing')
                if patches[high].source_key == patches[low].source_key and not allow_source_steps:
                    raise ValueError('Discontinuous compact source roof')
            direction = b-a
            normal = np.array([direction[1], -direction[0], 0.])
            # Polygon orientation determines which side of the edge is exposed.
            probe = (a+b)/2 + np.array([-direction[1], direction[0]])/np.linalg.norm(direction)*1e-6
            if not polygons[high].covers(Point(probe)):
                a, b = b, a
                top, bottom = top[::-1], bottom[::-1]
                normal = -normal
            ring = np.array([[*a, bottom[0]], [*b, bottom[1]], [*b, top[1]], [*a, top[0]]])
            ring = np.array([x for j, x in enumerate(ring) if np.linalg.norm(x-ring[j-1]) > TOLERANCE])
            walls.append(planar_face([ring], kind='wall', source=patches[high].source_key,
                                     markers=[high], normal=normal))
    return merge_planar(roofs + walls)


def geometric_graph(faces):
    """Node all face loops, then remove a vertex only from every incident loop.

    Tests use numerical 3D collinearity and incidence, with no detail tolerance.
    Equal XY at different elevations remains distinct.
    """
    vertices, lookup, boundaries = [], {}, {}
    def vertex(p):
        key = tuple(np.rint(p/TOLERANCE).astype(np.int64))
        # Inspect neighbouring hash cells to avoid grid-boundary mismatches.
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for i in lookup.get((key[0]+dx, key[1]+dy, key[2]+dz), []):
                        if np.linalg.norm(vertices[i]-p) <= TOLERANCE:
                            return i
        i = len(vertices)
        vertices.append(p)
        lookup.setdefault(key, []).append(i)
        return i
    for i, face in enumerate(faces):
        face.id = i
        rings = [[vertex(x) for x in r] for r in face_coordinates(face)]
        boundaries[i] = [(a, b) for r in rings for a, b in zip(r, r[1:]+r[:1])]
    v = np.asarray(vertices)
    # Node T junctions even when the endpoint belongs to a shorter vertical edge.
    for i, boundary in boundaries.items():
        edges = []
        for a, b in boundary:
            d = v[b]-v[a]
            t = (v-v[a]) @ d / (d @ d)
            residual = np.linalg.norm(v-(v[a]+t[:, None]*d), axis=1)
            ids = np.flatnonzero((t >= -1e-12) & (t <= 1+1e-12) & (residual <= TOLERANCE))
            ids = sorted(map(int, ids), key=lambda j: (t[j], j))
            edges.extend(zip(ids[:-1], ids[1:]))
        faces[i].rings = boundary_rings(edges)
        faces[i].rings.sort(key=lambda r: -Polygon((v[r]-faces[i].origin) @ faces[i].basis.T).area)
    raw_count = len(vertices)
    raw_per_face = {str(f.id):sum(map(len,f.rings)) for f in faces}
    removed = []
    while True:
        occurrences = defaultdict(list)
        for f in faces:
            for r in f.rings:
                for k, x in enumerate(r):
                    occurrences[x].append((f.id, r, k))
        discard = None
        for x, uses in sorted(occurrences.items()):
            neighbours = set()
            removable = True
            for _, r, k in uses:
                a, b = r[k-1], r[(k+1)%len(r)]
                neighbours.update((a, b))
                d = v[b]-v[a]
                t = (v[x]-v[a]) @ d / (d @ d)
                if (any(getattr(faces[use[0]], '_fixed_boundary', False) for use in uses) or
                    len(r) <= 3 or not 0 < t < 1 or
                    np.linalg.norm(v[x]-(v[a]+t*d)) > TOLERANCE):
                    removable = False
                    break
            if removable and len(neighbours) == 2:
                discard = x, uses
                break
        if discard is None:
            break
        x, uses = discard
        for _, r, _ in uses:
            r.remove(x)
        removed.append(x)
    reasons = []
    occurrences = defaultdict(list)
    for f in faces:
        for r in f.rings:
            for k, x in enumerate(r):
                occurrences[x].append((f, r[k-1], r[(k+1)%len(r)]))
    for x, uses in sorted(occurrences.items()):
        kinds = {f.kind for f, _, _ in uses}
        sources = {f.source_key for f, _, _ in uses if f.source_key is not None}
        corners, neighbours = [], set()
        for f, a, b in uses:
            d = v[b]-v[a]
            residual = np.linalg.norm(np.cross(v[x]-v[a], d))/np.linalg.norm(d)
            if residual > TOLERANCE:
                corners.append(f.id)
            neighbours.update((a, b))
        if 'terrain' in kinds and len(kinds) > 1:
            reason = 'necessary terrain contact'
        elif len(sources) > 1:
            reason = 'source interface'
        elif len(neighbours) > 2 or len(uses) > 2:
            reason = 'height crossing/junction'
        else:
            reason = 'geometric corner/crease'
        fixed = [f.id for f,_,_ in uses if getattr(f,'_fixed_boundary',False)]
        reasons.append(dict(vertex=x, xyz=v[x].tolist(), reason=reason,
                            incident_faces=[f.id for f, _, _ in uses],
                            corner_faces=corners, fixed_terrain_faces=fixed, neighbours=sorted(neighbours),
                            collinear_retained_because=('fixed outside-contact terrain triangulation' if fixed else
                                'junction changes edge incidence' if not corners else
                                'corner in incident faces '+str(corners))))
    for f in faces:
        f.boundary = [(a, b) for r in f.rings for a, b in zip(r, r[1:]+r[:1])]
    registry = SharedEdges(v, {f.id:f.boundary for f in faces}, tolerance=TOLERANCE)
    return registry, dict(raw_vertices=raw_count,raw_face_boundary_vertices=raw_per_face, removed_redundant_vertices=len(removed),
                          removed=[dict(vertex=i, xyz=v[i].tolist(), reason='jointly redundant collinear sample') for i in removed],
                          retained_vertices=reasons)


def terrain_domains(reference, extracted):
    """Permit retriangulation only in the original contact-vertex star.

    Fixed outside triangles constrain the edge graph, even across coplanar
    terrain. Thus a flat terrain cannot trigger a domain-wide remesh.
    """
    from .surface_faces import _make_face
    markers = np.asarray(reference.markers)
    contact = np.intersect1d(np.unique(reference.faces[markers < 0]),
                             np.unique(reference.faces[markers >= 0]))
    local = np.any(np.isin(reference.faces, contact), axis=1)
    result = []
    for face in extracted:
        if face.kind != 'terrain':
            continue
        eligible = face.original_faces[local[face.original_faces]]
        fixed = face.original_faces[~local[face.original_faces]]
        if not len(eligible) or not len(fixed):
            face._fixed_boundary = not len(eligible)
            result.append(face)
            continue
        for i in fixed:
            f = _make_face(reference, [i], face.origin, face.normal,
                           face_id=0, kind='terrain', source_key=None)
            f._fixed_boundary = True
            result.append(f)
        edges = defaultdict(list)
        for i in eligible:
            for a,b in zip(reference.faces[i],np.roll(reference.faces[i],-1)):
                edges[tuple(sorted((int(a),int(b))))].append(int(i))
        adjacent = defaultdict(set)
        for ids in edges.values():
            for i in ids:
                adjacent[i].update(ids)
        remaining = set(map(int,eligible))
        while remaining:
            root = min(remaining)
            group, todo = [], [root]
            remaining.remove(root)
            while todo:
                i = todo.pop()
                group.append(i)
                for j in sorted(adjacent[i]):
                    if j in remaining:
                        remaining.remove(j)
                        todo.append(j)
            result.append(_make_face(reference, group, face.origin, face.normal,
                                    face_id=0, kind='terrain', source_key=None))
    for f in result:
        f._geometric_rings = [reference.vertices[r].copy() for r in f.rings]
    return result, dict(contact_vertices=contact.tolist(),
                        outside_star_face_ids=np.flatnonzero(~local & (markers<0)).tolist(),
                        attachment_star_triangles=int(np.sum(local & (markers<0))),
                        outside_star_triangles=int(np.sum(~local & (markers<0))))
