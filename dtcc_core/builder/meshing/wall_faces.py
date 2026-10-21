"""Geometry-preserving triangulation of connected, source-owned planar walls.

The exposed panels and their orientation come from the shared native shell
builder. Only their tessellation is replaced here. All eave and terrain edges
are immutable; subdivision of vertical junction edges is shared by every
incident wall face. No roof-family information is used.
"""

from collections import defaultdict
import time

import numpy as np
import shapely
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import orient

from ...model import Field, Mesh, SemanticRegion
from .dtcc_mesher_backend import _load_dtcc_mesher
from .roof_interfaces import XY_PRECISION


QUALITY_TARGET = 0.02
_MAX_ROUNDS = 32
_MAX_FACE_POINTS = 4096
_MAX_ADDED_POINTS = 100000


def triangle_quality(points, triangles):
    """Return normalized triangle quality in [0, 1], with equilateral triangles at 1."""
    xyz = np.asarray(points)[np.asarray(triangles, dtype=np.int64)]
    a, b = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    if xyz.shape[2] == 2:
        area2 = abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])
    else:
        area2 = np.linalg.norm(np.cross(a, b), axis=1)
    lengths = sum(
        np.sum((xyz[:, i] - xyz[:, (i + 1) % 3]) ** 2, axis=1) for i in range(3)
    )
    return np.divide(
        2 * np.sqrt(3) * area2, lengths, out=np.zeros_like(lengths), where=lengths > 0
    )


def _quality_record(q):
    return dict(
        triangles=len(q),
        minimum=float(q.min()),
        below_002=int(np.sum(q < QUALITY_TARGET)),
        below_001=int(np.sum(q < 0.01)),
    )


def _wall_groups(mesh, owners):
    """Edge-connected components in one oriented plane and one source."""
    v, f, markers = mesh.vertices, mesh.faces, np.asarray(mesh.markers)
    cross = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    normals = cross / np.linalg.norm(cross, axis=1)[:, None]
    walls = np.flatnonzero((markers >= 0) & (markers < len(owners)))
    edges = defaultdict(list)
    for i in walls:
        for a, b in zip(f[i], np.roll(f[i], -1)):
            edges[tuple(sorted((int(a), int(b))))].append(int(i))
    adjacent = defaultdict(list)
    for incident in edges.values():
        if len(incident) == 2:
            i, j = incident
            if (
                owners[markers[i]] == owners[markers[j]]
                and np.dot(normals[i], normals[j]) > 1 - 1e-8
            ):
                adjacent[i].append(j)
                adjacent[j].append(i)
    seen, groups = set(), []
    for i in walls:
        if i in seen:
            continue
        origin, normal = v[f[i, 0]], normals[i]
        todo, selected = [int(i)], []
        seen.add(i)
        while todo:
            j = todo.pop()
            selected.append(j)
            for k in adjacent[j]:
                # A fixed root plane prevents transitive tolerance drift.
                if (
                    k not in seen
                    and np.max(abs((v[f[k]] - origin) @ normal)) <= XY_PRECISION
                ):
                    seen.add(k)
                    todo.append(k)
        groups.append((selected, origin, normal))
    return groups, edges


def _face_geometry(mesh, group, origin, normal):
    tangent = np.array([-normal[1], normal[0], 0.0])
    tangent /= np.linalg.norm(tangent)
    basis = np.array([tangent, [0.0, 0.0, 1.0]])
    ids = np.unique(mesh.faces[group])
    projection = (mesh.vertices[ids] - origin) @ basis.T
    lookup = {int(k): i for i, k in enumerate(ids)}
    local_faces = np.array([[lookup[int(k)] for k in t] for t in mesh.faces[group]])
    polygons = shapely.polygons(projection[local_faces])
    polygon = shapely.union_all(polygons)
    if polygon.geom_type != "Polygon" or not polygon.is_valid:
        raise ValueError("Wall union is not one valid planar polygon")
    edges = defaultdict(list)
    for t in mesh.faces[group]:
        for a, b in zip(t, np.roll(t, -1)):
            edges[tuple(sorted((int(a), int(b))))].append((int(a), int(b)))
    boundary = [items[0] for items in edges.values() if len(items) == 1]
    degree = defaultdict(int)
    for a, b in boundary:
        degree[a] += 1
        degree[b] += 1
    if any(n != 2 for n in degree.values()):
        raise ValueError("Wall boundary has a pinched or branching vertex")
    # An interior vertex used by another surface must not disappear.
    outside = np.ones(len(mesh.faces), dtype=bool)
    outside[group] = False
    required = np.intersect1d(ids, np.unique(mesh.faces[outside]))
    if set(required) - set(degree):
        raise ValueError("Wall interior contains a required surface junction")
    return polygon, boundary, basis


def _triangulate_face(
    points, segments, polygon, *, max_rounds=_MAX_ROUNDS, max_points=_MAX_FACE_POINTS
):
    """Refine a fixed-boundary CDT with interior off-centres.

    DTCC Mesher supplies the initial constrained triangulation. Subsequent
    SciPy Delaunay tessellations are accepted only when every fixed boundary
    segment survives, with exact polygon coverage (including holes). This
    avoids the backends' segment-splitting quality refinement. Unsupported
    constraints retain the valid panel triangulation through the caller.
    """
    from scipy.spatial import Delaunay

    dm = _load_dtcc_mesher()
    original = np.array(points, copy=True)
    boundary = np.array(segments, dtype=np.int64, copy=True)
    points, segments = original.copy(), boundary.copy()
    seed = np.array([polygon.representative_point().coords[0]])
    for iteration in range(max_rounds + 1):
        if len(points) > max_points:
            raise ValueError(f"Wall point budget exceeded ({max_points})")
        if iteration == 0:
            raw = dm.mesh(
                dm.CoverageGraph(points, boundary, seed, np.array([0], dtype=np.int32)),
                options=dm.MeshingOptions(refine=False, acute_protection="none"),
            )
            vertices, triangles = (
                np.asarray(raw.points),
                np.asarray(raw.triangles, dtype=np.int64),
            )
            if vertices.shape != points.shape or not np.array_equal(vertices, points):
                raise ValueError(
                    "Wall CDT attempted uncoordinated boundary subdivision"
                )
        else:
            vertices = points
            triangles = Delaunay(points).simplices
            centers = vertices[triangles].mean(axis=1)
            triangles = triangles[
                shapely.intersects_xy(polygon, centers[:, 0], centers[:, 1])
            ]
        _check_replacement(polygon, vertices, triangles, boundary)
        q = triangle_quality(vertices, triangles)
        bad = np.flatnonzero(q < QUALITY_TARGET * 1.25)
        if not len(bad):
            return vertices, triangles, iteration, []
        additions, immutable = [], []
        for ti in bad:
            t = triangles[ti]
            xyz = vertices[t]
            length = np.linalg.norm(np.roll(xyz, -1, axis=0) - xyz, axis=1)
            k = int(np.argmin(length))
            x, y, z = xyz[k], xyz[(k + 1) % 3], xyz[(k + 2) % 3]
            middle = (x + y) / 2
            inward = np.array([-(y - x)[1], (y - x)[0]]) / length[k]
            if np.dot(z - middle, inward) < 0:
                inward = -inward
            candidate = middle + 0.866025403784 * length[k] * inward
            # Off-centres may lie beyond this triangle, but stay in the face.
            if not polygon.contains(Point(candidate)):
                immutable.append(int(ti))
                continue
            if (
                polygon.boundary.distance(Point(candidate)) <= 1e-8
                or np.linalg.norm(vertices - candidate, axis=1).min() <= 1e-8
            ):
                immutable.append(int(ti))
                continue
            if (
                additions
                and np.linalg.norm(np.array(additions) - candidate, axis=1).min()
                <= 1e-8
            ):
                continue
            additions.append(candidate)
        if not additions:
            return vertices, triangles, iteration, immutable
        if iteration == max_rounds:
            raise ValueError(f"Wall refinement round budget exceeded ({max_rounds})")
        points = np.vstack([vertices, additions])
    raise AssertionError("Unreachable refinement state")


def _angle_constraint(polygon):
    angles = []
    oriented = orient(polygon, sign=1.0)
    for ring in [oriented.exterior, *oriented.interiors]:
        xy = np.asarray(ring.coords)[:-1]
        a, b = np.roll(xy, 1, axis=0) - xy, np.roll(xy, -1, axis=0) - xy
        cross = (-a[:, 0]) * b[:, 1] - (-a[:, 1]) * b[:, 0]
        cosine = np.clip(
            np.sum(a * b, axis=1)
            / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)),
            -1.0,
            1.0,
        )
        for i in np.flatnonzero(cross > 1e-14):
            theta = float(np.arccos(cosine[i]))
            angles.append((theta, xy[i]))
    if not angles:
        return {}
    theta, point = min(angles, key=lambda a: a[0])
    return dict(
        minimum_input_angle_degrees=float(np.degrees(theta)),
        corner_quality_upper_bound=float(
            np.sqrt(3) * np.sin(theta) / (2 - np.cos(theta))
        ),
        corner_local=point.tolist(),
    )


def _check_replacement(polygon, points, triangles, boundary):
    polygons = shapely.polygons(points[triangles])
    union = shapely.union_all(polygons)
    error = float(polygon.symmetric_difference(union).area)
    overlap = max(0.0, float(shapely.area(polygons).sum() - union.area))
    tolerance = max(1e-10, polygon.area * 1e-10)
    if error > tolerance or overlap > tolerance:
        raise ValueError(
            f"Wall coverage changed: loss={error:g}, overlap={overlap:g} m²"
        )
    all_edges = np.sort(
        np.concatenate(
            [triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]
        ),
        axis=1,
    )
    unique, counts = np.unique(all_edges, axis=0, return_counts=True)
    actual = {tuple(e) for e in unique[counts == 1]}
    if actual != {tuple(sorted(e)) for e in boundary}:
        raise ValueError("Wall CDT did not preserve every external boundary segment")
    if np.any(counts > 2) or triangle_quality(points, triangles).min() <= 1e-12:
        raise ValueError("Invalid replacement wall triangulation")
    return error, overlap


def refine_wall_faces(
    mesh,
    owners,
    *,
    max_rounds=_MAX_ROUNDS,
    max_face_points=_MAX_FACE_POINTS,
    max_added_points=_MAX_ADDED_POINTS,
):
    """Return a wall-only replacement and face/source/contributor diagnostics.

    Wall markers retain a representative patch of their source. The per-face
    ``catalogue_wall_face`` field and audit explicitly map spanning triangles
    to all contributing patches. Non-wall markers and triangles are untouched.
    Quality fallbacks keep the original valid geometry and never reject fits.
    """
    start = time.perf_counter()
    groups, wall_edges = _wall_groups(mesh, owners)
    old_q = triangle_quality(mesh.vertices, mesh.faces)
    records, geometry, requests = [], {}, defaultdict(set)
    face_group = {}
    for gi, (group, origin, normal) in enumerate(groups):
        for fi in group:
            face_group[fi] = gi
        markers = sorted(set(map(int, np.asarray(mesh.markers)[group])))
        record = dict(
            face_id=gi,
            source_key=owners[markers[0]],
            contributing_markers=markers,
            marker=markers[0],
            method="panels",
            original=_quality_record(old_q[group]),
            fallback_reason=None,
            refinement_rounds=0,
            added_vertices=0,
        )
        records.append(record)
        if not record["original"]["below_002"]:
            continue
        try:
            polygon, boundary, basis = _face_geometry(mesh, group, origin, normal)
            geometry[gi] = (polygon, boundary, basis)
            record["angle_constraint"] = _angle_constraint(polygon)
            if record["original"]["below_002"]:
                width = polygon.bounds[2] - polygon.bounds[0]
                # A genuinely narrow facade cannot meet the target with long
                # corner edges fixed. Split vertical wall junctions only.
                incident_lengths = defaultdict(list)
                for a, b in boundary:
                    va, vb = mesh.vertices[[a, b]]
                    if not np.array_equal(va[:2], vb[:2]):
                        length = float(np.linalg.norm(vb - va))
                        incident_lengths[a].append(length)
                        incident_lengths[b].append(length)
                for a, b in boundary:
                    va, vb = mesh.vertices[[a, b]]
                    if np.array_equal(va[:2], vb[:2]):
                        count = int(np.ceil(abs(vb[2] - va[2]) / (25 * width)))
                        if count > max_face_points:
                            raise ValueError(
                                f"Vertical junction requires {count} strips; wall point budget {max_face_points}"
                            )
                        key = tuple(sorted((a, b)))
                        if count > 1:
                            for k in range(1, count):
                                requests[key].add(
                                    k / count if a == key[0] else 1 - k / count
                                )
                        # A tiny eave/terrain edge at a corner also forces
                        # grading on the adjoining vertical boundary. Otherwise
                        # an interior off-centre is trapped against a long edge.
                        height = abs(vb[2] - va[2])
                        for endpoint in (a, b):
                            if not incident_lengths[endpoint]:
                                continue
                            spacing = 16 * min(incident_lengths[endpoint])
                            while spacing < height / 2:
                                t = spacing / height
                                requests[key].add(t if endpoint == key[0] else 1 - t)
                                spacing *= 2

        except (ValueError, RuntimeError) as exc:
            record["fallback_reason"] = str(exc)
    # A split is one global edge operation, shared across ownership and plane
    # boundaries. Vertical eave/terrain edges cannot occur in a valid shell.
    additions, split_edges, components = [], {}, {i: {i} for i in range(len(groups))}
    for edge, fractions in requests.items():
        incident = wall_edges[edge]
        all_incident = np.sum(np.isin(mesh.faces, edge).sum(axis=1) == 2)
        affected = {face_group[i] for i in incident}
        if all_incident != len(incident) or len(incident) != 2:
            for gi in affected:
                records[gi]["fallback_reason"] = (
                    "Vertical junction is not shared exclusively by two wall triangles"
                )
            continue
        cluster = set().union(*(components[i] for i in affected))
        for gi in cluster:
            components[gi] = cluster
        chain = [edge[0]]
        length = np.linalg.norm(mesh.vertices[edge[1]] - mesh.vertices[edge[0]])
        cuts = [0.0, *sorted(fractions), 1.0]
        # Proposals from two incident faces may nearly coincide. They are
        # optional points on the same straight curve, not geometric vertices.
        # Coalesce competing proposals to keep consecutive spacings graded.
        while len(cuts) > 2:
            gaps = np.diff(cuts)
            ratio = np.minimum(gaps[:-1], gaps[1:]) / np.maximum(gaps[:-1], gaps[1:])
            k = int(np.argmin(ratio))
            if ratio[k] >= 0.15:
                break
            del cuts[k + 1]
        if len(additions) + len(cuts) - 2 > max_added_points:
            for gi in cluster:
                records[gi]["fallback_reason"] = (
                    f"City wall vertex budget exceeded ({max_added_points})"
                )
            continue
        accepted = [0.0]
        for t in cuts[1:-1]:
            # These are optional new vertices on an unchanged straight edge.
            # Do not introduce numerically coincident junction subdivisions.
            if min(t - accepted[-1], 1 - t) * length <= 1e-6:
                continue
            accepted.append(t)
            chain.append(len(mesh.vertices) + len(additions))
            additions.append(
                mesh.vertices[edge[0]]
                + t * (mesh.vertices[edge[1]] - mesh.vertices[edge[0]])
            )
        chain.append(edge[1])
        split_edges[edge] = chain
    vertices = np.vstack([mesh.vertices, np.array(additions).reshape(-1, 3)])
    replacements = {}
    processed = set()
    for gi, cluster in components.items():
        if gi in processed:
            continue
        processed |= cluster
        pending, cluster_vertices = {}, []
        try:
            for ci in sorted(cluster):
                record = records[ci]
                if record["fallback_reason"]:
                    raise ValueError(record["fallback_reason"])
                group, origin, normal = groups[ci]
                if ci not in geometry:
                    if len(cluster) == 1:
                        continue
                    geometry[ci] = _face_geometry(mesh, group, origin, normal)
                    record["angle_constraint"] = _angle_constraint(geometry[ci][0])
                polygon, boundary, basis = geometry[ci]
                changed_boundary = any(
                    tuple(sorted(e)) in split_edges for e in boundary
                )
                if not changed_boundary and not record["original"]["below_002"]:
                    continue
                expanded = []
                for a, b in boundary:
                    key = tuple(sorted((a, b)))
                    chain = split_edges.get(key, [key[0], key[1]])
                    if a != key[0]:
                        chain = chain[::-1]
                    expanded.extend(zip(chain, chain[1:]))
                ids = np.unique(expanded)
                lookup = {int(k): i for i, k in enumerate(ids)}
                segments = np.array([[lookup[a], lookup[b]] for a, b in expanded])
                points = (vertices[ids] - origin) @ basis.T
                local, triangles, rounds, immutable = _triangulate_face(
                    points,
                    segments,
                    polygon,
                    max_rounds=max_rounds,
                    max_points=max_face_points,
                )
                error, overlap = _check_replacement(polygon, local, triangles, segments)
                interior = origin + local[len(ids) :] @ basis
                interior_ids = np.arange(
                    len(vertices) + len(cluster_vertices),
                    len(vertices) + len(cluster_vertices) + len(interior),
                )
                if (
                    len(vertices)
                    + len(cluster_vertices)
                    + len(interior)
                    - len(mesh.vertices)
                    > max_added_points
                ):
                    raise ValueError(
                        f"City wall vertex budget exceeded ({max_added_points})"
                    )
                cluster_vertices.extend(interior)
                mapped = np.r_[ids, interior_ids][triangles]
                # Align to the original outward orientation.
                if np.dot(np.cross(basis[0], basis[1]), normal) < 0:
                    mapped = mapped[:, [0, 2, 1]]
                q = triangle_quality(local, triangles)
                if np.any(q < QUALITY_TARGET):
                    record["rejected_quality"] = _quality_record(q)
                    record["immutable_triangles"] = immutable
                    record["reproducer"] = dict(
                        points=points.tolist(),
                        boundary=segments.tolist(),
                        polygon=shapely.to_wkt(polygon, rounding_precision=15),
                    )
                    raise ValueError(
                        "Fixed wall boundary or input angles prevent the quality target with interior refinement"
                    )
                # Coverage CDT is counterclockwise in its local plane.
                xy = local[triangles]
                reverse = (
                    (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1])
                    - (xy[:, 1, 1] - xy[:, 0, 1]) * (xy[:, 2, 0] - xy[:, 0, 0])
                ) < 0
                mapped[reverse] = mapped[reverse][:, [0, 2, 1]]
                pending[ci] = mapped
                record.update(
                    method="connected_planar_cdt",
                    transaction_accepted=True,
                    replacement=_quality_record(q),
                    refinement_rounds=rounds,
                    added_vertices=len(interior),
                    boundary_vertices_added=len(ids) - len(np.unique(boundary)),
                    coverage_error_m2=error,
                    overlap_m2=overlap,
                    plane_tolerance_m=XY_PRECISION,
                    maximum_original_plane_error_m=float(
                        np.max(
                            abs(
                                (mesh.vertices[np.unique(mesh.faces[group])] - origin)
                                @ normal
                            )
                        )
                    ),
                    immutable_triangles=immutable,
                )
            if (
                len(vertices) + len(cluster_vertices) - len(mesh.vertices)
                > max_added_points
            ):
                raise ValueError(
                    f"City wall vertex budget exceeded ({max_added_points})"
                )
        except (ValueError, RuntimeError) as exc:
            for ci in cluster:
                records[ci].update(
                    method="panels",
                    fallback_reason=str(exc),
                    transaction_accepted=False,
                )
            continue
        vertices = np.vstack([vertices, np.array(cluster_vertices).reshape(-1, 3)])
        replacements.update(pending)
    keep = np.ones(len(mesh.faces), dtype=bool)
    for gi in replacements:
        keep[groups[gi][0]] = False
    face_arrays = [np.asarray(mesh.faces[keep], dtype=np.int64)]
    marker_arrays = [np.asarray(mesh.markers)[keep]]
    face_ids = [np.full(int(keep.sum()), -1, dtype=np.int64)]
    for gi, triangles in replacements.items():
        face_arrays.append(triangles)
        marker_arrays.append(
            np.full(len(triangles), records[gi]["marker"], dtype=np.int64)
        )
        face_ids.append(np.full(len(triangles), gi, dtype=np.int64))
    # Discard unused proposed junction vertices after a transactional fallback.
    faces = np.concatenate(face_arrays)
    used_new = np.unique(faces[faces >= len(mesh.vertices)])
    retained = np.r_[np.arange(len(mesh.vertices)), used_new]
    inverse = np.full(len(vertices), -1, dtype=np.int64)
    inverse[retained] = np.arange(len(retained))
    result = Mesh(
        vertices=vertices[retained],
        faces=inverse[faces],
        markers=np.concatenate(marker_arrays),
    )
    result.fields.append(
        Field(
            name="catalogue_wall_face",
            values=np.concatenate(face_ids),
            association="face",
            description="Connected wall face ID; -1 denotes an unchanged triangle. See roof_catalogue.wall_triangulation.faces for source and contributing patch markers.",
        )
    )
    for gi in replacements:
        result.regions.append(
            SemanticRegion(
                "urn:dtcc:wall-face",
                id=str(gi),
                indices=np.flatnonzero(result.fields[0].values == gi),
                attributes=dict(
                    source_key=records[gi]["source_key"],
                    contributing_patch_markers=records[gi]["contributing_markers"],
                    representative_wall_marker=records[gi]["marker"],
                ),
            )
        )
    audit = dict(
        method="connected_planar_cdt",
        seconds=time.perf_counter() - start,
        faces=records,
        replaced_faces=len(replacements),
        quality_fallbacks=sum(r["fallback_reason"] is not None for r in records),
        added_vertices=len(retained) - len(mesh.vertices),
        limits=dict(
            rounds=max_rounds,
            face_points=max_face_points,
            added_points=max_added_points,
        ),
        boundary_policy="Original eave and terrain triangles fixed; vertical wall junction subdivisions shared globally",
    )
    return result, audit


def separate_wall_mesh(mesh, num_patches):
    """Restore terrain plus one mesh per patch, retaining spanning-wall metadata."""
    locations = np.where(mesh.markers < 0, 0, mesh.markers % num_patches + 1)
    result = []
    for index in range(num_patches + 1):
        selected = locations == index
        faces = mesh.faces[selected]
        ids, inverse = np.unique(faces, return_inverse=True)
        part = Mesh(
            vertices=mesh.vertices[ids],
            faces=inverse.reshape(-1, 3),
            markers=mesh.markers[selected],
        )
        for field in mesh.fields:
            if field.association == "face":
                part.fields.append(
                    Field(
                        name=field.name,
                        values=field.values[selected],
                        description=field.description,
                        association="face",
                    )
                )
        local_index = np.full(len(mesh.faces), -1, dtype=np.int64)
        local_index[selected] = np.arange(int(selected.sum()))
        for region in mesh.regions:
            membership = local_index[region.indices]
            membership = membership[membership >= 0]
            if len(membership):
                part.regions.append(
                    SemanticRegion(
                        region.semantic_type,
                        id=region.id,
                        indices=membership,
                        attributes=dict(region.attributes),
                    )
                )
        result.append(part)
    return result
