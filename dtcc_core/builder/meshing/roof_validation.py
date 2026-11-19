"""Partition, constrained triangulation and surface topology checks."""

import numpy as np
import shapely
from shapely.ops import unary_union

from ..geometry_policy import COINCIDENCE
from .roof_interfaces import XY_PRECISION, RoofAssemblyError


def validate_terrain_coverage(mesh, polygons):
    """Check outside terrain coverage before any building faces are attached."""
    target = unary_union(polygons)
    triangles = shapely.polygons(np.asarray(mesh.vertices)[np.asarray(mesh.faces), :2])
    union = shapely.union_all(triangles)
    error = float(target.symmetric_difference(union).area)
    overlap = max(0.0, float(shapely.area(triangles).sum() - union.area))
    tolerance = target.length * COINCIDENCE + 64 * np.finfo(float).eps * target.area
    if (
        error > tolerance
        or overlap > tolerance
        or np.any(np.asarray(mesh.markers) != -2)
    ):
        raise RoofAssemblyError(
            f"Terrain subdivision coverage failed: missing/different={error:g}, "
            f"overlap={overlap:g}, numerical budget={tolerance:g} m2"
        )
    return dict(
        coverage_error_m2=error,
        overlap_m2=overlap,
        numerical_area_budget_m2=float(tolerance),
    )


def validate_partition(patches, target):
    """Check complete, nonoverlapping roof coverage and return area diagnostics."""
    polygons = [p.polygon for p in patches]
    if any(not p.is_valid or p.is_empty or p.area <= 0 for p in polygons):
        raise RoofAssemblyError("Invalid or empty roof patch")
    union = shapely.union_all(polygons, grid_size=XY_PRECISION)
    missing = shapely.symmetric_difference(target, union, grid_size=XY_PRECISION).area
    overlap = sum(p.area for p in polygons) - union.area
    tolerance = max(1e-7, target.area * 1e-9)
    report = dict(coverage_error_m2=float(missing), overlap_m2=float(max(0.0, overlap)))
    report["coverage_tolerance_m2"] = float(tolerance + target.length * XY_PRECISION)
    if missing > report["coverage_tolerance_m2"] or overlap > tolerance:
        raise RoofAssemblyError(f"Roof partition fails exact coverage: {report}")
    return report


def validate_triangulation(mesh, patches):
    """Check full triangle containment, marker transfer and constrained edges.

    Containment plus complete coverage detects a lost crease, even when its
    endpoints survive. Backend Steiner vertices inherit their triangle's plane.
    """
    markers = np.asarray(mesh.markers)
    points = np.asarray(mesh.vertices)[np.asarray(mesh.faces), :2]
    triangles = shapely.polygons(np.concatenate([points, points[:, :1]], axis=1))
    for marker, patch in enumerate(patches):
        selected = triangles[markers == marker]
        if not len(selected):
            raise RoofAssemblyError(f"Backend lost patch {marker}")
        union = shapely.union_all(selected)
        tolerance = max(1e-7, patch.polygon.area * 1e-8)
        if union.symmetric_difference(patch.polygon).area > tolerance:
            raise RoofAssemblyError(f"Backend lost roof constraints in patch {marker}")
        # Area sums also exclude overlapping triangles within a patch.
        if abs(float(shapely.area(selected).sum()) - patch.polygon.area) > tolerance:
            raise RoofAssemblyError(f"Invalid triangulated coverage in patch {marker}")
    extra = set(map(int, markers[markers >= 0])) - set(range(len(patches)))
    if extra:
        raise RoofAssemblyError(
            f"Backend returned unknown patch markers {sorted(extra)}"
        )


def vertex_link_audit(mesh):
    """Each vertex link must be one cycle (interior) or one path (boundary).

    Edge pairing alone cannot detect two shells touching only at a vertex.
    """
    faces = np.asarray(mesh.faces, dtype=np.int64)
    links = {}
    for face in faces:
        for i, vertex in enumerate(face):
            links.setdefault(int(vertex), []).append(
                (int(face[(i + 1) % 3]), int(face[(i + 2) % 3]))
            )
    bad = []
    for vertex, edges in links.items():
        graph = {}
        for a, b in edges:
            graph.setdefault(a, []).append(b)
            graph.setdefault(b, []).append(a)
        degrees = [len(v) for v in graph.values()]
        seen, todo = set(), [next(iter(graph))]
        while todo:
            v = todo.pop()
            if v not in seen:
                seen.add(v)
                todo.extend(graph[v])
        if (
            len(seen) != len(graph)
            or max(degrees) > 2
            or degrees.count(1) not in (0, 2)
        ):
            bad.append(vertex)
    return dict(nonmanifold_vertices=len(bad), invalid_vertices=bad[:50])


def validate_shell(mesh, *, num_patches, closure_audit, base_z=None):
    """Audit manifoldness, winding and degeneracy with domain or base openings allowed."""
    closure = closure_audit(mesh, num_buildings=num_patches) if base_z is None else None
    links = vertex_link_audit(mesh)
    xyz = np.asarray(mesh.vertices)[np.asarray(mesh.faces)]
    cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    areas = np.linalg.norm(cross, axis=1) / 2
    lengths2 = sum(
        np.sum((xyz[:, i] - xyz[:, (i + 1) % 3]) ** 2, axis=1) for i in range(3)
    )
    quality = np.divide(
        4 * np.sqrt(3) * areas, lengths2, out=np.zeros_like(areas), where=lengths2 > 0
    )
    markers = np.asarray(mesh.markers)
    roofs = (markers >= num_patches) & (markers < 2 * num_patches)
    # Attribute topological failures to incident patch markers for block fallback.
    directed = np.concatenate(
        [mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]]
    )
    keys, inverse, counts = np.unique(
        np.sort(directed, axis=1), axis=0, return_inverse=True, return_counts=True
    )
    inverse = np.asarray(inverse).reshape(-1)
    forward = np.bincount(
        inverse, weights=directed[:, 0] < directed[:, 1], minlength=len(keys)
    )
    bad_edge = (counts > 2) | ((counts == 2) & (forward != 1))
    if base_z is None:
        # Naked edges on the city domain are intentional.
        xy = np.asarray(mesh.vertices)[:, :2]
        sides = np.concatenate(
            [abs(xy - xy.min(axis=0)) <= 1e-8, abs(xy - xy.max(axis=0)) <= 1e-8], axis=1
        )
        outline = np.any(sides[keys[:, 0]] & sides[keys[:, 1]], axis=1)
    else:
        outline = np.all(
            abs(np.asarray(mesh.vertices)[keys, 2] - base_z) <= 1e-8, axis=1
        )
        closure = dict(
            open_base_edges=int(np.sum((counts == 1) & outline)),
            unexpected_open_edges=int(np.sum((counts == 1) & ~outline)),
            nonmanifold_edges=int(np.sum(counts > 2)),
            inconsistent_edges=int(np.sum((counts == 2) & (forward != 1))),
            closed=bool(not np.any(counts == 1)),
        )
        closure["boundary_valid"] = not any(
            closure[k]
            for k in (
                "unexpected_open_edges",
                "nonmanifold_edges",
                "inconsistent_edges",
            )
        )
    bad_edge |= (counts == 1) & ~outline
    bad_face = bad_edge[inverse].reshape(3, -1).any(axis=0)
    bad_face |= np.isin(mesh.faces, links["invalid_vertices"]).any(axis=1)
    bad_face |= ~np.isfinite(areas) | (areas <= 1e-12) | (roofs & (cross[:, 2] <= 0))
    invalid_markers = sorted(
        set(
            int(m % num_patches)
            for m in markers[bad_face]
            if num_patches and 0 <= m < 2 * num_patches
        )
    )
    report = dict(
        closure=closure,
        **links,
        degenerate_faces=int(np.sum(~np.isfinite(areas) | (areas <= 1e-12))),
        inverted_roof_faces=int(np.sum(roofs & (cross[:, 2] <= 0))),
        minimum_quality=float(quality.min()) if len(quality) else None,
        quality_p01=float(np.quantile(quality, 0.01)) if len(quality) else None,
        triangles=len(areas),
    )
    report["invalid_patch_markers"] = invalid_markers
    report["valid"] = bool(
        closure.get("boundary_valid", closure["closed"])
        and not links["nonmanifold_vertices"]
        and not report["degenerate_faces"]
        and not report["inverted_roof_faces"]
    )
    return report
