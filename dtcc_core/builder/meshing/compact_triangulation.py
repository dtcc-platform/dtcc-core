"""Planar Mesher adapters and compact face coverage validation."""

import numpy as np
import shapely

from .wall_faces import _check_replacement, triangle_quality


def orient_triangles(points, triangles):
    """Return planar triangle indices with counterclockwise winding."""
    t = np.array(triangles, dtype=np.int64, copy=True)
    xyz = points[t]
    d, e = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    negative = d[:, 0] * e[:, 1] - d[:, 1] * e[:, 0] < 0
    t[negative] = t[negative][:, [0, 2, 1]]
    return t


def improve(
    points,
    triangles,
    segments,
    *,
    max_added_points=0,
    splittable_segments=(),
    allow_interior=True,
):
    """Borrow planar arrays; Mesher owns all optimization and insertion logic."""
    from .dtcc_mesher_backend import _load_dtcc_mesher

    dm = _load_dtcc_mesher()
    if not hasattr(dm, "improve"):
        raise RuntimeError(
            "Compact optimization needs dtcc-mesher improve; see docs/compact_mesher_interface.md"
        )
    mesh = dm.Mesh(
        np.asarray(points),
        orient_triangles(points, triangles),
        np.asarray(segments, dtype=np.uint32),
        summary=None,
    )
    return dm.improve(
        mesh,
        max_added_points=max_added_points,
        splittable_segments=splittable_segments,
        allow_interior=allow_interior,
    )


def triangulate(points, segments, polygon, *, backend="dtcc_mesher"):
    """Triangulate fixed planar constraints, optimize edge flips and validate coverage."""
    points = np.asarray(points, dtype=float)
    segments = np.asarray(segments, dtype=np.int64)
    if backend == "dtcc_mesher":
        from .dtcc_mesher_backend import _load_dtcc_mesher

        dm = _load_dtcc_mesher()
        holes = np.array(
            [
                shapely.Polygon(r).representative_point().coords[0]
                for r in polygon.interiors
            ],
            dtype=float,
        ).reshape(-1, 2)
        raw = dm.mesh(
            dm.Domain(points, segments, holes),
            options=dm.MeshingOptions(refine=False, acute_protection="none"),
        )
        output, triangles = (
            np.asarray(raw.points),
            np.asarray(raw.triangles, dtype=np.int64),
        )
        if {tuple(sorted(map(int, edge))) for edge in raw.segments} != {
            tuple(sorted(map(int, edge))) for edge in segments
        }:
            raise ValueError("Boundary-only triangulator changed a prescribed segment")
    elif backend == "triangle":
        from .. import _dtcc_builder
        from .surface_faces import boundary_rings

        rings = boundary_rings([tuple(map(int, e)) for e in segments])
        rings.sort(key=lambda r: -shapely.Polygon(points[r]).area)
        raw = _dtcc_builder.mesh_planar_polygon_bounded(
            points[rings[0]].tolist(),
            [points[r].tolist() for r in rings[1:]],
            -1.0,
            0.0,
            0,
        ).from_cpp()
        output = np.asarray(raw.vertices)[:, :2]
        distance = np.linalg.norm(output[:, None] - points[None, :], axis=2)
        if len(output) != len(points) or distance.min(axis=1).max() > 1e-8:
            raise ValueError(
                "Triangle added or changed vertices in zero-refinement mode"
            )
        mapping = distance.argmin(axis=1)
        if len(np.unique(mapping)) != len(points):
            raise ValueError("Triangle merged required vertices")
        triangles = mapping[np.asarray(raw.faces, dtype=np.int64)]
        output = points
    elif backend == "geos":
        raw = shapely.constrained_delaunay_triangles(polygon)
        coords = np.asarray([np.asarray(t.exterior.coords)[:3] for t in raw.geoms])
        distances = np.linalg.norm(coords[:, :, None] - points[None, None], axis=3)
        if distances.min(axis=2).max() > 1e-8:
            raise ValueError("GEOS introduced a vertex")
        triangles, output = distances.argmin(axis=2), points
    else:
        raise ValueError("Unknown compact triangulation backend")
    if output.shape != points.shape or not np.array_equal(output, points):
        raise ValueError("Boundary-only triangulator changed or added vertices")
    _check_replacement(polygon, points, triangles, segments)
    optimized = improve(points, triangles, segments)
    triangles, flips = (
        np.asarray(optimized.mesh.triangles, dtype=np.int64),
        optimized.flips,
    )
    _check_replacement(polygon, points, triangles, segments)
    boundary_count = len(np.unique(segments))
    expected = 2 * len(points) - boundary_count + 2 * len(polygon.interiors) - 2
    if len(triangles) != expected:
        raise ValueError(
            f"Boundary-only Euler count mismatch: {len(triangles)} != {expected}"
        )
    return triangles, flips


def point_trial(points, triangles, polygon, segments):
    """Adapt one native interior trial; all geometric validation stays explicit."""
    result = improve(points, triangles, segments, max_added_points=1)
    if result.inserted_point is None:
        return None
    output, faces = (
        result.mesh.points,
        np.asarray(result.mesh.triangles, dtype=np.int64),
    )
    _check_replacement(polygon, output, faces, segments)
    cost = len(faces) - len(triangles)
    return dict(
        points=output,
        triangles=faces,
        benefit_per_triangle=(result.before["penalty"] - result.after["penalty"])
        / cost,
        penalty_before=result.before["penalty"],
        penalty_after=result.after["penalty"],
        q_before=result.before["minimum_q"],
        q_after=result.after["minimum_q"],
        triangle_cost=cost,
        flips=result.flips,
    )
