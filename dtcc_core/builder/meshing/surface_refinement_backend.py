"""Adapters to maintained constrained quality meshers for planar face domains."""
from dataclasses import asdict

import numpy as np
from shapely.geometry import Polygon

from .dtcc_mesher_backend import _load_dtcc_mesher
from .surface_faces import boundary_rings


def mesh_face(points, segments, polygon, options, *, spacing, reuse_acute_corners=False):
    """Refine constraints and interiors; return all requested boundary splits.

    The DTCC API has no spatial size field. Its maintained circumcentre and
    encroached-segment refinement supplies local grading. A nominal area cap
    controls face interiors, while max_edge_length controls input segments.
    No custom point insertion or triangulation is performed here.
    """
    if options.backend == 'dtcc_mesher':
        dm = _load_dtcc_mesher()
        if 'reuse_acute_corners' not in dm.MeshingOptions.__dataclass_fields__:
            raise RuntimeError('Shared-edge pilot requires DTCC Mesher corner-reuse extension; see docs/shared_edge_pilot.md')
        holes = np.array([Polygon(r).representative_point().coords[0] for r in polygon.interiors], dtype=float).reshape(-1,2)
        raw = dm.mesh(dm.Domain(points, segments, holes), options=dm.MeshingOptions(
            min_angle=options.min_angle, max_edge_length=spacing,
            max_area=.5*spacing*spacing, refine=True, acute_protection='shell',protect_angle=options.min_angle,
            max_refine_steps=options.max_refine_steps, max_protection_levels=6,
            reuse_acute_corners=reuse_acute_corners))
        return np.asarray(raw.points), np.asarray(raw.triangles, dtype=np.int64), dict(
            backend='dtcc_mesher', refinement_enabled=True,
            quality_exempt_triangles=np.flatnonzero(raw.quality_exempt).tolist(), **asdict(raw.summary))
    if options.backend == 'triangle':
        from .. import _dtcc_builder
        rings = boundary_rings([tuple(map(int,e)) for e in segments])
        rings.sort(key=lambda r: -Polygon(points[r]).area)
        raw = _dtcc_builder.mesh_planar_polygon_bounded(
            points[rings[0]].tolist(), [points[r].tolist() for r in rings[1:]],
            spacing, options.min_angle, options.max_face_vertices)
        mesh = raw.from_cpp()
        return mesh.vertices[:,:2], np.asarray(mesh.faces,dtype=np.int64), dict(
            backend='triangle', refinement_enabled=True, max_steiner_points=options.max_face_vertices,
            protected_corner_count=None, exempt_triangle_count=None)
    raise ValueError(f'Unknown face-refinement backend {options.backend!r}')
