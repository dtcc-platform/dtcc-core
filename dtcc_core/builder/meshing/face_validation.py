"""Conditioned plane, boundary and area contracts for a triangulated face."""
import numpy as np
import shapely

from ..geometry_policy import COINCIDENCE,XY_PRECISION


def validate_face_geometry(face,vertices,triangles):
    """Bound geometry error in length first, then propagate it to area.

    Registry reconciliation and joint collinear simplification have a 10 nm
    length budget. Its possible area change scales with boundary length; an
    independent fixed square-metre cutoff is not a consistent length contract.
    Overlap receives only a local floating-point area allowance.
    """
    xyz = vertices[triangles]
    if not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError(f'{face.kind} triangulation is empty or nonfinite')
    residual = float(np.max(abs((xyz-face.origin)@face.normal)))
    if residual>XY_PRECISION:
        raise ValueError(f'{face.kind} plane changed by {residual:g} m')
    ids,inverse = np.unique(triangles,return_inverse=True)
    xy = ((vertices[ids]-face.origin)@face.basis.T)[inverse].reshape(-1,3,2)
    polygons = shapely.polygons(xy)
    union = shapely.union_all(polygons)
    if union.is_empty:
        raise ValueError(f'{face.kind} triangulation has no covered area')
    scale = max(1.,float(np.max(abs(xy))),max(map(abs,face.polygon.bounds)))
    eps = np.finfo(float).eps
    rounding_length = 64*eps*scale
    boundary = float(shapely.hausdorff_distance(face.polygon.boundary,union.boundary,densify=.25))
    if not np.isfinite(boundary) or boundary>COINCIDENCE+rounding_length:
        raise ValueError(f'{face.kind} boundary changed by {boundary:g} m '
                         f'(coincidence budget {COINCIDENCE:g} m)')
    corners = sum(len(r.coords)-1 for r in [face.polygon.exterior,*face.polygon.interiors])
    rounding_area = 64*eps*scale**2*(3*len(triangles)+corners)
    error = float(face.polygon.symmetric_difference(union).area)
    overlap = max(0.,float(shapely.area(polygons).sum()-union.area))
    # The boundary's delta-neighbourhood bounds the symmetric difference.
    displacement = boundary+rounding_length
    tolerance = (face.polygon.length+union.length)*displacement+2*np.pi*displacement**2+rounding_area
    if error>tolerance or overlap>rounding_area:
        raise ValueError(f'{face.kind} coverage changed: {error:g} m², overlap {overlap:g} m² '
                         f'(area budget {tolerance:g}, overlap budget {rounding_area:g})')
    return dict(plane_error_m=residual,coverage_error_m2=error,overlap_m2=overlap,
                boundary_deviation_m=boundary,boundary_tolerance_m=COINCIDENCE+rounding_length,
                coverage_tolerance_m2=float(tolerance),overlap_tolerance_m2=float(rounding_area),
                original_plane_error_m=face.original_plane_error)
