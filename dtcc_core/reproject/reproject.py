"""Horizontal coordinate reprojection with attached field values unchanged.

Coordinates must already be in the declared source frame. Fields retain their
values, associations, units and component conventions; vector fields warn on a
CRS change because components are not rotated or rescaled. Semantic regions,
normals, local affine transforms and DatasetContext still require explicit
transformation rules. Z values are retained, not reprojected.
"""

import warnings

import numpy as np
from pyproj import CRS, Transformer

from dtcc_core.model import (PointCloud, Object, Surface, Mesh, VolumeMesh, MultiSurface,
                             City, CityObject, Building, BuildingPart, Terrain, Landuse)


def _crs_pair(value, src_crs, target_crs, *, override=False):
    if not isinstance(target_crs, (str, CRS)) or not str(target_crs).strip():
        raise ValueError("Target CRS not provided")
    declared = getattr(getattr(value, 'transform', None), 'srs', '')
    source = src_crs or declared
    if not isinstance(source, (str, CRS)) or not str(source).strip():
        raise ValueError("Source CRS must be supplied or declared on the geometry")
    source, target = CRS.from_user_input(source), CRS.from_user_input(target_crs)
    if declared and not override and source != CRS.from_user_input(declared):
        raise ValueError("Source CRS conflicts with the geometry CRS")
    for crs in (source, target):
        if len(crs.axis_info) != 2 or not (crs.is_projected or crs.is_geographic):
            raise NotImplementedError("Reprojection supports horizontal two-axis CRSs only; Z is unchanged")
    return source, target


def _check_reprojection_state(value):
    if not np.array_equal(value.transform.affine, np.eye(4)):
        raise NotImplementedError("Apply the local affine transform explicitly before reprojection")
    if getattr(value, 'regions', None):
        raise NotImplementedError("Semantic regions require an explicit reprojection rule")
    if getattr(value, 'normal', np.array([])).size or getattr(value, 'normals', np.array([])).size:
        raise NotImplementedError("Stored normals require an explicit reprojection rule")
    if value.dataset_context is not None:
        raise NotImplementedError("DatasetContext requires an explicit reprojection/provenance update")
    if any(field.dim > 1 for field in getattr(value, 'fields', ())):
        warnings.warn(
            "Vector field components are preserved unchanged. Their orientation "
            "and units are not transformed to the target coordinate system.",
            UserWarning,
            stacklevel=3,
        )


def reproject_array(points: np.ndarray, src_crs: str, target_crs: str) -> np.ndarray:
    """Reproject finite XYZ coordinates horizontally, retaining Z values."""
    source, target = _crs_pair(None, src_crs, target_crs)
    if (not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3
            or points.dtype.kind not in 'fiu' or not np.isfinite(points).all()):
        raise ValueError("Reprojection requires finite real XYZ coordinates with shape (N, 3)")
    if source == target:
        return points
    if points.dtype.kind in 'iu' and any(int(float(z)) != int(z) for z in points[:, 2]):
        raise ValueError("Unchanged integer Z values must be exactly representable in floating-point output")
    transformer = Transformer.from_crs(source, target, always_xy=True)
    x, y = transformer.transform(points[:, 0], points[:, 1], errcheck=True)
    result = np.column_stack((x, y, points[:, 2]))
    if not np.isfinite(result).all():
        raise ValueError("Reprojection produced nonfinite coordinates")
    return result


def reproject_surface(surface: Surface, src_crs: str, target_crs: str) -> Surface:
    """Reproject a surface and its holes on a copy, keeping field values unchanged.

    Field units and associations are preserved. Vector components are not
    transformed; vector fields emit a UserWarning when the CRS changes.
    """
    source, target = _crs_pair(surface, src_crs, target_crs)
    if source == target:
        return surface
    _check_reprojection_state(surface)
    result = surface.copy()
    result.vertices = reproject_array(surface.vertices, source, target)
    result.holes = [reproject_array(hole, source, target) for hole in surface.holes]
    result.transform.srs = target.to_string()
    result._bounds = None
    return result


def reproject_mesh(mesh: Mesh | VolumeMesh, src_crs: str, target_crs: str) -> Mesh | VolumeMesh:
    """Reproject a triangle or volume mesh on a copy, keeping field values unchanged.

    Connectivity, markers, field units and associations are preserved. Vector
    components are not rotated or rescaled; vector fields emit a UserWarning
    when the CRS changes. Z is unchanged. This does not convert field components
    into the target coordinate basis or resample the field within elements.
    """
    source, target = _crs_pair(mesh, src_crs, target_crs)
    if source == target:
        return mesh
    _check_reprojection_state(mesh)
    result = mesh.copy()
    result.vertices = reproject_array(mesh.vertices, source, target)
    result.transform.srs = target.to_string()
    result._bounds = None
    return result


def reproject_multisurface(multisurface: MultiSurface, src_crs: str, target_crs: str) -> MultiSurface:
    """Reproject constituent surfaces, keeping parent and surface fields unchanged.

    Field units and associations are preserved. Vector components are not
    transformed; vector fields emit a UserWarning when the CRS changes.
    """
    source, target = _crs_pair(multisurface, src_crs, target_crs)
    if source == target:
        return multisurface
    _check_reprojection_state(multisurface)
    result = multisurface.copy()
    result.surfaces = [reproject_surface(surface, source, target) for surface in multisurface.surfaces]
    result.transform.srs = target.to_string()
    result._bounds = None
    return result


def reproject_object(obj: Object, src_crs: str, target_crs: str) -> Object:
    """Reproject one object's representations, keeping their field values unchanged.

    Field units and associations are preserved. Vector components are not
    transformed; vector fields emit a UserWarning when the CRS changes.
    Nested object frames remain unsupported.
    """
    source, target = _crs_pair(obj, src_crs, target_crs)
    if source == target:
        return obj
    _check_reprojection_state(obj)
    if type(obj) not in (Object, City, CityObject, Building, BuildingPart, Terrain, Landuse):
        raise NotImplementedError(f"{type(obj).__name__} object reprojection requires an explicit rule for intrinsic state")
    if any(obj.children.values()):
        raise NotImplementedError("Nested object reprojection requires explicit frame handling")
    handlers = {Surface: reproject_surface, MultiSurface: reproject_multisurface,
                Mesh: reproject_mesh, VolumeMesh: reproject_mesh,
                PointCloud: reproject_pointcloud}
    result = obj.copy()
    for record in result.geometry.values():
        handler = handlers.get(type(record.geometry))
        if handler is None:
            raise NotImplementedError(f"{type(record.geometry).__name__} reprojection is unsupported")
        record.geometry = handler(record.geometry, source, target)
    result.transform.srs = target.to_string()
    result._bounds = None
    return result


def reproject_pointcloud(pointcloud: PointCloud, src_crs: str | None, target_crs: str,
                         override_geometry_crs: bool = False) -> PointCloud:
    """Reproject a point cloud without modifying or sharing mutable source state.

    Field values, units and associations are preserved. Vector components are
    not transformed; vector fields emit a UserWarning when the CRS changes.
    Supply the source CRS or declare it in transform.srs. Conflicts fail unless
    override_geometry_crs=True explicitly selects the supplied source. Local
    affine transforms and metadata needing transformation rules are unsupported.
    """
    if type(override_geometry_crs) is not bool:
        raise TypeError("override_geometry_crs must be True or False")
    source, target = _crs_pair(pointcloud, src_crs, target_crs, override=override_geometry_crs)
    if source == target:
        return pointcloud
    _check_reprojection_state(pointcloud)
    result = pointcloud.copy()
    result.points = reproject_array(pointcloud.points, source, target)
    result.transform.srs = target.to_string()
    result._bounds = None
    return result
