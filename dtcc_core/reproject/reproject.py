import numpy as np
import pyproj
from pyproj import Transformer
from dtcc_core.model import PointCloud, Building, Object, Surface, Mesh, MultiSurface
from dtcc_core.common import warning


def reproject_array(points: np.ndarray, src_crs: str, target_crs: str) -> np.ndarray:
    """
    Reproject an Nx3 array of points from source CRS to target CRS.

    Args:
        points: Nx3 numpy array where columns are [X, Y, Z]
        src_crs: Source coordinate reference system (e.g., "EPSG:3006")
        target_crs: Target coordinate reference system (e.g., "EPSG:4326")

    Returns:
        Nx3 numpy array with reprojected X, Y coordinates and original Z
    """
    if src_crs == target_crs:
        return points

    transformer = Transformer.from_crs(src_crs, target_crs, always_xy=True)
    x_new, y_new = transformer.transform(points[:, 0], points[:, 1])
    reprojected = np.column_stack((x_new, y_new, points[:, 2]))

    return reprojected


def reproject_surface(surface: Surface, src_crs: str, target_crs: str) -> Surface:
    """
    Reproject a surface from source CRS to target CRS.
    Args:
        surface: Surface object to reproject
        src_crs: source coordinate reference system
        target_crs: target coordinate reference system

    Returns:
        reprojected surface
    """
    if src_crs == target_crs:
        return surface

    trans_verts = reproject_array(surface.vertices, src_crs, target_crs)
    trans_holes = [reproject_array(hole, src_crs, target_crs) for hole in surface.holes]
    return Surface(vertices=trans_verts, holes=trans_holes)


def reproject_mesh(mesh: Mesh, src_crs: str, target_crs: str) -> Mesh:
    trans_verts = reproject_array(mesh.vertices, src_crs, target_crs)
    new_mesh = mesh.copy()
    new_mesh.vertices = trans_verts
    return new_mesh


def reproject_multisurface(
    multisurface: MultiSurface, src_crs: str, target_crs: str
) -> MultiSurface:
    if src_crs == target_crs:
        return multisurface

    new_surfaces = []
    for surface in multisurface.surfaces:
        new_surface = reproject_surface(surface, src_crs, target_crs)
        new_surfaces.append(new_surface)
    return MultiSurface(surfaces=new_surfaces)


def reproject_object(obj: Object, src_crs: str, target_crs: str) -> Object:
    if src_crs == target_crs:
        return obj

    new_obj = obj.copy()

    for geom_type in new_obj.defined_geometries:
        geom = obj.geometry.get(geom_type, None)
        if geom is not None:
            if isinstance(geom, Surface):
                new_obj.geometry[geom_type] = reproject_surface(
                    geom, src_crs, target_crs
                )
            if isinstance(geom, MultiSurface):
                new_obj.geometry[geom_type] = reproject_multisurface(
                    geom, src_crs, target_crs
                )
            elif isinstance(geom, Mesh):
                new_obj.geometry[geom_type] = reproject_mesh(geom, src_crs, target_crs)
            elif isinstance(geom, PointCloud):
                new_obj.geometry[geom_type] = reproject_pointcloud(
                    geom, src_crs, target_crs
                )
            else:
                warning(
                    f"Geometry type {geom_type} not supported for reprojection. Leaving as is!"
                )
    return new_obj


def reproject_pointcloud(
    pointcloud: PointCloud,
    src_crs: str,
    target_crs: str,
    override_geometry_crs: bool = False,
) -> PointCloud:
    if src_crs == target_crs:
        return pointcloud

    if target_crs is None or len(target_crs) == 0:
        raise ValueError("Target CRS not provided")

    if not pointcloud.transform.is_identity:
        pts = pointcloud.transform(pointcloud.points)
    else:
        pts = pointcloud.points

    if src_crs is None or len(src_crs) == 0:
        src_crs = pointcloud.transform.srs
        if src_crs is None or len(src_crs) == 0:
            src_crs = "EPSG:3006"
            warning("Source CRS not provided, defaulting to EPSG:3006")

    else:
        if not override_geometry_crs:
            if pointcloud.transform.srs:
                src_crs = pointcloud.transform.srs
            else:
                warning("Pointcloud CRS not set, defaulting to EPSG:3006")
                src_crs = "EPSG:3006"

    pts = reproject_array(pts, src_crs, target_crs)
    pc = PointCloud(
        points=pts,
        classification=pointcloud.classification,
        intensity=pointcloud.intensity,
        return_number=pointcloud.return_number,
        num_returns=pointcloud.num_returns,
        transform=pointcloud.transform,
    )
    pc.calculate_bounds()
    pc.transform.srs = target_crs
    return pc
