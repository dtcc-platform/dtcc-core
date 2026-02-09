from ...model import (
    PointCloud,
    Terrain,
    Raster,
    Mesh,
    Surface,
    GeometryType,
    Bounds,
)

from ..model_conversion import (
    raster_to_builder_gridfield,
    builder_mesh_to_mesh,
    create_builder_polygon,
)

import numpy as np
from pypoints2grid import points2grid
from affine import Affine
from .. import _dtcc_builder
from typing import List, Union
from dtcc_core.common.progress import report_progress


def build_terrain_surface_mesh(
    data: Union[PointCloud, Raster],
    subdomains: list[Surface] = None,
    holes: list[Surface] = None,
    subdomain_resolution: Union[float, List[float]] = None,
    max_mesh_size=10,
    min_mesh_angle=20.7,
    smoothing=3,
    ground_points_only=True,
) -> Mesh:
    """
    Build a triangular surface mesh from terrain data.
    
    This function creates a triangular surface mesh representation of terrain from either
    point cloud or raster data, with optional subdomains for varying resolution.
    
    Parameters
    ----------
    data : Union[PointCloud, Raster]
        Input terrain data to mesh.
    subdomains : list[Surface], optional
        List of surface subdomains for varying mesh resolution.
    holes : list[Surface], optional
        Surfaces that should be treated as holes (excised from the mesh).
    subdomain_resolution : Union[float, List[float]], optional
        Resolution for each subdomain. If float, applies to all subdomains.
    max_mesh_size : float, default 10
        Maximum triangle size in meters.
    min_mesh_angle : float, default 20.7
        Minimum angle in degrees for mesh triangles (must be <= 33).
    smoothing : int, default 3
        Number of smoothing iterations to apply.
    ground_points_only : bool, default True
        Whether to use only ground-classified points from point cloud.
        
    Returns
    -------
    Mesh
        Triangular mesh representation of the terrain.
        
    Raises
    ------
    ValueError
        If min_mesh_angle > 33 degrees or data type is invalid.
    """

    """
    Build a surface mesh of the terrain from a point cloud.

    Parameters
    ----------
    data : PointCloud or Raster
        The point cloud or raster to build the terrain from.
    subdomains : list[Surface]
        The list of surface to use as subdomains.
    subdomain_resolution : Union[float, List[float]]
        The resolution of the subdomains. If a single value is given, it is
        used for all subdomains. If a list is given, it must have the same length
        as the subdomains list.
    max_mesh_size : float
        The maximum size of the triangles in the mesh.
    min_mesh_angle : float
        The minimum angle of the triangles in the mesh. Must be less than or equal
        to 33 degrees.
    smoothing : int
        The number of smoothing iterations to apply to the mesh.
    ground_points_only : bool
        If True, only ground points are used to build the terrain.

    Returns
    -------
    Mesh
        The surface mesh of the terrain.
    """
    if min_mesh_angle > 33:
        raise ValueError(
            "min_mesh_angle must be less than or equal to 33 degrees. "
            "This is a limitation of the meshing algorithm."
        )

    report_progress(percent=0, message="Preparing terrain data...")

    if isinstance(data, PointCloud):
        report_progress(percent=10, message="Building terrain raster from point cloud...")
        dem = build_terrain_raster(
            data, cell_size=max_mesh_size / 2, ground_only=ground_points_only,
            _report_progress=False,
        )
    elif isinstance(data, Raster):
        dem = data
    else:
        raise ValueError("data must be a PointCloud or a Raster.")

    report_progress(percent=30, message="Converting raster to grid field...")
    _builder_gridfield = raster_to_builder_gridfield(dem)

    if subdomains is None:
        subdomains = []
        subdomain_resolution = None
    else:
        subdomains = [create_builder_polygon(sub.to_polygon()) for sub in subdomains]
    if holes is None:
        hole_polygons: list = []
    else:
        hole_polygons = [create_builder_polygon(sub.to_polygon()) for sub in holes]
    if subdomain_resolution is None:
        subdomain_resolution = []
    elif isinstance(subdomain_resolution, (float, int)):
        subdomain_resolution = [subdomain_resolution] * len(subdomains)
    if (
        len(subdomains) > 0
        and len(subdomain_resolution) > 0
        and len(subdomains) != len(subdomain_resolution)
    ):
        raise ValueError(
            "subdomains and subdomain_resolution must have the same length."
        )

    subdomain_resolution = np.array(subdomain_resolution, dtype=np.float64)

    report_progress(percent=40, message="Building terrain surface mesh (this may take a while)...")

    terrain_mesh = _dtcc_builder.build_terrain_surface_mesh(
        subdomains,
        hole_polygons,
        subdomain_resolution,
        _builder_gridfield,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        False,
    )

    report_progress(percent=90, message="Converting mesh format...")
    terrain_mesh = builder_mesh_to_mesh(terrain_mesh)

    report_progress(percent=100, message="Terrain surface mesh complete")
    return terrain_mesh


def build_terrain_raster(
    pc: PointCloud, cell_size, bounds=None, window_size=3, radius=0, ground_only=True,
    _report_progress=True,
) -> Raster:
    """
    Rasterize a point cloud into a `Raster` object.

    Args:
        cell_size (float): The size of the raster cells in meters.
        bounds (Bounds): The bounds of the area to rasterize (default None, uses the bounds of the point cloud).
        window_size (int): The size of the window for the interpolation (default 3).
        radius (float): The radius of the search for the interpolation (default 0).
        ground_only (bool): Whether to only use ground points for the rasterization (default True).

    Returns:
        Raster: A `Raster` object representing the rasterized point cloud.
    """
    if _report_progress:
        report_progress(percent=0, message="Filtering ground points...")

    if (
        ground_only
        and (len(pc.classification) == len(pc.points))
        and 2 in pc.used_classifications()
    ):
        ground_point_idx = np.where(np.isin(pc.classification, [2, 9]))[0]
        ground_points = pc.points[ground_point_idx]
    else:
        ground_points = pc.points
    if bounds is None:
        if pc.bounds is None or pc.bounds.area == 0:
            pc.calculate_bounds()
        bounds = pc.bounds

    if _report_progress:
        report_progress(percent=30, message="Rasterizing points to grid...")

    dem = points2grid(
        ground_points, cell_size, bounds.tuple, window_size=window_size, radius=radius
    )

    if _report_progress:
        report_progress(percent=70, message="Creating raster object...")

    dem_raster = Raster()
    dem_raster.data = dem
    dem_raster.nodata = 0
    dem_raster.georef = Affine.translation(bounds.xmin, bounds.ymax) * Affine.scale(
        cell_size, -cell_size
    )

    if _report_progress:
        report_progress(percent=90, message="Filling holes in raster...")
    dem_raster = dem_raster.fill_holes()

    if _report_progress:
        report_progress(percent=100, message="Raster complete")
    return dem_raster


def flat_terrain(height, bounds: Bounds) -> Terrain:
    """
    Create a flat terrain.

    Args:
        height (float): The height of the terrain.
        bounds (Bounds): The bounds of the terrain.

    Returns:
        Terrain: A `Terrain` object representing the flat terrain.
    """
    terrain = Terrain()
    raster = Raster()
    raster.data = np.ones((1, 1)) * height
    raster.nodata = -9999
    raster.georef = Affine.translation(bounds.xmin, bounds.ymax) * Affine.scale(
        bounds.width, -bounds.height
    )
    terrain.add_geometry(raster, GeometryType.RASTER)
    return terrain
