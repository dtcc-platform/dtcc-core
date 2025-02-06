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


def build_terrain_mesh(
    data: Union[PointCloud, Raster],
    subdomains: list[Surface] = None,
    subdomain_resolution: Union[float, List[float]] = None,
    max_mesh_size=10,
    min_mesh_angle=25,
    smoothing=3,
    ground_points_only=True,
) -> Mesh:
    if isinstance(data, PointCloud):
        dem = build_terrain_raster(
            data, cell_size=max_mesh_size / 2, ground_only=ground_points_only
        )
    elif isinstance(data, Raster):
        dem = data
    else:
        raise ValueError("data must be a PointCloud or a Raster.")
    _builder_gridfield = raster_to_builder_gridfield(dem)
    if subdomains is None:
        subdomains = []
        subdomain_resolution = None
    else:
        subdomains = [create_builder_polygon(sub.to_polygon()) for sub in subdomains]
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

    terrain_mesh = _dtcc_builder.build_terrain_mesh(
        subdomains,
        subdomain_resolution,
        _builder_gridfield,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        False,
    )
    terrain_mesh = builder_mesh_to_mesh(terrain_mesh)
    return terrain_mesh


def build_terrain_raster(
    pc: PointCloud, cell_size, bounds=None, window_size=3, radius=0, ground_only=True
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
    dem = points2grid(
        ground_points, cell_size, bounds.tuple, window_size=window_size, radius=radius
    )
    dem_raster = Raster()
    dem_raster.data = dem
    dem_raster.nodata = 0
    dem_raster.georef = Affine.translation(bounds.xmin, bounds.ymax) * Affine.scale(
        cell_size, -cell_size
    )
    dem_raster = dem_raster.fill_holes()
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
