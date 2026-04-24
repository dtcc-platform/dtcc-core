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
    mesh_to_builder_mesh,
)

import numpy as np
from pypoints2grid import points2grid
from affine import Affine
from .. import _dtcc_builder
from typing import List, Optional, Union
from shapely.geometry import Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
from dtcc_core.common.progress import report_progress
from ..logging import info
from ..raster.interpolation import fill_holes as fill_raster_holes
from ..meshing.backends import resolve_2d_mesher
from ..meshing.flat_mesh_backends import build_city_flat_mesh_from_coverage


def _iter_polygon_components(geometry) -> list[Polygon]:
    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    if hasattr(geometry, "geoms"):
        polygons: list[Polygon] = []
        for child in geometry.geoms:
            polygons.extend(_iter_polygon_components(child))
        return polygons
    return []


def _terrain_mesh_regions(
    *,
    bounds: tuple[float, float, float, float],
    subdomains: list[Surface],
    holes: list[Surface],
    subdomain_resolution: list[float],
) -> tuple[list[Polygon], list[int], dict[int, float]]:
    subdomain_polygons: list[Polygon] = []
    for surface in subdomains:
        if surface is None:
            continue
        polygon = surface.to_polygon()
        if polygon is None or polygon.is_empty:
            continue
        subdomain_polygons.append(orient(polygon, sign=1.0))

    hole_polygons: list[Polygon] = []
    for surface in holes:
        if surface is None:
            continue
        polygon = surface.to_polygon()
        if polygon is None or polygon.is_empty:
            continue
        hole_polygons.append(orient(polygon, sign=1.0))

    ground_domain = box(*bounds)
    excluded_regions = [*subdomain_polygons, *hole_polygons]
    if excluded_regions:
        ground_domain = ground_domain.difference(unary_union(excluded_regions))

    ground_polygons = [
        orient(polygon, sign=1.0)
        for polygon in _iter_polygon_components(ground_domain)
        if polygon is not None and not polygon.is_empty
    ]

    region_polygons = [*ground_polygons, *subdomain_polygons]
    region_markers = [-2] * len(ground_polygons) + list(range(len(subdomain_polygons)))
    region_triangle_sizes = {
        index: float(resolution)
        for index, resolution in enumerate(subdomain_resolution[: len(subdomain_polygons)])
        if resolution is not None and float(resolution) > 0.0
    }
    return region_polygons, region_markers, region_triangle_sizes


def adaptive_terrain_mesh(
    data: Union[PointCloud, Raster], max_error: float, raster_size=1, smoothing=0
) -> Mesh:
    """builds an adaptive terrain mesh from a point cloud or raster data. The mesh generates the minimum number of
    triangles necessary to represent the terrain within a specified error margin, using the Zemlya algorithm."""

    if isinstance(data, PointCloud):
        dem = build_terrain_raster(data, cell_size=raster_size)
    elif isinstance(data, Raster):
        dem = data
    else:
        raise ValueError("data must be a PointCloud or a Raster.")
    _builder_gridfield = raster_to_builder_gridfield(dem)

    if max_error < 0:
        raise ValueError("max_error must be a positive number.")

    _builder_mesh = _dtcc_builder.build_terrain_mesh_zemlya(
        _builder_gridfield, max_error, smoothing
    )

    return builder_mesh_to_mesh(_builder_mesh)


def build_terrain_surface_mesh(
    data: Union[PointCloud, Raster],
    subdomains: list[Surface] = None,
    holes: list[Surface] = None,
    subdomain_resolution: Union[float, List[float]] = None,
    max_mesh_size=10,
    min_mesh_angle=20.7,
    smoothing=3,
    ground_points_only=True,
    report_mesh_quality=True,
    mesher: str | None = None,
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
    mesher : {"auto", "dtcc_mesher", "triangle", "spade"}, optional
        Select the 2D meshing backend for the ground triangulation. When
        omitted, the global default backend is used.
    Returns
    -------
    Mesh
        Triangular mesh representation of the terrain.

    Raises
    ------
    ValueError
        If min_mesh_angle > 33 degrees or data type is invalid.
    """
    if min_mesh_angle > 33:
        raise ValueError(
            "min_mesh_angle must be less than or equal to 33 degrees. "
            "This is a limitation of the meshing algorithm."
        )

    report_progress(percent=0, message="Preparing terrain data...")

    if isinstance(data, PointCloud):
        report_progress(
            percent=10, message="Building terrain raster from point cloud..."
        )
        dem = build_terrain_raster(
            data,
            cell_size=max_mesh_size / 2,
            ground_only=ground_points_only,
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
    if holes is None:
        holes = []
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
    active_mesher = resolve_2d_mesher(mesher)

    report_progress(
        percent=40, message="Building terrain surface mesh (this may take a while)..."
    )
    if active_mesher == "dtcc_mesher":
        if np.any(subdomain_resolution > 0.0):
            raise NotImplementedError(
                "Terrain subdomain-specific resolution is not yet supported with "
                "the dtcc_mesher default terrain mesher."
            )
        bounds = dem.bounds.tuple
        region_polygons, region_markers, region_triangle_sizes = _terrain_mesh_regions(
            bounds=bounds,
            subdomains=subdomains,
            holes=holes,
            subdomain_resolution=subdomain_resolution.tolist(),
        )
        ground_mesh = build_city_flat_mesh_from_coverage(
            region_polygons=region_polygons,
            region_markers=region_markers,
            bounds=bounds,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            backend=active_mesher,
            sort_triangles=False,
            region_triangle_sizes=region_triangle_sizes or None,
        )
        report_progress(percent=90, message="Converting mesh format...")
        terrain_mesh = _dtcc_builder.build_terrain_surface_mesh_from_ground_mesh(
            mesh_to_builder_mesh(ground_mesh),
            _builder_gridfield,
            smoothing,
        ).from_cpp()
    else:
        builder_subdomains = [create_builder_polygon(sub.to_polygon()) for sub in subdomains]
        builder_holes = [create_builder_polygon(sub.to_polygon()) for sub in holes]
        terrain_mesh = _dtcc_builder.build_terrain_surface_mesh(
            builder_subdomains,
            builder_holes,
            subdomain_resolution,
            _builder_gridfield,
            max_mesh_size,
            min_mesh_angle,
            smoothing,
            False,
        )

        report_progress(percent=90, message="Converting mesh format...")
        terrain_mesh = builder_mesh_to_mesh(terrain_mesh)

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import (
            triangle_mesh_quality,
            report_quality,
        )

        q = triangle_mesh_quality(terrain_mesh.vertices, terrain_mesh.faces)
        report_quality(q, log_fn=info)

    report_progress(percent=100, message="Terrain surface mesh complete")
    return terrain_mesh


def build_terrain_raster(
    pc: PointCloud,
    cell_size,
    bounds=None,
    window_size=3,
    radius=0,
    ground_only=True,
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
    dem_raster = fill_raster_holes(dem_raster)

    if _report_progress:
        report_progress(percent=100, message="Raster complete")
    return dem_raster


def flatten_terrain_raster(raster: Raster, height: Optional[float] = None) -> Raster:
    """
    Create a flat raster while preserving the source raster grid and georeferencing.

    Parameters
    ----------
    raster : Raster
        Source raster to flatten.
    height : float, optional
        Target ground level. If omitted, the minimum valid value in the source
        raster is used.

    Returns
    -------
    Raster
        Copy of the input raster with all valid cells set to one height.

    Raises
    ------
    ValueError
        If the raster does not contain any valid cells.
    """
    flat_raster = raster.copy()

    valid_mask = np.isfinite(flat_raster.data)
    if not np.isnan(flat_raster.nodata):
        valid_mask &= flat_raster.data != flat_raster.nodata

    if not np.any(valid_mask):
        raise ValueError("Cannot flatten a raster without any valid values.")

    if height is None:
        height = float(flat_raster.data[valid_mask].min())

    flat_raster.data[valid_mask] = height
    return flat_raster


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
