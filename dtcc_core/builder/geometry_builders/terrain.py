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
from numbers import Real
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
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
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
    *,
    hole_fill="nearest",
) -> Raster:
    """Rasterize elevations by inverse-distance weighting at cell centers.

    ``ground_only=True`` requires LAS class 2 points; water (class 9) and
    unclassified points are not ground. Pass False explicitly to use all points.
    Cell size and radius use the input coordinate units. Zero radius selects the
    backend default; window_size=0 disables its window interpolation. Remaining
    gaps use NaN, with optional nearest-neighbor filling (the historical default).
    No reprojection or vertical-unit conversion is performed. Nonidentity point
    transforms must be applied explicitly before calling this builder.
    """
    if not isinstance(pc, PointCloud):
        raise TypeError("pc must be a PointCloud")
    points = np.asarray(pc.points)
    if (points.ndim != 2 or points.shape[1] != 3 or len(points) == 0
            or points.dtype.kind not in "fiu" or not np.isfinite(points).all()):
        raise ValueError("PointCloud must contain finite real points with shape (N, 3)")
    if not np.array_equal(pc.transform.affine, np.eye(4)):
        raise ValueError("Apply the PointCloud transform before rasterization")
    if (isinstance(cell_size, (bool, np.bool_)) or not isinstance(cell_size, Real)
            or not np.isfinite(cell_size) or cell_size <= 0):
        raise ValueError("cell_size must be finite and positive")
    if (isinstance(radius, (bool, np.bool_)) or not isinstance(radius, Real)
            or not np.isfinite(radius) or radius < 0):
        raise ValueError("radius must be finite and nonnegative")
    if (isinstance(window_size, (bool, np.bool_))
            or not isinstance(window_size, (int, np.integer))
            or window_size < 0 or (window_size != 0 and window_size % 2 == 0)):
        raise ValueError("window_size must be zero or a positive odd integer")
    if not isinstance(ground_only, bool):
        raise ValueError("ground_only must be a boolean")
    if hole_fill not in {"nearest", "none"}:
        raise ValueError("hole_fill must be 'nearest' or 'none'")
    if _report_progress:
        report_progress(percent=0, message="Filtering ground points...")
    if ground_only:
        classifications = np.asarray(pc.classification)
        if classifications.shape != (len(points),) or classifications.dtype.kind not in "iu":
            raise ValueError("ground_only=True requires one integer LAS classification per point")
        points = points[classifications == 2]
        if len(points) == 0:
            raise ValueError("ground_only=True requires LAS class 2 ground points; use False explicitly for all points")
    if bounds is None:
        # Use the complete source footprint, even when ground selection leaves gaps.
        bounds = pc.calculate_bounds()
    if (not isinstance(bounds, Bounds) or not np.isfinite(bounds.tuple).all()
            or bounds.width <= 0 or bounds.height <= 0):
        raise ValueError("bounds must be finite with positive width and height")
    # The compiled grid backend indexes cells with signed integers. Reject a
    # numerically impossible request before narrowing dimensions or allocating.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        grid_shape = np.ceil(np.array([bounds.width, bounds.height]) / cell_size)
        cell_count = np.prod(grid_shape)
    if (not np.isfinite(cell_count) or np.any(grid_shape < 1)
            or cell_count > np.iinfo(np.int32).max):
        raise ValueError("Requested raster exceeds the backend signed cell-index limit")
    if _report_progress:
        report_progress(percent=30, message="Rasterizing points to grid...")
    dem = points2grid(points, cell_size, bounds.tuple, window_size=window_size, radius=radius)
    # The backend uses zero both for missing support and valid elevations. A
    # constant positive signal isolates support using identical interpolation;
    # shifting the elevations instead would lose precision near zero.
    if np.any(dem == 0):
        support_points = points.copy()
        support_points[:, 2] = 1
        support = points2grid(support_points, cell_size, bounds.tuple,
                              window_size=window_size, radius=radius)
        dem[support == 0] = np.nan
    if not np.isfinite(dem).any():
        raise ValueError("No selected points support the requested raster grid")
    if np.isinf(dem).any():
        raise ValueError("Raster interpolation produced nonfinite elevations")
    # points2grid rounds dimensions up from the lower-left bound. The upper
    # raster edge may therefore extend beyond the requested y maximum.
    dem_raster = Raster(
        data=dem, nodata=np.nan, crs=pc.transform.srs,
        georef=Affine.translation(bounds.xmin, bounds.ymin + dem.shape[0] * cell_size)
               * Affine.scale(cell_size, -cell_size),
    )
    if hole_fill == "nearest":
        if _report_progress:
            report_progress(percent=90, message="Filling holes in raster...")
        dem_raster = fill_raster_holes(dem_raster)
    if _report_progress:
        report_progress(percent=100, message="Raster complete")
    return dem_raster


def build_terrain_dem(
    pc: PointCloud, cell_size, *, unit: str, vertical_reference: str,
    bounds=None, ground_only=True, window_size=0, radius=0,
    hole_fill="none", source: str | None = None,
) -> Terrain:
    """Build a native Terrain with an explicitly interpreted elevation raster.

    The source must declare its coordinate reference system. ``unit`` describes
    its Z values; ``vertical_reference`` identifies their datum/reference without
    resolving or converting it. Values are point estimates at raster cell centers,
    not cell averages or direct measurements. No CityGML RasterRelief equivalence
    is claimed. Access the array with ``terrain.get_geometry(id="dem").data``.
    """
    if not isinstance(pc, PointCloud) or not isinstance(pc.transform.srs, str) or not pc.transform.srs.strip():
        raise ValueError("A DEM requires the source PointCloud CRS in transform.srs")
    from pyproj import CRS
    source_crs = CRS.from_user_input(pc.transform.srs)
    if not source_crs.is_projected:
        raise ValueError("A DEM requires a projected source CRS; reproject geographic coordinates first")
    vertical_axes = [axis for axis in source_crs.axis_info if axis.direction in {"up", "down"}]
    if vertical_axes:
        scale = {"m": 1.0, "cm": .01, "mm": .001, "km": 1000.0}.get(unit)
        if scale is not None and any(not np.isclose(axis.unit_conversion_factor, scale) for axis in vertical_axes):
            raise ValueError("DEM unit conflicts with the source CRS vertical axis unit")
        if any(axis.direction != "up" for axis in vertical_axes):
            raise ValueError("A DEM requires upward-positive elevations; transform depth coordinates first")
    raster = build_terrain_raster(
        pc, cell_size, bounds=bounds, ground_only=ground_only,
        window_size=window_size, radius=radius, hole_fill=hole_fill,
    )
    terrain = Terrain()
    from ...model._standard_schema import SEMANTIC_NAMESPACE
    terrain.semantic_type = SEMANTIC_NAMESPACE + "Terrain"
    terrain.transform.srs = raster.crs
    terrain.add_geometry(raster, id="dem", role="elevation")
    interpretation = {
        "geometry_id": "dem", "unit": unit, "vertical_reference": vertical_reference,
        "sampling": "cell_center", "interpolation": "inverse_distance_weighted",
        "hole_fill": hole_fill, "source_selection": "ground" if ground_only else "all",
        "window_size": int(window_size), "radius": float(radius),
    }
    if source is not None:
        interpretation["source"] = source
    terrain.attributes["elevation_rasters"] = [interpretation]
    # Use the standard authority for semantic metadata, without serializing or
    # rebuilding the numerical buffers merely to validate the new object.
    from ...model._standard_schema import SCHEMA_ID, DEFAULT_VERSION, validate_admitted
    validate_admitted(terrain, SCHEMA_ID, DEFAULT_VERSION)
    return terrain


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
