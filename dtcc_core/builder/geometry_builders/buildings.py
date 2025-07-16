from ...model import Building, GeometryType
from ...model import Surface, MultiSurface, PointCloud, Raster
from ..model_conversion import create_builder_polygon
from .terrain import build_terrain_raster
from .. import _dtcc_builder
from ..logging import debug, info, warning, error
from shapely.geometry import Polygon
import numpy as np

from typing import List, Tuple

from .surface import extrude_surface


def extrude_building(
    building: Building, default_ground_height=0, always_use_default=False
) -> MultiSurface:
    """
    Extrudes the LOD0 representation of a building from its height to the ground leve.

    Parameters
    ----------
    `building` : model.Building
        The building to extrude.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.
    `always_use_default_ground` : bool, optional
        Whether to always use the default ground height or use ground_height attribute, by default False.

    Returns
    -------
    `MultiSurface`
        The extruded building.
    """
    if always_use_default:
        ground_height = default_ground_height
    else:
        ground_height = building.attributes.get("ground_height", default_ground_height)

    geometry = building.lod0
    if geometry is None:
        error(f"Building {building.id} has no LOD0 geometry.")
        return None
    if isinstance(geometry, Surface):
        geometry = MultiSurface(surfaces=[geometry])
    if not isinstance(geometry, MultiSurface):
        error(f"Building {building.id} LOD0 geometry is not a (Multi)Surface.")
        return None

    extrusion = MultiSurface()

    for surface in geometry.surfaces:
        extrusion = extrusion.merge(extrude_surface(surface, ground_height))
    return extrusion


def compute_building_heights(
    buildings: List[Building],
    terrain: Raster,
    min_building_height=2.5,
    roof_percentile=0.9,
    overwrite=False,
) -> List[Building]:
    """
    Compute the heights of a list of buildings by finding the highest point in the point cloud
    for each building and subtracting the ground height from the terrain.

    Parameters
    ----------
    `buildings` : List[Building]
        The list of buildings to compute the heights for.
    `terrain` : Raster
        The terrain raster to use for computing ground heights.
    `min_building_height` : float, optional
        The minimum height to use for a building if the computed height is too low, by default 2.5.
    `roof_percentile` : float, optional
        The percentile to use for computing the roof height from the point cloud, by default 0.9.
    `overwrite` : bool, optional
        Whether to overwrite the existing height attribute of the building or not, by default False.

    Returns
    -------
    List[Building]
        The list of buildings with their heights computed.
    """
    info("Computing building heights...")
    for building in buildings:
        footprint = building.lod0
        if footprint is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        if len(footprint.vertices) < 3:
            warning(
                f"Building {building.id} has an invalid footprint with only {len(footprint.vertices)} vertices."
            )
            continue
        centroid = footprint.centroid
        if np.isnan(centroid[0]) or np.isnan(centroid[1]):
            warning(f"Building {building.id} has an invalid centroid.")
            continue
        ground_height = terrain.get_value(centroid[0], centroid[1])
        building.attributes["ground_height"] = ground_height
        if overwrite or footprint.zmax == 0:
            roof_points = building.point_cloud
            if roof_points is None or len(roof_points) == 0:
                warning(f"Building {building.id} has no roof points. using min height")
                footprint.set_z(ground_height + min_building_height)
                building.attributes["height"] = min_building_height
            else:
                z_values = roof_points.points[:, 2]
                roof_top = np.percentile(z_values, roof_percentile * 100)
                height = roof_top - ground_height
                if height < min_building_height:
                    height = min_building_height
                footprint.set_z(ground_height + height)
                building.attributes["height"] = height
    return buildings


def build_lod1_buildings(
    buildings: [Building],
    default_ground_height=0,
    always_use_default_ground=False,
    rebuild=True,
) -> List[Building]:
    """
    Build the LOD1 representation of the given buildings.

    Parameters
    ----------
    `buildings` : [model.Building]
        The buildings to build the LOD1 representation of.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.
    `always_use_default_ground` : bool, optional
        Whether to always use the default ground height or use groun_height attribute, by default False.
    `rebuild` : bool, optional
        Whether to rebuild the LOD1 representation if it already exists, by default True.
    Returns
    -------
    [model.Building]
        The buildings with the LOD1 representation built.
    """
    info(f"Building LOD1 representations of {len(buildings )} buildings...")
    for building in buildings:
        if building.lod1 is not None and not rebuild:
            continue
        if building.lod0 is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        geometry = extrude_building(
            building, default_ground_height, always_use_default_ground
        )
        if geometry is not None:
            building.add_geometry(geometry, GeometryType.LOD1)
        else:
            warning(f"Building {building.id} LOD1 geometry could not be built.")
    return buildings


def extract_roof_points(
    buildings: [Building],
    pointcloud: PointCloud,
    statistical_outlier_remover=True,
    roof_outlier_neighbors=5,
    roof_outlier_margin=1.5,
    ransac_outlier_remover=False,
    ransac_outlier_margin=3.0,
    ransac_iterations=250,
) -> List[Building]:
    """
    Extract roof points from a point cloud for a list of buildings.

    Parameters
    ----------
    buildings : list[Building]
        The list of buildings to extract roof points for.
    pointcloud : PointCloud
        The point cloud to extract roof points from.
    statistical_outlier_remover : bool, optional
        Whether to use a statistical outlier remover on the roof points. Default is True.
    roof_outlier_neighbors : int, optional
        The number of neighbors to consider for statistical outlier removal. Default is 5.
    roof_outlier_margin : float, optional
        The margin for statistical outlier removal. Default is 1.5.
    ransac_outlier_remover : bool, optional
        Whether to use a RANSAC outlier remover on the roof points. Default is False.
    ransac_outlier_margin : float, optional
        The margin for RANSAC outlier removal. Default is 3.0.
    ransac_iterations : int, optional
        The number of iterations for RANSAC outlier removal. Default is 250.

    Returns
    -------
    list[Building]
        The list of buildings with roof points extracted.
    """
    
    footprint_polygons = [b.get_footprint() for b in buildings]

    builder_polygon = [
        create_builder_polygon(p.to_polygon())
        for p in footprint_polygons
        if p is not None
    ]
    if len(pointcloud.points) == len(pointcloud.classification):
        ground_mask = np.logical_or(
            pointcloud.classification == 2, pointcloud.classification == 9
        )
        not_ground_mask = ~ground_mask
        points = pointcloud.points[not_ground_mask]
    else:
        points = pointcloud.points
    roof_points = _dtcc_builder.extract_building_points(
        builder_polygon,
        points,
        statistical_outlier_remover,
        roof_outlier_neighbors,
        roof_outlier_margin,
    )

    idx = 0
    # some buildings may not have a footprint, and thus not have roof points
    for fp in footprint_polygons:
        if fp is not None:
            pc = PointCloud(points=roof_points[idx])
            pc.calculate_bounds()
            buildings[idx].add_geometry(pc, GeometryType.POINT_CLOUD)
            idx += 1
    return buildings


def building_heights_from_pointcloud(
    buildings: [Building],
    pointcloud: PointCloud,
    terrain_raster: Raster = None,
    statistical_outlier_remover=True,
    roof_outlier_neighbors=5,
    roof_outlier_margin=1.5,
    overwrite=False,
    keep_roof_points=False,
) -> List[Building]:

    """
    Calculate building heights from a point cloud and optionally a terrain raster.

    This function extracts roof points from a point cloud, computes building heights,
    and optionally removes roof points from the building geometries.

    Parameters
    ----------
    buildings : list[Building]
        The buildings for which to calculate heights.
    pointcloud : PointCloud
        The point cloud used for extracting roof points.
    terrain_raster : Raster, optional
        A raster representing ground elevation. If not provided, it will be generated from the point cloud.
    statistical_outlier_remover : bool, optional
        Whether to use a statistical outlier remover on the roof points. Default is True.
    roof_outlier_neighbors : int, optional
        The number of neighbors to consider for statistical outlier removal. Default is 5.
    roof_outlier_margin : float, optional
        The margin for statistical outlier removal. Default is 1.5.
    overwrite : bool, optional
        Whether to overwrite existing building heights. Default is False.
    keep_roof_points : bool, optional
        Whether to retain the roof points in the building geometries. Default is False.

    Returns
    -------
    list[Building]
        The buildings with updated heights.
    """

    if terrain_raster is None:
        info("No terrain raster provided, building terrain raster from point cloud.")
        terrain_raster = build_terrain_raster(pointcloud, cell_size=2, ground_only=True)

    buildings = extract_roof_points(
        buildings,
        pointcloud,
        statistical_outlier_remover,
        roof_outlier_neighbors,
        roof_outlier_margin,
    )
    buildings = compute_building_heights(buildings, terrain_raster, overwrite=overwrite)
    if not keep_roof_points:
        for building in buildings:
            building.remove_geometry(GeometryType.POINT_CLOUD)
    return buildings
