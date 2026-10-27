import numpy as np
from shapely.geometry import Polygon

from dtcc_core.common.progress import report_progress

from ...model import Building, GeometryType, MultiSurface, PointCloud, Raster, Surface
from .. import _dtcc_builder
from ..cleaning.blocks import (
    invalidate_lod1,
    is_block,
    part_roof,
    record_lod1,
    set_block_roof,
    set_part_roof,
    source_footprint,
    update_walls,
)
from ..logging import error, info, warning
from ..model_conversion import create_builder_polygon
from .surface import extrude_surface
from .terrain import build_terrain_raster


def _source_vertices(part):
    """Exterior vertices of a part's source footprint, without the closing ones."""
    polygon = source_footprint(part)
    if polygon is None:
        return np.empty((0, 2))
    pieces = [polygon] if isinstance(polygon, Polygon) else list(polygon.geoms)
    return np.vstack([np.asarray(piece.exterior.coords)[:-1, :2] for piece in pieces])


def _roof_point_polygon(polygon):
    """The polygon that collects the roof points of a footprint, as for a
    building with this footprint (simplified as by ``get_footprint``)."""
    surface = Surface()
    surface.from_polygon(polygon, 0.0)
    footprint = Surface()
    footprint.from_polygon(surface.to_polygon(), 0.0)
    return footprint.to_polygon()


def _footprint_centroid(footprint):
    return _footprint_vertices(footprint).mean(axis=0)


def _footprint_vertices(footprint):
    if isinstance(footprint, MultiSurface):
        return np.vstack([surface.vertices for surface in footprint.surfaces])
    return footprint.vertices


def extrude_building(
    building: Building, default_ground_height=0, always_use_default=False
) -> MultiSurface:
    """
    Extrudes the LOD0 representation of a building from its height to the ground level.

    Parameters
    ----------
    `building` : model.Building
        The building to extrude.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.
    `always_use_default` : bool, optional
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


def set_building_heights_from_attribute(
    buildings: list[Building],
    terrain: Raster,
    height_attribute: str = "measured_height",
    default_ground_height: float = 0,
    always_use_default_ground: bool = False,
    min_building_height: float = 2.5,
    default_building_height: float = 10.0,
    ground_height_strategy: str = "centroid",
) -> list[Building]:
    """Calculates ground level and absolute height for each building footprint from a given height attribute.

    Each part of a cleaned block gets its height from its own attributes and
    its ground on the footprint of its source, also where it has lost its cell
    to a neighbour; the block's roof follows from its parts. Where the new
    roofs change which part wins an overlap, the walls of the block are
    recomputed (:func:`.cleaning.blocks.update_walls`).

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to set heights for.
    terrain : Raster
        Terrain raster used to determine ground elevation.
    height_attribute : str, default "measured_height"
        Source height attribute in metres. The modelling result is stored as
        estimated_height; the source measurement is preserved.
    default_ground_height : float, default 0
        Default ground height if not available from terrain.
    always_use_default_ground : bool, default False
        Whether to always use the default ground height.
    min_building_height : float, default 2.5
        Minimum height in meters for buildings.
    default_building_height : float, default 10.0
        Default building height if attribute is missing.
    ground_height_strategy : str, default "centroid"
        Strategy to determine ground height from terrain ("centroid", "vertex").
    """
    for building in buildings:
        if is_block(building):
            # A cleaned block: each part on the footprint of its source.
            for part in building.building_parts:
                roof = _height_from_attribute(
                    part,
                    _source_vertices(part),
                    terrain,
                    height_attribute=height_attribute,
                    default_ground_height=default_ground_height,
                    always_use_default_ground=always_use_default_ground,
                    min_building_height=min_building_height,
                    default_building_height=default_building_height,
                    ground_height_strategy=ground_height_strategy,
                )
                if roof is not None:
                    set_part_roof(part, roof)
            set_block_roof(building)
            update_walls(building)
            continue
        footprint = building.lod0
        if footprint is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        roof = _height_from_attribute(
            building,
            _footprint_vertices(footprint),
            terrain,
            height_attribute=height_attribute,
            default_ground_height=default_ground_height,
            always_use_default_ground=always_use_default_ground,
            min_building_height=min_building_height,
            default_building_height=default_building_height,
            ground_height_strategy=ground_height_strategy,
        )
        if roof is not None:
            footprint.set_z(roof)
            invalidate_lod1(building)
    return buildings


def _height_from_attribute(
    building,
    vertices,
    terrain,
    *,
    height_attribute,
    default_ground_height,
    always_use_default_ground,
    min_building_height,
    default_building_height,
    ground_height_strategy,
):
    """Set the ground and height of a building or part from an attribute,
    with the ground from the terrain at its footprint vertices. Returns the
    roof level, or None for an invalid footprint."""
    if len(vertices) < 3:
        warning(
            f"Building {building.id} has an invalid footprint with only {len(vertices)} vertices."
        )
        return None
    if always_use_default_ground:
        ground_height = default_ground_height
    else:
        ground_height_strategy = ground_height_strategy.lower()
        if ground_height_strategy == "centroid":
            centroid = vertices.mean(axis=0)
            if np.isnan(centroid[0]) or np.isnan(centroid[1]):
                warning(f"Building {building.id} has an invalid centroid.")
                ground_height = default_ground_height
            else:
                ground_height = terrain.get_value(centroid[0], centroid[1])
        elif ground_height_strategy == "vertex":
            z_values = []
            for vertex in vertices:
                z = terrain.get_value(vertex[0], vertex[1])
                if not np.isnan(z):
                    z_values.append(z)
            if len(z_values) == 0:
                warning(
                    f"Building {building.id} has no valid terrain values at its vertices."
                )
                ground_height = default_ground_height
            else:
                ground_height = min(z_values)
        else:
            warning(
                f"Unknown ground height strategy '{ground_height_strategy}'. Using default ground height."
            )
            ground_height = default_ground_height
    building.attributes["ground_height"] = ground_height

    height = building.attributes.get(height_attribute)
    if height is None:
        height = default_building_height
    if height < min_building_height:
        warning(
            f"Building {building.id} has a height of {height}, which is less than the minimum building height of {min_building_height}. Setting height to {min_building_height}."
        )
        height = min_building_height
    building.estimated_height = height
    return ground_height + height


def compute_building_heights(
    buildings: list[Building],
    terrain: Raster,
    min_building_height=2.5,
    roof_percentile=0.9,
    overwrite=False,
) -> list[Building]:
    """
    Compute building heights from roof points and terrain elevation.

    This function calculates building heights by determining the ground elevation
    from the terrain raster and the roof elevation from building roof points.
    Results (including minimum-height defaults) are stored as estimated_height;
    measured_height is preserved. Terrain and point coordinates must be in metres.
    Each part of a cleaned block gets its height from its own roof points and
    its ground on the footprint of its source, also where it has lost its cell
    to a neighbour; the block's roof follows from its parts. Where the new
    roofs change which part wins an overlap, the walls of the block are
    recomputed (:func:`.cleaning.blocks.update_walls`).

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to compute heights for.
    terrain : Raster
        Terrain raster used to determine ground elevation.
    min_building_height : float, default 2.5
        Minimum height in meters for buildings.
    roof_percentile : float, default 0.9
        Percentile of roof points to use for determining roof elevation.
    overwrite : bool, default False
        Whether to overwrite existing height values.

    Returns
    -------
    List[Building]
        List of buildings with computed heights.
    """
    total_buildings = len(buildings)
    info("Computing building heights...")
    report_progress(
        percent=0, message=f"Computing heights for {total_buildings} buildings..."
    )
    for i, building in enumerate(buildings):
        if is_block(building):
            # A cleaned block: each part on the footprint of its source, with
            # its own roof points.
            for part in building.building_parts:
                roof = _compute_height(
                    part,
                    _source_vertices(part),
                    part_roof(part),
                    terrain,
                    min_building_height,
                    roof_percentile,
                    overwrite,
                )
                if roof is not None:
                    set_part_roof(part, roof)
            set_block_roof(building)
            update_walls(building)
        else:
            footprint = building.lod0
            if footprint is None:
                warning(f"Building {building.id} has no LOD0 geometry.")
                continue
            roof = _compute_height(
                building,
                (
                    _source_vertices(building)
                    if "source_footprint" in building.attributes
                    else _footprint_vertices(footprint)
                ),
                footprint.zmax,
                terrain,
                min_building_height,
                roof_percentile,
                overwrite,
            )
            if roof is not None:
                footprint.set_z(roof)
                invalidate_lod1(building)
        if (i + 1) % 20 == 0 or i + 1 == total_buildings:
            report_progress(
                current=i + 1,
                total=total_buildings,
                message=f"Computing heights ({i + 1}/{total_buildings})...",
            )
    return buildings


def _compute_height(
    building, vertices, roof_z, terrain, min_building_height, roof_percentile, overwrite
):
    """Ground and height of one building or part from its roof points.

    The ground is the terrain at the mean of the footprint vertices. Returns
    the new roof level, or None if the roof is kept.
    """
    if len(vertices) < 3:
        warning(
            f"Building {building.id} has an invalid footprint with only {len(vertices)} vertices."
        )
        return None
    centroid = vertices.mean(axis=0)
    if np.isnan(centroid[0]) or np.isnan(centroid[1]):
        warning(f"Building {building.id} has an invalid centroid.")
        return None
    ground_height = terrain.get_value(centroid[0], centroid[1])
    building.attributes["ground_height"] = ground_height
    if not (overwrite or roof_z == 0):
        return None
    roof_points = building.point_cloud
    if roof_points is None or len(roof_points) == 0:
        warning(f"Building {building.id} has no roof points. using min height")
        building.estimated_height = min_building_height
        return ground_height + min_building_height
    z_values = roof_points.points[:, 2]
    roof_top = np.percentile(z_values, roof_percentile * 100)
    height = max(roof_top - ground_height, min_building_height)
    building.estimated_height = height
    return ground_height + height


def build_lod1_buildings(
    buildings: [Building],
    default_ground_height=0,
    always_use_default_ground=False,
    rebuild=True,
) -> list[Building]:
    """
    Build the LOD1 representation of the given buildings.

    Each building is extruded from its ground to its roof. A cleaned block
    with walls is extruded per part instead: its walls are first brought up
    to date with its region, sources and roofs
    (:func:`.cleaning.blocks.update_walls`, which removes the LOD1 geometry
    of a block whose walls change), then each part with a footprint gets its
    own LOD1 prism, from its own ground to its own roof, and the block itself
    and its parts without a footprint get none, as for building parts in
    CityGML. A block that keeps one footprint gets one prism and its parts
    none.
    Neighbouring prisms share faces up to the lower roof.
    ``building.flatten_geometry(GeometryType.LOD1)`` gives all prisms of a
    block; :func:`dtcc_core.builder.build_building_shell` gives one closed
    shell with step walls.

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
        Whether to rebuild the LOD1 representation if it already exists, by
        default True. Generated prisms invalidated by changed cells, roof or
        ground are rebuilt even when False. Existing unmarked or edited LOD1
        geometry is retained when False.
    Returns
    -------
    [model.Building]
        The buildings with the LOD1 representation built.
    """
    total_buildings = len(buildings)
    info(f"Building LOD1 representations of {total_buildings} buildings...")
    report_progress(
        percent=0, message=f"Building LOD1 for {total_buildings} buildings..."
    )

    for i, building in enumerate(buildings):
        parts = []
        if is_block(building):
            # The walls of a block follow its region, sources and roofs; a
            # change removes the prisms of the block and its parts.
            update_walls(building)
            parts = [
                part
                for part in building.building_parts
                if part.get_geometry(GeometryType.LOD0) is not None
            ]
        targets = parts if parts else [building]
        if is_block(building):
            # Only the parts with a footprint, or the block without walls,
            # have prisms.
            keep = {id(target) for target in targets}
            for obj in [building, *building.building_parts]:
                if id(obj) not in keep:
                    obj.remove_geometry(GeometryType.LOD1)
        for target in targets:
            invalidate_lod1(target)
            if target.lod1 is not None and not rebuild:
                continue
            if target.lod0 is None:
                warning(f"Building {target.id} has no LOD0 geometry.")
                continue
            geometry = extrude_building(
                target, default_ground_height, always_use_default_ground
            )
            if geometry is not None:
                target.add_geometry(geometry, GeometryType.LOD1)
                record_lod1(target)
            else:
                warning(f"Building {target.id} LOD1 geometry could not be built.")

        # Report progress every 10 buildings or on last building
        if (i + 1) % 10 == 0 or i + 1 == total_buildings:
            report_progress(
                current=i + 1,
                total=total_buildings,
                message=f"Building LOD1 ({i + 1}/{total_buildings})...",
            )

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
) -> list[Building]:
    """
    Extract roof points from a point cloud for a list of buildings.

    Each part of a cleaned block gets the roof points inside the footprint of
    its source.

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

    report_progress(percent=10, message="Preparing building footprints...")
    # Each building gets the roof points inside its footprint, each part of a
    # cleaned block those inside the footprint of its source.
    targets, polygons, owners = [], [], []
    for building in buildings:
        if is_block(building):
            for part in building.building_parts:
                polygon = source_footprint(part)
                if polygon is None:
                    continue
                for piece in (
                    [polygon] if isinstance(polygon, Polygon) else polygon.geoms
                ):
                    polygons.append(_roof_point_polygon(piece))
                    owners.append(len(targets))
                targets.append(part)
            continue
        footprint = building.get_footprint()
        if "source_footprint" in building.attributes:
            polygon = source_footprint(building)
            if polygon is not None:
                for piece in (
                    [polygon] if isinstance(polygon, Polygon) else polygon.geoms
                ):
                    polygons.append(_roof_point_polygon(piece))
                    owners.append(len(targets))
                targets.append(building)
        elif footprint is not None:
            polygons.append(footprint.to_polygon())
            owners.append(len(targets))
            targets.append(building)

    builder_polygon = [create_builder_polygon(p) for p in polygons]
    if len(pointcloud.points) == len(pointcloud.classification):
        ground_mask = np.logical_or(
            pointcloud.classification == 2, pointcloud.classification == 9
        )
        not_ground_mask = ~ground_mask
        points = pointcloud.points[not_ground_mask]
    else:
        points = pointcloud.points
    report_progress(percent=30, message="Extracting roof points (C++)...")
    roof_points = _dtcc_builder.extract_building_points(
        builder_polygon,
        points,
        statistical_outlier_remover,
        roof_outlier_neighbors,
        roof_outlier_margin,
    )

    report_progress(percent=80, message="Assigning roof points to buildings...")
    found = [[] for _ in targets]
    for owner, values in zip(owners, roof_points):
        found[owner].append(np.asarray(values).reshape(-1, 3))
    for target, values in zip(targets, found):
        pc = PointCloud(points=np.vstack(values) if values else np.empty((0, 3)))
        pc.calculate_bounds()
        target.add_geometry(pc, GeometryType.POINT_CLOUD)
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
) -> list[Building]:
    """
    Compute building heights from point cloud data.

    This function combines roof point extraction and height computation to determine
    building heights from LiDAR or similar point cloud data. It first extracts roof
    points, then computes heights using terrain elevation.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to compute heights for.
    pointcloud : PointCloud
        Point cloud containing building and terrain points.
    terrain_raster : Raster, optional
        Terrain raster for ground elevation. If None, creates one from point cloud.
    statistical_outlier_remover : bool, default True
        Whether to apply statistical outlier removal to roof points.
    roof_outlier_neighbors : int, default 5
        Number of neighbors for outlier detection.
    roof_outlier_margin : float, default 1.5
        Margin for statistical outlier removal.
    overwrite : bool, default False
        Whether to overwrite existing height values.
    keep_roof_points : bool, default False
        Whether to keep extracted roof points as building geometry.

    Returns
    -------
    List[Building]
        List of buildings with computed heights.
    """

    if terrain_raster is None:
        info("No terrain raster provided, building terrain raster from point cloud.")
        terrain_raster = build_terrain_raster(
            pointcloud, cell_size=2, ground_only=True, _report_progress=False
        )

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
            for part in building.building_parts:
                part.remove_geometry(GeometryType.POINT_CLOUD)
    return buildings
