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
from dtcc_core.common.progress import report_progress


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


def set_building_heights_from_attribute(
    buildings: List[Building],
    terrain: Raster,
    height_attribute: str = "height",
    default_ground_height: float = 0,
    always_use_default_ground: bool = False,
    min_building_height: float = 2.5,
    default_building_height: float = 10.0,
    ground_height_strategy: str = "centroid",
) -> List[Building]:
    """Calculates ground level and absolute height for each building footprint from a given height attribute.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to set heights for.
    terrain : Raster
        Terrain raster used to determine ground elevation.
    height_attribute : str, default "height"
        Attribute name to use for building height.
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
        footprint = building.lod0
        if footprint is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        if len(footprint.vertices) < 3:
            warning(
                f"Building {building.id} has an invalid footprint with only {len(footprint.vertices)} vertices."
            )
            continue
        if always_use_default_ground:
            ground_height = default_ground_height
        else:
            ground_height_strategy = ground_height_strategy.lower()
            if ground_height_strategy == "centroid":
                centroid = footprint.centroid
                if np.isnan(centroid[0]) or np.isnan(centroid[1]):
                    warning(f"Building {building.id} has an invalid centroid.")
                    ground_height = default_ground_height
                else:
                    ground_height = terrain.get_value(centroid[0], centroid[1])
            elif ground_height_strategy == "vertex":
                z_values = []
                for vertex in footprint.vertices:
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

        height = building.attributes.get(height_attribute, default_building_height)
        if height < min_building_height:
            warning(
                f"Building {building.id} has a height of {height}, which is less than the minimum building height of {min_building_height}. Setting height to {min_building_height}."
            )
            height = min_building_height
        footprint.set_z(ground_height + height)
        building.attributes["height"] = height
    return buildings


def compute_building_heights(
    buildings: List[Building],
    terrain: Raster,
    min_building_height=2.5,
    roof_percentile=0.9,
    overwrite=False,
) -> List[Building]:
    """
    Compute building heights from roof points and terrain elevation.

    This function calculates building heights by determining the ground elevation
    from the terrain raster and the roof elevation from building roof points.

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
    report_progress(percent=0, message=f"Computing heights for {total_buildings} buildings...")
    for i, building in enumerate(buildings):
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
        if (i + 1) % 20 == 0 or i + 1 == total_buildings:
            report_progress(
                current=i + 1,
                total=total_buildings,
                message=f"Computing heights ({i + 1}/{total_buildings})...",
            )
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
    total_buildings = len(buildings)
    info(f"Building LOD1 representations of {total_buildings} buildings...")
    report_progress(percent=0, message=f"Building LOD1 for {total_buildings} buildings...")

    for i, building in enumerate(buildings):
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

        # Report progress every 10 buildings or on last building
        if (i + 1) % 10 == 0 or i + 1 == total_buildings:
            report_progress(
                current=i + 1,
                total=total_buildings,
                message=f"Building LOD1 ({i + 1}/{total_buildings})..."
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

    report_progress(percent=10, message="Preparing building footprints...")
    footprint_polygons = [b.get_footprint() for b in buildings]

    # Track which buildings have footprints for correct index alignment
    valid_building_indices = [
        i for i, p in enumerate(footprint_polygons) if p is not None
    ]
    builder_polygon = [
        create_builder_polygon(footprint_polygons[i].to_polygon())
        for i in valid_building_indices
    ]

    has_classification = len(pointcloud.points) == len(pointcloud.classification)
    if has_classification:
        ground_mask = np.logical_or(
            pointcloud.classification == 2, pointcloud.classification == 9
        )
        not_ground_mask = ~ground_mask
        points = pointcloud.points[not_ground_mask]
        not_ground_classification = pointcloud.classification[not_ground_mask]
    else:
        points = pointcloud.points
        not_ground_classification = None

    report_progress(percent=30, message="Extracting roof points (C++)...")
    roof_points = _dtcc_builder.extract_building_points(
        builder_polygon,
        points,
        statistical_outlier_remover,
        roof_outlier_neighbors,
        roof_outlier_margin,
    )

    report_progress(percent=80, message="Assigning roof points to buildings...")
    # Reconstruct classification via nearest-neighbor lookup
    tree = None
    if not_ground_classification is not None:
        from scipy.spatial import cKDTree
        tree = cKDTree(points)

    for roof_idx, building_idx in enumerate(valid_building_indices):
        if roof_idx >= len(roof_points):
            break
        rp = roof_points[roof_idx]
        if rp is None or len(rp) == 0:
            continue
        cls = None
        if tree is not None and not_ground_classification is not None:
            _, original_indices = tree.query(rp, k=1)
            cls = not_ground_classification[original_indices]
        pc = PointCloud(points=rp, classification=cls if cls is not None else np.empty(0))
        pc.calculate_bounds()
        buildings[building_idx].add_geometry(pc, GeometryType.POINT_CLOUD)

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
        terrain_raster = build_terrain_raster(pointcloud, cell_size=2, ground_only=True, _report_progress=False)

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


def build_lod2_buildings(
    buildings,
    config=None,
    rebuild=True,
    fallback_to_flat=True,
    **kwargs,
):
    """Build LoD2 geometry for buildings with roof point clouds.

    Requires LoD1 to be already built. Detects roof planes from per-building
    point clouds, classifies roof type, and constructs LoD2 geometry with
    semantic surface labels (RoofSurface, WallSurface, GroundSurface).

    Parameters
    ----------
    buildings : list[Building]
        Buildings with LoD1 geometry and point clouds.
    config : RoofDetectionConfig, optional
        Full configuration object. If None, one is constructed from kwargs.
    rebuild : bool
        If False, skip buildings that already have LoD2 geometry.
    fallback_to_flat : bool
        If True, store flat fallback geometry when detection fails.
        If False, leave LoD2 slot empty on failure.
    **kwargs
        Passed to RoofDetectionConfig for convenience.
    """
    from dtcc_core.builder.geometry_builders.roof_config import RoofDetectionConfig
    from dtcc_core.builder.geometry_builders.roof_detection import (
        filter_roof_points,
        detect_roof_planes,
    )
    from dtcc_core.builder.geometry_builders.roof_classification import (
        merge_planes,
        classify_roof,
    )
    from dtcc_core.builder.geometry_builders.roof_geometry import (
        build_flat_geometry,
        build_gabled_geometry,
        build_hipped_geometry,
        build_fallback_geometry,
    )
    from dtcc_core.builder.geometry.shell_validation import validate_shell
    from dtcc_core.model.enums import RoofType

    if config is None:
        valid_keys = {f.name for f in RoofDetectionConfig.__dataclass_fields__.values()}
        config_kwargs = {k: v for k, v in kwargs.items() if k in valid_keys}
        config = RoofDetectionConfig(**config_kwargs)

    for building in buildings:
        if not rebuild and building.lod2 is not None:
            continue
        if building.lod1 is None:
            warning(f"Building {building.id} has no LoD1 geometry, skipping LoD2.")
            continue

        # Stage 1: Filter roof points
        pc = building.point_cloud
        if pc is None or len(pc.points) < config.min_roof_points:
            _set_lod2_fallback(
                building, "insufficient_points", 0.0,
                fallback_to_flat, build_fallback_geometry,
            )
            continue

        filtered_pc, normals = filter_roof_points(pc, config)
        if filtered_pc is None:
            _set_lod2_fallback(
                building, "insufficient_points", 0.0,
                fallback_to_flat, build_fallback_geometry,
            )
            continue

        # Stage 2: Detect planes and classify
        planes = detect_roof_planes(filtered_pc, normals, config)
        planes = merge_planes(planes, config)

        footprint = building.lod0
        footprint_area = 0.0
        if footprint is not None:
            try:
                poly = Polygon(footprint.vertices[:, :2])
                footprint_area = poly.area
            except Exception:
                footprint_area = 0.0

        fp_verts = footprint.vertices if footprint is not None else None
        result = classify_roof(
            planes, footprint_area, len(filtered_pc.points), config,
            footprint_vertices=fp_verts,
            all_points=filtered_pc.points,
        )

        building.attributes["roof_type"] = result.roof_type.name
        building.attributes["roof_confidence"] = result.confidence

        if result.roof_type == RoofType.UNKNOWN:
            _set_lod2_fallback(
                building, "classification_failed", result.confidence,
                fallback_to_flat, build_fallback_geometry,
            )
            continue

        # Stage 3: Geometry gate (for non-flat types)
        if result.roof_type != RoofType.FLAT and footprint is not None:
            try:
                poly = Polygon(footprint.vertices[:, :2])
                convex_hull = poly.convex_hull
                convexity = poly.area / max(convex_hull.area, 1e-10)
                mbr = poly.minimum_rotated_rectangle
                rectangularity = poly.area / max(mbr.area, 1e-10)
                has_holes = len(list(poly.interiors)) > 0

                if (
                    has_holes
                    or convexity < config.min_convexity_ratio
                    or rectangularity < config.min_rectangularity_ratio
                ):
                    _set_lod2_fallback(
                        building, "complex_footprint", result.confidence,
                        fallback_to_flat, build_fallback_geometry,
                    )
                    continue
            except Exception:
                _set_lod2_fallback(
                    building, "complex_footprint", result.confidence,
                    fallback_to_flat, build_fallback_geometry,
                )
                continue

        # Stage 4: Build geometry
        lod2 = None
        try:
            if result.roof_type == RoofType.FLAT:
                lod2 = build_flat_geometry(building.lod1)
            elif result.roof_type == RoofType.GABLED and result.ridge_line is not None:
                lod2 = build_gabled_geometry(
                    building.lod1, footprint, result.ridge_line,
                    result.eave_height, result.ridge_height,
                )
            elif result.roof_type == RoofType.HIPPED and result.ridge_line is not None:
                lod2 = build_hipped_geometry(
                    building.lod1, footprint, result.ridge_line,
                    result.eave_height, result.ridge_height,
                )
            else:
                lod2 = build_flat_geometry(building.lod1)
        except Exception as e:
            warning(f"Building {building.id}: geometry construction failed: {e}")
            _set_lod2_fallback(
                building, "geometry_construction_failed", result.confidence,
                fallback_to_flat, build_fallback_geometry,
            )
            continue

        # Stage 5: Validate (warn but don't fall back for non-flat types,
        # since gabled/hipped geometry may have edge gaps at wall-roof junctions)
        is_valid, issues = validate_shell(lod2, config.edge_snap_tolerance)
        if not is_valid:
            if result.roof_type == RoofType.FLAT:
                _set_lod2_fallback(
                    building, "validation_failed", result.confidence,
                    fallback_to_flat, build_fallback_geometry,
                )
                continue
            else:
                warning(
                    f"Building {building.id}: shell not watertight "
                    f"({len(issues)} issues), storing geometry anyway"
                )

        building.add_geometry(lod2, GeometryType.LOD2)

    return buildings


def _set_lod2_fallback(building, reason, confidence, fallback_to_flat, fallback_fn):
    """Apply LoD2 fallback: set attributes, optionally store flat geometry."""
    building.attributes["fallback_reason"] = reason
    if "roof_confidence" not in building.attributes:
        building.attributes["roof_confidence"] = confidence
    if "roof_type" not in building.attributes:
        building.attributes["roof_type"] = "UNKNOWN"

    if fallback_to_flat and building.lod1 is not None:
        lod2 = fallback_fn(building.lod1)
        building.add_geometry(lod2, GeometryType.LOD2)
