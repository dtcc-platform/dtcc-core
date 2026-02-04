import numpy as np

from dtcc_core import builder
from dtcc_core.builder.trees.utils import (
    detect_local_maxima,
    suppress_peaks_height_aware,
    extract_crown_polygons_from_labels,
    segment_tree_crowns,
)
from dtcc_core.model import Raster, PointCloud, Building, Tree
from scipy import ndimage as ndi

from shapely.geometry import Polygon

from dataclasses import dataclass


@dataclass
class TreeType:
    min_height: float
    smoothing_sigma: float
    peak_footprint: int
    min_radius: float


urban_trees = TreeType(2.5, 0.9, 3, 0.6)
mixed_trees = TreeType(3.0, 0.7, 3, 0.8)
dense_trees = TreeType(3.5, 0.6, 3, 1.0)
arrid_trees = TreeType(2.0, 0.2, 3, 0.4)

tree_types = {
    "urban": urban_trees,
    "mixed": mixed_trees,
    "dense": dense_trees,
    "arid": arrid_trees,
}


def tree_raster_from_pointcloud(
    pc: PointCloud,
    terrain_raster: Raster = None,
    buildings: [Building] = None,
    tree_type: str = "urban",
    cell_size: float = 0.5,
    smallest_cluster: float = 4,
    fill_hole_size: float = 2,
) -> Raster:
    """
    Generate a tree height raster from a point cloud.

    This function processes a point cloud to create a raster representing tree heights.
    It optionally uses a terrain raster to subtract ground elevation and applies various
    filters to clean and smooth the resulting raster.

    Parameters
    ----------
    pc : PointCloud
        The input point cloud containing vegetation points.
    terrain_raster : Raster, optional
        A raster representing ground elevation. If not provided, it will be generated
        from the point cloud.
    buildings : list[Building], optional
        If a list of buildings is provided, points within building footprints will be excluded
        from the tree rasterization.
    cell_size : float, optional
        The size of each raster cell in the output raster. Default is 0.5.
    shortest_tree : float, optional
        The minimum tree height to include in the raster. Values below this will be set
        to 0. Default is 2.0.
    smallest_cluster : float, optional
        The minimum size of tree clusters to retain (in pixels). Smaller clusters will be removed.
    fill_hole_size : float, optional
        The maximum size of holes to fill in the raster (in pixels). Default is 100.
    sigma : float, optional
        The standard deviation for the Gaussian filter applied to smooth the raster.
        Default is 1.0.

    Returns
    -------
    Raster
        A raster representing tree heights, with ground elevation subtracted and
        cleaned using the specified parameters.
    """

    if terrain_raster is None:
        terrain_raster = builder.build_terrain_raster(
            pc, cell_size=cell_size, ground_only=True
        )

    tree_settings = tree_types[tree_type]

    trees = pc.get_vegetation()
    sigma = tree_settings.smoothing_sigma
    shortest_tree = tree_settings.min_height

    if buildings is not None:
        building_footprints = [b.get_footprint() for b in buildings]
        building_footprints = [b for b in building_footprints if b is not None]
        trees = builder.pointcloud.filter.remove_points_in_polygons(
            trees, building_footprints
        )

    tree_raster: Raster = trees.rasterize(
        cell_size=abs(terrain_raster.cell_size[0]),
        bounds=terrain_raster.bounds,
        fill_holes=False,
        window_size=3,
    )

    tree_raster = tree_raster.erode_small_lines(neighborhood_size=5, nodata=0)
    if smallest_cluster > 0:
        tree_raster = tree_raster.remove_small_masks(
            min_size=smallest_cluster, nodata=0
        )
    if fill_hole_size > 0:
        tree_raster = tree_raster.fill_small_holes(hole_size=fill_hole_size, nodata=0)

    tree_raster_mask = np.zeros_like(tree_raster.data, dtype=bool)
    tree_raster_mask[tree_raster.data > 0] = True

    tree_raster.data -= terrain_raster.data
    tree_raster.data[tree_raster.data < shortest_tree] = 0

    tree_raster.data *= tree_raster_mask

    # Apply Gaussian filter to smooth the raster
    if sigma > 0:
        tree_raster.data = ndi.gaussian_filter(tree_raster.data, sigma=sigma)

    return tree_raster


def find_tree_tops(
    chm: Raster,
    tree_type: str = "urban",
    max_relative_drop: float = 0.3,
    max_absolute_drop: float = 8.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Detect dominant canopy tree tops from a CHM.
    """
    pixel_size = abs(chm.cell_size[0])
    tree_settings = tree_types.get(tree_type, None)
    if tree_settings is None:
        raise ValueError(
            f"Tree type {tree_type} not found. use one of: {list(tree_types.keys())}"
        )
    min_height = tree_settings.min_height
    min_radius = tree_settings.min_radius

    chm_data = chm.data
    mask = np.zeros_like(chm_data, dtype=bool)
    mask[chm_data >= min_height] = True

    coords, heights = detect_local_maxima(
        chm_data, mask, footprint_size=tree_settings.peak_footprint
    )
    min_radius = 0.6
    coords, heights, radius = suppress_peaks_height_aware(
        coords,
        heights,
        chm_data,
        pixel_size,
        min_radius,
        max_relative_drop,
        max_absolute_drop,
    )

    return coords, heights, radius


def tree_crown_polygons(
    tree_raster: Raster,
    tree_coords: np.ndarray,
    tree_type: str = "urban",
    min_area: float = 1.0,
) -> (list[Polygon], np.ndarray):
    """
    Create polygons for tree crowns from a raster and tree coordinates.
    """
    tree_params = tree_types.get(tree_type, None)
    if tree_params is None:
        raise ValueError(
            f"Tree type {tree_type} not found. use one of: {list(tree_types.keys())}"
        )
    min_height = tree_params.min_height
    labels = segment_tree_crowns(tree_raster.data, tree_coords, min_height)
    pixel_size = abs(tree_raster.cell_size[0])

    crown_polygons = extract_crown_polygons_from_labels(
        labels, tree_raster.georef, pixel_size * 2
    )
    valid_trees = np.array([v is not None for v in crown_polygons], dtype=bool)
    valid_polygons = [p for p in crown_polygons if p is not None]
    return valid_polygons, valid_trees


def trees_from_pointcloud(
    pc: PointCloud,
    terrain_raster: Raster = None,
    buildings: list[Building] = None,
    tree_type: str = "urban",
    cell_size: float = 0.5,
    smallest_cluster: float = 4,
    fill_hole_size: float = 2,
) -> list[Tree]:
    if terrain_raster is None:
        terrain_raster = builder.build_terrain_raster(
            pc, cell_size=cell_size, ground_only=True
        )
    tree_raster = tree_raster_from_pointcloud(
        pc,
        terrain_raster,
        buildings,
        tree_type,
        cell_size,
        smallest_cluster,
        fill_hole_size,
    )
    coords, height, radius = find_tree_tops(tree_raster, tree_type)
    trees = []
    for c, h, r in zip(coords, height, radius):
        tree_base = terrain_raster.data[c[0], c[1]]
        geo_ref_coord = terrain_raster.pixel_to_georef(c[0], c[1])
        trees.append(
            Tree(
                position=np.array([geo_ref_coord[0], geo_ref_coord[1], tree_base]),
                height=h,
                crown_radius=r,
            )
        )
    return trees
