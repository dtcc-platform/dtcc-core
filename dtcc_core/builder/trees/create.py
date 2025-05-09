from dtcc_core import io, builder
from dtcc_core.model import Tree, Raster, PointCloud

import scipy.ndimage


def tree_raster_from_pointcloud(
    pc: PointCloud,
    terrain_raster: Raster = None,
    cell_size: float = 0.5,
    shortest_tree: float = 2.0,
    smallest_cluster: float = 100,
    fill_hole_size: float = 100,
    sigma: float = 1.0,
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
    cell_size : float, optional
        The size of each raster cell in the output raster. Default is 0.5.
    shortest_tree : float, optional
        The minimum tree height to include in the raster. Values below this will be set
        to 0. Default is 2.0.
    smallest_cluster : float, optional
        The minimum size of tree clusters to retain (in pixels). Smaller clusters will be removed.
        Default is 100.
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

    trees = pc.get_vegetation()

    tree_raster: Raster = trees.rasterize(
        cell_size=abs(terrain_raster.cell_size[0]),
        bounds=terrain_raster.bounds,
        fill_holes=False,
        window_size=1,
    )

    tree_raster = tree_raster.erode_small_lines(neighborhood_size=None, nodata=0)
    if smallest_cluster > 0:
        tree_raster = tree_raster.remove_small_masks(
            min_size=smallest_cluster, nodata=0
        )
    if fill_hole_size > 0:
        tree_raster = tree_raster.fill_small_holes(hole_size=fill_hole_size, nodata=0)

    # Apply Gaussian filter to smooth the raster
    if sigma > 0:
        tree_raster.data = scipy.ndimage.gaussian_filter(tree_raster.data, sigma=1.0)

    tree_raster.data -= terrain_raster.data
    tree_raster.data[tree_raster.data < shortest_tree] = 0

    return tree_raster
