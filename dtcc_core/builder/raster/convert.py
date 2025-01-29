from ...model import Raster, PointCloud

from ..register import register_model_method

import numpy as np


@register_model_method
def to_pointcloud(raster: Raster, point_classification=1) -> PointCloud:
    """
    Convert a raster to a point cloud. The points will be in the center of each cell.

    Parameters
    ----------
    raster : Raster
        The raster to convert to a point cloud.

    Returns
    -------
    PointCloud
        The point cloud representation of the raster.

    """
    height = raster.height
    width = raster.width
    raster_affine = raster.georef

    xx, yy = np.meshgrid(np.arange(width), np.arange(height))

    xx = xx.flatten()
    yy = yy.flatten()
    zz = raster.data
    if len(zz.shape) != 2:
        ValueError("Only 2D rasters are supported. Got shape: ", zz.shape)
    zz = zz.flatten()
    transformed_points = raster_affine * (xx, yy)
    xx, yy = transformed_points
    cell_x, cell_y = raster.cell_size
    # offset to center of cell
    xx += cell_x / 2
    yy += cell_y / 2
    points = np.array([xx, yy, zz]).T
    classification = np.ones(len(points), dtype=int) * point_classification
    return PointCloud(points=points, classification=classification)
