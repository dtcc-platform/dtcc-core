import numpy as np
from typing import List

from ...model import PointCloud, Bounds
from ..logging import info, warning, error
from ..register import register_model_method

from .. import _dtcc_builder


@register_model_method
def remove_global_outliers(pc: PointCloud, margin: float):
    """
    Remove outliers from a `PointCloud` whose Z-values are more than `margin`
    standard deviations from the mean.

    Args:
        margin (float): The margin in standard deviations to consider a point an outlier.

    Returns:
        PointCloud: A new `PointCloud` object with the outliers removed.
    """

    z_pts = pc.points[:, 2]
    z_mean = np.mean(z_pts)
    z_std = np.std(z_pts)
    outliers = np.where(np.abs(z_pts - z_mean) > margin * z_std)[0]
    new_pc = pc.copy()
    new_pc.remove_points(outliers)
    return new_pc


@register_model_method
def statistical_outlier_filter(pc: PointCloud, neighbours, outlier_margin):
    """
    Remove statistical outliers from a `PointCloud` object.

    Args:
        pc (PointCloud): The `PointCloud` object to remove outliers from.
        neighbours (int): The number of neighbours to consider for the outlier detection.
        outlier_margin (float): The margin in standard deviations to consider a point an outlier.

    Returns:
        PointCloud: A new `PointCloud` object with the outliers removed.
    """

    outliers = _dtcc_builder.statistical_outlier_finder(
        pc.points, neighbours, outlier_margin
    )
    return pc.remove_points(outliers)


@register_model_method
def classification_filter(pc: PointCloud, classes: List[int], keep: bool = False):
    """
    Filter a `PointCloud` object based on its classification.

    Args:
        classes (List[int]): The classification values to keep or remove.
        keep (bool): Whether to keep the points with the specified classification values (default False, remove them).

    Returns:
        PointCloud: A new `PointCloud` object with the specified points removed.
    """
    if len(pc.points) != len(pc.classification):
        warning("Pointcloud not classified, returning original pointcloud.")
        return pc
    mask = np.isin(pc.classification, classes)
    if keep:
        mask = np.logical_not(mask)
    pc.remove_points(mask)
    return pc


@register_model_method
def z_range_filter(pc: PointCloud, min=None, max=None):
    """
    Filter a `PointCloud` object based on its Z-values.

    Args:
        min (float): The minimum Z-value to keep.
        max (float): The maximum Z-value to keep.

    Returns:
        PointCloud: A new `PointCloud` object with the specified points removed.
    """
    mask = np.ones(len(pc.points), dtype=bool)
    filtered = False
    if min is not None:
        mask = np.logical_and(mask, pc.points[:, 2] >= min)
        filtered = True
    if max is not None:
        mask = np.logical_and(mask, pc.points[:, 2] <= max)
        filtered = True
    if filtered:
        pc.remove_points(np.logical_not(mask))
    return pc


@register_model_method
def remove_vegetation(pc: PointCloud) -> PointCloud:
    """
    Return a pioint cloud with vegetation removed.

    Args:
        pc (PointCloud): The `PointCloud` object to remove vegetation from.

    Returns:
        PointCloud: A new `PointCloud` object with the vegetation removed.
    """
    new_pc = pc.copy()
    veg_indices = _find_vegetation(pc)
    new_pc.remove_points(veg_indices)
    return new_pc


def _find_vegetation(pc: PointCloud, filter_on_return_number=True):
    """Find the indices of points that belong to vegetation in a point cloud.

    Args:
        pc: A `PointCloud` object representing the point cloud to filter.
        filter_on_return_number: A boolean indicating whether to filter on return number (default True).

    Returns:
        A 1D NumPy array of indices of points that belong to vegetation.
    """

    has_classification = len(pc.classification) == len(pc.points)
    has_return_number = len(pc.return_number) == len(pc.points)
    if not has_classification and not has_return_number:
        warning(
            "Classification and return number are not set for all points. Ignoring vegetation filter."
        )
        return np.array([])
    if not has_classification:
        warning("Classification is not set for all points. Ignoring")

    if filter_on_return_number and not has_return_number:
        filter_on_return_number = False
        warning("Return number is not set for all points. Ignoring")

    classes_with_vegetation = set([3, 4, 5])
    used_classes = pc.used_classifications()
    veg_classes = classes_with_vegetation.intersection(used_classes)
    if len(veg_classes) == 0:
        has_classification = False
    else:
        veg_classes = np.array(list(veg_classes))
        filter_on_return_number = False

    vegetation_indices = np.array([])
    if has_classification:
        vegetation_indices = np.where(np.isin(pc.classification, veg_classes))[0]

    elif filter_on_return_number:
        is_veg = pc.return_number != pc.num_returns

        # only reclassify points that are not already classified
        if len(pc.classification) == len(pc.points):
            is_veg = np.logical_and(is_veg, pc.classification == 1)
        vegetation_indices = np.where(is_veg)[0]

    return vegetation_indices


@register_model_method
def crop(pc: PointCloud, bounds: Bounds, xy_only=True) -> PointCloud:
    """
    Crop a `PointCloud` object only include point inside given bounds object

    Args:
        bounds (Bounds): The bounds to keep.

    Returns:
        PointCloud: A new `PointCloud` object with all points inside the bounds.
    """
    x_min, x_max = bounds.xmin, bounds.xmax
    y_min, y_max = bounds.ymin, bounds.ymax

    x_keep_idx = np.where((pc.points[:, 0] >= x_min) & (pc.points[:, 0] <= x_max))[0]
    y_keep_idx = np.where((pc.points[:, 1] >= y_min) & (pc.points[:, 1] <= y_max))[0]
    keep_idx = np.intersect1d(x_keep_idx, y_keep_idx)
    if not xy_only:
        z_min, z_max = bounds.zmin, bounds.zmax
        z_keep_idx = np.where((pc.points[:, 2] >= z_min) & (pc.points[:, 2] <= z_max))[
            0
        ]
        keep_idx = np.intersect1d(keep_idx, z_keep_idx)

    new_pc = pc.copy()
    new_pc.keep_points(keep_idx)
    return new_pc
