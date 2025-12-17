import numpy as np
from typing import List, Union

from ...model import PointCloud, Bounds
from ..logging import info, warning, error
from ..register import register_model_method

from .. import _dtcc_builder


def find_global_outliers(pc: PointCloud, margin: float) -> np.ndarray:
    """
    Find points whose Z is farther than ``margin`` standard deviations from the mean.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to analyze.
    margin : float
        Standard deviation threshold for outlier detection.

    Returns
    -------
    numpy.ndarray
        Indices of outlier points.
    """

    z_pts = pc.points[:, 2]
    z_mean = np.mean(z_pts)
    z_std = np.std(z_pts)
    outliers = np.where(np.abs(z_pts - z_mean) > margin * z_std)[0]
    return outliers


def remove_global_outliers(pc: PointCloud, margin: float = 3.0) -> PointCloud:
    """
    Remove points whose Z is farther than ``margin`` standard deviations from the mean.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to filter.
    margin : float, optional
        Standard deviation threshold for outlier detection; default is 3.0.

    Returns
    -------
    PointCloud
        Point cloud with outliers removed.
    """
    outliers = find_global_outliers(pc, margin)
    new_pc = pc.copy()
    return new_pc.remove_points(outliers)


def find_statistical_outliers(
    pc: PointCloud, neighbours: int, outlier_margin: float
) -> np.ndarray:
    """
    Find statistical outliers in a point cloud using nearest neighbor analysis.
    
    This function identifies points that are statistical outliers based on their
    distance to neighboring points using a k-nearest neighbor approach.
    
    Parameters
    ----------
    pc : PointCloud
        The point cloud to analyze for outliers.
    neighbours : int
        Number of nearest neighbors to consider for each point.
    outlier_margin : float
        Standard deviation margin for outlier detection.
        
    Returns
    -------
    np.ndarray
        Array of indices of points identified as statistical outliers.
    """
    outliers = _dtcc_builder.statistical_outlier_finder(
        pc.points, neighbours, outlier_margin
    )
    return outliers


def statistical_outlier_filter(pc: PointCloud, neighbours, outlier_margin):
    """
    Remove statistical outliers from a point cloud.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to filter.
    neighbours : int
        Number of nearest neighbors used in detection.
    outlier_margin : float
        Standard deviation margin for outlier detection.

    Returns
    -------
    PointCloud
        Point cloud with statistical outliers removed.
    """

    outliers = find_statistical_outliers(pc, neighbours, outlier_margin)
    new_pc = pc.copy()
    return new_pc.remove_points(outliers)


def find_classification(pc: PointCloud, classes: Union[int, List[int]]) -> np.ndarray:
    """
    Find point indices that match given classification codes.

    Parameters
    ----------
    pc : PointCloud
        Point cloud with classification set.
    classes : int or list[int]
        Classification code(s) to match.

    Returns
    -------
    numpy.ndarray
        Indices of points with the requested classification.
    """
    if len(pc.points) != len(pc.classification):
        warning("Pointcloud not classified")
        return np.array([])
    if isinstance(classes, int):
        classes = [classes]

    cls_indices = np.where(np.isin(pc.classification, classes))[0]

    return cls_indices


def classification_filter(
    pc: PointCloud, classes: Union[int, List[int]], keep: bool = False
):
    """
    Filter a point cloud by classification.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to filter.
    classes : int or list[int]
        Classification code(s) to keep or remove.
    keep : bool, optional
        If True, keep matching points; otherwise remove them. Default is False.

    Returns
    -------
    PointCloud
        Filtered point cloud.
    """

    if len(pc.points) != len(pc.classification):
        error("Pointcloud not classified")
        return pc
    cls_indices = find_classification(pc, classes)

    new_pc = pc.copy()
    if keep:
        new_pc.keep_points(cls_indices)
    else:
        new_pc.remove_points(cls_indices)

    return new_pc


def z_range_filter(pc: PointCloud, min=None, max=None):
    """
    Filter points outside a Z-value range.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to filter in place.
    min : float, optional
        Minimum Z to keep.
    max : float, optional
        Maximum Z to keep.

    Returns
    -------
    PointCloud
        Point cloud with points outside the range removed.
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


def remove_vegetation(pc: PointCloud) -> PointCloud:
    """
    Remove vegetation points from a point cloud.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to filter.

    Returns
    -------
    PointCloud
        Point cloud with vegetation removed.
    """
    new_pc = pc.copy()
    veg_indices = _find_vegetation(pc)
    new_pc.remove_points(veg_indices)
    return new_pc


def get_vegetation(pc: PointCloud) -> PointCloud:
    """
    Extract vegetation points from a point cloud.
    
    This function creates a new point cloud containing only points classified
    as vegetation, either through LiDAR classification codes or return number analysis.
    
    Parameters
    ----------
    pc : PointCloud
        The point cloud to extract vegetation from.
        
    Returns
    -------
    PointCloud
        A new point cloud containing only vegetation points.
    """
    new_pc = pc.copy()
    veg_indices = _find_vegetation(pc)
    new_pc.keep_points(veg_indices)
    return new_pc


def _find_vegetation(pc: PointCloud, filter_on_return_number=True):
    """
    Find indices of points classified as vegetation.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to inspect.
    filter_on_return_number : bool, optional
        If True and classifications are absent, use return number to infer vegetation.

    Returns
    -------
    numpy.ndarray
        Indices of vegetation points.
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


def pts_in_bounds(pc: PointCloud, bounds: Bounds, xy_only=True) -> np.ndarray:
    """
    Find indices of points within specified bounds.
    
    This function identifies which points in a point cloud fall within the given
    bounds, with options to consider only XY coordinates or include Z dimension.
    
    Parameters
    ----------
    pc : PointCloud
        The point cloud to filter.
    bounds : Bounds
        The bounds object defining the spatial extent.
    xy_only : bool, default True
        Whether to consider only XY coordinates (True) or include Z dimension (False).
        
    Returns
    -------
    np.ndarray
        Array of indices of points within the specified bounds.
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
    return keep_idx


def crop(pc: PointCloud, bounds: Bounds, xy_only=True) -> PointCloud:
    """
    Crop a point cloud to the given bounds.

    Parameters
    ----------
    pc : PointCloud
        Point cloud to crop.
    bounds : Bounds
        Bounds to keep.
    xy_only : bool, optional
        If True, filter on XY only; if False, include Z. Default is True.

    Returns
    -------
    PointCloud
        Point cloud with points inside the bounds retained.
    """

    new_pc = pc.copy()
    indices = pts_in_bounds(new_pc, bounds, xy_only=xy_only)
    new_pc.keep_points(indices)
    return new_pc
