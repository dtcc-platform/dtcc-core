from ...model import Raster, PointCloud
from ..register import register_model_method

import numpy as np
import skimage as ski


@register_model_method
def remove_small_masks(raster: Raster, min_size=1, nodata=None) -> Raster:
    """Remove small objects from a binary mask.

    This function removes objects smaller than `min_size` from a raster. An object is defined as a connected
    pixels not equal to `nodata`.


    Args:
        raster (Raster): The input raster where you want to remove small objects.
        min_size (int, optional): The minimum size of objects to retain in the mask. Objects smaller than this will be removed.
                                  The default is 1, which means that objects with 1 pixel will be removed.

        nodata: (float): which value to consider nodata. Any other value will be considered valid. If None, the nodata
        value of the raster will be used.


    Returns:
        ndarray: The binary mask array with small objects removed.
    """

    if min_size <= 0:
        return raster
    if nodata is None:
        nodata = raster.nodata
    if nodata is None:
        raise ValueError("No nodata value provided and raster has no nodata value.")
    mask_data = raster.data != nodata
    # inv_mask = 1 - mask_data
    objects = ski.measure.label(mask_data)
    # inv_objects = ski.measure.label(inv_mask)

    mask_data = ski.morphology.remove_small_objects(objects, min_size=min_size)
    mask_data = mask_data > 0
    # inv_mask = ski.morphology.remove_small_objects(inv_objects, min_size=min_size)
    # inv_mask = inv_objects ^ inv_mask
    # return inv_mask.astype(bool)
    # inv_mask = 1 - inv_mask
    # mask_data = np.logical_or(mask_data, inv_mask)
    masked_raster = raster.copy(no_data=True)
    masked_raster.data = mask_data * raster.data
    masked_raster.nodata = nodata

    return masked_raster


@register_model_method
def erode_small_lines(raster: Raster, neighborhood_size=0, nodata=None) -> Raster:
    """
    Remove small linear features from a raster using morphological operations.
    
    This function applies binary erosion followed by dilation to remove thin
    linear features and small objects from a raster while preserving larger features.
    
    Parameters
    ----------
    raster : Raster
        Input raster to process.
    neighborhood_size : int, default 0
        Size of the morphological neighborhood. If 0, uses default 3x3 neighborhood.
    nodata : float, optional
        Value to consider as nodata. If None, uses raster's nodata value.
        
    Returns
    -------
    Raster
        Processed raster with small lines removed.
        
    Raises
    ------
    ValueError
        If no nodata value is provided and raster has no nodata value.
    """
    if nodata is None:
        nodata = raster.nodata
    if nodata is None:
        raise ValueError("No nodata value provided and raster has no nodata value.")
    mask_data = raster.data != nodata
    if neighborhood_size == 0 or neighborhood_size is None:
        neighborhood_size = None
    else:
        neighborhood_size = np.ones((neighborhood_size, neighborhood_size), dtype=bool)

    erosion_mask = ski.morphology.binary_erosion(mask_data, neighborhood_size)
    erosion_mask = ski.morphology.binary_dilation(erosion_mask, neighborhood_size)

    masked_raster = raster.copy(no_data=True)
    masked_raster.data = erosion_mask * raster.data
    masked_raster.nodata = nodata

    return masked_raster
