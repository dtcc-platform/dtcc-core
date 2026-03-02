import scipy.ndimage
import numpy as np
import rasterio
from ..register import register_model_method
from ...model import Raster
from logging import info, warning, error

import skimage as ski


@register_model_method
def fill_holes(raster: Raster) -> Raster:
    """
    Fill nodata holes in a raster using nearest-neighbor values.

    Parameters
    ----------
    raster : Raster
        Raster with nodata holes to fill.

    Returns
    -------
    Raster
        Copy of the raster with nodata pixels filled.
    """
    filled_raster = raster.copy()
    data = filled_raster.data
    nodata = filled_raster.nodata
    mask = data == nodata
    if np.any(mask):
        info(f"filling {mask.sum()} holes in raster")
        ind = scipy.ndimage.distance_transform_edt(
            mask, return_distances=False, return_indices=True
        )
        data = data[tuple(ind)]
    filled_raster.data = data
    return filled_raster


@register_model_method
def fill_small_holes(raster: Raster, hole_size=1, nodata=None) -> Raster:
    """
    Fill small holes in a raster using nearest-neighbor interpolation.

    Parameters
    ----------
    raster : Raster
        Raster with nodata holes to fill.
    hole_size : int, optional
        Largest hole area to fill; non-positive values skip filling. Default is 1.
    nodata : float, optional
        Value treated as nodata. If ``None``, uses the raster nodata value.

    Returns
    -------
    Raster
        Raster with small holes filled.

    Raises
    ------
    ValueError
        If no nodata value is provided and the raster has no nodata.
    """
    if hole_size <= 0:
        return raster

    if nodata is None:
        nodata = raster.nodata
    if nodata is None:
        raise ValueError("No nodata value provided and raster has no nodata value.")
    mask_data = raster.data != nodata

    # convert hole_size from area to pixel count
    hole_size = abs(hole_size / (raster.cell_size[0] * raster.cell_size[1]))

    filled_mask_data = ski.morphology.remove_small_holes(
        mask_data, area_threshold=hole_size, out=mask_data
    )

    filled_raster = fill_holes(raster)
    # filled_holes_mask = np.logical_xor(mask_data, filled_mask_data)

    filled_raster.data *= filled_mask_data

    return filled_raster


@register_model_method
def resample(raster: Raster, cell_size=None, scale=None, method="bilinear"):
    """
    Resample a raster to a new cell size or by a scale factor.

    Parameters
    ----------
    raster : Raster
        Raster to resample.
    cell_size : float, optional
        New cell size in meters; overrides ``scale`` if provided.
    scale : float, optional
        Multiplicative factor for cell size; required if ``cell_size`` is None.
    method : {"bilinear", "nearest", "cubic"}, optional
        Resampling method; default is "bilinear".

    Returns
    -------
    Raster
        Resampled raster.

    Raises
    ------
    ValueError
        If neither ``cell_size`` nor ``scale`` is provided, or if ``method`` is invalid.
    """
    sample_methods = {
        "bilinear": 1,
        "nearest": 0,
        "cubic": 3,
    }
    if cell_size is None and scale is None:
        raise ValueError("Either cell_size or scale must be specified")
    if not method in sample_methods:
        raise ValueError(
            f"Invalid resampling method, use one of {list(sample_methods.keys())}"
        )
    _raster = raster.copy()
    if cell_size is not None:
        scale = cell_size / _raster.cell_size[0]
    if scale == 1:
        return _raster
    _raster.data = scipy.ndimage.zoom(
        _raster.data,
        scale,
        order=sample_methods[method],
        mode="nearest",
        grid_mode=True,
    )
    _raster.georef *= _raster.georef.scale(scale, scale)

    return _raster
