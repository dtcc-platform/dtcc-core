import rasterio
from rasterio.features import rasterize

from ...model import Raster
from ..register import register_model_method
from shapely.geometry import Polygon


@register_model_method
def burn_polygons(
    raster: Raster, poly_values: tuple, all_touched: bool = True
) -> Raster:
    """
    Burn polygons into a raster with specified values.

    Parameters
    ----------
    raster : Raster
        Input raster to copy and burn values into.
    poly_values : tuple
        Iterable of ``(Polygon, value)`` pairs to burn.
    all_touched : bool, optional
        If True, rasterize all pixels touched by polygons; if False, only pixels
        whose center lies inside. Default is True.

    Returns
    -------
    Raster
        Copy of the input raster with polygons burned in.
    """
    out_raster = raster.copy()
    arr = out_raster.data

    transform = rasterio.transform.from_bounds(
        raster.bounds.west,
        raster.bounds.south,
        raster.bounds.east,
        raster.bounds.north,
        raster.width,
        raster.height,
    )
    for poly, value in poly_values:
        mask = rasterize(
            [poly],
            out_shape=arr.shape,
            transform=transform,
            fill=0,
            all_touched=all_touched,
        )
        arr[mask == 1] = value

    out_raster.data = arr
    return out_raster
