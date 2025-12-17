import rasterstats
from ...model import Raster
from shapely.geometry import Polygon
from typing import Union, List
from affine import Affine
from ..register import register_model_method


@register_model_method
def stats(raster: Raster, polygons: Union[Polygon, List[Polygon]], stats=["mean"]):
    """
    Compute zonal statistics for a raster within one or more polygons.

    Parameters
    ----------
    raster : Raster
        Raster to sample.
    polygons : Polygon or list[Polygon]
        Polygon(s) defining zones for statistics.
    stats : list[str], optional
        Statistics to compute. Supported: ``["count", "min", "max", "mean",
        "median", "majority", "minority", "unique", "sum", "std", "var",
        "percentile_X"]`` where ``X`` is between 0 and 100. Default is ``["mean"]``.

    Returns
    -------
    Any
        If one polygon and one stat: scalar; if one polygon and multiple stats:
        dict; if multiple polygons: list of dicts/scalars matching the request.
    """
    if isinstance(polygons, Polygon):
        polygons = [polygons]

    stats_str = " ".join(stats)
    rstats = rasterstats.zonal_stats(
        polygons, raster.data, affine=raster.georef, stats=stats_str
    )
    if len(stats) == 1:
        rstats = [s[stats_str] for s in rstats]
    if len(polygons) == 1:
        rstats = rstats[0]
    return rstats
