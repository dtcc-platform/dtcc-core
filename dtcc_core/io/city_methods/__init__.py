from skimage.measure import regionprops_table

from dtcc_core.model import City, Building, PointCloud, Bounds, Terrain, GeometryType


from dtcc_core.builder import register_model_method

import dtcc_core.io as io

from pathlib import Path
from typing import Union


@register_model_method
def load_footprints(city: City, path: Union[str, Path], user_city_bounds=False) -> City:
    """
    Load building footprints from a shapefile or geopackage into a City object.

    Args:
        city (City): The City object to load the footprints into.
        path (str): The path to the shapefile containing the building footprints.
        user_city_bounds (bool, optional): only load footprints within the set bounds

    Returns:
        City: The updated City object with the loaded footprints.
    """

    if user_city_bounds:
        if city.bounds is None or city.bounds.area == 0:
            raise ValueError(
                "City bounds must be set before loading footprints with user_city_bounds=True"
            )
        footprints = io.footprints.load(path, bounds=city.bounds)
    else:
        bounds = io.footprints.building_bounds(path, 2)
        footprints = io.footprints.load(path)
        if city.bounds is None:
            city.bounds = bounds
        else:
            city.bounds.union(bounds)
    city.add_buildings(footprints)
    return city


@register_model_method
def load_pointcloud(
    city: City,
    path: Union[Path, str],
    remove_global_outliers: float = 3.0,
    user_city_bounds=False,
) -> City:
    """
    Load pointcloud from las or csv file into city object
    """
    if user_city_bounds:
        if city.bounds is None or city.bounds.area == 0:
            raise ValueError(
                "City bounds must be set before loading footprints with user_city_bounds=True"
            )
        pc = io.pointcloud.load(path, bounds=city.bounds)
    else:
        pc = io.load_pointcloud(path)
        bounds = pc.bounds
        if city.bounds is None:
            city.bounds = bounds
        else:
            city.bounds.union(bounds)
    if remove_global_outliers > 0:
        pc = pc.remove_global_outliers(remove_global_outliers)
    city.add_geometry(pc, GeometryType.POINT_CLOUD)
    return city


@register_model_method
def load_terrain_raster(city: City, path: Union[Path, str]) -> City:
    raster = io.load_raster(path)
    city.add_terrain(raster)
    return City
