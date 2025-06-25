from skimage.measure import regionprops_table

from dtcc_core.model import City, Building, PointCloud, Bounds, Terrain, GeometryType


from dtcc_core.builder import register_model_method

from dtcc_core.logging import info, warning, error

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


def _get_download_bounds(city_bounds, user_bounds):
    """
    Determine the bounds to use for downloading data.
    """

    if user_bounds is not None and not isinstance(user_bounds, Bounds):
        raise TypeError("bounds must be a Bounds object")

    if city_bounds is None:
        if user_bounds is None:
            error("No bounds found, either set city.bounds or pass in bounds")
        else:
            return bounds
    else:
        if user_bounds is not None:
            bounds = city_bounds.intersect(user_bounds)
        else:
            bounds = city_bounds
    return bounds


@register_model_method
def download_footprints(
    city: City,
    bounds: Union[Bounds, None] = None,
) -> City:
    """
    Download building footprints from a URL and load them into a City object.

    Args:
        city (City): The City object to load the footprints into.
        url (str): The URL to download the building footprints from.
        user_city_bounds (bool, optional): only load footprints within the set bounds

    Returns:
        City: The updated City object with the loaded footprints.
    """

    download_bounds = _get_download_bounds(city.bounds, bounds)

    footprints = io.data.download_footprints(bounds=download_bounds)
    city.add_buildings(footprints)
    city.calculate_bounds()
    return city


@register_model_method
def download_pointcloud(city: City, bounds: Union[Bounds, None] = None) -> City:
    """
    Download pointcloud from a URL and load it into a City object.

    Args:
        city (City): The City object to load the pointcloud into.
        bounds (Bounds, optional): The bounds to filter the pointcloud.

    Returns:
        City: The updated City object with the loaded pointcloud.
    """

    download_bounds = _get_download_bounds(city.bounds, bounds)

    pc = io.data.download_pointcloud(bounds=download_bounds)
    city.add_point_cloud(pc)
    city.calculate_bounds()
    return city
