from collections import defaultdict
from functools import partial
from pathlib import Path

from shapely.geometry import Polygon

from ..model import City, GeometryType
from . import generic
from .cityjson import cityjson, write_cityjson
from .logging import error, info, warning
from .meshes import load_mesh_as_city
from .model import load_model, save_model

HAS_GEOPANDAS = False
try:
    import geopandas as gpd
    import pandas as pd

    HAS_GEOPANDAS = True
except ImportError:
    warning("Geopandas not found, some functionality may be disabled")


def _load_json(path, *, strict=False, extent_policy="validate", validate_schema=None):
    """Use the same CityJSON admission for public JSON and ZIP loading."""
    return cityjson.load(
        path,
        strict=strict,
        extent_policy=extent_policy,
        validate_schema=validate_schema,
    )


def _load_mesh_city(
    filename, lod=GeometryType.LOD1, merge_coplanar_surfaces=True
) -> City:
    return load_mesh_as_city(
        filename, lod=lod, merge_coplanar_surfaces=merge_coplanar_surfaces
    )


def load(path, **kwargs):
    """
    Load a City object from a file.

    Supports various formats including protobuf, CityJSON, and mesh formats.
    The format is automatically detected based on the file extension.
    Canonical .dtcc and strict=True CityJSON input validate the standard schema
    by default; validate_schema=False bypasses only semantic evaluation.
    Explicit True on CityJSON requires strict=True.

    Parameters
    ----------
    path : str or Path
        Path to the file to load.

    Returns
    -------
    City
        The loaded city object.
    """
    return generic.load(path, "city", City, _load_formats, **kwargs)


def save(city, path, **kwargs):
    """
    Save a City object to a file.

    Supports DTCC Protobuf (.dtcc) serialization.
    The format is automatically determined from the file extension.
    Canonical .dtcc and strict=True CityJSON output validate the standard schema
    by default; validate_schema=False bypasses only semantic evaluation.
    Explicit True on CityJSON requires strict=True.

    Parameters
    ----------
    city : City
        The city object to save.
    path : str or Path
        Path where the city will be saved.
    """
    return generic.save(city, path, "city", _save_formats, **kwargs)


def buildings_to_df(city: City, include_geometry=True, crs=None):
    """
    Convert city buildings to a pandas DataFrame or GeoDataFrame.

    Creates a tabular representation of building data with optional geometry
    information. Requires geopandas for geometric operations.

    Parameters
    ----------
    city : City
        The city object containing buildings to convert.
    include_geometry : bool, default=True
        If True, includes building footprint geometry in the DataFrame.
        Results in a GeoDataFrame if geopandas is available.
    crs : str or CRS object, optional
        Coordinate reference system for the geometry. Not currently used.

    Returns
    -------
    pandas.DataFrame or geopandas.GeoDataFrame or None
        DataFrame with building attributes and optionally geometry.
        Returns None if geopandas is not available when geometry is requested.
    """
    if not HAS_GEOPANDAS:
        warning("Geopandas not found, cannot convert buildings to dataframe")
        return None
    if include_geometry:
        try:
            import dtcc_core.builder
        except ImportError:
            warning("builder not found, cannot convert building geometry to dataframe")
            return None
    city_buildings = city.buildings

    building_attributes = city.get_building_attributes()
    if not include_geometry:
        return pd.DataFrame.from_dict(building_attributes)

    ## include geometry
    building_footprints = [b.get_footprint() for b in city_buildings]
    building_footprints = list(
        map(
            lambda x: x.to_polygon() if x is not None else Polygon(),
            building_footprints,
        )
    )

    df = gpd.GeoDataFrame(building_attributes, geometry=building_footprints)
    return df


def _save_cityjson(
    city: City,
    filename: str,
    scale: float = 0.001,
    *,
    strict=False,
    validate_schema=None,
):
    filename = Path(filename)
    write_cityjson.save(
        city, filename, scale=scale, strict=strict, validate_schema=validate_schema
    )


_load_formats = {
    City: {
        ".dtcc": partial(load_model, expected_type=City),
        ".json": _load_json,
        ".json.zip": _load_json,
        ".obj": _load_mesh_city,
        ".ply": _load_mesh_city,
        ".stl": _load_mesh_city,
        ".vtk": _load_mesh_city,
        ".vtu": _load_mesh_city,
    }
}

_save_formats = {
    City: {
        ".dtcc": save_model,
        ".json": _save_cityjson,
        ".json.zip": _save_cityjson,
    }
}
