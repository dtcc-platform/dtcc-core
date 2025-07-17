from ..model import City, GeometryType
from pathlib import Path
from .cityjson import cityjson
from .logging import info, warning, error
from .meshes import load_mesh_as_city
from . import generic
import json
import zipfile
from collections import defaultdict

from shapely.geometry import Polygon

HAS_GEOPANDAS = False
try:
    import geopandas as gpd
    import pandas as pd

    HAS_GEOPANDAS = True
except ImportError:
    warning("Geopandas not found, some functionality may be disabled")


def _load_json(path):
    """Load a city from a file.

    Args:
        path (str or Path): Path to the file.

    Returns:
        City: The loaded city.
    """
    path = Path(path)
    if path.suffix == ".json":
        with open(path, "r") as file:
            data = json.load(file)
    elif path.suffix == ".zip":
        with zipfile.ZipFile(path, "r") as z:
            files = z.namelist()
            if len(files) != 1 or not files[0].endswith("json"):
                raise ValueError("Invalid cityjson zip file")
            with z.open(files[0]) as f:
                data = json.load(f)
    else:
        raise ValueError(f"Unknown file format: {path.suffix}")
    if data.get("type") == "CityJSON":
        return cityjson.load(data)
    else:
        raise ValueError(f"{path} is not a CityJSON file")



def _load_proto_city(filename) -> City:
    with open(filename, "rb") as f:
        city = City()
        city.from_proto(f.read())
    return city


def _load_mesh_city(filename, lod=GeometryType.LOD1, merge_coplanar_surfaces=True) -> City:
    return load_mesh_as_city(filename, lod=lod, merge_coplanar_surfaces=merge_coplanar_surfaces)


def _save_proto_city(city: City, filename):
    with open(filename, "wb") as f:
        f.write(city.to_proto().SerializeToString())


def load(path):
    """
    Load a City object from a file.
    
    Supports various formats including protobuf, CityJSON, and mesh formats.
    The format is automatically detected based on the file extension.
    
    Parameters
    ----------
    path : str or Path
        Path to the file to load.
    
    Returns
    -------
    City
        The loaded city object.
    """
    return generic.load(path, "city", City, _load_formats)


def save(city, path):
    """
    Save a City object to a file.
    
    Supports protobuf format (.pb, .pb2) for binary serialization.
    The format is automatically determined from the file extension.
    
    Parameters
    ----------
    city : City
        The city object to save.
    path : str or Path
        Path where the city will be saved.
    """
    return generic.save(city, path, "city", _save_formats)


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
            warning(
                "builder not found, cannot convert building geometry to dataframe"
            )
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


_load_formats = {
    City: {".pb": _load_proto_city,
           ".pb2": _load_proto_city,
           ".json": _load_json,
           ".json.zip": _load_json,
           ".obj": _load_mesh_city,
           ".ply": _load_mesh_city,
           ".stl": _load_mesh_city,
           ".vtk": _load_mesh_city,
           ".vtu": _load_mesh_city
           }}

_save_formats = {City: {".pb": _save_proto_city, ".pb2": _save_proto_city}}
