from ..model import City
from pathlib import Path
from .cityjson import cityjson
from .logging import info, warning, error
from . import generic
import json
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
        if data.get("type") == "CityJSON":
            return cityjson.load(path)
        else:
            raise ValueError(f"{path} is not a CityJSON file")
    else:
        raise ValueError(f"Unknown file format: {path.suffix}")


def _load_proto_city(filename) -> City:
    with open(filename, "rb") as f:
        city = City()
        city.from_proto(f.read())
    return city


def _save_proto_city(city: City, filename):
    with open(filename, "wb") as f:
        f.write(city.to_proto().SerializeToString())


def load(path):
    return generic.load(path, "city", City, _load_formats)


def save(city, path):
    return generic.save(city, path, "city", _save_formats)


def buildings_to_df(city: City, include_geometry=True, crs=None):
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
    City: {".pb": _load_proto_city, ".pb2": _load_proto_city, ".json": _load_json}
}

_save_formats = {City: {".pb": _save_proto_city, ".pb2": _save_proto_city}}
