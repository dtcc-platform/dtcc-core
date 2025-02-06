import fiona

from ..model.object import RoadNetwork, GeometryType
from ..model.geometry import LineString, MultiLineString
from . import generic

import shapely.geometry
import shapely.ops
import shapely.affinity

from pathlib import Path

import numpy as np

from .logging import info, warning, error

HAS_GEOPANDAS = False
try:
    import geopandas as gpd
    import pandas as pd

    HAS_GEOPANDAS = True
except ImportError:
    warning("Geopandas not found, some functionality may be disabled")


def _load_fiona(filename, id_field="id", round_coordinates=2, load_geometry=True, bounds=None):
    road_network = RoadNetwork()
    filename = Path(filename)
    if not filename.is_file():
        raise FileNotFoundError(f"File {filename} not found")
    bounds_filter = None
    if bounds is not None:
        bounds_filter = shapely.geometry.box(*bounds.tuple)

    with fiona.open(filename) as src:
        attr_keys = src.schema["properties"].keys()
        features = [f for f in src if bounds_filter is None or bounds_filter.intersects(shapely.geometry.shape(f["geometry"]))]
        attrs = [dict(f["properties"]) for f in features]

        shapely_geom = [shapely.geometry.shape(f["geometry"]) for f in features]
        coords = [(r.coords[0], r.coords[-1]) for r in shapely_geom]
        lengths = [r.length for r in shapely_geom]


    if round_coordinates is not None:
        rounded_coords = [
            (
                (round(c[0][0], round_coordinates), round(c[0][1], round_coordinates)),
                (round(c[1][0], round_coordinates), round(c[1][1], round_coordinates)),
            )
            for c in coords
        ]
    else:
        rounded_coords = coords

    coord_set = set([c[0] for c in rounded_coords] + [c[1] for c in rounded_coords])
    coord_lookup = {c: i for i, c in enumerate(coord_set)}
    vertices = np.array(list(coord_set))

    edges = []
    for c_start, c_end in rounded_coords:
        edges.append(
            (
                coord_lookup[c_start],
                coord_lookup[c_end],
            )
        )
    road_network.vertices = vertices
    road_network.edges = np.array(edges)
    road_network.length = np.array(lengths)

    road_network.attributes = {}
    for key in attr_keys:
        road_network.attributes[key] = [a[key] for a in attrs]

    rn_geom = MultiLineString()
    if load_geometry:
        for g in shapely_geom:
            ls = LineString()
            ls.from_shapely(g)
            rn_geom.linestrings.append(ls)

        road_network.geometry[GeometryType.MULTILINESTRING] = rn_geom

    return road_network


def load(
    filename, id_field="id", round_coordinates=2, load_geometry=True, bounds=None
) -> RoadNetwork:
    filename = Path(filename)
    if not filename.is_file():
        raise FileNotFoundError(f"File {filename} not found")
    return generic.load(
        filename,
        "road_network",
        RoadNetwork,
        _load_formats,
        id_field=id_field,
        round_coordinates=round_coordinates,
        load_geometry=load_geometry,
        bounds=bounds,
    )


def to_dataframe(road_network: RoadNetwork, crs=None):
    if HAS_GEOPANDAS is False:
        warning("Geopandas not found, unable to convert to dataframe")
        return None
    df = gpd.GeoDataFrame.from_dict(road_network.attributes)
    road_geometry = [
        linestring.to_shapely()
        for linestring in road_network.geometry[
            GeometryType.MULTILINESTRING
        ].linestrings
    ]
    df.set_geometry(road_geometry, inplace=True, crs=crs)
    return df


_load_formats = {
    RoadNetwork: {
        ".json": _load_fiona,
        ".shp": _load_fiona,
        ".geojson": _load_fiona,
        ".gpkg": _load_fiona,
    }
}
