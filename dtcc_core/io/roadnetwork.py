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
from .vector_utils import (
    validate_vector_file,
    create_bounds_filter,
    determine_io_crs,
    safe_reproject_geometry,
    set_geometry_crs,
)
from .utils import get_epsg

HAS_GEOPANDAS = False
try:
    import geopandas as gpd
    import pandas as pd

    HAS_GEOPANDAS = True
except ImportError:
    warning("Geopandas not found, some functionality may be disabled")


def _load_fiona(filename, id_field="id", round_coordinates=2, load_geometry=True, bounds=None, target_crs=None):
    road_network = RoadNetwork()
    filename = validate_vector_file(filename)
    bounds_filter = create_bounds_filter(bounds, strategy='intersects')

    with fiona.open(filename) as src:
        # Read source CRS
        source_crs = get_epsg(src.crs)
        target_crs = determine_io_crs(source_crs, target_crs, context="road network")

        attr_keys = src.schema["properties"].keys()
        features = [f for f in src if not bounds_filter or bounds_filter['strategy'](bounds_filter['geometry'], shapely.geometry.shape(f["geometry"]))]
        attrs = [dict(f["properties"]) for f in features]

        shapely_geom = [shapely.geometry.shape(f["geometry"]) for f in features]
        shapely_geom = safe_reproject_geometry(
            shapely_geom, source_crs, target_crs, error_context="road network"
        )

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

    # Set CRS on road network
    set_geometry_crs(road_network, target_crs)

    return road_network


def load(
    filename, id_field="id", round_coordinates=2, load_geometry=True, bounds=None, target_crs=None
) -> RoadNetwork:
    """
    Load a road network from a supported vector file.

    Parameters
    ----------
    filename : str or Path
        Path to the input file.
    id_field : str, default "id"
        Attribute field containing edge identifiers.
    round_coordinates : int, default 2
        Number of decimal places to round coordinates to when building topology.
    load_geometry : bool, default True
        Whether to populate geometry on the resulting road network.
    bounds : Bounds, optional
        Bounding box filter; only features intersecting the bounds are loaded.
    target_crs : str, optional
        Target coordinate reference system (e.g., "EPSG:4326").
        If specified and different from source CRS, geometries will be
        automatically reprojected. If None, uses the file's native CRS.

    Returns
    -------
    RoadNetwork
        Parsed road network.

    Raises
    ------
    FileNotFoundError
        If ``filename`` does not exist.
    """
    filename = validate_vector_file(filename)
    return generic.load(
        filename,
        "road_network",
        RoadNetwork,
        _load_formats,
        id_field=id_field,
        round_coordinates=round_coordinates,
        load_geometry=load_geometry,
        bounds=bounds,
        target_crs=target_crs,
    )


def to_dataframe(road_network: RoadNetwork, crs=None):
    """
    Convert a road network to a GeoPandas DataFrame.

    Parameters
    ----------
    road_network : RoadNetwork
        Network to convert.
    crs : str or CRS, optional
        Coordinate reference system to set on the resulting GeoDataFrame.

    Returns
    -------
    geopandas.GeoDataFrame or None
        DataFrame with road geometry and attributes, or ``None`` if GeoPandas is unavailable.
    """
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
