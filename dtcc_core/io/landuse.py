from ..model import Landuse, LanduseClasses, GeometryType
from ..model.geometry import Surface, MultiSurface

from pathlib import Path
import fiona
import shapely.geometry
from .logging import info, warning, error
from . import generic
from .vector_utils import (
    validate_vector_file,
    determine_io_crs,
    safe_reproject_geometry,
    set_geometry_crs,
)
from .utils import get_epsg


LM_landuse_map = {
    "VATTEN": LanduseClasses.WATER,
    "SKOGLÖV": LanduseClasses.FOREST,
    "SKOGSMARK": LanduseClasses.FOREST,
    "SKOGBARR": LanduseClasses.FOREST,
    "ODLÅKER": LanduseClasses.FARMLAND,
    "ÅKERMARK": LanduseClasses.FARMLAND,
    "ODLFRUKT": LanduseClasses.FARMLAND,
    "ÖPMARK": LanduseClasses.GRASS,
    "BEBLÅG": LanduseClasses.LIGHT_URBAN,
    "BEBHÖG": LanduseClasses.HEAVY_URBAN,
    "BEBSLUT": LanduseClasses.URBAN,
    "BEBIND": LanduseClasses.INDUSTRIAL,
}

landuse_mappings = {"LM": LM_landuse_map}


def _get_landuse_class(properties, key, lookup_map):
    value = properties.get(key, "")
    if value:
        landuse_code = lookup_map.get(value, LanduseClasses.UNKNOWN)
    else:
        landuse_code = LanduseClasses.UNKNOWN
    if landuse_code == LanduseClasses.UNKNOWN:
        warning(f"Unknown landuse code {value}")
    return landuse_code


def _load_fiona(filename, landuse_field="DETALJTYP", landuse_datasource="LM", target_crs=None):
    landuse = Landuse()
    landuse_surfaces = MultiSurface()
    landuse_map = landuse_mappings.get(landuse_datasource, {})
    filename = validate_vector_file(filename)
    with fiona.open(filename) as src:
        # Read source CRS
        source_crs = get_epsg(src.crs)
        target_crs = determine_io_crs(source_crs, target_crs, context="landuse")

        for f in src:
            geom = shapely.geometry.shape(f["geometry"])
            geom = safe_reproject_geometry(geom, source_crs, target_crs)

            if geom.geom_type == "Polygon":
                surface = Surface()
                landuse_surfaces.surfaces.append(surface.from_polygon(geom, 0))
                landuse.landuses.append(
                    _get_landuse_class(f["properties"], landuse_field, LM_landuse_map)
                )
            elif geom.geom_type == "MultiPolygon":
                landuse_class = _get_landuse_class(
                    f["properties"], landuse_field, LM_landuse_map
                )
                for poly in list(geom.geoms):
                    surface = Surface()
                    landuse_surfaces.surfaces.append(surface.from_polygon(poly, 0))
                    landuse.landuses.append(landuse_class)
            else:
                warning(f"Unsupported geometry type {geom.geom_type}")

    landuse.add_geometry(landuse_surfaces, GeometryType.MULTISURFACE)
    # Set CRS on landuse object
    set_geometry_crs(landuse, target_crs)
    return landuse


def load(filename, landuse_field="DETALJTYP", landuse_datasource="LM", target_crs=None):
    """
    Load land use data from supported vector formats.

    Parameters
    ----------
    filename : str or Path
        Path to a GIS file (e.g., SHP, GeoJSON, GPKG).
    landuse_field : str, default "DETALJTYP"
        Attribute field used to classify land use.
    landuse_datasource : str, default "LM"
        Key selecting the land use mapping dictionary.
    target_crs : str, optional
        Target coordinate reference system (e.g., "EPSG:4326").
        If specified and different from source CRS, geometries will be
        automatically reprojected. If None, uses the file's native CRS.

    Returns
    -------
    Landuse
        Parsed land use dataset with geometries attached.

    Raises
    ------
    FileNotFoundError
        If ``filename`` does not exist.
    """
    filename = validate_vector_file(filename)
    return generic.load(
        filename,
        "landuse",
        Landuse,
        _load_formats,
        landuse_field="DETALJTYP",
        landuse_datasource="LM",
        target_crs=target_crs,
    )


_load_formats = {
    Landuse: {
        ".json": _load_fiona,
        ".shp": _load_fiona,
        ".geojson": _load_fiona,
        ".gpkg": _load_fiona,
    }
}
