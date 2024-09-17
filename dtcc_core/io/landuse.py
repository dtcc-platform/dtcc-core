from ..model import Landuse, LanduseClasses, GeometryType
from ..model.geometry import Surface, MultiSurface

from pathlib import Path
import fiona
import shapely.geometry
from .logging import info, warning, error
from . import generic


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


def _load_fiona(filename, landuse_field="DETALJTYP", landuse_datasource="LM"):
    landuse = Landuse()
    landuse_surfaces = MultiSurface()
    landuse_map = landuse_mappings.get(landuse_datasource, {})
    filename = Path(filename)
    if not filename.is_file():
        raise FileNotFoundError(f"File {filename} not found")
    with fiona.open(filename) as src:
        for f in src:
            geom = shapely.geometry.shape(f["geometry"])
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
    return landuse


def load(filename, landuse_field="DETALJTYP", landuse_datasource="LM"):
    filename = Path(filename)
    if not filename.is_file():
        raise FileNotFoundError(f"File {filename} not found")
    return generic.load(
        filename,
        "landuse",
        Landuse,
        _load_formats,
        landuse_field="DETALJTYP",
        landuse_datasource="LM",
    )


_load_formats = {
    Landuse: {
        ".json": _load_fiona,
        ".shp": _load_fiona,
        ".geojson": _load_fiona,
        ".gpkg": _load_fiona,
    }
}
