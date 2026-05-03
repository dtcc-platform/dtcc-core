import os
from pathlib import Path
from typing import Any

import numpy as np
import requests
from shapely.geometry import box
from shapely.validation import make_valid

from dtcc_core.model import Bounds, DeSO, GeometryType, MultiSurface, Object, Surface
from .cache import cache_dir
from .logging import info, warning

try:
    import geopandas as gpd
    import pandas as pd

    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    gpd = None
    pd = None


SCB_DESO_WFS_URL = "https://geodata.scb.se/geoserver/stat/wfs"
SCB_DESO_SUPPORTED_YEARS = (2018, 2025)
_REQUEST_TIMEOUT_SECONDS = 60


def download_deso(
    bounds: Bounds,
    year: int = 2025,
    source: str = "SCB",
    use_cache: bool = True,
) -> DeSO:
    """Download DeSO polygons and return a DTCC DeSO object."""
    gdf = download_deso_geodataframe(
        bounds=bounds,
        year=year,
        source=source,
        use_cache=use_cache,
    )
    return deso_from_geodataframe(gdf, bounds=bounds, year=year, source=source)


def download_deso_geodataframe(
    bounds: Bounds,
    year: int = 2025,
    source: str = "SCB",
    use_cache: bool = True,
):
    """Download DeSO polygons from SCB WFS as a GeoDataFrame."""
    if not HAS_GEOPANDAS:
        raise ImportError("GeoPandas is required to download DeSO data.")
    if source != "SCB":
        raise ValueError("Only source='SCB' is supported for DeSO v1.")
    if year not in SCB_DESO_SUPPORTED_YEARS:
        raise ValueError(
            f"Unsupported DeSO year {year}. "
            f"Supported years: {SCB_DESO_SUPPORTED_YEARS}."
        )

    bounds = _ensure_bounds(bounds)
    cache_path = _cache_path(bounds, year)
    if use_cache and cache_path.exists():
        info(f"Using cached DeSO {year} data from {cache_path}")
        gdf = gpd.read_file(cache_path)
    else:
        url = _wfs_url(bounds, year)
        info(f"Downloading DeSO {year} polygons from SCB")
        response = requests.get(url, timeout=_REQUEST_TIMEOUT_SECONDS)
        response.raise_for_status()
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = cache_path.with_suffix(cache_path.suffix + ".part")
        tmp_path.write_bytes(response.content)
        os.replace(tmp_path, cache_path)
        gdf = gpd.read_file(cache_path)

    if gdf.crs is None:
        warning("SCB DeSO GeoPackage has no CRS metadata; assuming EPSG:3006.")
        gdf = gdf.set_crs("EPSG:3006")
    elif str(gdf.crs).upper() != "EPSG:3006":
        gdf = gdf.to_crs("EPSG:3006")

    return filter_deso_geodataframe(gdf, bounds)


def filter_deso_geodataframe(gdf, bounds: Bounds):
    """Filter a DeSO GeoDataFrame to features intersecting bounds."""
    bounds = _ensure_bounds(bounds)
    bbox = box(bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
    return gdf[gdf.geometry.intersects(bbox)].copy()


def deso_from_geodataframe(
    gdf,
    bounds: Bounds | None = None,
    year: int = 2025,
    source: str = "SCB",
) -> DeSO:
    """Convert a DeSO GeoDataFrame to a DTCC DeSO object."""
    deso = DeSO()
    deso.attributes = {
        "dataset": "deso",
        "source": source,
        "year": year,
        "area_count": int(len(gdf)),
        "source_url": _wfs_url(bounds, year) if bounds is not None else SCB_DESO_WFS_URL,
    }
    deso.transform.srs = "EPSG:3006"
    if bounds is not None:
        deso.bounds = _ensure_bounds(bounds)

    for _, row in gdf.iterrows():
        geometry = _geometry_to_multisurface(row.geometry)
        if geometry is None or len(geometry.surfaces) == 0:
            continue

        attributes = _row_attributes(row)
        area = Object()
        area.id = str(
            attributes.get("desokod")
            or attributes.get("objektidentitet")
            or attributes.get("objectid")
            or area.id
        )
        area.attributes = attributes
        area.transform.srs = "EPSG:3006"
        area.add_geometry(geometry, GeometryType.LOD0)
        deso.add_child(area)

    deso.attributes["area_count"] = len(deso)
    return deso


def _geometry_to_multisurface(geometry) -> MultiSurface | None:
    if geometry is None or geometry.is_empty:
        return None
    geometry = make_valid(geometry)
    polygons = []
    if geometry.geom_type == "Polygon":
        polygons = [geometry]
    elif geometry.geom_type == "MultiPolygon":
        polygons = list(geometry.geoms)
    else:
        warning(f"Skipping unsupported DeSO geometry type: {geometry.geom_type}")
        return None

    surfaces = []
    for polygon in polygons:
        if polygon.is_empty:
            continue
        surface = Surface()
        surface.from_polygon(polygon, height=0.0)
        surfaces.append(surface)
    return MultiSurface(surfaces=surfaces)


def _row_attributes(row) -> dict[str, Any]:
    attributes = {}
    for key, value in row.items():
        if key == "geometry":
            continue
        attributes[key] = _to_python_value(value)
    return attributes


def _to_python_value(value):
    if pd is not None:
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            pass
    if isinstance(value, np.generic):
        return value.item()
    if hasattr(value, "isoformat"):
        return value.isoformat()
    return value


def _ensure_bounds(bounds: Bounds | tuple | list) -> Bounds:
    if isinstance(bounds, Bounds):
        return bounds
    if isinstance(bounds, (tuple, list)):
        return Bounds(xmin=bounds[0], ymin=bounds[1], xmax=bounds[2], ymax=bounds[3])
    raise TypeError("bounds must be a dtcc_core.model.Bounds or a bbox tuple/list.")


def _wfs_url(bounds: Bounds | None, year: int) -> str:
    params = (
        f"service=WFS&version=1.1.0&request=GetFeature"
        f"&typeName=stat%3ADeSO_{year}"
        f"&outputFormat=geopackage&srsName=EPSG%3A3006"
    )
    if bounds is not None:
        bounds = _ensure_bounds(bounds)
        bbox = f"{bounds.xmin},{bounds.ymin},{bounds.xmax},{bounds.ymax},EPSG:3006"
        params += f"&bbox={bbox}"
    return f"{SCB_DESO_WFS_URL}?{params}"


def _cache_path(bounds: Bounds, year: int) -> Path:
    bounds = _ensure_bounds(bounds)
    xmin, ymin, xmax, ymax = (round(value) for value in bounds.tuple)
    return (
        Path(cache_dir)
        / "downloaded_deso"
        / f"deso_{year}_{xmin}_{ymin}_{xmax}_{ymax}.gpkg"
    )


__all__ = [
    "download_deso",
    "download_deso_geodataframe",
    "filter_deso_geodataframe",
    "deso_from_geodataframe",
    "SCB_DESO_WFS_URL",
    "SCB_DESO_SUPPORTED_YEARS",
]
