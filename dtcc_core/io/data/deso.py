import os
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import requests
from shapely.geometry import box
from shapely.validation import make_valid

from dtcc_core.model import (
    Bounds,
    DeSO,
    Field,
    GeometryType,
    MultiSurface,
    Object,
    Surface,
)
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
SCB_PXWEB_API_BASE_URL = "https://api.scb.se/OV0104/v1/doris/sv/ssd/START"
SCB_DESO_SUPPORTED_YEARS = (2018, 2025)
SCB_DESO_STATISTICS = ("population", "households", "cars")
_REQUEST_TIMEOUT_SECONDS = 60


_STATISTIC_QUERIES = {
    "population": {
        "url": f"{SCB_PXWEB_API_BASE_URL}/BE/BE0101/BE0101Y/FolkmDesoAldKon",
        "year_values": (2024, 2025),
        "fields": [
            {
                "name": "population_total",
                "unit": "persons",
                "description": "Total population by DeSO area from SCB.",
                "selections": {
                    "Alder": "totalt",
                    "Kon": "1+2",
                    "ContentsCode": "000007Y7",
                },
            }
        ],
    },
    "households": {
        "url": f"{SCB_PXWEB_API_BASE_URL}/BE/BE0101/BE0101Y/HushallDesoTyp",
        "year_values": (2024, 2025),
        "fields": [
            {
                "name": "households_total",
                "unit": "households",
                "description": "Total number of households by DeSO area from SCB.",
                "selections": {
                    "Hushallstyp": "TOTALT",
                    "ContentsCode": "000007Y1",
                },
            }
        ],
    },
    "cars": {
        "url": f"{SCB_PXWEB_API_BASE_URL}/TK/TK1001/TK1001Z/PersBilarDesoN",
        "year_values": (2024, 2025),
        "fields": [
            {
                "name": "cars_total",
                "unit": "cars",
                "description": (
                    "Total passenger cars registered to residents by DeSO area from SCB."
                ),
                "selections": {
                    "Bestand": "TOT",
                    "ContentsCode": "000007ZL",
                },
            },
            {
                "name": "cars_in_traffic",
                "unit": "cars",
                "description": (
                    "Passenger cars in traffic registered to residents by DeSO area from SCB."
                ),
                "selections": {
                    "Bestand": "ITRAF",
                    "ContentsCode": "000007ZL",
                },
            },
        ],
    },
}


def download_deso(
    bounds: Bounds,
    year: int = 2025,
    source: str = "SCB",
    statistics: Sequence[str] | None = None,
    statistics_year: int | None = None,
    use_cache: bool = True,
) -> DeSO:
    """Download DeSO polygons and return a DTCC DeSO object."""
    gdf = download_deso_geodataframe(
        bounds=bounds,
        year=year,
        source=source,
        use_cache=use_cache,
    )
    deso = deso_from_geodataframe(gdf, bounds=bounds, year=year, source=source)
    if statistics:
        attach_deso_statistics(
            deso,
            statistics=statistics,
            year=statistics_year,
            source=source,
        )
    return deso


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


def attach_deso_statistics(
    deso: DeSO,
    statistics: Sequence[str],
    year: int | None = None,
    source: str = "SCB",
) -> DeSO:
    """Attach SCB DeSO statistics as area-aligned fields."""
    if source != "SCB":
        raise ValueError("Only source='SCB' is supported for DeSO statistics v1.")
    if not statistics:
        return deso

    topics = _normalize_statistics(statistics)
    statistics_year = year or _default_statistics_year(deso)
    fields = download_deso_statistics(
        codes=deso.codes,
        statistics=topics,
        year=statistics_year,
        source=source,
    )
    for field in fields:
        deso.attach_field(field)

    deso.attributes["statistics"] = topics
    deso.attributes["statistics_year"] = statistics_year
    deso.attributes["statistics_source"] = "SCB Statistikdatabasen"
    return deso


def download_deso_statistics(
    codes: Sequence[str],
    statistics: Sequence[str],
    year: int | None = None,
    source: str = "SCB",
) -> list[Field]:
    """Download selected SCB statistics for DeSO codes as DTCC Fields."""
    if source != "SCB":
        raise ValueError("Only source='SCB' is supported for DeSO statistics v1.")
    topics = _normalize_statistics(statistics)
    if not codes:
        return []

    statistics_year = year or 2025
    query_codes = [_scb_deso_region_code(code, statistics_year) for code in codes]
    fields = []
    for topic in topics:
        query = _STATISTIC_QUERIES[topic]
        if statistics_year not in query["year_values"]:
            supported = ", ".join(str(value) for value in query["year_values"])
            raise ValueError(
                f"Statistic '{topic}' does not support year {statistics_year}. "
                f"Supported years: {supported}."
            )

        for field_spec in query["fields"]:
            values = _query_scb_statistic(
                url=query["url"],
                region_codes=query_codes,
                year=statistics_year,
                selections=field_spec["selections"],
            )
            fields.append(
                Field(
                    name=field_spec["name"],
                    unit=field_spec["unit"],
                    description=field_spec["description"],
                    values=np.asarray(values, dtype=float).reshape((-1, 1)),
                    dim=1,
                )
            )
    return fields


def _query_scb_statistic(
    url: str,
    region_codes: Sequence[str],
    year: int,
    selections: dict[str, str],
) -> np.ndarray:
    query = [
        {
            "code": "Region",
            "selection": {"filter": "item", "values": list(region_codes)},
        }
    ]
    for code, value in selections.items():
        query.append(
            {
                "code": code,
                "selection": {"filter": "item", "values": [value]},
            }
        )
    query.append(
        {
            "code": "Tid",
            "selection": {"filter": "item", "values": [str(year)]},
        }
    )
    payload = {"query": query, "response": {"format": "JSON"}}
    response = requests.post(url, json=payload, timeout=_REQUEST_TIMEOUT_SECONDS)
    response.raise_for_status()
    return _parse_scb_statistic_response(response.json(), region_codes)


def _parse_scb_statistic_response(
    payload: dict[str, Any],
    region_codes: Sequence[str],
) -> np.ndarray:
    values_by_region = {code: np.nan for code in region_codes}
    for row in payload.get("data", []):
        key = row.get("key", [])
        values = row.get("values", [])
        if not key or not values:
            continue
        values_by_region[key[0]] = _to_float(values[0])
    return np.asarray([values_by_region[code] for code in region_codes], dtype=float)


def _normalize_statistics(statistics: Sequence[str]) -> list[str]:
    topics = []
    for topic in statistics:
        normalized = str(topic).lower()
        if normalized not in SCB_DESO_STATISTICS:
            supported = ", ".join(SCB_DESO_STATISTICS)
            raise ValueError(
                f"Unsupported DeSO statistic '{topic}'. Supported statistics: {supported}."
            )
        if normalized not in topics:
            topics.append(normalized)
    return topics


def _default_statistics_year(deso: DeSO) -> int:
    year = deso.attributes.get("year", 2025)
    if year == 2025:
        return 2025
    raise ValueError("DeSO statistics v1 currently supports DeSO 2025 only.")


def _scb_deso_region_code(code: str, year: int) -> str:
    code = str(code)
    if year >= 2024 and "_DeSO2025" not in code:
        return f"{code}_DeSO2025"
    return code


def _to_float(value):
    if value in (None, "", "..", ".", "-"):
        return np.nan
    if isinstance(value, str):
        value = value.replace(" ", "").replace(",", ".")
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


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
    "download_deso_statistics",
    "attach_deso_statistics",
    "filter_deso_geodataframe",
    "deso_from_geodataframe",
    "SCB_DESO_WFS_URL",
    "SCB_DESO_SUPPORTED_YEARS",
    "SCB_DESO_STATISTICS",
]
