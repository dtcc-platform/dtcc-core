"""Download the SCB 1 km population grid from the SCB WFS."""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import requests
import shapely.geometry
from affine import Affine
from shapely.geometry import box

from ...model import Bounds, Raster
from .cache import cache_dir
from .deso import SCB_DESO_WFS_URL, SCB_WFS_REQUEST_HEADERS, _scb_request_timeout
from .logging import info, warning

try:
    import geopandas as gpd

    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False


POPULATION_GRID_YEARS = tuple(range(2015, 2026))
POPULATION_GRID_DEFAULT_YEAR = 2025
POPULATION_GRID_CELL_SIZE = 1000.0
POPULATION_GRID_PAGE_SIZE = 10000
# SCB cell corners are up to a few millimetres off whole kilometres, so a cell
# that only touches the bounds can overlap them by a sliver. Cells are 1 km².
_MIN_OVERLAP_M2 = 1.0

# dtcc-core variable names mapped to SCB WFS fields.
POPULATION_GRID_VARIABLES = {
    "total": "beftotalt",
    "men": "man",
    "women": "kvinna",
    **{f"age_{lo}_{lo + 4}": f"ald{lo}_{lo + 4}" for lo in range(0, 100, 5)},
    "age_100_plus": "ald100w",
}


def download_population_grid(
    bounds: Bounds,
    year: int = POPULATION_GRID_DEFAULT_YEAR,
    variable: str = "total",
    use_cache: bool = True,
) -> Raster:
    """Download one SCB 1 km population variable for ``bounds`` as a :class:`Raster`.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006.
    year : int, default 2025
        Reference year from 2015 to 2025 (population on 31 December).
    variable : str, default "total"
        ``"total"``, ``"men"``, ``"women"``, or a five-year age band such as
        ``"age_20_24"`` or ``"age_100_plus"``.
    use_cache : bool, default True
        Reuse a cached download for the same bounds and year.

    Returns
    -------
    Raster
        Persons per 1 km cell on the SCB grid, extended outward to whole cells
        so it covers ``bounds``. SCB omits cells without residents, so those
        cells, and cells outside Sweden, are 0.
    """
    if variable not in POPULATION_GRID_VARIABLES:
        raise ValueError(
            f"Unsupported population grid variable {variable!r}. "
            f"Supported: {tuple(POPULATION_GRID_VARIABLES)}."
        )
    gdf = download_population_grid_geodataframe(bounds, year=year, use_cache=use_cache)
    return population_grid_raster(gdf, _ensure_bounds(bounds), variable)


def download_population_grid_geodataframe(
    bounds: Bounds, year: int = POPULATION_GRID_DEFAULT_YEAR, use_cache: bool = True
):
    """Download SCB 1 km population cells overlapping ``bounds`` as a GeoDataFrame.

    The result is in EPSG:3006 with one row per inhabited cell, a ``cell_id``
    (SCB ruta id), a ``reference_date`` and one column per variable in
    :data:`POPULATION_GRID_VARIABLES`.
    """
    if not HAS_GEOPANDAS:
        raise ImportError("GeoPandas is required to download the population grid.")
    if year not in POPULATION_GRID_YEARS:
        raise ValueError(
            f"Unsupported population grid year {year}. "
            f"Supported years: {POPULATION_GRID_YEARS}."
        )
    bounds = _ensure_bounds(bounds)

    cache_path = _cache_path(bounds, year)
    if use_cache and cache_path.exists():
        info(f"Using cached SCB population grid {year} from {cache_path}")
        return gpd.read_file(cache_path)

    info(f"Downloading SCB 1 km population grid {year}")
    gdf = _cells_overlapping(_query_cells(bounds, year), bounds)
    if len(gdf) == 0:
        warning(
            f"No SCB population grid {year} cells in bounds {bounds.tuple}. "
            "The grid covers inhabited cells in Sweden only."
        )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_name(cache_path.name + ".part.gpkg")
    gdf.to_file(tmp_path, layer="population_grid", driver="GPKG")
    os.replace(tmp_path, cache_path)
    return gdf


def population_grid_raster(gdf, bounds: Bounds, variable: str = "total") -> Raster:
    """Rasterize one variable of a population grid GeoDataFrame onto the 1 km grid."""
    cell = POPULATION_GRID_CELL_SIZE
    xmin = math.floor(bounds.xmin / cell) * cell
    ymin = math.floor(bounds.ymin / cell) * cell
    xmax = max(math.ceil(bounds.xmax / cell) * cell, xmin + cell)
    ymax = max(math.ceil(bounds.ymax / cell) * cell, ymin + cell)
    width, height = round((xmax - xmin) / cell), round((ymax - ymin) / cell)
    data = np.zeros((height, width), dtype=np.float32)
    for geometry, value in zip(gdf.geometry, gdf[variable]):
        # Rounding absorbs the millimetre offsets of SCB cell corners.
        cell_xmin, _, _, cell_ymax = geometry.bounds
        col = round((cell_xmin - xmin) / cell)
        row = round((ymax - cell_ymax) / cell)
        if 0 <= row < height and 0 <= col < width and value is not None:
            data[row, col] = value
    transform = Affine(cell, 0, xmin, 0, -cell, ymax)
    return Raster(data=data, georef=transform, nodata=np.nan, crs="EPSG:3006")


def _query_cells(bounds: Bounds, year: int) -> list:
    params = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeNames": f"stat:befolkning_1km_{year}",
        "outputFormat": "application/json",
        "srsName": "EPSG:3006",
        "bbox": f"{bounds.xmin},{bounds.ymin},{bounds.xmax},{bounds.ymax},EPSG:3006",
        "sortBy": "objectid",
        "count": POPULATION_GRID_PAGE_SIZE,
    }
    features = []
    while True:
        params["startIndex"] = len(features)
        try:
            response = requests.get(
                SCB_DESO_WFS_URL,
                params=params,
                headers=SCB_WFS_REQUEST_HEADERS,
                timeout=_scb_request_timeout(),
            )
            response.raise_for_status()
            payload = response.json()
        except (requests.exceptions.RequestException, ValueError) as exc:
            raise RuntimeError(
                f"Failed to download SCB population grid {year}: {exc}"
            ) from exc
        page = payload.get("features", [])
        features.extend(page)
        if len(page) < POPULATION_GRID_PAGE_SIZE:
            return features


def _cells_overlapping(features: list, bounds: Bounds):
    """Keep cells whose area overlaps ``bounds``, not those only touching an edge."""
    area = box(*bounds.tuple)
    rows = []
    for feature in features:
        props = feature.get("properties") or {}
        geometry = feature.get("geometry")
        if not geometry:
            continue
        shape = shapely.geometry.shape(geometry)
        if shape.intersection(area).area < _MIN_OVERLAP_M2:
            continue
        row = {
            "cell_id": props.get("rutid_scb"),
            "reference_date": props.get("referenstid"),
            "geometry": shape,
        }
        for variable, field in POPULATION_GRID_VARIABLES.items():
            row[variable] = props.get(field)
        rows.append(row)
    columns = ["cell_id", "reference_date", *POPULATION_GRID_VARIABLES, "geometry"]
    return gpd.GeoDataFrame(rows, columns=columns, crs="EPSG:3006")


def _ensure_bounds(bounds) -> Bounds:
    if isinstance(bounds, Bounds):
        return bounds
    return Bounds(*bounds)


def _cache_path(bounds: Bounds, year: int) -> Path:
    xmin, ymin, xmax, ymax = (round(value) for value in bounds.tuple)
    return (
        Path(cache_dir)
        / "downloaded_population_grid"
        / f"population_grid_{year}_{xmin}_{ymin}_{xmax}_{ymax}.gpkg"
    )


__all__ = [
    "POPULATION_GRID_DEFAULT_YEAR",
    "POPULATION_GRID_VARIABLES",
    "POPULATION_GRID_YEARS",
    "download_population_grid",
    "download_population_grid_geodataframe",
    "population_grid_raster",
]
