"""Download the SCB 1 km population grid from the SCB WFS."""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import requests
import shapely.geometry
from affine import Affine
from shapely.geometry import box

from ...model import Bounds, PopulationGrid
from .cache import cache_dir
from .deso import (
    SCB_DESO_WFS_URL,
    SCB_WFS_REQUEST_HEADERS,
    _scb_request_timeout,
    download_deso_geodataframe,
)
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
POPULATION_GRID_CACHE_VERSION = 3
# DeSO areas tile all of Sweden, including its territorial sea, so they define
# which cells lie in Sweden.
COVERAGE_DESO_YEAR = 2025
# A cell overlapping Sweden by less than this (in m², of 1 km²) counts as
# outside it, so a DeSO boundary running along a cell edge adds no sliver.
_MIN_OVERLAP_M2 = 1.0
# SCB cell corners are up to a few millimetres off whole kilometres.
_CELL_SIZE_TOLERANCE_M = 0.1

# dtcc-core variable names mapped to SCB WFS fields.
POPULATION_GRID_VARIABLES = {
    "total": "beftotalt",
    "men": "man",
    "women": "kvinna",
    **{f"age_{lo}_{lo + 4}": f"ald{lo}_{lo + 4}" for lo in range(0, 100, 5)},
    "age_100_plus": "ald100w",
}


def download_population_grid(
    bounds: Bounds | Sequence[float],
    year: int = POPULATION_GRID_DEFAULT_YEAR,
    variable: str = "total",
    use_cache: bool = True,
) -> PopulationGrid:
    """Download one SCB 1 km population variable for ``bounds`` as a :class:`Raster`.

    Parameters
    ----------
    bounds : Bounds or sequence of float
        Area of interest in EPSG:3006: a ``Bounds``, ``(xmin, ymin, xmax, ymax)``
        or ``(xmin, ymin, zmin, xmax, ymax, zmax)``.
    year : int, default 2025
        Reference year from 2015 to 2025 (population on 31 December).
    variable : str, default "total"
        ``"total"``, ``"men"``, ``"women"``, or a five-year age band such as
        ``"age_20_24"`` or ``"age_100_plus"``.
    use_cache : bool, default True
        Reuse a cached download for the same bounds and year.

    Returns
    -------
    PopulationGrid
        Persons per 1 km cell on the SCB grid, extended outward to whole cells
        so it covers ``bounds``. SCB omits cells without residents, so those
        cells are 0. Cells with no part in Sweden are nodata (NaN).

    Raises
    ------
    RuntimeError
        If SCB cannot be reached or returns an incomplete or malformed
        response, for the population cells or the DeSO areas used to find
        which cells lie in Sweden.
    """
    if variable not in POPULATION_GRID_VARIABLES:
        raise ValueError(
            f"Unsupported population grid variable {variable!r}. "
            f"Supported: {tuple(POPULATION_GRID_VARIABLES)}."
        )
    bounds = _ensure_bounds(bounds)
    gdf = download_population_grid_geodataframe(bounds, year=year, use_cache=use_cache)
    areas = _sweden_areas(grid_bounds(bounds), use_cache=use_cache)
    return population_grid_raster(gdf, bounds, variable, year=year, sweden=areas)


def download_population_grid_geodataframe(
    bounds: Bounds | Sequence[float],
    year: int = POPULATION_GRID_DEFAULT_YEAR,
    use_cache: bool = True,
):
    """Download the SCB 1 km population cells covering ``bounds`` as a GeoDataFrame.

    The result is in EPSG:3006 with one row per inhabited cell, a ``cell_id``
    (SCB ruta id), a ``reference_date`` and one column per variable in
    :data:`POPULATION_GRID_VARIABLES`. It holds every inhabited cell of the
    whole 1 km cells covering ``bounds`` (see :func:`grid_bounds`), which is
    also the cache key, so nearby bounds share a cache entry only when they
    cover the same cells. Only complete, valid responses are cached.
    """
    if not HAS_GEOPANDAS:
        raise ImportError("GeoPandas is required to download the population grid.")
    if year not in POPULATION_GRID_YEARS:
        raise ValueError(
            f"Unsupported population grid year {year}. "
            f"Supported years: {POPULATION_GRID_YEARS}."
        )
    bounds = _ensure_bounds(bounds)
    cells = grid_bounds(bounds)

    cache_path = _cache_path(cells, year)
    if use_cache and cache_path.exists():
        info(f"Using cached SCB population grid {year} from {cache_path}")
        return gpd.read_file(cache_path)

    info(f"Downloading SCB 1 km population grid {year}")
    gdf = _cells_inside(_query_cells(cells, year), cells, year)
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


def grid_bounds(bounds: Bounds) -> Bounds:
    """Return ``bounds`` extended outward to whole 1 km SCB grid cells."""
    cell = POPULATION_GRID_CELL_SIZE
    xmin = math.floor(bounds.xmin / cell) * cell
    ymin = math.floor(bounds.ymin / cell) * cell
    xmax = max(math.ceil(bounds.xmax / cell) * cell, xmin + cell)
    ymax = max(math.ceil(bounds.ymax / cell) * cell, ymin + cell)
    return Bounds(xmin, ymin, xmax, ymax)


def population_grid_raster(
    gdf, bounds: Bounds, variable: str = "total", year: int | None = None, sweden=None
) -> PopulationGrid:
    """Rasterize one variable of a population grid GeoDataFrame onto the 1 km grid.

    ``sweden`` is a GeoDataFrame or sequence of polygons covering Sweden, such
    as DeSO areas. Uninhabited cells with no part in it become nodata (NaN);
    without it, every uninhabited cell is 0.
    """
    cell = POPULATION_GRID_CELL_SIZE
    cells = grid_bounds(bounds)
    width = round((cells.xmax - cells.xmin) / cell)
    height = round((cells.ymax - cells.ymin) / cell)
    data = np.zeros((height, width), dtype=np.float32)
    inhabited = np.zeros((height, width), dtype=bool)
    for geometry, value in zip(gdf.geometry, gdf[variable], strict=True):
        index = _grid_index(geometry, cells)
        if index is not None:
            data[index] = value
            inhabited[index] = True

    boundaries = None
    if sweden is not None:
        polygons = list(getattr(sweden, "geometry", sweden))
        coverage = (
            shapely.union_all(polygons) if polygons else shapely.geometry.Polygon()
        )
        data[_outside(coverage, cells, height, width) & ~inhabited] = np.nan
        extent = box(*cells.tuple)
        boundaries = (
            shapely.union_all([p.boundary for p in polygons]).intersection(extent)
            if polygons
            else None
        )

    transform = Affine(cell, 0, cells.xmin, 0, -cell, cells.ymax)
    return PopulationGrid(
        data=data,
        georef=transform,
        nodata=np.nan,
        crs="EPSG:3006",
        year=year,
        variable=variable,
        requested_bounds=Bounds(*bounds.tuple),
        area_boundaries=boundaries,
    )


def _grid_index(geometry, cells: Bounds) -> tuple[int, int] | None:
    """Return the (row, col) of a cell polygon in the whole-cell bounds ``cells``.

    Rows run north to south. Rounding absorbs the millimetre offsets of SCB
    cell corners. Returns None for a cell outside ``cells``.
    """
    cell = POPULATION_GRID_CELL_SIZE
    xmin, _, _, ymax = geometry.bounds
    col = round((xmin - cells.xmin) / cell)
    row = round((cells.ymax - ymax) / cell)
    width = round((cells.xmax - cells.xmin) / cell)
    height = round((cells.ymax - cells.ymin) / cell)
    if 0 <= row < height and 0 <= col < width:
        return row, col
    return None


def _outside(coverage, cells: Bounds, height: int, width: int) -> np.ndarray:
    """Return a (height, width) mask of grid cells with no part in ``coverage``."""
    cell = POPULATION_GRID_CELL_SIZE
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    x = cells.xmin + cols * cell
    y = cells.ymax - (rows + 1) * cell
    boxes = shapely.box(x, y, x + cell, y + cell)
    shapely.prepare(coverage)
    outside = ~shapely.intersects(coverage, boxes)
    # Only cells crossing the coverage boundary need an exact overlap area.
    edge = ~outside & ~shapely.contains_properly(coverage, boxes)
    overlap = shapely.area(shapely.intersection(coverage, boxes[edge]))
    outside[edge] = overlap < _MIN_OVERLAP_M2
    return outside


def _sweden_areas(cells: Bounds, use_cache: bool = True):
    try:
        return download_deso_geodataframe(
            cells, year=COVERAGE_DESO_YEAR, use_cache=use_cache
        )
    except RuntimeError as exc:
        raise RuntimeError(
            f"Could not download SCB DeSO areas to find which population grid cells "
            f"in {cells.tuple} lie in Sweden: {exc}. Retry later."
        ) from exc


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
    seen = set()
    total = None
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
        if not isinstance(payload, dict) or not isinstance(
            payload.get("features"), list
        ):
            raise _invalid_response(year, "the response has no 'features' list")
        page = payload["features"]
        # A server that ignores startIndex repeats cells, so stop at the first repeat.
        for cell_id in _cell_ids(page):
            if cell_id in seen:
                raise _invalid_response(year, f"cell {cell_id} appears more than once")
            seen.add(cell_id)
        features.extend(page)
        # The server may cap pages below the requested count, so a short page
        # does not mark the end: paging stops at SCB's reported total.
        matched = payload.get("numberMatched")
        if isinstance(matched, bool) or not isinstance(matched, int):
            raise _invalid_response(
                year, f"the response has no cell count (numberMatched={matched!r})"
            )
        if total is None:
            total = matched
        elif matched != total:
            raise _invalid_response(
                year, f"SCB's cell count changed from {total} to {matched} while paging"
            )
        if page and len(features) < total:
            continue
        if len(features) != total:
            raise _invalid_response(
                year, f"SCB reported {total} cells but returned {len(features)}"
            )
        return features


def _cell_ids(page: list) -> list[str]:
    """Return the cell ids of the features in ``page`` that have one."""
    ids = []
    for feature in page:
        props = feature.get("properties") if isinstance(feature, dict) else None
        if isinstance(props, dict) and props.get("rutid_scb"):
            ids.append(str(props["rutid_scb"]))
    return ids


def _invalid_response(year: int, problem: str) -> RuntimeError:
    return RuntimeError(
        f"Invalid SCB population grid {year} response: {problem}. "
        "Nothing was cached; retry later, and report it if it persists."
    )


def _cells_inside(features: list, cells: Bounds, year: int):
    """Keep the cells whose north-west corner, rounded to whole kilometres,
    lies inside the whole-cell bounds ``cells``.

    Neighbouring cells that only touch ``cells`` are dropped.

    Raises RuntimeError if a feature lacks a cell id, a 1 km square polygon
    geometry or a finite, non-negative count for any variable.
    """
    size = POPULATION_GRID_CELL_SIZE
    rows = []
    for index, feature in enumerate(features):
        if not isinstance(feature, dict):
            raise _invalid_response(year, f"feature {index} is not an object")
        props = feature.get("properties")
        if not isinstance(props, dict) or not props.get("rutid_scb"):
            raise _invalid_response(year, f"feature {index} has no cell id (rutid_scb)")
        cell_id = props["rutid_scb"]
        try:
            shape = shapely.geometry.shape(feature.get("geometry"))
        except Exception:
            shape = None
        if shape is None or shape.geom_type != "Polygon" or shape.is_empty:
            raise _invalid_response(year, f"cell {cell_id} has no polygon geometry")
        xmin, ymin, xmax, ymax = shape.bounds
        if (
            abs(xmax - xmin - size) > _CELL_SIZE_TOLERANCE_M
            or abs(ymax - ymin - size) > _CELL_SIZE_TOLERANCE_M
        ):
            raise _invalid_response(
                year,
                f"cell {cell_id} is not a 1 km cell "
                f"({xmax - xmin:.3f} m x {ymax - ymin:.3f} m)",
            )
        values = {}
        for variable, field in POPULATION_GRID_VARIABLES.items():
            value = props.get(field)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not 0 <= value < math.inf
            ):
                raise _invalid_response(
                    year, f"cell {cell_id} has no valid count for {field!r} ({value!r})"
                )
            values[variable] = value
        if _grid_index(shape, cells) is None:
            continue
        row = {
            "cell_id": cell_id,
            "reference_date": props.get("referenstid"),
            "geometry": shape,
            **values,
        }
        rows.append(row)
    columns = ["cell_id", "reference_date", *POPULATION_GRID_VARIABLES, "geometry"]
    return gpd.GeoDataFrame(rows, columns=columns, crs="EPSG:3006")


def _ensure_bounds(bounds) -> Bounds:
    """Bounds from a ``Bounds`` or a 4- or 6-value sequence in dataset order.

    Six values are ``(xmin, ymin, zmin, xmax, ymax, zmax)``, as the datasets
    take them. Values that are not finite, or a minimum above its maximum,
    raise ``ValueError``.
    """
    if not isinstance(bounds, Bounds):
        values = [float(v) for v in bounds]
        if len(values) == 4:
            bounds = Bounds(*values)
        elif len(values) == 6:
            xmin, ymin, zmin, xmax, ymax, zmax = values
            bounds = Bounds(xmin, ymin, xmax, ymax, zmin, zmax)
        else:
            raise ValueError(f"Bounds must have 4 or 6 values, got {len(values)}.")
    # A box in another order shows up as a minimum above its maximum.
    corners = (
        bounds.xmin,
        bounds.ymin,
        bounds.zmin,
        bounds.xmax,
        bounds.ymax,
        bounds.zmax,
    )
    if not all(math.isfinite(value) for value in corners):
        raise ValueError(f"Bounds {corners} must be finite numbers.")
    if (
        bounds.xmin > bounds.xmax
        or bounds.ymin > bounds.ymax
        or bounds.zmin > bounds.zmax
    ):
        raise ValueError(
            f"Bounds {bounds.tuple} must have each minimum at or below its maximum; "
            "six values go in the order (xmin, ymin, zmin, xmax, ymax, zmax)."
        )
    return bounds


def _cache_path(cells: Bounds, year: int) -> Path:
    """Cache path for the whole-cell bounds ``cells``, whose corners are whole kilometres."""
    xmin, ymin, xmax, ymax = (int(value) for value in cells.tuple)
    return (
        Path(cache_dir)
        / "downloaded_population_grid"
        / f"population_grid_v{POPULATION_GRID_CACHE_VERSION}_{year}_{xmin}_{ymin}_{xmax}_{ymax}.gpkg"
    )


__all__ = [
    "POPULATION_GRID_DEFAULT_YEAR",
    "POPULATION_GRID_VARIABLES",
    "POPULATION_GRID_YEARS",
    "download_population_grid",
    "download_population_grid_geodataframe",
    "grid_bounds",
    "population_grid_raster",
]
