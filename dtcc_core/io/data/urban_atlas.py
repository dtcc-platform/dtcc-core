"""Download Copernicus Urban Atlas land use/land cover from the EEA REST service."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import requests
import shapely.geometry
from shapely.geometry import box
from shapely.validation import make_valid

from ...model import Bounds, GeometryType, Landuse, LanduseClasses
from ...model.geometry import MultiSurface, Surface
from ..vector_utils import set_geometry_crs
from .cache import cache_dir
from .logging import info, warning

try:
    import geopandas as gpd

    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False


URBAN_ATLAS_SUPPORTED_YEARS = (2018,)
URBAN_ATLAS_SERVICE_URL = (
    "https://image.discomap.eea.europa.eu/arcgis/rest/services/UrbanAtlas/"
    "UA_UrbanAtlas_{year}/MapServer/2/query"
)
URBAN_ATLAS_DOI = {2018: "https://doi.org/10.2909/fb4dffa1-6ceb-4cc0-8372-1ed354c285e6"}
URBAN_ATLAS_ATTRIBUTION = "Generated using European Union's Copernicus Land Monitoring Service information; {doi}"
URBAN_ATLAS_PAGE_SIZE = 1000
# Stops paging if the service keeps reporting more results (a million features).
URBAN_ATLAS_MAX_PAGES = 1000
URBAN_ATLAS_CACHE_VERSION = 2
# "Other roads" is one polygon per Functional Urban Area (Göteborg: 8,284 ha,
# 1.1 M vertices). Unsimplified bbox queries that include it fail with HTTP 500,
# so it is fetched separately with server-side simplification.
URBAN_ATLAS_ROAD_CODE = "12220"
URBAN_ATLAS_ROAD_SIMPLIFY_M = 0.5
URBAN_ATLAS_CONNECT_TIMEOUT_SECONDS = float(
    os.environ.get("DTCC_URBAN_ATLAS_CONNECT_TIMEOUT", "10")
)
URBAN_ATLAS_READ_TIMEOUT_SECONDS = float(
    os.environ.get("DTCC_URBAN_ATLAS_READ_TIMEOUT", "120")
)
_MIN_PIECE_AREA_M2 = 0.01
URBAN_ATLAS_COLUMNS = ["ua_code", "ua_class", "fua_name", "geometry"]

# The 27 Urban Atlas 2018 classes mapped onto the coarser LanduseClasses.
# The original code and label are kept as surface-aligned attributes.
URBAN_ATLAS_LANDUSE_MAP = {
    "11100": LanduseClasses.HEAVY_URBAN,  # Continuous urban fabric (>80%)
    "11210": LanduseClasses.URBAN,  # Discontinuous dense (50-80%)
    "11220": LanduseClasses.URBAN,  # Discontinuous medium density (30-50%)
    "11230": LanduseClasses.LIGHT_URBAN,  # Discontinuous low density (10-30%)
    "11240": LanduseClasses.LIGHT_URBAN,  # Discontinuous very low density (<10%)
    "11300": LanduseClasses.LIGHT_URBAN,  # Isolated structures
    "12100": LanduseClasses.INDUSTRIAL,  # Industrial, commercial, public, military
    "12210": LanduseClasses.ROAD,  # Fast transit roads
    "12220": LanduseClasses.ROAD,  # Other roads
    "12230": LanduseClasses.RAIL,  # Railways
    "12300": LanduseClasses.INDUSTRIAL,  # Port areas
    "12400": LanduseClasses.INDUSTRIAL,  # Airports
    "13100": LanduseClasses.INDUSTRIAL,  # Mineral extraction and dump sites
    "13300": LanduseClasses.INDUSTRIAL,  # Construction sites
    "13400": LanduseClasses.GRASS,  # Land without current use
    "14100": LanduseClasses.GRASS,  # Green urban areas
    "14200": LanduseClasses.GRASS,  # Sports and leisure facilities
    "21000": LanduseClasses.FARMLAND,  # Arable land
    "22000": LanduseClasses.FARMLAND,  # Permanent crops
    "23000": LanduseClasses.FARMLAND,  # Pastures
    "24000": LanduseClasses.FARMLAND,  # Complex and mixed cultivation
    "25000": LanduseClasses.FARMLAND,  # Orchards
    "31000": LanduseClasses.FOREST,  # Forests
    "32000": LanduseClasses.GRASS,  # Herbaceous vegetation associations
    "33000": LanduseClasses.GRASS,  # Open spaces with little or no vegetation
    "40000": LanduseClasses.WATER,  # Wetlands
    "50000": LanduseClasses.WATER,  # Water
}


def download_urban_atlas(
    bounds: Bounds, year: int = 2018, use_cache: bool = True
) -> Landuse:
    """Download Urban Atlas land use for ``bounds`` as a :class:`Landuse`.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006.
    year : int, default 2018
        Urban Atlas reference year. Only 2018 is served anonymously.
    use_cache : bool, default True
        Reuse a cached download for the same bounds and year.

    Returns
    -------
    Landuse
        One surface per Urban Atlas polygon piece, clipped to ``bounds``, with
        surface-aligned ``ua_code``, ``ua_class`` and ``fua_name`` attributes
        and a ``requested_bounds`` attribute.

    Raises
    ------
    DatasetUpstreamError
        If the EEA service cannot be reached or returns an error or an
        invalid response. Its ``failure_class`` tells transient failures
        (connection, timeout, HTTP 5xx) from others.
    """
    bounds = _ensure_bounds(bounds)
    gdf = download_urban_atlas_geodataframe(bounds, year=year, use_cache=use_cache)
    landuse = landuse_from_urban_atlas(gdf, year=year)
    landuse.attributes["requested_bounds"] = list(bounds.tuple)
    return landuse


def download_urban_atlas_geodataframe(
    bounds: Bounds, year: int = 2018, use_cache: bool = True
):
    """Download Urban Atlas polygons clipped to ``bounds`` as a GeoDataFrame.

    The result is in EPSG:3006 with columns ``ua_code``, ``ua_class`` and
    ``fua_name``, one row per polygon piece. Only complete downloads are
    cached, keyed by the exact bounds. Raises :class:`DatasetUpstreamError`
    as :func:`download_urban_atlas` does.
    """
    if not HAS_GEOPANDAS:
        raise ImportError("GeoPandas is required to download Urban Atlas data.")
    if year not in URBAN_ATLAS_SUPPORTED_YEARS:
        raise ValueError(
            f"Unsupported Urban Atlas year {year}. "
            f"Supported years: {URBAN_ATLAS_SUPPORTED_YEARS}."
        )
    bounds = _ensure_bounds(bounds)

    cache_path = _cache_path(bounds, year)
    if use_cache and cache_path.exists():
        info(f"Using cached Urban Atlas {year} data from {cache_path}")
        return gpd.read_file(cache_path)

    info(f"Downloading Urban Atlas {year} land use from the EEA service")
    features = _query_features(
        bounds, year, f"code_{year} <> '{URBAN_ATLAS_ROAD_CODE}'"
    )
    features += _query_features(
        bounds,
        year,
        f"code_{year} = '{URBAN_ATLAS_ROAD_CODE}'",
        simplify=URBAN_ATLAS_ROAD_SIMPLIFY_M,
    )
    gdf = _clip_features(features, bounds, year)
    if len(gdf) == 0:
        warning(
            f"No Urban Atlas {year} data in bounds {bounds.tuple}. Urban Atlas only "
            "covers Functional Urban Areas."
        )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_name(cache_path.name + ".part.gpkg")
    gdf.to_file(tmp_path, layer="urban_atlas", driver="GPKG")
    os.replace(tmp_path, cache_path)
    return gdf


def landuse_from_urban_atlas(gdf, year: int = 2018) -> Landuse:
    """Convert an Urban Atlas GeoDataFrame to a :class:`Landuse`."""
    landuse = Landuse()
    surfaces = MultiSurface()
    codes, classes, fua_names = [], [], []
    for row in gdf.itertuples(index=False):
        code = str(row.ua_code)
        surfaces.surfaces.append(Surface().from_polygon(row.geometry, 0))
        landuse.landuses.append(
            URBAN_ATLAS_LANDUSE_MAP.get(code, LanduseClasses.UNKNOWN)
        )
        codes.append(code)
        classes.append(row.ua_class)
        fua_names.append(row.fua_name)
    landuse.add_geometry(surfaces, GeometryType.MULTISURFACE)
    set_geometry_crs(landuse, "EPSG:3006")
    doi = URBAN_ATLAS_DOI.get(year, "")
    landuse.attributes = {
        "ua_code": codes,
        "ua_class": classes,
        "fua_name": fua_names,
        "year": year,
        "source": f"Copernicus Urban Atlas {year}",
        "doi": doi,
        "attribution": URBAN_ATLAS_ATTRIBUTION.format(doi=doi),
    }
    return landuse


def _upstream_error(
    year: int, url: str, failure_class: str, message: str, status_code=None
):
    from ...datasets.dataset import DatasetUpstreamError

    return DatasetUpstreamError(
        dataset="urban_atlas",
        operation="query_features",
        target=url,
        failure_class=failure_class,
        status_code=status_code,
        message=f"Urban Atlas {year} query failed: {message}",
    )


def _query_features(
    bounds: Bounds, year: int, where: str, simplify: float | None = None
):
    from ...datasets.dataset import DatasetDescriptor

    url = URBAN_ATLAS_SERVICE_URL.format(year=year)
    params = {
        "geometry": ",".join(str(v) for v in bounds.tuple),
        "geometryType": "esriGeometryEnvelope",
        "inSR": 3006,
        "outSR": 3006,
        "spatialRel": "esriSpatialRelIntersects",
        "where": where,
        "outFields": f"code_{year},class_{year},fua_name",
        "orderByFields": "OBJECTID",
        "resultRecordCount": URBAN_ATLAS_PAGE_SIZE,
        "f": "geojson",
    }
    if simplify is not None:
        params["maxAllowableOffset"] = simplify

    features = []
    for _ in range(URBAN_ATLAS_MAX_PAGES):
        params["resultOffset"] = len(features)
        try:
            response = requests.get(
                url,
                params=params,
                timeout=(
                    URBAN_ATLAS_CONNECT_TIMEOUT_SECONDS,
                    URBAN_ATLAS_READ_TIMEOUT_SECONDS,
                ),
            )
            response.raise_for_status()
        except requests.exceptions.RequestException as exc:
            raise DatasetDescriptor.build_upstream_error(
                "urban_atlas", "query_features", url, exc
            ) from exc
        try:
            payload = response.json()
        except ValueError as exc:
            raise _upstream_error(
                year, url, "invalid_json", f"invalid JSON: {exc}"
            ) from exc
        page = _checked_page(payload, year, url)
        features.extend(page)
        exceeded = payload.get("exceededTransferLimit") or (
            payload.get("properties") or {}
        ).get("exceededTransferLimit")
        if not exceeded or not page:
            return features
    raise _upstream_error(
        year,
        url,
        "invalid_payload",
        f"the service still reported more results after {URBAN_ATLAS_MAX_PAGES} pages",
    )


def _checked_page(payload, year: int, url: str) -> list:
    """Return the features of one response page, or raise for an error or invalid page."""
    if not isinstance(payload, dict):
        raise _upstream_error(
            year, url, "invalid_payload", "the response is not a JSON object"
        )
    if "error" in payload:
        error = payload["error"]
        code = error.get("code") if isinstance(error, dict) else None
        failure_class = (
            "http_5xx" if isinstance(code, int) and code >= 500 else "invalid_payload"
        )
        raise _upstream_error(
            year,
            url,
            failure_class,
            f"service error: {error}",
            status_code=code if isinstance(code, int) else None,
        )
    page = payload.get("features")
    if not isinstance(page, list):
        raise _upstream_error(
            year, url, "invalid_payload", "the response has no 'features' list"
        )
    for index, feature in enumerate(page):
        props = feature.get("properties") if isinstance(feature, dict) else None
        if not isinstance(props, dict) or not props.get(f"code_{year}"):
            raise _upstream_error(
                year, url, "invalid_payload", f"feature {index} has no code_{year}"
            )
        if not feature.get("geometry"):
            raise _upstream_error(
                year, url, "invalid_payload", f"feature {index} has no geometry"
            )
    return page


def _clip_features(features, bounds: Bounds, year: int):
    area = box(*bounds.tuple)
    rows = []
    for feature in features:
        props = feature.get("properties") or {}
        geometry = feature.get("geometry")
        if not geometry:
            continue
        clipped = make_valid(shapely.geometry.shape(geometry)).intersection(area)
        for piece in getattr(clipped, "geoms", [clipped]):
            if piece.geom_type != "Polygon" or piece.area < _MIN_PIECE_AREA_M2:
                continue
            rows.append(
                {
                    "ua_code": str(props.get(f"code_{year}")),
                    "ua_class": props.get(f"class_{year}"),
                    "fua_name": props.get("fua_name"),
                    "geometry": piece,
                }
            )
    return gpd.GeoDataFrame(rows, columns=URBAN_ATLAS_COLUMNS, crs="EPSG:3006")


def empty_urban_atlas_geodataframe():
    """Return an empty Urban Atlas GeoDataFrame with the usual columns, in EPSG:3006."""
    return gpd.GeoDataFrame(
        columns=URBAN_ATLAS_COLUMNS, geometry="geometry", crs="EPSG:3006"
    )


def _ensure_bounds(bounds) -> Bounds:
    """Return ``bounds`` as Bounds; sequences are [minx, miny, maxx, maxy] or 3D."""
    if isinstance(bounds, Bounds):
        return bounds
    values = list(bounds)
    if len(values) == 4:
        return Bounds(*values)
    if len(values) == 6:
        minx, miny, _, maxx, maxy, _ = values
        return Bounds(minx, miny, maxx, maxy)
    raise ValueError(
        "bounds must be [minx, miny, maxx, maxy] or [minx, miny, minz, maxx, maxy, maxz], "
        f"got {len(values)} values"
    )


def _cache_path(bounds: Bounds, year: int) -> Path:
    exact = ",".join(repr(float(value)) for value in bounds.tuple)
    digest = hashlib.sha256(exact.encode("ascii")).hexdigest()[:16]
    return (
        Path(cache_dir)
        / "downloaded_urban_atlas"
        / f"urban_atlas_v{URBAN_ATLAS_CACHE_VERSION}_{year}_{digest}.gpkg"
    )


__all__ = [
    "URBAN_ATLAS_LANDUSE_MAP",
    "URBAN_ATLAS_SUPPORTED_YEARS",
    "download_urban_atlas",
    "download_urban_atlas_geodataframe",
    "empty_urban_atlas_geodataframe",
    "landuse_from_urban_atlas",
]
