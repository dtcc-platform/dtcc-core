import hashlib
import json
import os
from collections.abc import Sequence
from math import prod
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

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

from .bbox import as_bounds
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
SCB_MAX_VALUES = 110_000
DESO_AGE_BANDS = ("-4", *(f"{age}-{age + 4}" for age in range(5, 80, 5)), "80-")
DESO_HOUSEHOLD_TYPES = ("ESUB", "ESMB", "SBUB", "SBMB", "OVRIGA")
SCB_DESO_STATISTICS = (
    "population",
    "households",
    "cars",
    "employment",
    "age",
    "household_types",
)
SCB_CONNECT_TIMEOUT_SECONDS = float(os.environ.get("DTCC_SCB_CONNECT_TIMEOUT", "8"))
SCB_READ_TIMEOUT_SECONDS = float(os.environ.get("DTCC_SCB_READ_TIMEOUT", "30"))
SCB_USER_AGENT = os.environ.get(
    "DTCC_SCB_USER_AGENT",
    "dtcc-core/0.9.8 (https://github.com/dtcc-platform/dtcc-core)",
)
SCB_WFS_REQUEST_HEADERS = {"User-Agent": SCB_USER_AGENT}
SCB_API_REQUEST_HEADERS = {
    "User-Agent": SCB_USER_AGENT,
    "Accept": "application/json",
}


_STATISTIC_QUERIES = {
    "population": {
        "url": f"{SCB_PXWEB_API_BASE_URL}/BE/BE0101/BE0101Y/FolkmDesoAldKon",
        "year_values": tuple(range(2018, 2026)),
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
        "year_values": tuple(range(2018, 2026)),
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
    "employment": {
        "url": f"{SCB_PXWEB_API_BASE_URL}/AM/AM0210/AM0210G/ArRegDesoStatusN",
        "year_values": (2020, 2021, 2022, 2023, 2024),
        "fields": [
            {
                "name": "employed_residents_total",
                "unit": "persons",
                "description": (
                    "Employed residents aged 15-74 by DeSO area from SCB BAS."
                ),
                "selections": {
                    "Kon": "1+2",
                    "Alder": "15-74",
                    "ContentsCode": "0000089X",
                },
            }
        ],
    },
}


def _category_fields(topic, base_field, dimension, codes):
    fields = []
    for code in codes:
        if topic == "age":
            label = {"-4": "0–4 år", "80-": "80– år"}.get(
                code, code.replace("-", "–") + " år"
            )
            suffix = {"-4": "0_4", "80-": "80_plus"}.get(code, code.replace("-", "_"))
            name, description = f"age_{suffix}", f"SCB DeSO population aged {label}."
        else:
            name = f"{topic}_{code}"
            description = f"SCB DeSO household count for category {code}."
        fields.append(
            {
                **base_field,
                "name": name,
                "description": description,
                "selections": {**base_field["selections"], dimension: code},
            }
        )
    return fields


for _topic, _base, _dimension, _codes in (
    ("age", "population", "Alder", DESO_AGE_BANDS),
    ("household_types", "households", "Hushallstyp", DESO_HOUSEHOLD_TYPES),
):
    _definition = _STATISTIC_QUERIES[_base]
    _STATISTIC_QUERIES[_topic] = {
        **_definition,
        "fields": _category_fields(
            _topic, _definition["fields"][0], _dimension, _codes
        ),
    }


def fetch_deso_control_table(codes: Sequence[str], year: int, topic: str) -> dict:
    """Return the SCB request and response for population or household controls.

    Large selections retain each raw request and response in `batches`, alongside
    the combined selection and rows in `request` and `response`.
    """
    categories = {"population": "age", "households": "household_types"}
    if topic not in categories:
        raise ValueError("DeSO control topic must be 'population' or 'households'.")
    definition = _STATISTIC_QUERIES[topic]
    columns = _selection_columns(
        definition["fields"] + _STATISTIC_QUERIES[categories[topic]]["fields"]
    )
    statistic_year = _statistic_year(topic, year)
    url = definition["url"]
    regions = list(
        dict.fromkeys(_scb_deso_region_code(code, statistic_year) for code in codes)
    )
    batches = list(_query_scb_batches(url, regions, statistic_year, columns, topic))
    if len(batches) == 1:
        return {"url": url, **batches[0]}
    response = dict(batches[0]["response"]) if batches else {}
    response["data"] = [row for batch in batches for row in batch["response"]["data"]]
    archive = {
        "url": url,
        "request": _scb_payload(regions, statistic_year, columns),
        "response": response,
    }
    if batches:
        archive["batches"] = batches
    return archive


def _selection_columns(fields):
    return {
        name: list(dict.fromkeys(field["selections"][name] for field in fields))
        for name in fields[0]["selections"]
    }


def _scb_payload(regions, year, columns):
    selections = {"Region": regions, **columns, "Tid": [str(year)]}
    return {
        "query": [
            {"code": code, "selection": {"filter": "item", "values": values}}
            for code, values in selections.items()
        ],
        "response": {"format": "JSON"},
    }


def _query_scb_batches(url, regions, year, columns, name):
    batch_size = SCB_MAX_VALUES // prod(map(len, columns.values()))
    for start in range(0, len(regions), batch_size):
        payload = _scb_payload(regions[start : start + batch_size], year, columns)
        yield {
            "request": payload,
            "response": _post_scb_statistic(url, payload, name, year),
        }


def _post_scb_statistic(url, payload, name, year):
    try:
        response = requests.post(
            url,
            json=payload,
            headers=SCB_API_REQUEST_HEADERS,
            timeout=_scb_request_timeout(),
        )
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(
            f"Failed to query SCB DeSO statistic {name} ({year})"
        ) from exc
    if not isinstance(data, dict) or not isinstance(data.get("data"), list):
        raise ValueError(f"The SCB DeSO {name} response has no data rows")
    return data


def download_deso(
    bounds: Bounds,
    year: int = 2025,
    source: str = "SCB",
    statistics: Sequence[str] | None = None,
    statistics_year: int | None = None,
    use_cache: bool = True,
) -> DeSO:
    """Download DeSO polygons and return a DTCC DeSO object."""
    if statistics:
        statistics = _normalize_statistics(statistics)
        _statistics_years(statistics, statistics_year, geometry_year=year)
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

    bounds = as_bounds(bounds)
    cache_path = _cache_path(bounds, year)
    if use_cache and cache_path.exists():
        info(f"Using cached DeSO {year} data from {cache_path}")
        gdf = gpd.read_file(cache_path)
    else:
        url = _wfs_url(bounds, year)
        info(f"Downloading DeSO {year} polygons from SCB")
        try:
            response = requests.get(
                url,
                headers=SCB_WFS_REQUEST_HEADERS,
                timeout=_scb_request_timeout(),
            )
            response.raise_for_status()
        except requests.exceptions.RequestException as exc:
            raise RuntimeError(
                f"Failed to download DeSO {year} polygons from SCB"
            ) from exc
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
    bounds = as_bounds(bounds)
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
        "source_url": _wfs_url(bounds, year)
        if bounds is not None
        else SCB_DESO_WFS_URL,
    }
    deso.transform.srs = "EPSG:3006"
    if bounds is not None:
        deso.bounds = as_bounds(bounds)

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
    statistics_years = _statistics_years(
        topics, year, geometry_year=deso.attributes.get("year", 2025)
    )
    unique_years = set(statistics_years.values())
    common_year = next(iter(unique_years)) if len(unique_years) == 1 else None
    fields = download_deso_statistics(
        codes=deso.codes,
        statistics=topics,
        year=common_year,
        source=source,
    )
    for field in fields:
        deso.attach_field(field)

    deso.attributes["statistics"] = topics
    deso.attributes["statistics_years"] = statistics_years
    if common_year is not None:
        deso.attributes["statistics_year"] = common_year
    else:
        deso.attributes.pop("statistics_year", None)
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
    if len(codes) == 0:
        return []

    fields = []
    statistics_years = _statistics_years(topics, year)
    topic_labels = ", ".join(f"{topic} ({statistics_years[topic]})" for topic in topics)
    info(f"Downloading DeSO statistics from SCB for {len(codes)} areas: {topic_labels}")
    for topic in topics:
        query = _STATISTIC_QUERIES[topic]
        statistics_year = statistics_years[topic]
        query_codes = [_scb_deso_region_code(code, statistics_year) for code in codes]
        specs = query["fields"]
        cache_paths = {
            spec["name"]: _statistics_cache_path(
                name=spec["name"],
                region_codes=query_codes,
                year=statistics_year,
                selections=spec["selections"],
            )
            for spec in specs
        }
        values = {
            name: _read_statistics_cache(path, len(query_codes))
            for name, path in cache_paths.items()
        }
        missing = [spec for spec in specs if values[spec["name"]] is None]
        for name, cached in values.items():
            if cached is not None:
                info(f"Using cached SCB DeSO statistic {name} ({statistics_year})")
        if missing:
            names = ", ".join(spec["name"] for spec in missing)
            info(f"Querying SCB DeSO statistics {names} ({statistics_year})")
            downloaded = _query_scb_statistics(
                query["url"], query_codes, statistics_year, missing
            )
            for name, counts in downloaded.items():
                _write_statistics_cache(cache_paths[name], counts)
            values.update(downloaded)
        for spec in specs:
            name = spec["name"]
            fields.append(
                Field(
                    name=name,
                    unit=spec["unit"],
                    description=spec["description"],
                    values=values[name].reshape((-1, 1)),
                    dim=1,
                )
            )
    return fields


def _query_scb_statistics(url, region_codes, year, fields):
    columns = _selection_columns(fields)
    dimensions = [name for name in columns if name != "ContentsCode"]
    varying = [name for name in dimensions if len(columns[name]) > 1]
    names = {
        tuple(spec["selections"][name] for name in varying): spec["name"]
        for spec in fields
    }
    data = {name: {} for name in names.values()}
    regions = list(dict.fromkeys(region_codes))
    for batch in _query_scb_batches(
        url, regions, year, columns, ", ".join(names.values())
    ):
        payload = batch["response"]
        response_keys = [
            item["code"] for item in payload.get("columns", []) if item["type"] == "d"
        ] or ["Region", *dimensions, "Tid"]
        region_idx = response_keys.index("Region")
        indices = [response_keys.index(name) for name in varying]
        last_index = max([region_idx, *indices])
        for row in payload["data"]:
            key, values = row.get("key", []), row.get("values", [])
            if not values or len(key) <= last_index:
                continue
            name = names.get(tuple(key[i] for i in indices))
            if name is not None:
                data[name][key[region_idx]] = _to_float(values[0])
    return {
        name: np.array([counts.get(code, np.nan) for code in region_codes], dtype=float)
        for name, counts in data.items()
    }


def _scb_request_timeout() -> tuple[float, float]:
    return (SCB_CONNECT_TIMEOUT_SECONDS, SCB_READ_TIMEOUT_SECONDS)


def _statistics_cache_path(
    name: str,
    region_codes: Sequence[str],
    year: int,
    selections: dict[str, str],
) -> Path:
    key = {
        "name": name,
        "region_codes": list(region_codes),
        "selections": selections,
        "year": year,
    }
    digest = hashlib.sha256(
        json.dumps(key, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    return (
        Path(cache_dir)
        / "downloaded_deso"
        / "statistics"
        / f"{name}_{year}_{digest}.npy"
    )


def _read_statistics_cache(path: Path, expected_len: int) -> np.ndarray | None:
    if not path.exists():
        return None
    try:
        values = np.load(path, allow_pickle=False)
    except (OSError, ValueError) as exc:
        warning(f"Ignoring unreadable cached SCB DeSO statistic {path}: {exc}")
        return None
    if len(values) != expected_len:
        warning(f"Ignoring cached SCB DeSO statistic with unexpected length: {path}")
        return None
    return np.asarray(values, dtype=float)


def _write_statistics_cache(path: Path, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = NamedTemporaryFile(dir=path.parent, suffix=".part", delete=False)
    tmp_path = Path(handle.name)
    try:
        with handle:
            np.save(handle, np.asarray(values, dtype=float))
        os.replace(tmp_path, path)
    finally:
        tmp_path.unlink(missing_ok=True)


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


def _statistics_years(
    topics: Sequence[str], year: int | None = None, *, geometry_year: int | None = None
) -> dict[str, int]:
    return {
        topic: _statistic_year(topic, year, geometry_year=geometry_year)
        for topic in topics
    }


def _statistic_year(
    topic: str, year: int | None = None, *, geometry_year: int | None = None
) -> int:
    supported_years = _STATISTIC_QUERIES[topic]["year_values"]
    if year is not None and year not in supported_years:
        supported = ", ".join(str(value) for value in supported_years)
        raise ValueError(
            f"Statistic '{topic}' does not support year {year}. "
            f"Supported years: {supported}."
        )
    candidates = supported_years if year is None else (year,)
    if geometry_year is not None:
        candidates = [
            value for value in candidates if _geometry_vintage(value) == geometry_year
        ]
        if not candidates:
            raise ValueError(
                f"Statistic '{topic}' has no compatible year for DeSO {geometry_year} geometry. "
                "Statistics through 2023 require DeSO 2018 geometry. "
                "Statistics from 2024 require DeSO 2025 geometry."
            )
    return max(candidates)


def _geometry_vintage(year: int) -> int:
    return 2018 if year <= 2023 else 2025


def _scb_deso_region_code(code: str, year: int) -> str:
    code = str(code)
    suffixed = code.endswith("_DeSO2025")
    if _geometry_vintage(year) == 2018:
        if suffixed:
            raise ValueError(
                f"DeSO 2025 code {code!r} cannot be used for statistics in {year}."
            )
        return code
    return code if suffixed else f"{code}_DeSO2025"


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


def _wfs_url(bounds: Bounds | None, year: int) -> str:
    params = (
        f"service=WFS&version=1.1.0&request=GetFeature"
        f"&typeName=stat%3ADeSO_{year}"
        f"&outputFormat=geopackage&srsName=EPSG%3A3006"
    )
    if bounds is not None:
        bounds = as_bounds(bounds)
        bbox = f"{bounds.xmin},{bounds.ymin},{bounds.xmax},{bounds.ymax},EPSG:3006"
        params += f"&bbox={bbox}"
    return f"{SCB_DESO_WFS_URL}?{params}"


def _cache_path(bounds: Bounds, year: int) -> Path:
    bounds = as_bounds(bounds)
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
    "fetch_deso_control_table",
    "attach_deso_statistics",
    "filter_deso_geodataframe",
    "deso_from_geodataframe",
    "SCB_DESO_WFS_URL",
    "SCB_DESO_SUPPORTED_YEARS",
    "SCB_DESO_STATISTICS",
]
