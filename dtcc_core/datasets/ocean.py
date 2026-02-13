# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

"""Ocean dataset from SMHI Oceanographic Observations (OcObs).

This module provides a dataset interface to fetch the latest-hour snapshot
of oceanographic observations from the SMHI open-data OcObs download API
(opendata-download-ocobs.smhi.se).  The dataset returns a SensorCollection
containing measurement stations/platforms as child Objects.

Each station has:
- Point geometry for location (z = 0 for ocean platforms)
- One Field per requested oceanographic parameter (e.g. sea_temperature,
  sea_level, salinity)
- Metadata attributes (station ID, name, quality codes, etc.)

The data is fetched as one bulk CSV per parameter (station-set/all mode) and
then filtered to the requested bounding box on the client side.
"""

from typing import Optional, Literal, Tuple, List, Dict, Any, Sequence
from pydantic import Field
import re
import numpy as np
from datetime import datetime, timezone

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ..model.object import Object, SensorCollection
from ..model.geometry import Point
from ..model.values import Field as DtccField
from ..common import info
from ..reproject.reproject import reproject_array


# ── Parameter name mapping ───────────────────────────────────────────────

# SMHI OcObs parameter ID → stable DTCC-friendly English field name.
#
# Available parameters (as of 2026):
#   1  Våghöjd, signifikant 30 min      – significant wave height [m]
#   2  Strömriktning                      – current direction [°]
#   3  Strömhastighet                     – current speed [cm/s]
#   4  Salthalt                           – salinity [PSU]
#   5  Havstemperatur                     – sea temperature [°C]
#   6  Havsvattenstånd                    – sea level [cm]
#   7  Vågriktning, medelvärde 30 min    – mean wave direction [°]
#   8  Vågriktning vid Tp                 – wave direction at peak [°]
#   9  Vågperiod, peakvärde 30 min       – peak wave period [s]
#  10  Vågperiod, medelvärde 30 min      – mean wave period [s]
#  11  Våghöjd, maximal 30 min           – max wave height [m]
#  12  Havsvattenstånd, RW timvärde      – sea level RW hourly [cm]
#  13  Havsvattenstånd, minutvärde       – sea level minute [cm]
#  14  Havsvattenstånd, RW minutvärde    – sea level RW minute [cm]
#  15  Syrehalt                           – dissolved oxygen [ml/l]

PARAMETER_NAMES: Dict[int, str] = {
    1: "wave_height_significant",
    2: "current_direction",
    3: "current_speed",
    4: "salinity",
    5: "sea_temperature",
    6: "sea_level",
    7: "wave_direction_mean",
    8: "wave_direction_peak",
    9: "wave_period_peak",
    10: "wave_period_mean",
    11: "wave_height_max",
    12: "sea_level_rw_hourly",
    13: "sea_level_minute",
    14: "sea_level_rw_minute",
    15: "dissolved_oxygen",
}

# Default parameters to fetch if none specified
DEFAULT_PARAMETERS = [5, 6]  # sea temperature + sea level

# Reverse lookup: name → parameter ID.  Accepts the canonical DTCC name
# as well as common short aliases.
_NAME_TO_ID: Dict[str, int] = {}
for _id, _name in PARAMETER_NAMES.items():
    _NAME_TO_ID[_name] = _id

# Short aliases
_NAME_TO_ID.update(
    {
        "temperature": 5,
        "sea_temp": 5,
        "water_temperature": 5,
        "water_temp": 5,
        "level": 6,
        "water_level": 6,
        "salinity": 4,
        "salt": 4,
        "waves": 1,
        "wave_height": 1,
        "current": 3,
        "oxygen": 15,
    }
)


def _resolve_parameter(p) -> int:
    """Resolve a parameter specification to an integer ID.

    Accepts an int, a stringified int, or a name / alias.

    Raises
    ------
    ValueError
        If the name cannot be resolved.
    """
    if isinstance(p, int):
        return p
    s = str(p).strip()
    # Numeric string
    try:
        return int(s)
    except ValueError:
        pass
    key = s.lower().replace(" ", "_").replace("-", "_")
    if key in _NAME_TO_ID:
        return _NAME_TO_ID[key]
    raise ValueError(
        f"Unknown ocean parameter '{p}'. "
        f"Use an integer ID or one of: {', '.join(sorted(_NAME_TO_ID.keys()))}"
    )


# ── HTTP helpers ─────────────────────────────────────────────────────────


def _get_text(url: str, timeout_s: float = 10.0) -> str:
    """Fetch text content from *url* with error handling.

    Parameters
    ----------
    url : str
        The URL to fetch.
    timeout_s : float
        Request timeout in seconds.

    Returns
    -------
    str
        Response body as text.

    Raises
    ------
    RuntimeError
        If the request fails.
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError(
            "requests library required for ocean dataset. "
            "Install with: pip install requests"
        )

    try:
        response = requests.get(url, timeout=timeout_s)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to fetch {url}: {e}")


# ── Coordinate helpers ───────────────────────────────────────────────────


def _transform_bounds_to_wgs84(
    bounds: Tuple[float, float, float, float], crs: str
) -> Tuple[float, float, float, float]:
    """Transform bounds from *crs* to WGS84 (lon, lat).

    Parameters
    ----------
    bounds : tuple
        (xmin, ymin, xmax, ymax) in source CRS.
    crs : str
        Source coordinate reference system (e.g. ``"EPSG:3006"``).

    Returns
    -------
    tuple
        (lon_min, lat_min, lon_max, lat_max) in WGS84.
    """
    if crs.upper() in ("CRS84", "EPSG:4326", "WGS84"):
        return bounds

    xmin, ymin, xmax, ymax = bounds
    corners = np.array([[xmin, ymin, 0], [xmax, ymax, 0]])
    transformed = reproject_array(corners, crs, "EPSG:4326")
    return (
        transformed[0, 0],
        transformed[0, 1],
        transformed[1, 0],
        transformed[1, 1],
    )


# ── CSV parser ───────────────────────────────────────────────────────────


def _parse_latest_hour_csv(text: str):
    """Parse a SMHI OcObs latest-hour station-set/all CSV payload.

    The CSV structure (observed for OcObs parameters) is::

        Parameternamn;Beskrivning;Enhet
        Havstemperatur;null;°C

        StationsId;Stationsnamn;Latitude;Longitude;Havstemperatur;Kvalitet;;...
        2541;UDDEVALLA;58.3475;11.8948;0.9;O;;...
        ...

    Unlike metobs, there is no Height column, and column 4 contains the
    measurement value directly (with the parameter name as header).
    Metadata comments are appended after ``;;`` on some station rows.

    Parameters
    ----------
    text : str
        Full CSV response body.

    Returns
    -------
    meta : dict
        ``parameter_name``, ``description``, ``unit`` fields.
    records : list[dict]
        One dict per station with keys ``station_id`` (int), ``station_name``,
        ``lat``, ``lon``, ``value`` (float or NaN), ``quality``.
    """
    meta: Dict[str, Any] = {}
    records: List[Dict[str, Any]] = []

    lines = text.splitlines()

    # ── Pass 1: extract parameter metadata (first two non-empty lines) ──
    non_empty = [ln for ln in lines if ln.strip()]
    if len(non_empty) >= 2:
        meta_header = non_empty[0].split(";")
        meta_values = non_empty[1].split(";")
        if len(meta_values) >= 3:
            meta["parameter_name"] = meta_values[0].strip()
            meta["description"] = meta_values[1].strip()
            meta["unit"] = meta_values[2].strip()

    # ── Pass 2: find the column header row (starts with StationsId) ──────
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("StationsId"):
            header_idx = i
            break

    if header_idx is None:
        return meta, records

    # Parse column headers — structure is:
    # StationsId;Stationsnamn;Latitude;Longitude;<ParamName>;Kvalitet;;...
    # So value is at index 4, quality at index 5.
    val_idx = 4
    qual_idx = 5

    # ── Pass 3: parse station rows ───────────────────────────────────
    for line in lines[header_idx + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split(";")
        if len(parts) < qual_idx + 1:
            continue

        # Station ID must be numeric
        try:
            station_id = int(parts[0])
        except (ValueError, IndexError):
            continue

        try:
            lat = float(parts[2])
            lon = float(parts[3])
        except (ValueError, IndexError):
            continue

        # Value may be empty or non-numeric
        raw_val = parts[val_idx].strip() if val_idx < len(parts) else ""
        try:
            value = float(raw_val)
        except (ValueError, TypeError):
            value = float("nan")

        quality = parts[qual_idx].strip() if qual_idx < len(parts) else ""

        records.append(
            {
                "station_id": station_id,
                "station_name": parts[1].strip(),
                "lat": lat,
                "lon": lon,
                "value": value,
                "quality": quality,
            }
        )

    return meta, records


# ── Dataset definition ───────────────────────────────────────────────────


class OceanDatasetArgs(DatasetBaseArgs):
    """Arguments for the ocean dataset.

    Attributes
    ----------
    bounds : Sequence[float]
        Geographic bounds ``(xmin, ymin, xmax, ymax)`` or 6-element variant.
    crs : str
        Coordinate reference system for bounds and output coordinates.
    format : str, optional
        Output format (``"pb"`` for protobuf bytes, ``None`` for Python object).
    timeout_s : float
        HTTP request timeout in seconds.
    max_stations : int
        Maximum number of stations to include (safety limit).
    drop_missing : bool
        If ``True``, skip a station entirely when it has *no* valid value for
        any of the requested parameters.
    base_url : str
        SMHI OcObs API base URL (override for testing).
    version : str
        API version path element.
    period : Literal["latest-hour"]
        Observation period.  Only ``"latest-hour"`` is supported in v1.
    parameters : list
        SMHI parameter IDs to fetch (default: sea temperature + sea level).
    field_name_style : str
        ``"dtcc"`` for stable English identifiers, ``"smhi"`` for the Swedish
        parameter name returned by the API.
    """

    crs: str = Field("EPSG:3006", description="Coordinate reference system")
    format: Optional[Literal["pb"]] = Field(
        None, description='Output format ("pb" for protobuf)'
    )
    timeout_s: float = Field(10.0, description="HTTP timeout in seconds", gt=0)
    max_stations: int = Field(2000, description="Maximum stations to include", gt=0)
    drop_missing: bool = Field(
        True,
        description="Skip stations that have no valid value for any parameter",
    )
    base_url: str = Field(
        "https://opendata-download-ocobs.smhi.se/api",
        description="SMHI OcObs API base URL",
    )
    version: str = Field("latest", description="API version path element")
    period: Literal["latest-hour"] = Field(
        "latest-hour", description="Observation period (only latest-hour in v1)"
    )
    parameters: List = Field(
        default=[5, 6],
        description=(
            "Ocean parameters to fetch.  Each element can be an integer "
            "SMHI parameter ID or a name/alias string such as "
            "'sea_temperature', 'temperature', 'sea_level', 'level', "
            "'salinity', 'wave_height', 'current_speed', 'dissolved_oxygen'."
        ),
    )
    field_name_style: Literal["dtcc", "smhi"] = Field(
        "dtcc",
        description="Field naming style: dtcc (English) or smhi (Swedish API name)",
    )


class OceanDataset(DatasetDescriptor):
    """Oceanographic observations from SMHI OcObs (latest-hour snapshot).

    This dataset fetches the most recent hourly observations from the SMHI
    open-data oceanographic observations API.  It downloads one bulk CSV
    per parameter (Sweden-wide, station-set/all), filters to the requested
    bounding box, and returns a ``SensorCollection`` where each station is
    an ``Object`` with a ``Point`` geometry carrying one ``Field`` per
    parameter.

    Default parameters: sea temperature (5) and sea level (6).

    Example
    -------
    >>> import dtcc_core.datasets as datasets
    >>> sensors = datasets.ocean(
    ...     bounds=(674000, 6580000, 676000, 6582000),
    ...     parameters=["sea_temperature", "sea_level"],
    ... )
    >>> print(f"Found {len(sensors.stations())} stations")
    """

    name = "ocean"
    description = (
        "Oceanographic observations from SMHI OcObs (latest-hour snapshot) "
        "as a SensorCollection with one Field per parameter per station "
        "within the specified bounds."
    )
    ArgsModel = OceanDatasetArgs

    def build(self, args: OceanDatasetArgs):
        """Build the ocean dataset.

        Parameters
        ----------
        args : OceanDatasetArgs
            Validated dataset arguments.

        Returns
        -------
        SensorCollection or bytes
            ``SensorCollection`` object, or protobuf bytes if
            ``format="pb"``.
        """
        # ── Resolve parameter names to IDs ────────────────────────
        param_ids = [_resolve_parameter(p) for p in args.parameters]

        # ── Parse bounds ─────────────────────────────────────────────
        bounds = self.parse_bounds(args.bounds)
        bounds_tuple = (bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
        wgs84_bounds = _transform_bounds_to_wgs84(bounds_tuple, args.crs)
        lon_min, lat_min, lon_max, lat_max = wgs84_bounds

        info(
            f"Ocean dataset: fetching {len(param_ids)} parameter(s) "
            f"for bbox ({lon_min:.4f}, {lat_min:.4f}) – "
            f"({lon_max:.4f}, {lat_max:.4f}) WGS84"
        )

        # ── Fetch CSV per parameter and merge ────────────────────────
        # station_id -> { "station_name", "lat", "lon",
        #                  "fields": { field_name: (value, unit, quality) } }
        station_map: Dict[int, Dict[str, Any]] = {}
        param_meta: Dict[int, Dict[str, Any]] = {}

        for pid in param_ids:
            url = (
                f"{args.base_url}/version/{args.version}/parameter/{pid}"
                f"/station-set/all/period/{args.period}/data.csv"
            )
            info(f"  Fetching parameter {pid} …")

            try:
                text = _get_text(url, timeout_s=args.timeout_s)
            except RuntimeError as exc:
                info(f"  Warning: failed to fetch parameter {pid}: {exc}")
                continue

            meta, records = _parse_latest_hour_csv(text)
            param_meta[pid] = meta

            # Determine field name
            if args.field_name_style == "dtcc":
                field_name = PARAMETER_NAMES.get(pid, f"parameter_{pid}")
            else:
                field_name = meta.get("parameter_name", f"parameter_{pid}")

            unit = meta.get("unit", "")

            # Filter to bbox and merge
            for rec in records:
                if not (
                    lon_min <= rec["lon"] <= lon_max
                    and lat_min <= rec["lat"] <= lat_max
                ):
                    continue

                sid = rec["station_id"]
                if sid not in station_map:
                    station_map[sid] = {
                        "station_name": rec["station_name"],
                        "lat": rec["lat"],
                        "lon": rec["lon"],
                        "fields": {},
                    }

                station_map[sid]["fields"][field_name] = (
                    rec["value"],
                    unit,
                    rec["quality"],
                )

        info(f"  {len(station_map)} station(s) within bounds")

        # ── Optionally drop stations with all-missing values ─────────
        if args.drop_missing:
            station_map = {
                sid: data
                for sid, data in station_map.items()
                if any(not np.isnan(v) for v, _u, _q in data["fields"].values())
            }
            info(f"  {len(station_map)} station(s) after dropping all-missing")

        # ── Coordinate transform for output ──────────────────────────
        need_reproject = args.crs.upper() not in ("CRS84", "EPSG:4326", "WGS84")

        if need_reproject and station_map:
            sids = list(station_map.keys())
            lonlat = np.array(
                [[station_map[s]["lon"], station_map[s]["lat"], 0.0] for s in sids]
            )
            projected = reproject_array(lonlat, "EPSG:4326", args.crs)
            for i, sid in enumerate(sids):
                station_map[sid]["x"] = projected[i, 0]
                station_map[sid]["y"] = projected[i, 1]
        else:
            for sid, data in station_map.items():
                data["x"] = data["lon"]
                data["y"] = data["lat"]

        # ── Build SensorCollection ───────────────────────────────────
        sensor_collection = SensorCollection()
        sensor_collection.attributes = {
            "dtcc_type": "sensor_collection",
            "source": "smhi_ocobs",
            "dataset": "ocean",
            "crs": args.crs,
            "bounds": (
                f"({bounds_tuple[0]}, {bounds_tuple[1]}, "
                f"{bounds_tuple[2]}, {bounds_tuple[3]})"
            ),
            "retrieval_time": datetime.now(timezone.utc).isoformat(),
            "period": args.period,
            "parameters": param_ids,
        }

        # Collect field→unit mapping for collection-level metadata
        parameter_fields: Dict[str, str] = {}
        first_field_name = None
        first_field_unit = ""
        for pid, pmeta in param_meta.items():
            if args.field_name_style == "dtcc":
                fname = PARAMETER_NAMES.get(pid, f"parameter_{pid}")
            else:
                fname = pmeta.get("parameter_name", f"parameter_{pid}")
            parameter_fields[fname] = pmeta.get("unit", "")
            if first_field_name is None:
                first_field_name = fname
                first_field_unit = pmeta.get("unit", "")
        sensor_collection.attributes["parameter_fields"] = parameter_fields

        # Apply max_stations limit
        sorted_sids = sorted(station_map.keys())
        if len(sorted_sids) > args.max_stations:
            sorted_sids = sorted_sids[: args.max_stations]

        for sid in sorted_sids:
            data = station_map[sid]

            station = Object()
            station.attributes = {
                "station_id": sid,
                "station_name": data["station_name"],
            }

            # Set 'value' and 'unit' from the first field so that
            # SensorCollection.__str__ can show it in the summary.
            if first_field_name and first_field_name in data["fields"]:
                fv, fu, _fq = data["fields"][first_field_name]
                if not np.isnan(fv):
                    station.attributes["value"] = float(fv)
                station.attributes["unit"] = fu

            # Add per-parameter quality codes
            for fname, (_val, _unit, qual) in data["fields"].items():
                station.attributes[f"q_{fname}"] = qual

            point = Point(x=data["x"], y=data["y"], z=0.0)
            point.fields = []

            for fname, (val, unit, _qual) in data["fields"].items():
                field = DtccField()
                field.name = fname
                field.unit = unit
                field.dim = 1
                field.values = np.array([val], dtype=np.float32)
                point.fields.append(field)

            station.geometry["location"] = point
            sensor_collection.add_station(station)

        sensor_collection.calculate_bounds()

        info(
            f"Ocean dataset: returning {len(sensor_collection.stations())} "
            f"station(s) with {len(parameter_fields)} field(s)"
        )

        if args.format == "pb":
            return sensor_collection.to_proto().SerializeToString()
        return sensor_collection


__all__ = ["OceanDatasetArgs", "OceanDataset"]
