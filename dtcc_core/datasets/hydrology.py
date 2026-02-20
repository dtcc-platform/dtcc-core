# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

"""Hydrology dataset from SMHI Hydrological Observations (HydroObs).

This module provides a dataset interface to fetch the latest-day snapshot
of hydrological observations from the SMHI open-data HydroObs download API
(opendata-download-hydroobs.smhi.se).  The dataset returns a SensorCollection
containing measurement stations as child Objects.

Each station has:
- Point geometry for location
- One Field per requested hydrological parameter (e.g. discharge, water_level)
- Metadata attributes (station ID, name, catchment, quality codes, etc.)

Unlike the metobs API used by the weather dataset, HydroObs has no bulk
``station-set/all`` endpoint.  Instead, this module:

1. Fetches the station list for each parameter (JSON).
2. Filters stations by bounding box and active status.
3. Fetches ``latest-day`` data per station (JSON).
4. Merges results and builds a ``SensorCollection``.
"""

from typing import Optional, Literal, Tuple, List, Dict, Any, Sequence
from pydantic import Field
import numpy as np
from datetime import datetime, timezone

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ..model.object import Object, SensorCollection
from ..model.geometry import Point
from ..model.values import Field as DtccField
from ..common import info
from ..reproject.reproject import reproject_array


# ── Parameter name mapping ───────────────────────────────────────────────

# SMHI HydroObs parameter ID → stable DTCC-friendly English field name.
#
# Available parameters (as of 2026):
#   1  Vattenföring (Dygn)        – daily discharge [m³/s]
#   2  Vattenföring (15 min)      – 15-minute discharge [m³/s]
#   3  Vattenstånd                – water level [cm]
#   4  Vattendragstemperatur      – stream temperature [°C]
#   5  Isläggning                 – ice formation [null]
#   6  Islossning                 – ice breakup [null]
#   7  Istjocklek                 – ice thickness [cm]
#   8  Snödensitet                – snow density [g/cm³]
#   9  Vatteninnehåll             – water equivalent [mm]
#  10  Vattenföring (Månad)       – monthly discharge [m³/s]

PARAMETER_NAMES: Dict[int, str] = {
    1: "discharge_daily",
    2: "discharge_15min",
    3: "water_level",
    4: "water_temperature",
    5: "ice_formation",
    6: "ice_breakup",
    7: "ice_thickness",
    8: "snow_density",
    9: "water_equivalent",
    10: "discharge_monthly",
}

# Reverse lookup: name → parameter ID.  Accepts the canonical DTCC name
# as well as common short aliases so users can write e.g.
# parameters=["discharge", "water_level"] instead of [1, 3].
_NAME_TO_ID: Dict[str, int] = {}
for _id, _name in PARAMETER_NAMES.items():
    _NAME_TO_ID[_name] = _id

# Short aliases
_NAME_TO_ID.update(
    {
        "discharge": 1,
        "flow": 1,
        "flow_daily": 1,
        "flow_15min": 2,
        "level": 3,
        "water_temp": 4,
        "temperature": 4,
        "ice": 7,
        "snow": 8,
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
        f"Unknown hydrology parameter '{p}'. "
        f"Use an integer ID or one of: {', '.join(sorted(_NAME_TO_ID.keys()))}"
    )


# ── HTTP helpers ─────────────────────────────────────────────────────────


def _get_json(url: str, timeout_s: float = 10.0) -> dict:
    """Fetch JSON content from *url* with error handling.

    Parameters
    ----------
    url : str
        The URL to fetch.
    timeout_s : float
        Request timeout in seconds.

    Returns
    -------
    dict
        Parsed JSON response.

    Raises
    ------
    RuntimeError
        If the request fails.
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError(
            "requests library required for hydrology dataset. "
            "Install with: pip install requests"
        )

    try:
        response = requests.get(url, timeout=timeout_s)
        response.raise_for_status()
        return response.json()
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


# ── Station list + data fetchers ─────────────────────────────────────────


def _fetch_station_list(
    base_url: str, version: str, param_id: int, timeout_s: float
) -> List[Dict[str, Any]]:
    """Fetch the station list for a given parameter.

    Returns
    -------
    list[dict]
        Each dict has keys: ``key`` (station ID as str), ``name``,
        ``latitude``, ``longitude``, ``active``, ``catchmentName``,
        ``catchmentNumber``, ``catchmentSize``, and other metadata.
    """
    url = f"{base_url}/version/{version}/parameter/{param_id}.json"
    data = _get_json(url, timeout_s=timeout_s)
    return data.get("station", [])


def _fetch_latest_day(
    base_url: str, version: str, param_id: int, station_key: str, timeout_s: float
) -> Dict[str, Any]:
    """Fetch latest-day data for a single station and parameter.

    Returns
    -------
    dict
        The parsed JSON response with keys ``value``, ``parameter``,
        ``station``, ``period``, ``position``, etc.  The ``value`` array
        contains dicts with ``date``, ``value``, and ``quality``.
    """
    url = (
        f"{base_url}/version/{version}/parameter/{param_id}"
        f"/station/{station_key}/period/latest-day/data.json"
    )
    return _get_json(url, timeout_s=timeout_s)


# ── Dataset definition ───────────────────────────────────────────────────


class HydrologyDatasetArgs(DatasetBaseArgs):
    """Arguments for the hydrology dataset.

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
    active_only : bool
        If ``True``, only include stations marked as active.
    base_url : str
        SMHI HydroObs API base URL (override for testing).
    version : str
        API version path element.
    parameters : list
        SMHI parameter IDs to fetch (default: daily discharge and water level).
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
    active_only: bool = Field(
        True,
        description="Only include stations marked as currently active",
    )
    base_url: str = Field(
        "https://opendata-download-hydroobs.smhi.se/api",
        description="SMHI HydroObs API base URL",
    )
    version: str = Field("latest", description="API version path element")
    parameters: List = Field(
        default=[1, 3],
        description=(
            "Hydrology parameters to fetch.  Each element can be an integer "
            "SMHI parameter ID or a name/alias string such as "
            "'discharge_daily', 'discharge', 'flow', 'water_level', 'level', "
            "'water_temperature', 'ice_thickness', 'snow_density'."
        ),
    )
    field_name_style: Literal["dtcc", "smhi"] = Field(
        "dtcc",
        description="Field naming style: dtcc (English) or smhi (Swedish API name)",
    )


class HydrologyDataset(DatasetDescriptor):
    """Hydrological observations from SMHI HydroObs (latest-day snapshot).

    This dataset fetches the most recent daily observations from the SMHI
    open-data hydrological observations API.  It first retrieves the
    station list per parameter (filtering by bounding box and active
    status), then fetches ``latest-day`` JSON data for each matching
    station, merges fields, and returns a ``SensorCollection``.

    Default parameters: daily discharge (1) and water level (3).

    Example
    -------
    >>> import dtcc_core.datasets as datasets
    >>> sensors = datasets.hydrology(
    ...     bounds=(674000, 6580000, 676000, 6582000),
    ...     parameters=[1, 3],
    ... )
    >>> print(f"Found {len(sensors.stations())} stations")
    """

    name = "hydrology"
    description = (
        "Hydrological observations from SMHI HydroObs (latest-day snapshot) "
        "as a SensorCollection with one Field per parameter per station "
        "within the specified bounds."
    )
    ArgsModel = HydrologyDatasetArgs

    def build(self, args: HydrologyDatasetArgs):
        """Build the hydrology dataset.

        Parameters
        ----------
        args : HydrologyDatasetArgs
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
            f"Hydrology dataset: fetching {len(param_ids)} parameter(s) "
            f"for bbox ({lon_min:.4f}, {lat_min:.4f}) – "
            f"({lon_max:.4f}, {lat_max:.4f}) WGS84"
        )

        # ── Per-parameter: list stations → filter → fetch data ──────
        # station_key -> { "station_name", "lat", "lon",
        #                   "catchment_name", "catchment_number",
        #                   "catchment_size",
        #                   "fields": { field_name: (value, unit, quality) } }
        station_map: Dict[str, Dict[str, Any]] = {}
        param_meta: Dict[int, Dict[str, Any]] = {}

        for pid in param_ids:
            info(f"  Fetching station list for parameter {pid} …")

            try:
                stations = _fetch_station_list(
                    args.base_url, args.version, pid, args.timeout_s
                )
            except RuntimeError as exc:
                info(
                    f"  Warning: failed to fetch station list for parameter {pid}: {exc}"
                )
                continue

            # Filter by active status and bounding box
            candidates = []
            for st in stations:
                if args.active_only and not st.get("active", False):
                    continue
                lat = st.get("latitude")
                lon = st.get("longitude")
                if lat is None or lon is None:
                    continue
                if not (lon_min <= lon <= lon_max and lat_min <= lat <= lat_max):
                    continue
                candidates.append(st)

            info(f"  {len(candidates)} station(s) in bbox for parameter {pid}")

            # Fetch latest-day data for each candidate
            for st in candidates:
                skey = str(st["key"])
                info(f"    Fetching station {skey} ({st.get('name', '?')}) …")

                try:
                    data = _fetch_latest_day(
                        args.base_url, args.version, pid, skey, args.timeout_s
                    )
                except RuntimeError as exc:
                    info(f"    Warning: skipping station {skey}: {exc}")
                    continue

                # Extract the most recent value
                values = data.get("value", [])
                if not values:
                    value = float("nan")
                    quality = ""
                else:
                    # Take the last (most recent) entry
                    latest = values[-1]
                    raw_val = latest.get("value")
                    try:
                        value = float(raw_val)
                    except (ValueError, TypeError):
                        value = float("nan")
                    quality = latest.get("quality", "")

                # Extract parameter metadata
                param_info = data.get("parameter", {})
                param_name_smhi = param_info.get("name", f"parameter_{pid}")
                unit = param_info.get("unit", "")

                param_meta[pid] = {
                    "parameter_name": param_name_smhi,
                    "unit": unit,
                }

                # Determine field name
                if args.field_name_style == "dtcc":
                    field_name = PARAMETER_NAMES.get(pid, f"parameter_{pid}")
                else:
                    field_name = param_name_smhi

                # Merge into station map
                if skey not in station_map:
                    station_map[skey] = {
                        "station_name": st.get("name", ""),
                        "lat": st["latitude"],
                        "lon": st["longitude"],
                        "catchment_name": st.get("catchmentName", ""),
                        "catchment_number": st.get("catchmentNumber", ""),
                        "catchment_size": st.get("catchmentSize", ""),
                        "fields": {},
                    }

                station_map[skey]["fields"][field_name] = (value, unit, quality)

        info(f"  {len(station_map)} station(s) within bounds (all parameters)")

        # ── Optionally drop stations with all-missing values ─────────
        if args.drop_missing:
            station_map = {
                skey: data
                for skey, data in station_map.items()
                if any(not np.isnan(v) for v, _u, _q in data["fields"].values())
            }
            info(f"  {len(station_map)} station(s) after dropping all-missing")

        # ── Coordinate transform for output ──────────────────────────
        need_reproject = args.crs.upper() not in ("CRS84", "EPSG:4326", "WGS84")

        if need_reproject and station_map:
            skeys = list(station_map.keys())
            lonlat = np.array(
                [[station_map[s]["lon"], station_map[s]["lat"], 0.0] for s in skeys]
            )
            projected = reproject_array(lonlat, "EPSG:4326", args.crs)
            for i, skey in enumerate(skeys):
                station_map[skey]["x"] = projected[i, 0]
                station_map[skey]["y"] = projected[i, 1]
        else:
            for skey, data in station_map.items():
                data["x"] = data["lon"]
                data["y"] = data["lat"]

        # ── Build SensorCollection ───────────────────────────────────
        sensor_collection = SensorCollection()
        sensor_collection.attributes = {
            "dtcc_type": "sensor_collection",
            "source": "smhi_hydroobs",
            "dataset": "hydrology",
            "crs": args.crs,
            "bounds": (
                f"({bounds_tuple[0]}, {bounds_tuple[1]}, "
                f"{bounds_tuple[2]}, {bounds_tuple[3]})"
            ),
            "retrieval_time": datetime.now(timezone.utc).isoformat(),
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
        sorted_skeys = sorted(station_map.keys())
        if len(sorted_skeys) > args.max_stations:
            sorted_skeys = sorted_skeys[: args.max_stations]

        for skey in sorted_skeys:
            data = station_map[skey]

            station = Object()
            station.attributes = {
                "station_id": skey,
                "station_name": data["station_name"],
                "catchment_name": data["catchment_name"],
                "catchment_number": data["catchment_number"],
                "catchment_size": data["catchment_size"],
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

            point = Point(x=data["x"], y=data["y"])
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
            f"Hydrology dataset: returning {len(sensor_collection.stations())} "
            f"station(s) with {len(parameter_fields)} field(s)"
        )

        if args.format == "pb":
            return sensor_collection.to_proto().SerializeToString()
        return sensor_collection


__all__ = ["HydrologyDatasetArgs", "HydrologyDataset"]
