# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

"""Air quality dataset from datavardluft.smhi.se.

This module provides a dataset interface to fetch air quality measurements
from the Swedish SMHI air quality API (datavardluft.smhi.se). The dataset
returns a SensorCollection containing measurement stations as child objects.

Each station has:
- Point geometry for location
- Field with the latest measured value
- Metadata attributes (ID, name, timestamp, etc.)
"""

from typing import Optional, Literal, Tuple, List, Dict, Any
from pydantic import Field
import numpy as np
import json
from datetime import datetime, timezone

from .dataset import DatasetBaseArgs, DatasetDescriptor, DatasetUpstreamError
from .geospatial import bounds_to_wgs84, is_wgs84_crs
from .providers import provider_entry
from ..model.object import Object, SensorCollection
from ..model.geometry import Point
from ..model.values import Field as DtccField
from ..common import info
from ..common.progress import ProgressTracker, report_progress
from ..reproject.reproject import reproject_array


STALE_VALUE_DAYS = 7

# Fixture-verified phenomenon aliases. Other labels are resolved through the
# provider's /phenomena endpoint so stale hard-coded IDs do not silently win.
_FIXTURE_VERIFIED_PHENOMENA: Dict[str, Dict[str, str]] = {
    "NO2": {"id": "8", "label": "NO2"},
    "PM10": {"id": "5", "label": "PM10"},
    "O3": {"id": "7", "label": "O3"},
}


# API helper functions (kept private to this module)


def _get_json(
    url: str, params: Dict[str, Any] = None, timeout_s: float = 10.0
) -> Dict[str, Any]:
    """Fetch JSON from URL with error handling.

    Parameters
    ----------
    url : str
        The URL to fetch
    params : dict, optional
        Query parameters
    timeout_s : float
        Request timeout in seconds

    Returns
    -------
    dict
        Parsed JSON response

    Raises
    ------
    RuntimeError
        If request fails or response is not JSON
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError(
            "requests library required for air quality dataset. "
            "Install with: pip install requests"
        )

    try:
        response = requests.get(url, params=params, timeout=timeout_s)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        raise DatasetDescriptor.build_upstream_error(
            "air_quality", "fetch_json", url, e
        ) from e
    except json.JSONDecodeError as e:
        raise DatasetUpstreamError(
            dataset="air_quality",
            operation="decode_json",
            target=url,
            failure_class="invalid_json",
            message=f"air_quality decode_json failed for {url}: {e}",
        ) from e


def _resolve_phenomenon(
    base_url: str,
    phenomenon_str: str,
    timeout_s: float,
    strict_live: bool = False,
    upstream_errors: Optional[List[DatasetUpstreamError]] = None,
) -> Tuple[str, Dict[str, Any]]:
    """Resolve phenomenon name/ID to phenomenon ID.

    The API uses numeric phenomenon IDs. This function accepts either numeric
    IDs or labels. Only fixture-verified labels are resolved through a static
    fallback; all other labels are resolved through the provider endpoint.

    Parameters
    ----------
    base_url : str
        Base API URL
    phenomenon_str : str
        Phenomenon name or ID (e.g., "NO2", "PM10", "1", "5")
    timeout_s : float
        Request timeout

    Returns
    -------
    tuple
        ``(phenomenon_id, metadata)``.
    """
    requested = str(phenomenon_str).strip()
    if not requested:
        raise ValueError("Air quality phenomenon must be a non-empty string or ID.")

    metadata: Dict[str, Any] = {
        "requested": requested,
        "phenomenon_id": "",
        "label": requested,
        "resolution_source": "",
        "static_mapping_status": "not_used",
        "verified_static_aliases": sorted(_FIXTURE_VERIFIED_PHENOMENA),
    }

    # If it's already numeric, return it.
    if requested.isdigit():
        metadata.update(
            {
                "phenomenon_id": requested,
                "resolution_source": "numeric_id",
                "static_mapping_status": "not_applicable",
            }
        )
        info(f"Using numeric phenomenon ID: {requested}")
        return requested, metadata

    phenomenon_key = requested.upper()
    if phenomenon_key in _FIXTURE_VERIFIED_PHENOMENA:
        match = _FIXTURE_VERIFIED_PHENOMENA[phenomenon_key]
        phenomenon_id = match["id"]
        metadata.update(
            {
                "phenomenon_id": phenomenon_id,
                "label": match["label"],
                "resolution_source": "fixture_verified_alias",
                "static_mapping_status": "fixture_verified",
            }
        )
        info(f"Mapped phenomenon '{requested}' to fixture-verified ID {phenomenon_id}")
        return phenomenon_id, metadata

    # Try to fetch from API
    try:
        url = f"{base_url}/phenomena"
        data = _get_json(url, timeout_s=timeout_s)
        if not isinstance(data, list):
            raise DatasetUpstreamError(
                dataset="air_quality",
                operation="resolve_phenomenon",
                target=url,
                failure_class="invalid_payload",
                message=(
                    "air_quality resolve_phenomenon expected /phenomena to return "
                    f"a list, got {type(data).__name__}"
                ),
            )
        for phen in data:
            if not isinstance(phen, dict):
                continue
            if str(phen.get("label", "")).strip().upper() == phenomenon_key:
                if "id" not in phen:
                    raise ValueError(
                        "Air quality phenomenon lookup matched "
                        f"'{requested}' but the provider payload has no 'id'."
                    )
                phenomenon_id = str(phen["id"])
                label = str(phen.get("label", requested))
                metadata.update(
                    {
                        "phenomenon_id": phenomenon_id,
                        "label": label,
                        "resolution_source": "provider_lookup",
                        "static_mapping_status": "not_applicable",
                    }
                )
                return phenomenon_id, metadata
    except DatasetUpstreamError as exc:
        if upstream_errors is not None:
            upstream_errors.append(exc)
        if strict_live or exc.failure_class == "invalid_payload":
            raise
        metadata.update(
            {
                "phenomenon_id": requested,
                "resolution_source": "raw_unresolved_due_to_lookup_failure",
                "static_mapping_status": "unverified_raw",
            }
        )
    else:
        raise ValueError(f"Unknown air quality phenomenon: {requested}")

    # Upstream unreachable (non-strict mode): treat input as a raw ID
    return requested, metadata


def _resolve_phenomenon_id(
    base_url: str,
    phenomenon_str: str,
    timeout_s: float,
    strict_live: bool = False,
    upstream_errors: Optional[List[DatasetUpstreamError]] = None,
) -> str:
    """Resolve phenomenon name/ID to phenomenon ID."""
    phenomenon_id, _metadata = _resolve_phenomenon(
        base_url,
        phenomenon_str,
        timeout_s,
        strict_live=strict_live,
        upstream_errors=upstream_errors,
    )
    return phenomenon_id


def _parse_measurement_datetime(timestamp: str) -> Optional[datetime]:
    """Parse a provider measurement timestamp as UTC when possible."""
    if not timestamp:
        return None
    try:
        parsed = datetime.fromisoformat(str(timestamp).replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _measurement_staleness(
    timestamp: str,
    reference_time: datetime,
    *,
    stale_after_days: int = STALE_VALUE_DAYS,
) -> Tuple[bool, Optional[float]]:
    """Return whether a timestamp is stale and its age in days."""
    parsed = _parse_measurement_datetime(timestamp)
    if parsed is None:
        return False, None
    age = reference_time.astimezone(timezone.utc) - parsed
    age_days = age.total_seconds() / 86400.0
    return age_days > stale_after_days, age_days


def _fetch_stations(
    base_url: str,
    bounds: Tuple[float, float, float, float],
    crs: str,
    timeout_s: float,
    strict_live: bool = False,
    upstream_errors: Optional[List[DatasetUpstreamError]] = None,
) -> List[Dict[str, Any]]:
    """Fetch stations within bounding box.

    Parameters
    ----------
    base_url : str
        Base API URL
    bounds : tuple
        (xmin, ymin, xmax, ymax) bounding box
    crs : str
        Coordinate reference system
    timeout_s : float
        Request timeout

    Returns
    -------
    list
        List of station dictionaries
    """
    url = f"{base_url}/stations"
    info(f"Querying air quality API for stations...")

    # Transform bounds to WGS84 if necessary
    wgs84_bounds = bounds_to_wgs84(bounds, crs)

    # API expects bbox as: xmin,ymin,xmax,ymax in WGS84/CRS84
    bbox_str = (
        f"{wgs84_bounds[0]},{wgs84_bounds[1]},{wgs84_bounds[2]},{wgs84_bounds[3]}"
    )

    params = {
        "bbox": bbox_str,
        "crs": "CRS84",  # API always uses CRS84 for bbox
    }

    try:
        data = _get_json(url, params=params, timeout_s=timeout_s)
        return data if isinstance(data, list) else []
    except DatasetUpstreamError as e:
        if upstream_errors is not None:
            upstream_errors.append(e)
        if strict_live:
            raise
        info(f"Warning: Failed to fetch stations: {e}")
        return []


def _fetch_timeseries_for_station(
    base_url: str,
    station_id: str,
    phenomenon_id: str,
    timeout_s: float,
    strict_live: bool = False,
    upstream_errors: Optional[List[DatasetUpstreamError]] = None,
) -> List[Dict[str, Any]]:
    """Fetch timeseries metadata for a station and phenomenon.

    The API structure is:
    1. GET /stations/{id} returns station with embedded timeseries IDs
    2. GET /timeseries/{id} returns timeseries metadata including lastValue

    Parameters
    ----------
    base_url : str
        Base API URL
    station_id : str
        Station ID
    phenomenon_id : str
        Phenomenon ID (e.g., "8" for NO2)
    timeout_s : float
        Request timeout

    Returns
    -------
    list
        List of timeseries dictionaries with lastValue
    """
    try:
        # Step 1: Get station details which include timeseries IDs
        station_url = f"{base_url}/stations/{station_id}"
        station_data = _get_json(station_url, timeout_s=timeout_s)

        # Step 2: Extract timeseries from station properties
        timeseries_dict = station_data.get("properties", {}).get("timeseries", {})

        # Step 3: Filter by phenomenon and fetch full timeseries metadata
        result = []
        for ts_id, ts_info in timeseries_dict.items():
            # Check if this timeseries matches our phenomenon
            phenomenon = ts_info.get("phenomenon", {})
            if str(phenomenon.get("id")) == str(phenomenon_id):
                # Fetch the full timeseries metadata (includes lastValue)
                ts_url = f"{base_url}/timeseries/{ts_id}"
                try:
                    ts_data = _get_json(ts_url, timeout_s=timeout_s)
                    result.append(ts_data)
                except DatasetUpstreamError as exc:
                    if upstream_errors is not None:
                        upstream_errors.append(exc)
                    if strict_live:
                        raise
                    # Skip this timeseries if we can't fetch it
                    pass

        return result

    except DatasetUpstreamError as exc:
        if upstream_errors is not None:
            upstream_errors.append(exc)
        if strict_live:
            raise
        return []


def _extract_latest_value(
    timeseries_dict: Dict[str, Any],
) -> Optional[Tuple[float, str, str]]:
    """Extract latest value from timeseries metadata.

    The API returns lastValue as a dict with 'timestamp' (milliseconds since epoch)
    and 'value' (float).

    Parameters
    ----------
    timeseries_dict : dict
        Timeseries metadata dictionary from /timeseries/{id}

    Returns
    -------
    tuple or None
        (value, timestamp_iso, unit) if available, else None
    """
    # Check for lastValue in metadata
    last_value_dict = timeseries_dict.get("lastValue")
    if last_value_dict and isinstance(last_value_dict, dict):
        try:
            value = float(last_value_dict.get("value"))
            timestamp_ms = last_value_dict.get("timestamp")

            # Convert timestamp from milliseconds to UTC ISO format.
            timestamp_iso = ""
            if timestamp_ms is not None:
                dt = datetime.fromtimestamp(float(timestamp_ms) / 1000.0, timezone.utc)
                timestamp_iso = dt.isoformat()
                info(
                    f"Extracted value {value} with timestamp {timestamp_iso} "
                    f"(raw: {timestamp_ms} ms)"
                )

            unit = timeseries_dict.get("uom", "")
            return (value, timestamp_iso, unit)
        except (ValueError, TypeError, KeyError):
            pass

    return None


def _fallback_get_latest_from_getData(
    base_url: str,
    timeseries_id: str,
    timeout_s: float,
    strict_live: bool = False,
    upstream_errors: Optional[List[DatasetUpstreamError]] = None,
) -> Optional[Tuple[float, str, str]]:
    """Fallback: fetch latest value using getData endpoint.

    Parameters
    ----------
    base_url : str
        Base API URL
    timeseries_id : str
        Timeseries ID
    timeout_s : float
        Request timeout

    Returns
    -------
    tuple or None
        (value, timestamp, unit) if available, else None
    """
    url = f"{base_url}/timeseries/{timeseries_id}/getData"

    # Request recent data - try last 7 days to get more recent measurements
    params = {"timespan": "P7D"}  # Last 7 days

    try:
        info(
            "Fallback: Fetching recent data from getData endpoint for "
            f"timeseries {timeseries_id}"
        )
        data = _get_json(url, params=params, timeout_s=timeout_s)
        values = data.get("values", [])

        if values:
            # Take the last value
            last = values[-1]
            value = float(last["value"])
            timestamp = last.get("timestamp", "")
            unit = data.get("uom", "")
            info(f"Fallback: Got value {value} with timestamp {timestamp} from getData")
            return (value, timestamp, unit)
        else:
            info(f"Fallback: No values returned from getData endpoint")
    except DatasetUpstreamError as e:
        if upstream_errors is not None:
            upstream_errors.append(e)
        if strict_live:
            raise
        info(f"Fallback: Failed to get data from getData endpoint: {e}")

    return None


class AirQualityDatasetArgs(DatasetBaseArgs):
    """Arguments for air quality dataset.

    Attributes
    ----------
    bounds : tuple
        Geographic bounds (xmin, ymin, xmax, ymax)
    phenomenon : str
        Phenomenon to fetch (e.g., "NO2", "PM10", "PM2.5", or numeric ID)
    crs : str
        Coordinate reference system for bounds and station coordinates
    format : str, optional
        Output format ("pb" for protobuf bytes, None for Python object)
    timeout_s : float
        HTTP request timeout in seconds
    max_stations : int
        Maximum number of stations to process (safety limit)
    drop_missing : bool
        Skip stations with missing values
    base_url : str
        API base URL (for testing/mocking)
    """

    phenomenon: str = Field(
        "NO2", description="Phenomenon name or ID (NO2, PM10, PM2.5, etc.)"
    )
    crs: str = Field("EPSG:3006", description="Coordinate reference system")
    format: Optional[Literal["pb"]] = Field(
        None, description="Output format (pb for protobuf)"
    )
    timeout_s: float = Field(10.0, description="HTTP timeout in seconds", gt=0)
    max_stations: int = Field(250, description="Maximum stations to process", gt=0)
    drop_missing: bool = Field(True, description="Skip stations with missing values")
    base_url: str = Field(
        "https://datavardluft.smhi.se/52North/api", description="API base URL"
    )


class AirQualityDataset(DatasetDescriptor):
    """Air quality sensor data from datavardluft.smhi.se.

    This dataset fetches snapshot air quality measurements from the Swedish SMHI
    air quality monitoring network. It returns a SensorCollection where each
    station is represented as an Object with:

    - Point geometry for station location
    - Field with the latest measured value
    - Metadata attributes (station ID, name, timestamp, etc.)

    The dataset queries stations within the specified geographic bounds and
    retrieves the most recent measurement for the selected phenomenon (pollutant).

    Fixture-verified aliases include NO2, O3, and PM10. Other phenomena can
    be requested with a numeric ID or resolved from the provider /phenomena
    endpoint when available.

    Example
    -------
    >>> import dtcc_core.datasets as datasets
    >>> # Fetch NO2 measurements in Stockholm area (SWEREF99 TM)
    >>> sensors = datasets.air_quality(
    ...     bounds=(674000, 6580000, 676000, 6582000),
    ...     phenomenon="NO2"
    ... )
    >>> print(f"Found {len(sensors.stations())} stations")

    >>> # Get data as arrays
    >>> points, values = sensors.to_arrays()
    >>> print(f"NO2 range: {values.min():.1f} - {values.max():.1f}")
    """

    name = "air_quality"
    title = "Air Quality Observations"
    description = (
        "Air-quality station measurements from the SMHI datavardluft API "
        "as a SensorCollection with latest available readings for one "
        "requested phenomenon within the specified bounds."
    )
    ArgsModel = AirQualityDatasetArgs
    data_category = "raw"
    result_kind = "sensor_collection"
    python_return_type = "dtcc_core.model.SensorCollection"
    provider = [provider_entry("smhi", role="source_provider")]
    source = [
        {
            "name": "SMHI datavardluft air-quality API",
            "role": "source_provider",
            "service": "datavardluft",
            "base_url": "https://datavardluft.smhi.se/52North/api",
            "endpoint_patterns": [
                "/phenomena",
                "/stations?bbox={bbox}&crs=CRS84",
                "/stations/{station_id}",
                "/timeseries/{timeseries_id}",
                "/timeseries/{timeseries_id}/getData?timespan=P7D",
            ],
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "Requires review: verify SMHI datavardluft source terms before "
        "redistribution."
    )
    collection_period = (
        "Latest available station snapshot. The provider /timeseries metadata "
        "may expose lastValue; stale or missing lastValue readings trigger a "
        "P7D getData lookup when a timeseries ID is available."
    )
    default_crs = "EPSG:3006"
    data_types = [
        "sensor_collection",
        "air_quality_observations",
        "points",
        "time_series_snapshot",
    ]
    geographic_coverage = "Sweden, constrained by station coverage and requested bounds"
    update_frequency = "latest available observation snapshot"
    processing_steps = [
        (
            "Resolve the requested air-quality phenomenon from a numeric ID, "
            "fixture-verified alias, or provider /phenomena lookup"
        ),
        "Transform requested bounds to WGS84/CRS84 for provider station filtering",
        (
            "Fetch station metadata and matching timeseries IDs for the "
            "requested phenomenon"
        ),
        "Read latest values from /timeseries/{id} lastValue metadata",
        (
            "When lastValue is missing or stale, check "
            "/timeseries/{id}/getData over the recent P7D window"
        ),
        (
            "Reproject station points to the requested output CRS and preserve "
            "station/operator metadata"
        ),
        (
            "Record phenomenon metadata, value timestamps, units, source "
            "endpoint, stale counts, fallback counts, and upstream errors"
        ),
        "Return a SensorCollection or serialize protobuf bytes when format='pb'",
    ]
    presentation_headline = "Latest SMHI Air-Quality Stations"
    presentation_summary = (
        "Point observations from SMHI air-quality monitoring stations, filtered "
        "to the requested bounds for one pollutant or phenomenon."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each point is a monitoring station. The attached field is the "
                "latest available value for the requested phenomenon, such as "
                "NO2, PM10, or a provider phenomenon ID."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Measurements are station observations, not an interpolated air "
                "quality surface. Station attributes preserve operator, unit, "
                "phenomenon ID, timestamp, value source, and staleness metadata."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Station coverage and update cadence vary by pollutant and "
                "operator. Empty results may mean no station or no current "
                "measurement in the requested area."
            ),
        },
    ]
    key_points = [
        "Fixture-verified aliases are NO2, PM10, and O3",
        (
            "Unverified phenomenon labels are resolved through the provider "
            "/phenomena endpoint"
        ),
        (
            "Station attributes preserve operator, unit, timestamp, timeseries "
            "ID, value source, and stale-value flags"
        ),
        (
            "Stale lastValue readings trigger a recent getData lookup before "
            "the stale value is retained"
        ),
        (
            "The dataset is table-visible as station-level context and upstream "
            "lineage for derived air-quality fields"
        ),
    ]
    presentation_legend = {
        "title": "Air-quality station field",
        "entries": [
            {
                "label": "NO2",
                "meaning": "Nitrogen dioxide, fixture-verified phenomenon ID 8",
            },
            {
                "label": "PM10",
                "meaning": "Particulate matter PM10, fixture-verified ID 5",
            },
            {
                "label": "O3",
                "meaning": "Ozone, fixture-verified phenomenon ID 7",
            },
            {
                "label": "value_source",
                "meaning": "Station attribute identifying lastValue or getData",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "points",
        "default_crs": "EPSG:3006",
        "table_role": "station_context",
        "quality_attribute": "is_stale",
    }
    presentation_warnings = [
        (
            "Source and license terms require review before redistributing "
            "generated packages."
        ),
        (
            "Sparse monitoring stations should not be interpreted as continuous "
            "area coverage."
        ),
        (
            "Stale readings are retained only with explicit station attributes; "
            "inspect is_stale and value_age_days before analysis."
        ),
        (
            "Interpolation from station points can misrepresent street-canyon, "
            "source-proximity, and meteorological effects."
        ),
        (
            "A non-strict run may return partial results when station or "
            "timeseries requests fail; inspect result health attributes."
        ),
    ]
    presentation_limitations = [
        "Only one phenomenon is fetched per dataset call.",
        "Station and phenomenon availability vary by monitoring network and operator.",
        (
            "The provider's lastValue metadata can be stale or absent; getData "
            "fallback checks only the recent P7D window."
        ),
        (
            "Quality assurance beyond unit/timestamp/source preservation "
            "remains a domain review task."
        ),
    ]

    def build(self, args: AirQualityDatasetArgs):
        """Build the air quality dataset.

        Parameters
        ----------
        args : AirQualityDatasetArgs
            Dataset arguments

        Returns
        -------
        SensorCollection or bytes
            SensorCollection object, or protobuf bytes if format="pb"
        """
        progress_phases = {
            "resolve_api": 0.05,
            "fetch_stations": 0.15,
            "fetch_measurements": 0.65,
            "build_collection": 0.10,
            "export": 0.05,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            bounds_tuple = (bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
            upstream_errors: List[DatasetUpstreamError] = []

            with progress.phase(
                "resolve_api", f"Resolving phenomenon '{args.phenomenon}'..."
            ):
                phenomenon_id, phenomenon_metadata = _resolve_phenomenon(
                    args.base_url,
                    args.phenomenon,
                    args.timeout_s,
                    strict_live=args.strict_live,
                    upstream_errors=upstream_errors,
                )
                info(f"Resolved phenomenon {args.phenomenon} to ID {phenomenon_id}")

            retrieval_time = datetime.now(timezone.utc)
            with progress.phase("fetch_stations", "Fetching stations within bounds..."):
                stations_data = _fetch_stations(
                    args.base_url,
                    bounds_tuple,
                    args.crs,
                    args.timeout_s,
                    strict_live=args.strict_live,
                    upstream_errors=upstream_errors,
                )
                report_progress(
                    percent=100,
                    message=f"Found {len(stations_data)} stations",
                )

            with progress.phase(
                "fetch_measurements", "Fetching measurements from stations..."
            ):
                sensor_collection = SensorCollection()
                sensor_collection.attributes = {
                    "dtcc_type": "sensor_collection",
                    "source": "datavardluft.smhi.se",
                    "phenomenon": args.phenomenon,
                    "phenomenon_id": phenomenon_id,
                    "phenomenon_metadata": phenomenon_metadata,
                    "crs": args.crs,
                    "bounds": f"{bounds_tuple}",
                    "retrieval_time": retrieval_time.isoformat(),
                    "total_stations_found": len(stations_data),
                    "stale_after_days": STALE_VALUE_DAYS,
                }

                stations_used = 0
                stations_skipped_upstream = 0
                stations_skipped_no_coords = 0
                stations_skipped_no_timeseries = 0
                stations_skipped_no_value = 0
                stale_value_count = 0
                fallback_getdata_attempts = 0
                fallback_getdata_successes = 0

                stations_to_process = stations_data[: args.max_stations]
                total_stations = len(stations_to_process)
                need_reproject = not is_wgs84_crs(args.crs)

                for i, station_data in enumerate(stations_to_process):
                    report_progress(
                        current=i,
                        total=total_stations,
                        message=f"Processing station {i + 1}/{total_stations}...",
                    )
                    station_id = str(station_data.get("id", ""))

                    # Extract coordinates
                    geometry = station_data.get("geometry", {})
                    coordinates = geometry.get("coordinates", [])

                    if len(coordinates) < 2:
                        stations_skipped_no_coords += 1
                        continue

                    x, y = coordinates[0], coordinates[1]
                    z = coordinates[2] if len(coordinates) > 2 else 0.0

                    if need_reproject:
                        projected = reproject_array(
                            np.array([[x, y, z]], dtype=float), "EPSG:4326", args.crs
                        )
                        x_out = projected[0, 0]
                        y_out = projected[0, 1]
                        z_out = projected[0, 2] if projected.shape[1] > 2 else z
                    else:
                        x_out, y_out, z_out = x, y, z

                    if not self.point_within_bounds(x_out, y_out, bounds):
                        continue

                    station_error_count_before = len(upstream_errors)

                    # Fetch timeseries for this station
                    timeseries_list = _fetch_timeseries_for_station(
                        args.base_url,
                        station_id,
                        phenomenon_id,
                        args.timeout_s,
                        strict_live=args.strict_live,
                        upstream_errors=upstream_errors,
                    )

                    if not timeseries_list:
                        if len(upstream_errors) > station_error_count_before:
                            stations_skipped_upstream += 1
                        else:
                            stations_skipped_no_timeseries += 1
                        continue

                    # Try to get latest value
                    value = None
                    timestamp = ""
                    unit = ""
                    timeseries_id = ""
                    value_source = ""
                    fallback_reason = ""
                    fallback_status = ""
                    metadata_last_value_timestamp = ""

                    for ts in timeseries_list:
                        ts_id = str(ts.get("id", ""))
                        if ts.get("uom") and not unit:
                            unit = str(ts.get("uom", ""))

                        # First try to extract from metadata
                        result = _extract_latest_value(ts)
                        if result:
                            value, timestamp, unit = result
                            timeseries_id = ts_id
                            value_source = "metadata.lastValue"
                            metadata_last_value_timestamp = timestamp

                            is_stale, _age_days = _measurement_staleness(
                                timestamp, retrieval_time
                            )
                            if is_stale and ts_id:
                                fallback_reason = "stale_lastValue"
                                fallback_getdata_attempts += 1
                                fallback_error_count_before = len(upstream_errors)
                                fallback = _fallback_get_latest_from_getData(
                                    args.base_url,
                                    ts_id,
                                    args.timeout_s,
                                    strict_live=args.strict_live,
                                    upstream_errors=upstream_errors,
                                )
                                if fallback:
                                    (
                                        fallback_value,
                                        fallback_timestamp,
                                        fallback_unit,
                                    ) = fallback
                                    metadata_dt = _parse_measurement_datetime(timestamp)
                                    fallback_dt = _parse_measurement_datetime(
                                        fallback_timestamp
                                    )
                                    if (
                                        metadata_dt is None
                                        or fallback_dt is None
                                        or fallback_dt >= metadata_dt
                                    ):
                                        value = fallback_value
                                        timestamp = fallback_timestamp
                                        unit = fallback_unit or unit
                                        value_source = "getData"
                                        fallback_status = "used"
                                        fallback_getdata_successes += 1
                                    else:
                                        fallback_status = "older_than_lastValue"
                                elif len(upstream_errors) > fallback_error_count_before:
                                    fallback_status = "upstream_error"
                                else:
                                    fallback_status = "no_recent_values"
                            break

                        # Fallback: fetch from getData
                        if ts_id:
                            timeseries_id = ts_id
                            fallback_reason = "missing_lastValue"
                            fallback_getdata_attempts += 1
                            fallback_error_count_before = len(upstream_errors)
                            result = _fallback_get_latest_from_getData(
                                args.base_url,
                                ts_id,
                                args.timeout_s,
                                strict_live=args.strict_live,
                                upstream_errors=upstream_errors,
                            )
                            if result:
                                value, timestamp, unit = result
                                value_source = "getData"
                                fallback_status = "used"
                                fallback_getdata_successes += 1
                                break
                            if len(upstream_errors) > fallback_error_count_before:
                                fallback_status = "upstream_error"
                            else:
                                fallback_status = "no_recent_values"

                    if value is None:
                        if args.drop_missing:
                            if len(upstream_errors) > station_error_count_before:
                                stations_skipped_upstream += 1
                            else:
                                stations_skipped_no_value += 1
                            continue
                        else:
                            value = np.nan
                            value_source = "missing"

                    is_stale, value_age_days = _measurement_staleness(
                        timestamp, retrieval_time
                    )
                    if is_stale:
                        stale_value_count += 1
                        rounded_age = int(value_age_days) if value_age_days else 0
                        info(
                            f"Warning: Station {station_id} has stale air-quality data "
                            f"(age: {rounded_age} days)"
                        )

                    # Create station object
                    station = Object()
                    station_props = station_data.get("properties", {})
                    station.attributes = {
                        "station_id": station_id,
                        "station_name": station_props.get("label", ""),
                        "operator": station_props.get("operator", ""),
                        "phenomenon": args.phenomenon,
                        "phenomenon_id": phenomenon_id,
                        "unit": unit,
                        "timestamp": timestamp,
                        "value": value,
                        "timeseries_id": timeseries_id,
                        "value_source": value_source,
                        "fallback_reason": fallback_reason,
                        "fallback_status": fallback_status,
                        "metadata_last_value_timestamp": metadata_last_value_timestamp,
                        "is_stale": is_stale,
                        "value_age_days": (
                            round(value_age_days, 3)
                            if value_age_days is not None
                            else None
                        ),
                    }

                    # Create Point geometry
                    point = Point(x=x_out, y=y_out, z=z_out)

                    # Create Field with value
                    field = DtccField()
                    field.name = args.phenomenon
                    field.unit = unit
                    field.dim = 1
                    field.values = np.array([value], dtype=np.float32)

                    # Attach field to geometry
                    point.fields = [field]

                    # Add geometry to station
                    station.geometry["location"] = point

                    # Add station to collection
                    sensor_collection.add_station(station)
                    stations_used += 1

                report_progress(
                    percent=100,
                    message=f"Processed {stations_used} stations with measurements",
                )

            with progress.phase("build_collection", "Building sensor collection..."):
                sensor_collection.attributes.update(
                    {
                        "stations_used": stations_used,
                        "stations_skipped_upstream": stations_skipped_upstream,
                        "stations_skipped_no_coords": stations_skipped_no_coords,
                        "stations_skipped_no_timeseries": (
                            stations_skipped_no_timeseries
                        ),
                        "stations_skipped_no_value": stations_skipped_no_value,
                        "stale_value_count": stale_value_count,
                        "fallback_getdata_attempts": fallback_getdata_attempts,
                        "fallback_getdata_successes": fallback_getdata_successes,
                    }
                )
                self.apply_result_health_metadata(
                    sensor_collection.attributes,
                    upstream_errors=upstream_errors,
                    stations_skipped_upstream=stations_skipped_upstream,
                )

                info(
                    f"Successfully retrieved {stations_used} stations with measurements"
                )
                if stations_skipped_no_value > 0:
                    info(
                        f"Skipped {stations_skipped_no_value} stations with no "
                        "current measurements"
                    )
                report_progress(percent=50, message="Calculating bounds...")
                sensor_collection.calculate_bounds()

            with progress.phase(
                "export",
                (
                    "Exporting to protobuf..."
                    if args.format == "pb"
                    else "Preparing sensor collection..."
                ),
            ):
                if args.format == "pb":
                    return sensor_collection.to_proto().SerializeToString()
                else:
                    return sensor_collection


__all__ = ["AirQualityDatasetArgs", "AirQualityDataset"]
