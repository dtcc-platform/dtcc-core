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

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ..model.object import Object, SensorCollection
from ..model.geometry import Point
from ..model.values import Field as DtccField
from ..common import info
from ..reproject.reproject import reproject_array


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
            "requests library required for air quality dataset. Install with: pip install requests"
        )

    try:
        response = requests.get(url, params=params, timeout=timeout_s)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to fetch {url}: {e}")
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Invalid JSON response from {url}: {e}")


def _resolve_phenomenon_id(base_url: str, phenomenon_str: str, timeout_s: float) -> str:
    """Resolve phenomenon name/ID to phenomenon ID.

    The API uses phenomenon IDs like "1" for NO2, "5" for PM10, etc.
    This function accepts either numeric IDs or common names.

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
    str
        Phenomenon ID
    """
    # If it's already numeric, return it
    if phenomenon_str.isdigit():
        info(f"Using numeric phenomenon ID: {phenomenon_str}")
        return phenomenon_str

    # Common mappings (verified against API 2024)
    common_mappings = {
        "NO2": "8",
        "O3": "7",
        "PM10": "5",
        "PM2.5": "6001",
        "SO2": "1",
        "CO": "2064",
    }

    if phenomenon_str.upper() in common_mappings:
        result = common_mappings[phenomenon_str.upper()]
        info(f"Mapped phenomenon '{phenomenon_str}' to ID {result}")
        return result

    # Try to fetch from API
    try:
        url = f"{base_url}/phenomena"
        data = _get_json(url, timeout_s=timeout_s)
        for phen in data:
            if phen.get("label", "").upper() == phenomenon_str.upper():
                return str(phen["id"])
    except Exception:
        pass

    # Default fallback - treat as ID
    return phenomenon_str


def _transform_bounds_to_wgs84(
    bounds: Tuple[float, float, float, float], crs: str
) -> Tuple[float, float, float, float]:
    """Transform bounds to WGS84 (CRS84) for API requests.

    The SMHI API expects coordinates in WGS84/CRS84 (longitude, latitude).
    This function transforms bounds from the specified CRS to WGS84 using
    the existing reproject_array function.

    Parameters
    ----------
    bounds : tuple
        (xmin, ymin, xmax, ymax) bounding box in source CRS
    crs : str
        Source coordinate reference system (e.g., "EPSG:3006", "CRS84")

    Returns
    -------
    tuple
        (lon_min, lat_min, lon_max, lat_max) in WGS84
    """
    # If already in CRS84/WGS84, return as-is
    if crs.upper() in ["CRS84", "EPSG:4326", "WGS84"]:
        return bounds

    # Transform corner points using existing reproject functionality
    xmin, ymin, xmax, ymax = bounds
    corners = np.array([[xmin, ymin, 0], [xmax, ymax, 0]])
    transformed = reproject_array(corners, crs, "EPSG:4326")

    return (transformed[0, 0], transformed[0, 1], transformed[1, 0], transformed[1, 1])


def _fetch_stations(
    base_url: str, bounds: Tuple[float, float, float, float], crs: str, timeout_s: float
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
    wgs84_bounds = _transform_bounds_to_wgs84(bounds, crs)

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
    except Exception as e:
        info(f"Warning: Failed to fetch stations: {e}")
        return []


def _fetch_timeseries_for_station(
    base_url: str, station_id: str, phenomenon_id: str, timeout_s: float
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
                except Exception:
                    # Skip this timeseries if we can't fetch it
                    pass

        return result

    except Exception:
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

            # Convert timestamp from milliseconds to ISO format
            timestamp_iso = ""
            if timestamp_ms:
                from datetime import datetime

                dt = datetime.fromtimestamp(timestamp_ms / 1000.0)
                timestamp_iso = dt.isoformat()
                info(
                    f"Extracted value {value} with timestamp {timestamp_iso} (raw: {timestamp_ms} ms)"
                )

            unit = timeseries_dict.get("uom", "")
            return (value, timestamp_iso, unit)
        except (ValueError, TypeError, KeyError):
            pass

    return None


def _fallback_get_latest_from_getData(
    base_url: str, timeseries_id: str, timeout_s: float
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
            f"Fallback: Fetching recent data from getData endpoint for timeseries {timeseries_id}"
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
    except Exception as e:
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

    Supported phenomena include NO2, O3, PM10, PM2.5, SO2, CO, and others.

    Example
    -------
    >>> import dtcc_core.datasets as datasets
    >>> # Fetch NO2 measurements in Stockholm area (SWEREF99 TM)
    >>> sensors = datasets.airquality(
    ...     bounds=(674000, 6580000, 676000, 6582000),
    ...     phenomenon="NO2"
    ... )
    >>> print(f"Found {len(sensors.stations())} stations")

    >>> # Get data as arrays
    >>> points, values = sensors.to_arrays()
    >>> print(f"NO2 range: {values.min():.1f} - {values.max():.1f}")
    """

    name = "airquality"
    description = (
        "Air quality sensor measurements from datavardluft.smhi.se. "
        "Returns a SensorCollection with latest snapshot readings from stations "
        "within the specified bounds. Each station Object has Point geometry and "
        "Field with measured value plus metadata attributes."
    )
    ArgsModel = AirQualityDatasetArgs

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
        # Parse bounds
        bounds = self.parse_bounds(args.bounds)
        bounds_tuple = (bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
        info(
            f"Fetching air quality data for {args.phenomenon} within bounds {bounds_tuple}"
        )

        # Resolve phenomenon
        phenomenon_id = _resolve_phenomenon_id(
            args.base_url, args.phenomenon, args.timeout_s
        )
        info(f"Resolved phenomenon {args.phenomenon} to ID {phenomenon_id}")

        # Fetch stations
        stations_data = _fetch_stations(
            args.base_url, bounds_tuple, args.crs, args.timeout_s
        )
        info(f"Found {len(stations_data)} stations in the specified area")

        # Create sensor collection
        sensor_collection = SensorCollection()
        sensor_collection.attributes = {
            "dtcc_type": "sensor_collection",
            "source": "datavardluft.smhi.se",
            "phenomenon": args.phenomenon,
            "phenomenon_id": phenomenon_id,
            "crs": args.crs,
            "bounds": f"{bounds_tuple}",
            "retrieval_time": datetime.now(timezone.utc).isoformat(),
            "total_stations_found": len(stations_data),
        }

        stations_used = 0
        stations_skipped_no_coords = 0
        stations_skipped_no_timeseries = 0
        stations_skipped_no_value = 0

        processing_logged = False
        for station_data in stations_data[: args.max_stations]:
            if not processing_logged:
                info(f"Processing stations and fetching measurements...")
                processing_logged = True
            station_id = str(station_data.get("id", ""))

            # Extract coordinates
            geometry = station_data.get("geometry", {})
            coordinates = geometry.get("coordinates", [])

            if len(coordinates) < 2:
                stations_skipped_no_coords += 1
                continue

            x, y = coordinates[0], coordinates[1]
            z = coordinates[2] if len(coordinates) > 2 else 0.0

            # Fetch timeseries for this station
            timeseries_list = _fetch_timeseries_for_station(
                args.base_url, station_id, phenomenon_id, args.timeout_s
            )
            # Progress logged per station would be too verbose

            if not timeseries_list:
                stations_skipped_no_timeseries += 1
                continue

            # Try to get latest value
            value = None
            timestamp = ""
            unit = ""

            for ts in timeseries_list:
                # First try to extract from metadata
                result = _extract_latest_value(ts)
                if result:
                    value, timestamp, unit = result
                    # Check if data is old (more than 7 days)
                    if timestamp:
                        from datetime import timedelta

                        try:
                            ts_dt = datetime.fromisoformat(
                                timestamp.replace("Z", "+00:00")
                            )
                            age = datetime.now(timezone.utc) - ts_dt.replace(
                                tzinfo=timezone.utc
                            )
                            if age > timedelta(days=7):
                                info(
                                    f"Warning: Station {station_id} has old data (age: {age.days} days)"
                                )
                        except Exception:
                            pass
                    break

                # Fallback: fetch from getData
                ts_id = str(ts.get("id", ""))
                if ts_id:
                    result = _fallback_get_latest_from_getData(
                        args.base_url, ts_id, args.timeout_s
                    )
                    if result:
                        value, timestamp, unit = result
                        break

            if value is None:
                if args.drop_missing:
                    stations_skipped_no_value += 1
                    continue
                else:
                    value = np.nan

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
            }

            # Create Point geometry
            point = Point(x=x, y=y, z=z)

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

        # Update collection attributes with statistics
        sensor_collection.attributes.update(
            {
                "stations_used": stations_used,
                "stations_skipped_no_coords": stations_skipped_no_coords,
                "stations_skipped_no_timeseries": stations_skipped_no_timeseries,
                "stations_skipped_no_value": stations_skipped_no_value,
            }
        )

        info(f"Successfully retrieved {stations_used} stations with measurements")
        if stations_skipped_no_value > 0:
            info(
                f"Skipped {stations_skipped_no_value} stations with no current measurements"
            )
        # Calculate bounds for collection
        sensor_collection.calculate_bounds()

        # Return based on format
        if args.format == "pb":
            return self.export_to_bytes(sensor_collection, "pb")
        else:
            return sensor_collection


__all__ = ["AirQualityDatasetArgs", "AirQualityDataset"]
