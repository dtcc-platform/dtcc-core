"""ERA5 histories sampled through Open-Meteo's archive API."""

import numpy as np

from ._weather_history import _guard, _ns, _output_point, _store
from .dataset import DatasetUpstreamError

_VARIABLES = {
    "air_temperature": ("temperature_2m", "°C", "celsius", 2),
    "wind_speed": ("wind_speed_10m", "m/s", "m/s", 10),
    "wind_direction": ("wind_direction_10m", "°", "degrees", 10),
}
_ENDPOINT = "https://archive-api.open-meteo.com/v1/archive"
_CUSTOMER = "https://customer-archive-api.open-meteo.com/v1/archive"


def _build_open_meteo(args, result, client, lon, lat, errors):
    endpoint = _CUSTOMER if args.api_key else _ENDPOINT
    location_id = "open_meteo:location"
    records = {name: {} for name in args.parameter_names}
    successful = False
    returned_location = None
    left, end = _ns(args.start), _ns(args.end)
    while left < end:
        right = min(end, left + args.chunk_days * 86_400_000_000_000)
        query = {
            "latitude": lat,
            "longitude": lon,
            "start_date": str(np.datetime64(left, "ns").astype("datetime64[D]")),
            "end_date": str(np.datetime64(right - 1, "ns").astype("datetime64[D]")),
            "hourly": ",".join(_VARIABLES[name][0] for name in args.parameter_names),
            "models": "era5",
            "temperature_unit": "celsius",
            "wind_speed_unit": "ms",
            "timeformat": "unixtime",
            "timezone": "GMT",
        }
        if args.api_key:
            query["apikey"] = args.api_key
        try:
            payload = client.json(endpoint, query)
            with _guard(endpoint):
                parsed, location = _parse(
                    payload, args.parameter_names, left, right, client
                )
                if returned_location is not None and returned_location != location:
                    raise ValueError("Open-Meteo changed grid location between chunks")
                returned_location = location
                for name, readings in parsed.items():
                    if records[name].keys() & readings.keys():
                        raise ValueError("Overlapping Open-Meteo chunks")
                    records[name].update(readings)
                successful = True
        except DatasetUpstreamError as exc:
            if args.strict_live:
                raise
            errors.append(exc)
        left = right
    if not successful:
        return
    grid_lon, grid_lat, elevation = returned_location
    result.locations[location_id] = {
        "kind": "model_sample",
        "crs": args.crs,
        "requested_point": result.attributes["requested_point"],
        "point": _output_point(grid_lon, grid_lat, elevation, args.crs),
        "latitude": grid_lat,
        "longitude": grid_lon,
        "elevation": elevation,
        "requested_wgs84": {"longitude": lon, "latitude": lat},
        "elevation_policy": "provider_default_downscaling",
    }
    for name, readings in records.items():
        variable, source_unit, unit, height = _VARIABLES[name]
        ordered = sorted(readings)
        _store(
            result,
            (location_id, name),
            np.array(ordered, dtype="datetime64[ns]"),
            [readings[time] for time in ordered],
            unit,
            "reanalysis",
            {
                "source": "open_meteo",
                "data_kind": "reanalysis",
                "model": "era5",
                "model_requested": "era5",
                "variable": variable,
                "source_unit": source_unit,
                "height_above_ground_m": height,
                "temporal_semantics": "instantaneous",
                "description": f"ERA5 {variable} sampled through Open-Meteo",
                "quality_available": False,
                "endpoint": endpoint,
            },
            {},
            args,
        )


def _parse(payload, parameters, left, right, client):
    hourly, units = payload["hourly"], payload["hourly_units"]
    if payload["utc_offset_seconds"] != 0 or units["time"] != "unixtime":
        raise ValueError("Expected UTC Unix seconds")
    seconds = np.asarray(hourly["time"])
    if seconds.ndim != 1 or (seconds.size and seconds.dtype.kind not in "iu"):
        raise ValueError("Expected integer Unix seconds")
    if (
        seconds.size
        and seconds.dtype.kind == "u"
        and seconds.max() > np.iinfo(np.int64).max
    ):
        raise ValueError("Unix seconds exceed the signed integer range")
    times = seconds.astype("datetime64[s]")
    # TimeSeries' converter checks datetime range before converting to nanoseconds.
    from ..model.values.time_series import _to_times

    times = _to_times(times)[0]
    if len(times) > 1 and not (times[1:] > times[:-1]).all():
        raise ValueError("Unsorted or duplicated timestamps")
    ticks = times.astype(np.int64)
    mask = (ticks >= left) & (ticks < right)
    parsed = {}
    for name in parameters:
        variable, expected_unit, _unit, _height = _VARIABLES[name]
        if units[variable] != expected_unit:
            raise ValueError(f"Unexpected unit for {variable}")
        values = np.asarray(hourly[variable], dtype=float)
        if values.shape != times.shape or np.isinf(values).any():
            raise ValueError("Invalid observation array")
        parsed[name] = dict(
            zip(ticks[mask].tolist(), values[mask].tolist(), strict=True)
        )
        for _ in parsed[name]:
            client.count_observation()
    longitude, latitude, elevation = (
        float(payload[key]) for key in ("longitude", "latitude", "elevation")
    )
    if (
        not all(np.isfinite(value) for value in (longitude, latitude, elevation))
        or not -180 <= longitude <= 180
        or not -90 <= latitude <= 90
    ):
        raise ValueError("Invalid returned grid location")
    return parsed, (longitude, latitude, elevation)
