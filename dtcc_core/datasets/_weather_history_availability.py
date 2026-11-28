"""Point-based advertised coverage discovery without downloading observations."""

from contextlib import contextmanager
from datetime import UTC, datetime

import numpy as np
from pyproj import Geod

from ._weather_history import _Client, _guard, _link, _output_point
from ._weather_history_open_meteo import _VARIABLES
from ._weather_history_smhi import _PARAMETERS, _ROOT, _consistent_positions, _interval
from .dataset import DatasetDescriptor, DatasetUpstreamError

_ERA5_DOCUMENTATION = "https://open-meteo.com/en/docs/historical-weather-api"


@contextmanager
def _capture(args, result, parameter, target):
    try:
        with _guard(target):
            yield
    except DatasetUpstreamError as exc:
        if args.strict_live:
            raise
        result["complete"] = False
        if parameter is None:
            for item in result["parameters"].values():
                item["complete"] = False
        else:
            parameter["complete"] = False
        result["upstream_errors"].append(
            DatasetDescriptor.serialize_upstream_error(exc)
        )


def _bounds(resource):
    for key in ("from", "to"):
        if type(resource[key]) is not int:
            raise ValueError("Expected integer epoch-millisecond bounds")
    left, right = _interval(resource)
    limits = np.iinfo(np.int64)
    if not limits.min < left < right <= limits.max:
        raise ValueError("Invalid advertised interval")
    return left, right


def _iso(value):
    return str(np.datetime64(value, "ns"))


def _coordinates(resource):
    lon, lat = float(resource["longitude"]), float(resource["latitude"])
    height = resource.get("height")
    height = None if height is None else float(height)
    if (
        not np.isfinite([lon, lat]).all()
        or not -180 <= lon <= 180
        or not -90 <= lat <= 90
        or (height is not None and not np.isfinite(height))
    ):
        raise ValueError("Invalid station coordinates")
    return lon, lat, height


def _location(lon, lat, height, crs):
    point = _output_point(lon, lat, height or 0, crs)
    if not np.isfinite([point.x, point.y]).all():
        raise ValueError("Station coordinates cannot be represented in requested CRS")
    return dict(
        x=point.x, y=point.y, crs=crs, longitude=lon, latitude=lat, height_m=height
    )


def _near_positions(station, entry, args, geod):
    first, last = _bounds(entry)
    positions = []
    lon, lat = args._wgs84
    if not station["position"]:
        raise ValueError("Missing station position history")
    for position in station["position"]:
        left, right = _bounds(position)
        left, right = max(left, first), min(right, last)
        if left >= right:
            continue
        x, y, height = _coordinates(position)
        distance = float(geod.inv(lon, lat, x, y)[2])
        positions.append(
            dict(
                start=left,
                end=right,
                longitude=x,
                latitude=y,
                height=height,
                distance_m=distance,
            )
        )
    _consistent_positions(positions, positions, geod)
    return [p for p in positions if p["distance_m"] <= args.max_distance_m]


def _advertised_periods(period, name, url, positions, crs):
    first, last = _bounds(period)
    result = []
    for position in positions:
        left, right = max(first, position["start"]), min(last, position["end"])
        if left >= right:
            continue
        result.append(
            dict(
                start=_iso(left),
                end_exclusive=_iso(right),
                archive=name,
                source_url=url,
                coverage_kind="advertised",
                observations_verified=False,
                gap_free=None,
                distance_m=position["distance_m"],
                location=_location(
                    position["longitude"], position["latitude"], position["height"], crs
                ),
            )
        )
    return result


def _station(args, client, result, parameter, entry, geod):
    url = _link(entry, "station", "application/json")
    with _capture(args, result, parameter, url):
        station = client.json(url)
        positions = _near_positions(station, entry, args, geod)
        if not positions:
            return
        # Validate the full descriptor list before making nested requests.
        descriptors = [
            (item["key"], _link(item, "period", "application/json"))
            for item in station["period"]
            if item["key"] in {"latest-months", "corrected-archive"}
        ]
        periods = []
        for name, period_url in descriptors:
            with _capture(args, result, parameter, period_url):
                period = client.json(period_url)
                periods.extend(
                    _advertised_periods(period, name, period_url, positions, args.crs)
                )
        if periods:
            x, y, height = _coordinates(entry)
            parameter["candidates"].append(
                dict(
                    station_id=int(entry["key"]),
                    station_name=entry.get("name") or station["title"],
                    source_url=url,
                    catalogue_location=_location(x, y, height, args.crs),
                    distance_m=min(p["distance_m"] for p in periods),
                    distance_reference="minimum_over_advertised_periods",
                    periods=sorted(
                        periods,
                        key=lambda p: (p["start"], p["end_exclusive"], p["archive"]),
                    ),
                )
            )


def _parameter(args, client, result, name, resources, geod, historical):
    parameter = result["parameters"][name]
    pid, source_unit, _unit = _PARAMETERS[name]
    if pid not in resources:
        return
    with _capture(args, result, parameter, _ROOT):
        url = _link(resources[pid], "parameter", "application/json")
        catalog = client.json(url)
        if catalog["valueType"] != "SAMPLING" or catalog["unit"] != source_unit:
            raise ValueError("Unsupported SMHI measurement semantics or unit")
        for entry in catalog["station"]:
            with _capture(args, result, parameter, url):
                x, y, _height = _coordinates(entry)
                if historical or geod.inv(*args._wgs84, x, y)[2] <= args.max_distance_m:
                    _station(args, client, result, parameter, entry, geod)
    parameter["candidates"].sort(key=lambda c: (c["distance_m"], c["station_id"]))


def _smhi(args, client, result, historical):
    resources = {}
    with _capture(args, result, None, _ROOT):
        version = client.json(_ROOT)
        resources = {int(item["key"]): item for item in version["resource"]}
    if not result["complete"]:
        return
    geod = Geod(ellps="WGS84")
    for name in args.parameter_names:
        _parameter(args, client, result, name, resources, geod, historical)


def _era5(args, result):
    # Documented model-wide coverage; it is not a point-specific availability check.
    result["notes"].append(
        "ERA5 coverage follows provider documentation. The latest published date "
        "and the actual grid point are unknown; updates normally lag by five days."
    )
    for name in args.parameter_names:
        result["parameters"][name]["candidates"].append(
            dict(
                model="era5",
                spatial_coverage="global",
                grid_point_verified=False,
                source_url=_ERA5_DOCUMENTATION,
                periods=[
                    dict(
                        start="1940-01-01T00:00:00.000000000",
                        end_exclusive=None,
                        end_description="latest published ERA5 data",
                        typical_delay_days=5,
                        coverage_kind="advertised",
                        observations_verified=False,
                        gap_free=None,
                    )
                ],
            )
        )


def _status(parameter):
    if not parameter["complete"]:
        return "partial" if parameter["candidates"] else "unknown"
    return "advertised" if parameter["candidates"] else "not_advertised"


def _availability(args, historical):
    lon, lat = args._wgs84
    result = dict(
        source=args.source,
        query=dict(
            x=args.location[0],
            y=args.location[1],
            crs=args.crs,
            longitude=lon,
            latitude=lat,
            max_distance_m=args.max_distance_m,
        ),
        queried_at=datetime.now(UTC).isoformat(),
        coverage_kind="advertised",
        observations_verified=False,
        gap_free=None,
        complete=True,
        spatial_scope=("position_history" if historical else "catalogue_positions")
        if args.source == "smhi"
        else "global_model_documentation",
        notes=[
            "Only parameters supported by weather_history are listed. "
            "Advertised periods may contain missing readings; observations are not downloaded."
        ],
        parameters={
            name: dict(
                unit=_VARIABLES[name][2],
                data_kind="observation" if args.source == "smhi" else "reanalysis",
                candidates=[],
                complete=True,
            )
            for name in args.parameter_names
        },
        upstream_errors=[],
        requests=0,
        download_bytes=0,
    )
    if args.source == "smhi":
        if not historical:
            result["notes"].append(
                "Candidates are filtered by catalogue coordinates. Stations now "
                "outside the radius may have older nearby observations; use "
                "include_historical_locations=True to scan all position histories."
            )
        with _Client(args) as client:
            _smhi(args, client, result, historical)
            result.update(requests=client.requests, download_bytes=client.bytes)
    else:
        _era5(args, result)
    for parameter in result["parameters"].values():
        parameter["status"] = _status(parameter)
    return result
