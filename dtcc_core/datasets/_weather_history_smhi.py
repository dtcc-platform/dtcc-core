"""SMHI station selection and historical metobs CSV observations."""

import csv
import io

import numpy as np
from pyproj import Geod

from ._weather_history import (
    NoMatchingStation,
    _guard,
    _link,
    _ns,
    _output_point,
    _stamp,
    _store,
)
from .dataset import DatasetUpstreamError

_ROOT = "https://opendata-download-metobs.smhi.se/api/version/latest.json"
_PARAMETERS = {
    "air_temperature": (1, "celsius", "celsius"),
    "wind_speed": (4, "meter per sekund", "m/s"),
    "wind_direction": (3, "grader", "degrees"),
}
_MS = 1_000_000


def _interval(resource):
    return int(resource["from"]) * _MS, (int(resource["to"]) + 1) * _MS


def _covers(intervals, start, end):
    cursor = start
    for left, right in sorted(intervals):
        if left > cursor:
            break
        cursor = max(cursor, right)
        if cursor >= end:
            return True
    return False


def _positions(station, start, end, lon, lat, geod):
    segments = []
    for position in station["position"]:
        left, right = _interval(position)
        if left >= end or right <= start:
            continue
        x, y, z = (float(position[key]) for key in ("longitude", "latitude", "height"))
        if (
            not np.isfinite([x, y, z]).all()
            or not -180 <= x <= 180
            or not -90 <= y <= 90
        ):
            raise ValueError("Invalid station position")
        distance = float(geod.inv(lon, lat, x, y)[2])
        segments.append(
            {
                "start": max(left, start),
                "end": min(right, end),
                "longitude": x,
                "latitude": y,
                "height": z,
                "distance_m": distance,
            }
        )
    if not _covers([(p["start"], p["end"]) for p in segments], start, end):
        return None
    segments.sort(key=lambda p: p["start"])
    return segments


def _consistent_positions(first, second, geod):
    for a in first:
        for b in second:
            if a["start"] < b["end"] and b["start"] < a["end"]:
                if (
                    geod.inv(
                        a["longitude"], a["latitude"], b["longitude"], b["latitude"]
                    )[2]
                    > 1
                ):
                    raise ValueError("Inconsistent station positions across parameters")


def _periods(station, start, end, client):
    result = []
    available = {item["key"]: item for item in station["period"]}
    for name in ("latest-months", "corrected-archive"):
        if name not in available:
            continue
        period = client.json(_link(available[name], "period", "application/json"))
        left, right = _interval(period)
        if left < end and right > start:
            result.append((name, period))
    return result


def _select(args, client, lon, lat, excluded):
    start, end = _ns(args.start), _ns(args.end)
    geod = Geod(ellps="WGS84")
    version = client.json(_ROOT)
    resources = {int(item["key"]): item for item in version["resource"]}
    catalogs, entries = {}, {}
    for name in args.parameter_names:
        pid, expected_unit, _unit = _PARAMETERS[name]
        catalog = client.json(_link(resources[pid], "parameter", "application/json"))
        if catalog["valueType"] != "SAMPLING" or catalog["unit"] != expected_unit:
            raise ValueError("Unsupported SMHI measurement semantics or unit")
        catalogs[name] = catalog
        entries[name] = {
            int(item["key"]): item
            for item in catalog["station"]
            if _covers([_interval(item)], start, end)
        }
    common = set.intersection(*(set(items) for items in entries.values()))
    ranked = []
    for station_id in sorted(common):
        documents, positions = {}, None
        for name in args.parameter_names:
            station = client.json(
                _link(entries[name][station_id], "station", "application/json")
            )
            current = _positions(station, start, end, lon, lat, geod)
            if current is None:
                excluded.append(
                    {
                        "station_id": station_id,
                        "parameter": name,
                        "reason": "position_coverage",
                    }
                )
                break
            _consistent_positions(current, current, geod)
            if positions is not None:
                _consistent_positions(positions, current, geod)
            else:
                positions = current
            documents[name] = station
        if len(documents) != len(args.parameter_names):
            continue
        distance = max(p["distance_m"] for p in positions)
        if distance <= args.max_distance_m:
            ranked.append((distance, station_id, documents, positions))
    for distance, station_id, documents, positions in sorted(
        ranked, key=lambda row: row[:2]
    ):
        periods = {
            name: _periods(station, start, end, client)
            for name, station in documents.items()
        }
        if all(
            _covers([_interval(period) for _key, period in items], start, end)
            for items in periods.values()
        ):
            return station_id, distance, positions, documents, periods, catalogs
    return None


def _csv_header(reader, unit):
    description = ""
    height = None
    for row in reader:
        if row[:1] == ["Stationsnamn"]:
            station = next(reader)
            if len(station) > 3 and station[3].strip():
                height = float(station[3])
                if not np.isfinite(height):
                    raise ValueError("Invalid measurement height")
        if row[:1] == ["Parameternamn"]:
            parameter = next(reader)
            if parameter[2].strip() != unit:
                raise ValueError("SMHI CSV unit differs from catalogue")
            description = parameter[1]
        if row[:2] == ["Datum", "Tid (UTC)"]:
            if not description or row[3] != "Kvalitet":
                raise ValueError("Invalid SMHI history CSV header")
            return description, height
    raise ValueError("Missing SMHI history CSV header")


def _parse_csv(text, unit, args, client):
    """Read station histories, preserving per-reading quality and missing values."""
    reader = csv.reader(io.StringIO(text.lstrip("\ufeff")), delimiter=";")
    description, height = _csv_header(reader, unit)
    records = {}
    start, end = _stamp(args.start), _stamp(args.end)
    first_day, last_day = str(start)[:10], str(end)[:10]
    for row in reader:
        if not row or not any(row):
            continue
        if len(row) < 4:
            raise ValueError("Incomplete SMHI observation")
        # Full archives can span decades; convert times only in the requested days.
        if not first_day <= row[0] <= last_day:
            continue
        stamp = _stamp(f"{row[0]}T{row[1]}Z")
        if not start <= stamp < end:
            continue
        key = str(stamp)
        if key in records:
            raise ValueError("Duplicate timestamp within a SMHI period")
        value = float(row[2]) if row[2].strip() else float("nan")
        if np.isinf(value):
            raise ValueError("Infinite SMHI observation")
        client.count_observation()
        records[key] = {"value": value, "quality": row[3].strip()}
    return records, description, height


def _build_smhi(args, result, client, lon, lat, errors):
    selection = {
        "status": "no_matching_station",
        "method": "nearest_common_station",
        "requested_point": result.attributes["requested_point"],
        "crs": args.crs,
        "requested_wgs84": {"longitude": lon, "latitude": lat},
        "max_distance_m": args.max_distance_m,
        "distance_reference": "maximum_over_requested_interval",
        "excluded_stations": [],
    }
    result.attributes["selection"] = selection
    try:
        with _guard(_ROOT):
            selected = _select(args, client, lon, lat, selection["excluded_stations"])
    except DatasetUpstreamError as exc:
        selection["status"] = "discovery_failed"
        if args.strict_live:
            raise
        errors.append(exc)
        return
    if selected is None:
        if args.strict_live:
            raise NoMatchingStation(
                f"No station offers {args.parameter_names} for [{args.start}, {args.end}) "
                f"within {args.max_distance_m:g} m"
            )
        return
    station_id, distance, positions, documents, periods, catalogs = selected
    station = documents[args.parameter_names[0]]
    position = positions[0]
    point = _output_point(
        position["longitude"], position["latitude"], position["height"], args.crs
    )
    history = [
        {
            "from": str(np.datetime64(p["start"], "ns")),
            "to_exclusive": str(np.datetime64(p["end"], "ns")),
            **{
                key: p[key] for key in ("longitude", "latitude", "height", "distance_m")
            },
            "point": _output_point(
                p["longitude"], p["latitude"], p["height"], args.crs
            ),
        }
        for p in positions
    ]
    selection.update(
        status="matched",
        station_id=station_id,
        station_name=station["title"],
        station_location=point,
        distance_m=distance,
        longitude=position["longitude"],
        latitude=position["latitude"],
        position_history=history,
    )
    location_id = f"smhi:{station_id}"
    result.locations[location_id] = {
        **selection,
        "point": point,
        "kind": "station",
        "elevation": position["height"],
    }
    for name in args.parameter_names:
        pid, source_unit, unit = _PARAMETERS[name]
        records, sources = {}, []
        description, measurement_height = catalogs[name]["summary"], None
        for period_name, period in periods[name]:
            target = _ROOT
            try:
                with _guard(_ROOT):
                    target = next(
                        _link(data, "data", "text/plain")
                        for data in period["data"]
                        if any(
                            link["rel"] == "data" and link["type"] == "text/plain"
                            for link in data["link"]
                        )
                    )
                with _guard(target):
                    text = client.text(target, cache_token=period["updated"])
                    incoming, description, measurement_height = _parse_csv(
                        text, source_unit, args, client
                    )
                    for reading in incoming.values():
                        reading.update(source_period=period_name, source_url=target)
                    # Corrected records replace provisional values and ancillary data together.
                    records.update(incoming)
                    sources.append({"url": target, "updated": period["updated"]})
            except DatasetUpstreamError as exc:
                if args.strict_live:
                    raise
                errors.append(exc)
        if sources:
            ordered = sorted(records)
            _store(
                result,
                (location_id, name),
                np.array(ordered, dtype="datetime64[ns]"),
                [records[t]["value"] for t in ordered],
                unit,
                "observed",
                {
                    "source": "smhi",
                    "data_kind": "observation",
                    "parameter_id": pid,
                    "source_unit": source_unit,
                    "description": description,
                    "sources": sources,
                    "height_above_ground_m": measurement_height,
                    "temporal_semantics": catalogs[name]["summary"],
                    "quality_available": True,
                },
                {
                    t: {k: v for k, v in records[t].items() if k != "value"}
                    for t in ordered
                },
                args,
            )
