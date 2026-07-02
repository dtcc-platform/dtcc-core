"""Synthetic calibration grid dataset for table-projector alignment."""

from __future__ import annotations

import json
from typing import Any, Literal, Optional

import numpy as np
from pydantic import Field

from dtcc_core.model import Bounds, CalibrationGrid

from .dataset import DatasetBaseArgs, DatasetDescriptor


class CalibrationGridArgs(DatasetBaseArgs):
    """Arguments for the synthetic calibration grid dataset."""

    divisions: int = Field(
        40,
        ge=1,
        le=1000,
        description=(
            "Number of grid cells per axis. The grid draws divisions + 1 "
            "lines per axis, including both bounds edges. The default 40 "
            "puts lines 1 cm apart on the 40 cm printed table model when "
            "used with the 500 m table bounds."
        ),
    )
    crs: Optional[str] = Field(
        "EPSG:3006",
        description=(
            "Coordinate reference system declared in GeoJSON outputs for GIS "
            "readers. Set to None to omit the legacy GeoJSON CRS member."
        ),
    )
    format: Optional[Literal["geojson"]] = Field(
        None,
        description=(
            "Serialized output format. If omitted, the dataset returns a "
            "CalibrationGrid model."
        ),
    )


class CalibrationGridDataset(DatasetDescriptor):
    name = "calibration_grid"
    description = (
        "Synthetic alignment grid of evenly spaced lines spanning the "
        "requested bounds, for checking table-projector calibration."
    )
    ArgsModel = CalibrationGridArgs
    data_category = "derived"
    result_kind = "calibration_grid"
    python_return_type = "dtcc_core.model.CalibrationGrid"
    timeout_hint = 2
    provider = [{"name": "DTCC Platform", "role": "generator"}]
    source = ["Synthetic grid generated from requested bounds"]
    license = "MIT"
    geographic_coverage = "requested synthetic bounds"
    update_frequency = "generated on demand"
    processing_steps = ["Generate evenly spaced grid lines for requested bounds"]
    presentation_summary = (
        "Synthetic alignment grid for checking table-projector calibration."
    )
    key_points = [
        "Generated locally from the request",
        "Useful for table and projector alignment checks",
    ]
    view_hints = {"preferred_geometry": "lines", "default_style": "calibration_grid"}

    def build(self, args: CalibrationGridArgs):
        bounds = self.parse_bounds(args.bounds)
        geojson = _grid_geojson(bounds, args)
        if args.format is None:
            return geojson
        return json.dumps(geojson, separators=(",", ":")).encode("utf-8")

    def prepare_result(self, result, validated_args: CalibrationGridArgs):
        if validated_args.format is None and isinstance(result, dict):
            return CalibrationGrid.from_geojson(result)
        return result


def _grid_geojson(bounds: Bounds, args: CalibrationGridArgs) -> dict[str, Any]:
    xs = np.linspace(bounds.xmin, bounds.xmax, args.divisions + 1)
    ys = np.linspace(bounds.ymin, bounds.ymax, args.divisions + 1)

    features = [
        _line_feature("vertical", index, x, bounds.ymin, bounds.ymax)
        for index, x in enumerate(xs)
    ]
    features.extend(
        _line_feature("horizontal", index, y, bounds.xmin, bounds.xmax)
        for index, y in enumerate(ys)
    )

    collection = {
        "type": "FeatureCollection",
        "name": "calibration_grid",
        "features": features,
        "metadata": {
            "dataset": "calibration_grid",
            "divisions": args.divisions,
            "line_count": len(features),
            "spacing": [
                bounds.width / args.divisions,
                bounds.height / args.divisions,
            ],
            "bounds": [bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax],
        },
    }
    if args.crs is not None:
        collection["crs"] = {
            "type": "name",
            "properties": {"name": args.crs},
        }
        collection["metadata"]["crs"] = args.crs
    return collection


def _line_feature(
    orientation: str,
    index: int,
    position: float,
    start: float,
    stop: float,
) -> dict[str, Any]:
    position = float(position)
    if orientation == "vertical":
        coordinates = [[position, start], [position, stop]]
    else:
        coordinates = [[start, position], [stop, position]]
    return {
        "type": "Feature",
        "geometry": {"type": "LineString", "coordinates": coordinates},
        "properties": {
            "orientation": orientation,
            "index": index,
            "position": position,
        },
    }
