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
    title = "Calibration Grid"
    description = (
        "Deterministic synthetic line grid spanning the requested bounds, "
        "used to check tangible-table projector alignment against known map "
        "coordinates."
    )
    ArgsModel = CalibrationGridArgs
    data_category = "synthetic"
    result_kind = "calibration_grid"
    python_return_type = "dtcc_core.model.CalibrationGrid"
    timeout_hint = 2
    provider = [{"name": "DTCC Platform", "role": "synthetic_generator"}]
    source = [
        {
            "name": "dtcc-core calibration grid generator",
            "role": "synthetic_generator",
            "url": "https://github.com/dtcc-platform/dtcc-core",
        }
    ]
    license = "MIT"
    collection_period = (
        "Not applicable: deterministic synthetic geometry generated from "
        "the request bounds and divisions."
    )
    data_types = ["vector", "line", "calibration_grid"]
    geographic_coverage = "requested synthetic bounds"
    update_frequency = "generated on demand"
    generated_at = "Computed at request time by dtcc-core."
    processing_steps = [
        "Validate requested bounds, CRS, and positive grid division count",
        "Sample divisions + 1 x positions and divisions + 1 y positions including both bounds edges",
        "Create vertical and horizontal LineString features that span the opposite axis",
        "Attach deterministic grid metadata including spacing, line count, bounds, divisions, and CRS",
    ]
    presentation_headline = "Table Calibration Grid"
    presentation_summary = (
        "A deterministic coordinate grid for checking whether the projected "
        "table overlay lines up with the physical 500 m model bounds."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Vertical and horizontal lines cover the requested bounds. "
                "The default 40 divisions create 41 lines in each direction, "
                "including the four outer edges."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "On the Gothenburg 500 m table model, 40 divisions across "
                "500 m corresponds to 12.5 m in map space or 1 cm on the "
                "1:1250 physical model."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "This is synthetic alignment geometry only. It does not "
                "confirm projector focus, physical fabrication accuracy, or "
                "whether any real-world dataset is spatially correct."
            ),
        },
    ]
    key_points = [
        "Generated locally from request bounds and division count",
        "Default spacing is 12.5 m for the canonical 500 m table bounds",
        "Every boundary edge is represented by a grid line",
        "No provider, network, or external data dependency",
    ]
    presentation_legend = {
        "title": "Calibration grid",
        "unit": "m",
        "entries": [
            {
                "label": "Grid line",
                "meaning": "constant x or y coordinate in the declared CRS",
            },
            {
                "label": "Outer line",
                "meaning": "requested table bounds edge",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "lines",
        "default_style": "calibration_grid",
        "table_role": "alignment",
        "default_crs": "EPSG:3006",
        "legend_required": False,
    }
    presentation_warnings = [
        "Use only with table/model bounds that match the physical table profile."
    ]
    presentation_limitations = [
        "Synthetic alignment helper; not a surveyed control network.",
        "A correct grid overlay does not validate real dataset source accuracy.",
        "The 1 cm physical spacing claim applies to the 500 m by 500 m, 1:1250 table profile.",
    ]

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
