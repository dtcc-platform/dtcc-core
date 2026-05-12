"""Generate the seven-case smoke dataset table pack.

Companion to demos/smoke.py. Organizes the same family of products as
one (product, format) case per subdirectory under output/smoke/table_cases/
so that each case has a clean artifact + manifest pair.

When DTCC_UPLOAD_URL and DTCC_UPLOAD_TOKEN are set in the environment,
each case is also published to a dtcc-upload catalog under a
table-smoke-<name> dataset_key. Otherwise the script writes local files
only and prints an INFO line.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

BOUNDS: list[int] = [319720, 6397660, 320220, 6398160]
OUTPUT_DIR: Path = Path("output/smoke/table_cases")

CASES: list[dict[str, Any]] = [
    {
        "name": "field_vtu",
        "filename": "smoke_field.vtu",
        "dataset_key": "table-smoke-field-vtu",
        "params": {"product": "field", "resolution": 17},
    },
    {
        "name": "field_pb",
        "filename": "smoke_field.pb",
        "dataset_key": "table-smoke-field-pb",
        "params": {"product": "field", "resolution": 17},
    },
    {
        "name": "slice_geojson",
        "filename": "smoke_slice.geojson",
        "dataset_key": "table-smoke-slice-geojson",
        "params": {
            "product": "slice",
            "resolution": 31,
            "slice_axis": "z",
            "slice_position": 0.5,
        },
    },
    {
        "name": "streamlines_geojson",
        "filename": "smoke_streamlines.geojson",
        "dataset_key": "table-smoke-streamlines-geojson",
        "params": {
            "product": "streamlines",
            "streamline_count": 25,
            "streamline_steps": 120,
        },
    },
    {
        "name": "slice_png",
        "filename": "smoke_slice.png",
        "dataset_key": "table-smoke-slice-png",
        "params": {
            "product": "slice",
            "resolution": 128,
            "slice_axis": "z",
            "slice_position": 0.5,
            "profile": "table",
            "width": 1920,
            "height": 1920,
        },
    },
    {
        "name": "streamlines_png",
        "filename": "smoke_streamlines.png",
        "dataset_key": "table-smoke-streamlines-png",
        "params": {
            "product": "streamlines",
            "streamline_count": 64,
            "streamline_steps": 180,
            "profile": "table",
            "width": 1920,
            "height": 1920,
        },
    },
    {
        "name": "streamlines_mp4",
        "filename": "smoke_streamlines.mp4",
        "dataset_key": "table-smoke-streamlines-mp4",
        "params": {
            "product": "streamlines",
            "streamline_count": 48,
            "streamline_steps": 150,
            "profile": "table",
            "format": "mp4",
            "width": 1920,
            "height": 1920,
            "fps": 30,
            "duration": 8.0,
            "period": 8.0,
        },
    },
]
