"""Tests for the synthetic calibration grid dataset."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.calibration_grid import (
    CalibrationGridArgs,
    CalibrationGridDataset,
)

# The DTCC table bounds: 500 m x 500 m mapped onto the 40 cm printed model.
BOUNDS = [319720.0, 6397660.0, 320220.0, 6398160.0]


def test_calibration_grid_registered_name():
    """The dataset class should expose the correct registration name."""
    assert CalibrationGridDataset().name == "calibration_grid"


def test_calibration_grid_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("calibration_grid")
    assert ds is not None
    assert ds.name == "calibration_grid"


def test_calibration_grid_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "calibration_grid")
    assert callable(datasets.calibration_grid)


def test_default_build_returns_feature_collection_with_sweref_crs():
    """Without format, the dataset returns a GeoJSON dict declaring EPSG:3006."""
    result = datasets.calibration_grid(bounds=BOUNDS)

    assert isinstance(result, dict)
    assert result["type"] == "FeatureCollection"
    assert result["crs"]["properties"]["name"] == "EPSG:3006"


def test_default_divisions_40_gives_41_lines_per_axis():
    """divisions cells need divisions + 1 lines per axis, both edges included.

    The printed table model is 40 cm x 40 cm over the 500 m bounds, so the
    default 40 cells put lines exactly 1 cm apart on the model.
    """
    result = datasets.calibration_grid(bounds=BOUNDS)

    features = result["features"]
    vertical = [f for f in features if f["properties"]["orientation"] == "vertical"]
    horizontal = [
        f for f in features if f["properties"]["orientation"] == "horizontal"
    ]
    assert len(vertical) == 41
    assert len(horizontal) == 41
    assert len(features) == 82


def test_grid_lines_span_bounds_with_exact_spacing():
    """Lines must start at the min edge, end at the max edge, and be evenly spaced."""
    result = datasets.calibration_grid(bounds=BOUNDS, divisions=40)

    spacing = 12.5  # 500 m / 40 cells = exactly 1 cm on the 1:1250 model
    vertical = sorted(
        f["properties"]["position"]
        for f in result["features"]
        if f["properties"]["orientation"] == "vertical"
    )
    assert vertical[0] == pytest.approx(319720.0)
    assert vertical[-1] == pytest.approx(320220.0)
    for left, right in zip(vertical, vertical[1:]):
        assert right - left == pytest.approx(spacing)

    horizontal = sorted(
        f["properties"]["position"]
        for f in result["features"]
        if f["properties"]["orientation"] == "horizontal"
    )
    assert horizontal[0] == pytest.approx(6397660.0)
    assert horizontal[-1] == pytest.approx(6398160.0)


def test_grid_line_geometry_spans_the_opposite_axis():
    """A vertical line runs from ymin to ymax at its position; 2D coordinates only."""
    result = datasets.calibration_grid(bounds=BOUNDS)

    for feature in result["features"]:
        geometry = feature["geometry"]
        assert geometry["type"] == "LineString"
        coordinates = geometry["coordinates"]
        assert all(len(point) == 2 for point in coordinates)

        position = feature["properties"]["position"]
        if feature["properties"]["orientation"] == "vertical":
            assert [point[0] for point in coordinates] == [position, position]
            assert [point[1] for point in coordinates] == [6397660.0, 6398160.0]
        else:
            assert [point[1] for point in coordinates] == [position, position]
            assert [point[0] for point in coordinates] == [319720.0, 320220.0]


def test_geojson_format_returns_bytes_matching_dict():
    """format='geojson' should serialize the same FeatureCollection to bytes."""
    payload = datasets.calibration_grid(bounds=BOUNDS, format="geojson")

    assert isinstance(payload, bytes)
    assert json.loads(payload) == datasets.calibration_grid(bounds=BOUNDS)


def test_metadata_describes_the_grid():
    """The embedded metadata block ships in the published artifact."""
    result = datasets.calibration_grid(bounds=BOUNDS)

    metadata = result["metadata"]
    assert metadata["dataset"] == "calibration_grid"
    assert metadata["divisions"] == 40
    assert metadata["line_count"] == 82
    assert metadata["spacing"] == [12.5, 12.5]
    assert metadata["bounds"] == BOUNDS
    assert metadata["crs"] == "EPSG:3006"


def test_crs_none_omits_crs_member():
    """crs=None should omit the legacy GeoJSON crs member, mirroring smoke."""
    result = datasets.calibration_grid(bounds=BOUNDS, crs=None)
    assert "crs" not in result
    assert "crs" not in result["metadata"]


def test_export_writes_payload_and_manifest(tmp_path):
    """export() infers geojson from the suffix and writes artifact + manifest."""
    target = tmp_path / "calibration_grid.geojson"
    datasets.calibration_grid.export(target, bounds=BOUNDS, divisions=40)

    data = json.loads(target.read_text())
    assert data["type"] == "FeatureCollection"
    assert len(data["features"]) == 82

    manifest = json.loads(target.with_suffix(".manifest.json").read_text())
    assert manifest["format"] == "geojson"
    assert manifest["bounds"] == BOUNDS
    assert manifest["parameters"]["divisions"] == 40


def test_divisions_one_draws_only_the_frame():
    """divisions=1 is the degenerate grid: just the four boundary lines."""
    result = datasets.calibration_grid(bounds=BOUNDS, divisions=1)
    assert len(result["features"]) == 4


def test_divisions_zero_rejected():
    with pytest.raises(ValidationError):
        datasets.calibration_grid(bounds=BOUNDS, divisions=0)


def test_args_model_defaults():
    args = CalibrationGridArgs(bounds=BOUNDS)
    assert args.divisions == 40
    assert args.crs == "EPSG:3006"
    assert args.format is None
