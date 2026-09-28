"""Tests for the satellite imagery dataset wrapper."""

from __future__ import annotations

import typing
from unittest.mock import Mock, patch

import numpy as np
import pytest
from affine import Affine
from pydantic import ValidationError

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.satellite_imagery import (
    ImageryStyle,
    SatelliteImageryArgs,
    SatelliteImageryDataset,
)
from dtcc_core.io.data.digitalearth import VALID_STYLES
from dtcc_core.model import Raster


BOUNDS = (319000.0, 6397000.0, 321000.0, 6399000.0)


def _rgba_raster():
    raster = Raster()
    raster.data = np.full((2, 3, 4), 255, dtype=np.uint8)
    raster.georef = Affine(10.0, 0.0, BOUNDS[0], 0.0, -10.0, BOUNDS[3])
    raster.crs = "EPSG:3006"
    return raster


def test_satellite_imagery_registered_name():
    assert SatelliteImageryDataset().name == "satellite_imagery"


def test_satellite_imagery_get_dataset():
    ds = get_dataset("satellite_imagery")
    assert ds is not None
    assert ds.name == "satellite_imagery"


def test_satellite_imagery_module_attribute():
    assert hasattr(datasets, "satellite_imagery")
    assert callable(datasets.satellite_imagery)


def test_style_choices_match_download_module():
    assert list(typing.get_args(ImageryStyle)) == VALID_STYLES


def test_rejects_malformed_date():
    with pytest.raises(ValidationError, match="YYYY-MM-DD"):
        SatelliteImageryArgs(bounds=BOUNDS, date="22/09/2026")


def test_rejects_unknown_style():
    with pytest.raises(ValidationError):
        SatelliteImageryArgs(bounds=BOUNDS, style="infrared")


def test_rejects_out_of_range_cloud_cover():
    with pytest.raises(ValidationError):
        SatelliteImageryArgs(bounds=BOUNDS, max_cloud=150)


@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_default_build_forwards_options_and_returns_raster(mock_download):
    downloaded = Mock(name="raster")
    mock_download.return_value = downloaded

    result = SatelliteImageryDataset().build(
        SatelliteImageryArgs(bounds=BOUNDS, date="2026-09-22", style="ndvi")
    )

    assert result is downloaded
    kwargs = mock_download.call_args.kwargs
    assert kwargs["bounds"].tuple == BOUNDS
    assert kwargs["provider"] == "DES"
    assert kwargs["epsg"] == "3006"
    assert kwargs["date"] == "2026-09-22"
    assert kwargs["style"] == "ndvi"
    assert kwargs["max_cloud"] == 10.0
    assert kwargs["resolution"] == 10.0
    assert kwargs["search_days"] == 120


@pytest.mark.parametrize(
    ("format", "magic"),
    [("tif", (b"II*\x00", b"MM\x00*")), ("png", (b"\x89PNG",))],
)
@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_build_serializes_requested_format(mock_download, format, magic):
    mock_download.return_value = _rgba_raster()

    result = SatelliteImageryDataset().build(
        SatelliteImageryArgs(bounds=BOUNDS, format=format)
    )

    assert isinstance(result, bytes)
    assert result.startswith(magic)


@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_call_attaches_epsg3006_context(mock_download):
    mock_download.return_value = _rgba_raster()

    result = datasets.satellite_imagery(bounds=BOUNDS)

    assert "EPSG:3006" in result.dataset_context.metadata.crs
    assert result.dataset_context.identity.name == "satellite_imagery"
