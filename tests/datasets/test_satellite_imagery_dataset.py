"""Tests for the satellite imagery dataset wrapper."""

from __future__ import annotations

import datetime
import typing
from unittest.mock import Mock, patch

import numpy as np
import pytest
import requests
from affine import Affine
from pydantic import ValidationError

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.satellite_imagery import (
    ImageryStyle,
    SatelliteImageryArgs,
    SatelliteImageryDataset,
)
from dtcc_core.io.data import digitalearth
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


@pytest.mark.parametrize("date", ["22/09/2026", "2026-02-30"])
def test_rejects_invalid_date(date):
    with pytest.raises(ValidationError, match="valid date"):
        SatelliteImageryArgs(bounds=BOUNDS, date=date)


@pytest.mark.parametrize("date", ["2026-09-22", datetime.date(2026, 9, 22)])
def test_accepts_date_string_or_object(date):
    args = SatelliteImageryArgs(bounds=BOUNDS, date=date)
    assert args.date == datetime.date(2026, 9, 22)
    assert args.model_dump(mode="json")["date"] == "2026-09-22"


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
    assert kwargs["date"] == datetime.date(2026, 9, 22)
    assert kwargs["style"] == "ndvi"
    assert kwargs["max_cloud"] == 10.0
    assert kwargs["resolution"] == 10.0
    assert kwargs["search_days"] == 120
    assert kwargs["tile_errors"] == []


@pytest.mark.parametrize(("strict_live", "expected"), [(False, []), (True, None)])
@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_strict_live_decides_whether_failed_tiles_raise(
    mock_download, strict_live, expected
):
    mock_download.return_value = _rgba_raster()
    SatelliteImageryDataset().build(
        SatelliteImageryArgs(bounds=BOUNDS, strict_live=strict_live)
    )
    assert mock_download.call_args.kwargs["tile_errors"] == expected


def _failed_request(exc_type, url, status_code=None):
    request = requests.Request("GET", url).prepare()
    if status_code is None:
        return exc_type("boom", request=request)
    response = requests.Response()
    response.status_code = status_code
    response.request = request
    return exc_type("boom", request=request, response=response)


@pytest.mark.parametrize(
    ("exc", "operation", "failure_class"),
    [
        (
            _failed_request(requests.ConnectionError, digitalearth.WMS_URL + "/"),
            "wms_getmap",
            "connection",
        ),
        (
            _failed_request(requests.Timeout, digitalearth.STAC_URL + "/search"),
            "stac_search",
            "timeout",
        ),
        (
            _failed_request(requests.HTTPError, digitalearth.WMS_URL + "/", 503),
            "wms_getmap",
            "http_5xx",
        ),
    ],
)
@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_build_raises_typed_upstream_error(
    mock_download, exc, operation, failure_class
):
    mock_download.side_effect = exc
    with pytest.raises(DatasetUpstreamError) as raised:
        SatelliteImageryDataset().build(SatelliteImageryArgs(bounds=BOUNDS))
    assert raised.value.dataset == "satellite_imagery"
    assert raised.value.operation == operation
    assert raised.value.failure_class == failure_class
    assert raised.value.__cause__ is exc


@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_build_keeps_bad_input_errors_untyped(mock_download):
    mock_download.side_effect = ValueError("outside the date range")
    with pytest.raises(ValueError, match="outside the date range"):
        SatelliteImageryDataset().build(SatelliteImageryArgs(bounds=BOUNDS))


def test_build_reports_progress_per_tile(monkeypatch, tmp_path):
    from dtcc_core.common import progress as progress_module

    monkeypatch.setattr(digitalearth, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(digitalearth, "create_session", lambda: object())
    monkeypatch.setattr(digitalearth, "_imagery_coverage", lambda *a, **k: 1.0)
    monkeypatch.setattr(digitalearth, "get_map", lambda *a, **k: b"")
    monkeypatch.setattr(
        digitalearth,
        "_tile_to_array",
        lambda png, width, height: np.full((4, height, width), 255, dtype=np.uint8),
    )
    messages = []
    real_report = progress_module.report_progress

    def recording_report(*args, **kwargs):
        messages.append(kwargs.get("message"))
        return real_report(*args, **kwargs)

    monkeypatch.setattr(digitalearth, "report_progress", recording_report)

    SatelliteImageryDataset().build(
        SatelliteImageryArgs(bounds=(0.0, 0.0, 6000.0, 6000.0), date="2026-09-22")
    )
    tile_messages = [m for m in messages if m and m.startswith("Downloading tile")]
    assert tile_messages == [
        f"Downloading tile {i}/4 for 2026-09-22" for i in (1, 2, 3, 4)
    ]


def test_presentation_warns_that_imagery_is_not_quantitative():
    warnings = " ".join(SatelliteImageryDataset.presentation_warnings)
    assert "not surface reflectance or index values" in warnings


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


@patch("dtcc_core.datasets.satellite_imagery.dtcc_core.io.data.download_imagery")
def test_build_types_service_error_reports_as_invalid_payload(mock_download):
    exc = _failed_request(digitalearth.ServiceError, digitalearth.WMS_URL + "/")
    mock_download.side_effect = exc
    with pytest.raises(DatasetUpstreamError) as raised:
        SatelliteImageryDataset().build(SatelliteImageryArgs(bounds=BOUNDS))
    assert raised.value.operation == "wms_getmap"
    assert raised.value.failure_class == "invalid_payload"
    assert not raised.value.is_transient
    assert raised.value.__cause__ is exc
