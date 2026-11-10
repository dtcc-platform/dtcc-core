import numpy as np
import pytest
import rasterio
from affine import Affine
from pydantic import ValidationError

import dtcc_core.datasets as datasets
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.model import Bounds, Raster


def _raster():
    return Raster(
        data=np.array([[0.0, 50.0], [100.0, np.nan]], dtype=np.float32),
        georef=Affine(10.0, 0, 0.0, 0, -10.0, 20.0),
        nodata=np.nan,
        crs="EPSG:3006",
    )


def _failing_download(bounds, year):
    raise DatasetUpstreamError(
        dataset="imperviousness",
        operation="download_image",
        target="https://example.test",
        failure_class="http_5xx",
        status_code=503,
        message="Imperviousness 2024 download failed: 503",
    )


def test_imperviousness_dataset_registered():
    assert "imperviousness" in datasets.list()
    assert callable(datasets.imperviousness)


def test_imperviousness_dataset_returns_raster(monkeypatch):
    calls = []

    def fake_download(bounds, year):
        calls.append((bounds, year))
        return _raster()

    monkeypatch.setattr("dtcc_core.io.data.download_imperviousness", fake_download)

    raster = datasets.imperviousness(bounds=(0.0, 0.0, 20.0, 20.0))

    assert isinstance(raster, Raster)
    assert isinstance(calls[0][0], Bounds)
    assert calls[0][1] == 2024
    assert raster.dataset_context.health["partial_result"] is False
    assert raster.dataset_context.health["status"] == "complete"


def test_imperviousness_dataset_exports_tif(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_imperviousness", lambda bounds, year: _raster()
    )

    payload = datasets.imperviousness(
        bounds=(0.0, 0.0, 20.0, 20.0), year=2018, format="tif"
    )
    path = tmp_path / "imperviousness.tif"
    path.write_bytes(payload)

    with rasterio.open(path) as src:
        assert src.crs.to_epsg() == 3006
        assert src.read(1)[0].tolist() == [0.0, 50.0]


def test_upstream_failure_returns_empty_partial_raster(monkeypatch):
    monkeypatch.setattr("dtcc_core.io.data.download_imperviousness", _failing_download)

    raster = datasets.imperviousness(bounds=(0.0, 0.0, 20.0, 20.0))

    assert raster.data.shape == (2, 2)
    assert np.isnan(raster.data).all()
    health = raster.dataset_context.health
    assert health["partial_result"] is True
    assert health["upstream_errors"][0]["failure_class"] == "http_5xx"
    assert health["status"] == "failed"
    presentation_warnings = raster.dataset_context.presentation.warnings
    assert presentation_warnings[0] == "Imperviousness 2024 download failed: 503"


def test_upstream_failure_raises_with_strict_live(monkeypatch):
    monkeypatch.setattr("dtcc_core.io.data.download_imperviousness", _failing_download)

    with pytest.raises(DatasetUpstreamError) as info:
        datasets.imperviousness(bounds=(0.0, 0.0, 20.0, 20.0), strict_live=True)
    assert info.value.failure_class == "http_5xx"


def test_imperviousness_rejects_unsupported_year():
    with pytest.raises(ValidationError):
        datasets.imperviousness(bounds=(0.0, 0.0, 20.0, 20.0), year=2015)
