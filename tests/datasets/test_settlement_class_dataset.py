import numpy as np
import pytest
import rasterio
from affine import Affine

import dtcc_core.datasets as datasets
from dtcc_core.io.data.ghsl import SETTLEMENT_CLASS_NODATA
from dtcc_core.model import Bounds, Raster


def _raster():
    return Raster(
        data=np.array([[30, 21], [11, SETTLEMENT_CLASS_NODATA]], dtype=np.int16),
        georef=Affine(1000.0, 0, 0.0, 0, -1000.0, 2000.0),
        nodata=SETTLEMENT_CLASS_NODATA,
        crs="EPSG:3006",
    )


def test_settlement_class_dataset_registered():
    assert "settlement_class" in datasets.list()
    assert callable(datasets.settlement_class)


def test_settlement_class_dataset_returns_raster(monkeypatch):
    calls = []

    def fake_download(bounds, epoch):
        calls.append((bounds, epoch))
        return _raster()

    monkeypatch.setattr("dtcc_core.io.data.download_settlement_class", fake_download)

    raster = datasets.settlement_class(bounds=(0.0, 0.0, 2000.0, 2000.0), epoch=1975)

    assert isinstance(raster, Raster)
    assert isinstance(calls[0][0], Bounds)
    assert calls[0][1] == 1975


def test_settlement_class_dataset_defaults_to_2020(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "dtcc_core.io.data.download_settlement_class",
        lambda bounds, epoch: calls.append(epoch) or _raster(),
    )

    datasets.settlement_class(bounds=(0.0, 0.0, 2000.0, 2000.0))

    assert calls == [2020]


def test_settlement_class_dataset_exports_tif_with_codes(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_settlement_class", lambda bounds, epoch: _raster()
    )

    payload = datasets.settlement_class(bounds=(0.0, 0.0, 2000.0, 2000.0), format="tif")
    path = tmp_path / "settlement_class.tif"
    path.write_bytes(payload)

    with rasterio.open(path) as src:
        assert src.dtypes[0] == "int16"
        assert src.nodata == SETTLEMENT_CLASS_NODATA
        assert src.read(1).tolist() == [[30, 21], [11, SETTLEMENT_CLASS_NODATA]]


def test_settlement_class_rejects_unsupported_epoch():
    with pytest.raises(Exception):
        datasets.settlement_class(bounds=(0.0, 0.0, 2000.0, 2000.0), epoch=2021)
