import numpy as np
import pytest
import rasterio
from affine import Affine

import dtcc_core.datasets as datasets
from dtcc_core.model import Bounds, Raster


def _raster():
    return Raster(
        data=np.array([[0.0, 2500.0], [10000.0, np.nan]], dtype=np.float32),
        georef=Affine(100.0, 0, 0.0, 0, -100.0, 200.0),
        nodata=np.nan,
        crs="EPSG:3006",
    )


def test_built_up_dataset_registered():
    assert "built_up" in datasets.list()
    assert callable(datasets.built_up)


def test_built_up_dataset_returns_raster(monkeypatch):
    calls = []

    def fake_download(bounds, epoch, component):
        calls.append((bounds, epoch, component))
        return _raster()

    monkeypatch.setattr("dtcc_core.io.data.download_built_up", fake_download)

    raster = datasets.built_up(bounds=(0.0, 0.0, 200.0, 200.0))

    assert isinstance(raster, Raster)
    assert isinstance(calls[0][0], Bounds)
    assert calls[0][1:] == (2020, "total")


def test_built_up_dataset_passes_epoch_and_component(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "dtcc_core.io.data.download_built_up",
        lambda bounds, epoch, component: calls.append((epoch, component)) or _raster(),
    )

    datasets.built_up(bounds=(0.0, 0.0, 200.0, 200.0), epoch=1975, component="non_residential")

    assert calls == [(1975, "non_residential")]


def test_built_up_dataset_exports_tif(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_built_up", lambda bounds, epoch, component: _raster()
    )

    payload = datasets.built_up(bounds=(0.0, 0.0, 200.0, 200.0), format="tif")
    path = tmp_path / "built_up.tif"
    path.write_bytes(payload)

    with rasterio.open(path) as src:
        assert src.read(1)[0].tolist() == [0.0, 2500.0]


def test_built_up_rejects_unsupported_epoch():
    with pytest.raises(Exception):
        datasets.built_up(bounds=(0.0, 0.0, 200.0, 200.0), epoch=2022)
