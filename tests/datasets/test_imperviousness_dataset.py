import numpy as np
import pytest
import rasterio
from affine import Affine

import dtcc_core.datasets as datasets
from dtcc_core.model import Bounds, Raster


def _raster():
    return Raster(
        data=np.array([[0.0, 50.0], [100.0, np.nan]], dtype=np.float32),
        georef=Affine(10.0, 0, 0.0, 0, -10.0, 20.0),
        nodata=np.nan,
        crs="EPSG:3006",
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


def test_imperviousness_dataset_exports_tif(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_imperviousness", lambda bounds, year: _raster()
    )

    payload = datasets.imperviousness(bounds=(0.0, 0.0, 20.0, 20.0), year=2018, format="tif")
    path = tmp_path / "imperviousness.tif"
    path.write_bytes(payload)

    with rasterio.open(path) as src:
        assert src.crs.to_epsg() == 3006
        assert src.read(1)[0].tolist() == [0.0, 50.0]


def test_imperviousness_rejects_unsupported_year():
    with pytest.raises(Exception):
        datasets.imperviousness(bounds=(0.0, 0.0, 20.0, 20.0), year=2015)
