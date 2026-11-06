import numpy as np
import pytest
import rasterio
from affine import Affine
from pydantic import ValidationError

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
    with pytest.raises(ValidationError, match="epoch"):
        datasets.settlement_class(bounds=(0.0, 0.0, 2000.0, 2000.0), epoch=2021)


@pytest.fixture
def settlement_1975(monkeypatch):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_settlement_class", lambda bounds, epoch: _raster()
    )
    return datasets.settlement_class(bounds=(0.0, 0.0, 2000.0, 2000.0), epoch=1975)


def test_every_settlement_class_has_a_fixed_colour():
    from dtcc_core.io.data.ghsl import SETTLEMENT_CLASSES

    hints = datasets.SettlementClassDataset.view_hints
    assert sorted(map(int, hints["class_colors"])) == sorted(SETTLEMENT_CLASSES)
    assert len(set(hints["class_colors"].values())) == len(SETTLEMENT_CLASSES)


def test_plot_draws_a_class_map_with_fixed_colours(settlement_1975):
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    colors = datasets.SettlementClassDataset.view_hints["class_colors"]

    ax = settlement_1975.plot(show=False)

    try:
        assert ax.name == "rectilinear"
        image = ax.images[0].get_array()
        np.testing.assert_allclose(image[0, 0], to_rgba(colors["30"]))
        np.testing.assert_allclose(image[0, 1], to_rgba(colors["21"]))
        np.testing.assert_allclose(image[1, 0], to_rgba(colors["11"]))
        assert image[1, 1][3] == 0  # no data is transparent
        assert ax.get_xlim() == pytest.approx((0.0, 2000.0))
        assert ax.get_title() == "Settlement Class, 1975"
    finally:
        plt.close("all")


def test_plot_legend_names_every_class_and_credits_the_source(settlement_1975):
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    from dtcc_core.io.data.ghsl import SETTLEMENT_CLASSES

    ax = settlement_1975.plot(show=False)

    try:
        legend = ax.get_legend()
        assert legend.get_title().get_text() == "Settlement class"
        labels = [text.get_text() for text in legend.get_texts()]
        assert sorted(labels) == sorted(SETTLEMENT_CLASSES.values())
        assert any("European Union" in text.get_text() for text in ax.texts)
    finally:
        plt.close("all")
