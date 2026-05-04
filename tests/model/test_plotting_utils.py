import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg", force=True)
plt = pytest.importorskip("matplotlib.pyplot")
to_rgba = pytest.importorskip("matplotlib.colors").to_rgba

from dtcc_core import plotting


def test_plotting_exposes_dtcc_brand_defaults():
    assert plotting.DTCC_COLORS["yellow"] == "#FADA36"
    assert plotting.DTCC_COLORS["teal"] == "#78C8BE"
    assert plotting.DTCC_COLORS["orange"] == "#E35A1D"
    assert np.allclose(
        plotting.as_numeric_values([1, "2", None]),
        np.array([1.0, 2.0, np.nan]),
        equal_nan=True,
    )
    assert plotting.as_numeric_values([1, "primary"]) is None


def test_apply_dtcc_style_sets_axes_colors():
    _, ax = plt.subplots()

    plotting.apply_dtcc_style(ax, xlabel="x", ylabel="y", grid=True)

    assert ax.figure.get_facecolor() == to_rgba(plotting.DTCC_THEMES["dark"]["figure"])
    assert ax.xaxis.label.get_color() == plotting.DTCC_THEMES["dark"]["text"]
    assert ax.yaxis.label.get_color() == plotting.DTCC_THEMES["dark"]["text"]
    plt.close(ax.figure)


def test_add_plot_context_formats_metadata_and_bounds():
    _, ax = plt.subplots()

    plotting.add_plot_context(
        ax,
        title="DTCC Test",
        metadata={"Items": 2},
        bounds=(0.0, 1.0, 2.0, 3.0),
    )

    assert ax.get_title() == "DTCC Test"
    assert len(ax.texts) == 1
    assert "Items: 2" in ax.texts[0].get_text()
    assert "Bounds:" in ax.texts[0].get_text()
    assert "x 0.00 -> 2.00" in ax.texts[0].get_text()
    plt.close(ax.figure)
