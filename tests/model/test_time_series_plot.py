"""Temporal previews preserve the selected data and compose on existing axes."""

import datetime as dt
from unittest.mock import patch

import numpy as np
import pytest

from dtcc_core.model import TimeSeries
from dtcc_core.plotting.style import DTCC_THEMES

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg", force=True)
plt = pytest.importorskip("matplotlib.pyplot")


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def series(values, **kwargs):
    times = np.datetime64("2026-01-01", "ns") + np.arange(
        len(values)
    ) ** 2 * np.timedelta64(1, "h")
    return TimeSeries(
        name="Temperature", unit="C", times=times, values=values, **kwargs
    )


def test_scalar_dates_missing_values_labels_and_no_mutation():
    ts = series([1.0, np.nan, np.inf, 4.0], role="observed")
    before = ts.values.copy()
    with patch.object(plt, "show") as show:
        ax = ts.plot(show=False)
        show.assert_not_called()
    ax.figure.canvas.draw()
    assert ax.name == "rectilinear"
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), ts.times)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), [1, np.nan, np.nan, 4])
    np.testing.assert_array_equal(ts.values, before)
    assert ax.get_xlabel() == "Time (UTC)"
    assert ax.get_ylabel() == "Temperature [C]"
    assert ax.get_title() == "Temperature"
    assert ax.get_legend().get_texts()[0].get_text() == "Temperature (observed)"
    assert ax.xaxis.get_major_formatter()._tz == dt.UTC
    assert ax.xaxis.get_offset_text().get_color() == DTCC_THEMES["dark"]["muted"]
    assert ax.get_facecolor() == matplotlib.colors.to_rgba(DTCC_THEMES["dark"]["axes"])


def test_element_selection_preserves_order_and_indices():
    ts = series(np.arange(12).reshape(4, 3))
    ax = ts.plot(elements=[2, 0], show=False)
    assert len(ax.lines) == 2
    for line, index in zip(ax.lines, [2, 0], strict=True):
        np.testing.assert_array_equal(line.get_ydata(), ts.values[:, index])
        assert f"element {index}" in line.get_label()


def test_default_scalar_elements():
    ts = series([[1, 2], [3, 4]])
    ax = ts.plot(show=False)
    assert len(ax.lines) == 2
    np.testing.assert_array_equal(ax.lines[1].get_ydata(), [2, 4])


@pytest.mark.parametrize("component", [None, np.int64(1), "magnitude"])
def test_vector_components_and_magnitude(component):
    ts = series(np.arange(24).reshape(4, 2, 3), dim=3)
    ax = ts.plot(elements=np.int64(1), component=component, show=False)
    expected = ts.values[:, 1]
    if component == "magnitude":
        expected = np.linalg.norm(expected, axis=-1)[:, None]
        assert "magnitude" in ax.get_ylabel()
    elif component is not None:
        expected = expected[:, [component]]
    assert len(ax.lines) == expected.shape[1]
    for index, line in enumerate(ax.lines):
        np.testing.assert_allclose(line.get_ydata(), expected[:, index])
        assert "element 1" in line.get_label()
        assert (
            "magnitude" if component == "magnitude" else "component"
        ) in line.get_label()


def test_magnitude_avoids_squared_overflow_and_keeps_gaps():
    ts = series([[[3e200, 4e200]], [[np.nan, 1]]], dim=2)
    ax = ts.plot(component="magnitude", show=False)
    np.testing.assert_allclose(ax.lines[0].get_ydata(), [5e200, np.nan])
    integer = series([np.iinfo(np.int64).min])
    ax = integer.plot(component="magnitude", show=False)
    assert ax.lines[0].get_ydata()[0] > 0


def test_axes_overlay_style_kwargs_and_show():
    observed = series([1, 2], role="observed")
    predicted = series([3, 4], role="predicted")
    ax = observed.plot(show=False)
    with patch.object(plt, "show") as show:
        result = predicted.plot(
            ax=ax,
            theme="light",
            title="Comparison",
            label="Forecast",
            c="red",
            linestyle="--",
            marker="s",
        )
        show.assert_called_once_with()
    assert result is ax
    assert len(ax.lines) == 2
    assert ax.get_title() == "Comparison"
    assert ax.get_facecolor() == matplotlib.colors.to_rgba(DTCC_THEMES["light"]["axes"])
    assert ax.lines[1].get_color() == "red"
    assert ax.lines[1].get_linestyle() == "--"
    assert ax.lines[1].get_marker() == "s"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [
        "Temperature (observed)",
        "Forecast",
    ]


def test_single_timestamp_and_scalar_component():
    ts = series([-4])
    ax = ts.plot(component=0, title="", show=False)
    assert ax.lines[0].get_marker() == "o"
    assert ax.get_title() == ""
    ax.figure.canvas.draw()


def test_timezone_aware_input_is_plotted_in_utc_despite_matplotlib_defaults():
    ts = TimeSeries(times=["2026-01-01T02:00+02:00"], values=[1])
    with matplotlib.rc_context({"timezone": "US/Pacific"}):
        ax = ts.plot(show=False)
        ax.figure.canvas.draw()
    np.testing.assert_array_equal(
        ax.lines[0].get_xdata(), np.array(["2026-01-01T00:00"], dtype="datetime64[ns]")
    )
    assert ax.xaxis.get_major_formatter()._tz == dt.UTC


@pytest.mark.parametrize(
    "elements", [True, -1, 2, [0, 0], [0.0], [[0]], "0", 0.5, [False]]
)
def test_invalid_elements(elements):
    with pytest.raises(ValueError, match="element"):
        series([[1, 2]]).plot(elements=elements, show=False)
    assert not plt.get_fignums()


@pytest.mark.parametrize("component", [True, -1, 2, "norm", [0], 0.5])
def test_invalid_components(component):
    with pytest.raises(ValueError, match="component"):
        series([[[1, 2]]], dim=2).plot(component=component, show=False)
    assert not plt.get_fignums()


@pytest.mark.parametrize(
    "ts,elements",
    [
        (TimeSeries(), None),
        (series([1]), []),
        (series([np.nan, np.inf]), None),
        (series([[1, np.nan], [2, np.nan]]), 1),
    ],
)
def test_empty_or_nonfinite_selection(ts, elements):
    with pytest.raises(ValueError, match="No .*values"):
        ts.plot(elements=elements, show=False)
    assert not plt.get_fignums()


def test_limits_apply_before_rendering_and_selection_can_reduce_work():
    ts = series(np.ones((2, 21)))
    with pytest.raises(ValueError, match="20 curves"):
        ts.plot(show=False)
    assert not plt.get_fignums()
    assert len(ts.plot(elements=20, show=False).lines) == 1
    vector = series(np.ones((2, 7, 3)), dim=3)
    with pytest.raises(ValueError, match="20 curves"):
        vector.plot(show=False)
    assert len(vector.plot(component="magnitude", show=False).lines) == 7
    large = TimeSeries(
        times=np.datetime64("2026-01-01") + np.arange(50_001) * np.timedelta64(1, "s"),
        values=np.zeros((50_001, 20)),
    )
    with pytest.raises(ValueError, match="1,000,000 values"):
        large.plot(show=False)


@pytest.mark.parametrize("projection", ["3d", "polar"])
def test_rejects_incompatible_axes(projection):
    ax = plt.figure().add_subplot(projection=projection)
    with pytest.raises(ValueError, match="rectilinear"):
        series([1]).plot(ax=ax, show=False)
    assert not ax.lines


def test_time_filtering_and_resampling_compose_with_plot():
    ts = series([1, 2, 3, 4])
    selected = ts.between("2026-01-01T01:00", "2026-01-01T09:00").resample("1h")
    ax = selected.plot(show=False)
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), selected.times)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), selected.values[:, 0])
