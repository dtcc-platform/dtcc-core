"""Tests for smoke dataset plotting modes."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "dtcc-matplotlib-tests"),
)

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import dtcc_core.datasets as datasets


BOUNDS = (0.0, 0.0, 10.0, 20.0)


def _small_kwargs() -> dict:
    return {
        "bounds": BOUNDS,
        "resolution": 16,
        "streamline_count": 4,
        "streamline_steps": 20,
        "width": 640,
        "height": 360,
    }


def _figure_text(figure: Figure) -> str:
    texts = [text.get_text() for text in figure.texts]
    for ax in figure.axes:
        texts.extend(text.get_text() for text in ax.texts)
    return "\n".join(texts)


def _close(result) -> None:
    if isinstance(result, Figure):
        plt.close(result)
    else:
        plt.close(result.figure)


def test_smoke_plot_default_preview_returns_figure():
    fig = datasets.smoke.plot(bounds=BOUNDS, show=False)

    try:
        assert isinstance(fig, Figure)
        assert len(fig.axes) >= 3
        text = _figure_text(fig)
        assert "Synthetic Urban Smoke Flow" in text
        assert "What you are seeing" in text
        assert "How to read it" in text
        assert "Fast corridor" in text
        assert "Recirculation" in text
        assert "Key facts" in text
        assert "not a validated CFD simulation" in text
    finally:
        _close(fig)


def test_smoke_plot_explicit_preview_returns_figure():
    fig = datasets.smoke.plot(**_small_kwargs(), mode="preview", show=False)

    try:
        assert isinstance(fig, Figure)
        text = _figure_text(fig)
        assert "Synthetic Urban Smoke Flow" in text
        assert "Smoke speed" in text
        assert "PNG preview" in text
    finally:
        _close(fig)


def test_smoke_plot_preview_rejects_ax():
    fig, ax = plt.subplots()

    try:
        with pytest.raises(ValueError, match="full figure layout"):
            datasets.smoke.plot(**_small_kwargs(), mode="preview", ax=ax, show=False)
    finally:
        plt.close(fig)


def test_smoke_plot_with_ax_defaults_to_plot_mode():
    fig, ax = plt.subplots()

    try:
        returned = datasets.smoke.plot(**_small_kwargs(), ax=ax, show=False)

        assert returned is ax
        assert returned.axison
        assert returned.get_xlabel()
        assert returned.get_ylabel()
        assert returned.get_title() == "DTCC Smoke Slice"
    finally:
        plt.close(fig)


def test_smoke_field_slice_object_plot_defaults_to_preview():
    field_slice = datasets.smoke(**_small_kwargs(), product="slice")

    fig = field_slice.plot(show=False)

    try:
        assert isinstance(fig, Figure)
        text = _figure_text(fig)
        assert "Synthetic Urban Smoke Flow" in text
        assert "What you are seeing" in text
        assert "Fast corridor" in text
    finally:
        _close(fig)


def test_smoke_field_slice_object_plot_with_ax_defaults_to_plot_mode():
    field_slice = datasets.smoke(**_small_kwargs(), product="slice")
    fig, ax = plt.subplots()

    try:
        returned = field_slice.plot(ax=ax, show=False)

        assert returned is ax
        assert returned.axison
        assert returned.get_title() == "DTCC Smoke Slice"
    finally:
        plt.close(fig)


def test_smoke_streamline_object_plot_defaults_to_preview():
    streamlines = datasets.smoke(**_small_kwargs(), product="streamlines")

    fig = streamlines.plot(show=False)

    try:
        assert isinstance(fig, Figure)
        text = _figure_text(fig)
        assert "Synthetic Urban Smoke Flow" in text
        assert "Smoke speed" in text
        assert "Key facts" in text
    finally:
        _close(fig)


def test_smoke_plot_mode_plot_returns_axes():
    ax = datasets.smoke.plot(**_small_kwargs(), mode="plot", show=False)

    try:
        assert isinstance(ax, Axes)
        assert ax.axison
        assert ax.get_xlabel() == "x"
        assert ax.get_ylabel() == "y"
        assert ax.get_title() == "DTCC Smoke Slice"
        assert len(ax.figure.axes) >= 2
    finally:
        _close(ax)


def test_smoke_plot_mode_artifact_returns_axes_without_axis():
    ax = datasets.smoke.plot(**_small_kwargs(), mode="artifact", show=False)

    try:
        assert isinstance(ax, Axes)
        assert not ax.axison
        assert ax.get_title() == ""
        assert ax.get_xlabel() == ""
        assert ax.get_ylabel() == ""
        assert len(ax.figure.axes) == 1
        text = _figure_text(ax.figure)
        assert "Synthetic Urban Smoke Flow" not in text
        assert "Fast corridor" not in text
    finally:
        _close(ax)


def test_smoke_plot_invalid_mode_fails():
    with pytest.raises(ValueError, match="artifact, preview, plot"):
        datasets.smoke.plot(**_small_kwargs(), mode="dashboard", show=False)


def test_smoke_plot_profile_alias_table_maps_to_artifact():
    ax = datasets.smoke.plot(**_small_kwargs(), profile="table", show=False)

    try:
        assert isinstance(ax, Axes)
        assert not ax.axison
        assert len(ax.figure.axes) == 1
    finally:
        _close(ax)


def test_smoke_plot_profile_alias_python_maps_to_plot():
    ax = datasets.smoke.plot(**_small_kwargs(), profile="python", show=False)

    try:
        assert isinstance(ax, Axes)
        assert ax.axison
        assert ax.get_title() == "DTCC Smoke Slice"
    finally:
        _close(ax)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"mode": "artifact", "profile": "python"},
        {"mode": "plot", "profile": "table"},
        {"mode": "preview", "profile": "table"},
    ],
)
def test_smoke_plot_conflicting_mode_profile_fails(kwargs):
    with pytest.raises(ValueError, match="conflicts|not compatible"):
        datasets.smoke.plot(**_small_kwargs(), **kwargs, show=False)


def test_smoke_describe_includes_structured_presentation():
    presentation = datasets.smoke.describe()["presentation"]

    assert presentation["headline"] == "Synthetic Urban Smoke Flow"
    assert presentation["summary"]
    assert len(presentation["narrative"]) >= 2
    assert presentation["legend"]["title"] == "Smoke speed"
    assert len(presentation["annotations"]) >= 2
    assert presentation["key_points"]
    assert presentation["limitations"]
    for annotation in presentation["annotations"]:
        x, y = annotation["position"]
        assert 0.0 <= x <= 1.0
        assert 0.0 <= y <= 1.0


def test_smoke_dataset_context_includes_structured_presentation():
    mesh = datasets.smoke(bounds=BOUNDS, resolution=4)
    presentation = mesh.dataset_context.presentation

    assert presentation.headline == "Synthetic Urban Smoke Flow"
    assert presentation.summary.startswith("A deterministic smoke-test dataset")
    assert presentation.legend["title"] == "Smoke speed"
    assert presentation.annotations[0]["label"] == "Fast corridor"
    assert presentation.view_hints["default_plot_mode"] == "preview"
