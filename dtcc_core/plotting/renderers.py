"""Matplotlib renderers for DTCC plotting products."""

from __future__ import annotations

from io import BytesIO
import os
from pathlib import Path
import tempfile
from typing import Any

import numpy as np

from .style import (
    apply_dtcc_style,
    get_theme,
    require_matplotlib,
    resolve_colormap,
    show_plot,
    style_colorbar,
    style_plot_extras,
)

from .options import RasterRenderOptions
from .products import SliceProduct, StreamlineProduct


def render_product_png(product: Any, options: RasterRenderOptions) -> bytes:
    """Render a visualization product as PNG bytes."""
    _ensure_mpl_config_dir()
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    theme = get_theme(options.theme)
    background = _background_color(options, theme)

    fig = Figure(
        figsize=options.figsize,
        dpi=options.dpi,
        facecolor="none" if options.transparent else background,
    )
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("none" if options.transparent else background)

    draw_product(ax, product, options)
    _configure_axes(ax, product, options)

    buffer = BytesIO()
    fig.savefig(
        buffer,
        format="png",
        dpi=options.dpi,
        transparent=options.transparent,
        facecolor=fig.get_facecolor(),
        edgecolor="none",
    )
    return buffer.getvalue()


def plot_product(
    product: Any,
    options: RasterRenderOptions,
    *,
    ax=None,
    show: bool = True,
):
    """Plot a visualization product into a live Matplotlib axes."""
    plt = require_matplotlib("visualization plotting")
    if ax is None:
        _, ax = plt.subplots(figsize=options.figsize)
    theme = get_theme(options.theme)
    background = _background_color(options, theme)
    ax.figure.patch.set_facecolor("none" if options.transparent else background)
    ax.set_facecolor("none" if options.transparent else background)
    draw_product(ax, product, options)
    _configure_axes(ax, product, options)
    show_plot(show)
    return ax


def draw_product(ax, product: Any, options: RasterRenderOptions):
    """Draw a visualization product into an existing Matplotlib axes."""
    if isinstance(product, SliceProduct):
        return _draw_slice(ax, product, options)
    if isinstance(product, StreamlineProduct):
        return _draw_streamlines(ax, product, options)
    raise TypeError(f"Unsupported visualization product: {type(product).__name__}")


def _draw_slice(ax, product: SliceProduct, options: RasterRenderOptions):
    image = product.image
    artist = ax.imshow(
        image,
        extent=product.extent,
        origin="lower",
        cmap=resolve_colormap(options.cmap),
        vmin=options.vmin,
        vmax=options.vmax,
        interpolation=options.interpolation,
        aspect="auto",
    )
    if options.legend:
        label = product.field_name
        if product.field_unit:
            label = f"{label} [{product.field_unit}]"
        colorbar = ax.figure.colorbar(artist, ax=ax, label=label)
        style_colorbar(colorbar, theme=options.theme)
    return artist


def _draw_streamlines(ax, product: StreamlineProduct, options: RasterRenderOptions):
    from matplotlib.collections import LineCollection

    segments: list[np.ndarray] = []
    values: list[float] = []
    for line_index, line in enumerate(product.lines):
        if len(line) < 2:
            continue
        xy = np.asarray(line[:, product.axes], dtype=float)
        line_segments = np.stack((xy[:-1], xy[1:]), axis=1)
        segments.extend(line_segments)

        if product.line_values is not None and line_index < len(product.line_values):
            line_values = np.asarray(product.line_values[line_index], dtype=float)
            if len(line_values) == len(line):
                values.extend(((line_values[:-1] + line_values[1:]) * 0.5).tolist())

    if not segments:
        return None

    if options.glow:
        glow = LineCollection(
            segments,
            colors=options.line_color,
            linewidths=options.line_width * 4.0,
            alpha=0.16,
            capstyle="round",
            joinstyle="round",
            zorder=2,
        )
        ax.add_collection(glow)

    if options.streamline_color_by == product.value_name and len(values) == len(segments):
        collection = LineCollection(
            segments,
            linewidths=options.line_width,
            cmap=resolve_colormap(options.cmap),
            alpha=options.line_alpha,
            capstyle="round",
            joinstyle="round",
            zorder=3,
        )
        collection.set_array(np.asarray(values, dtype=float))
        if options.vmin is not None or options.vmax is not None:
            collection.set_clim(options.vmin, options.vmax)
        ax.add_collection(collection)
        if options.legend:
            label = product.value_name
            if product.value_unit:
                label = f"{label} [{product.value_unit}]"
            colorbar = ax.figure.colorbar(collection, ax=ax, label=label)
            style_colorbar(colorbar, theme=options.theme)
        return collection

    collection = LineCollection(
        segments,
        colors=options.line_color,
        linewidths=options.line_width,
        alpha=options.line_alpha,
        capstyle="round",
        joinstyle="round",
        zorder=3,
    )
    ax.add_collection(collection)
    return collection


def _configure_axes(ax, product: Any, options: RasterRenderOptions) -> None:
    extent = product.extent
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    table_profile = options.profile == "table"
    theme = get_theme(options.theme)
    background = _background_color(options, theme)
    figure_background = "none" if options.transparent else background
    apply_dtcc_style(
        ax,
        theme=options.theme,
        axis="off" if table_profile else "on",
        xlabel=None if table_profile else product.axes_names[0],
        ylabel=None if table_profile else product.axes_names[1],
        title=None if table_profile else options.title,
        facecolor="none" if options.transparent else background,
        figure_facecolor=figure_background,
    )
    ax.set_aspect("equal" if options.preserve_aspect else "auto", adjustable="box")

    if table_profile:
        return
    style_plot_extras(ax, theme=options.theme)


def _background_color(options: RasterRenderOptions, theme: dict[str, str]) -> str:
    if options.background is not None:
        return options.background
    if options.profile == "table":
        return theme["axes"]
    return theme["figure"]


def _ensure_mpl_config_dir() -> None:
    if os.environ.get("MPLCONFIGDIR"):
        return
    path = Path(tempfile.gettempdir()) / "dtcc-matplotlib"
    path.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(path)
