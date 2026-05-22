"""Matplotlib helpers for inspecting footprint conditioning results."""

from __future__ import annotations

from collections.abc import Sequence

from shapely.geometry.base import BaseGeometry

from dtcc_core.plotting.style import (
    add_metadata_box,
    apply_dtcc_style,
    geometry_bounds,
    get_theme,
    plot_polygon_geometries,
    require_matplotlib,
    set_axes_extent,
)


def plot_footprint_cleaning_comparison(
    raw_polygons: Sequence[BaseGeometry],
    cleaned_polygons: Sequence[BaseGeometry],
    *,
    title: str = "Footprint cleaning",
    raw_title: str = "Input footprints",
    cleaned_title: str = "Conditioned footprints",
    show: bool = True,
    block: bool = True,
    theme: str = "dark",
):
    """Plot raw and conditioned footprint coverages side by side.

    Parameters
    ----------
    raw_polygons : sequence of BaseGeometry
        Input polygonal footprint coverage before conditioning.
    cleaned_polygons : sequence of BaseGeometry
        Output polygonal footprint coverage after conditioning.
    title : str, optional
        Figure title, by default ``"Footprint cleaning"``.
    raw_title : str, optional
        Title of the raw-input panel.
    cleaned_title : str, optional
        Title of the conditioned-output panel.
    show : bool, optional
        Whether to call ``matplotlib.pyplot.show()`` before returning.
    block : bool, optional
        Passed through to ``plt.show(block=...)`` when ``show=True``.

    Returns
    -------
    tuple
        ``(fig, axes)`` from Matplotlib.
    """

    plt = require_matplotlib("footprint cleaning plots")
    theme_values = get_theme(theme)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    for ax in axes:
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=True,
            axis="off",
        )

    raw_parts = plot_polygon_geometries(
        axes[0],
        raw_polygons,
        edgecolor=theme_values["text"],
        linewidth=0.8,
        alpha=0.9,
    )
    cleaned_parts = plot_polygon_geometries(
        axes[1],
        cleaned_polygons,
        edgecolor=theme_values["text"],
        linewidth=0.8,
        alpha=0.92,
    )

    axes[0].set_title(raw_title)
    axes[1].set_title(cleaned_title)
    for ax in axes:
        ax.title.set_color(theme_values["text"])
        ax.title.set_fontweight(600)

    all_parts = [*raw_parts, *cleaned_parts]
    set_axes_extent(axes, all_parts)

    for ax, parts in zip(axes, (raw_parts, cleaned_parts)):
        add_metadata_box(
            ax,
            {"Polygons": len(parts)},
            bounds=geometry_bounds(parts),
            theme=theme,
            loc="lower left",
        )

    fig.suptitle(title, color=theme_values["text"], fontweight=600)
    if show:
        plt.show(block=block)
    return fig, axes
