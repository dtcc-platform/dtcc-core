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
    show_changes: bool = False,
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
    show_changes : bool, optional
        Add a third panel showing added and removed coverage. Invalid input
        polygons use ``make_valid`` only for this difference calculation.

    Returns
    -------
    tuple
        ``(fig, axes)`` from Matplotlib.
    """

    plt = require_matplotlib("footprint cleaning plots")
    theme_values = get_theme(theme)
    columns = 3 if show_changes else 2
    fig, axes = plt.subplots(
        1,
        columns,
        figsize=(6 * columns, 6),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )
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
    if show_changes:
        from matplotlib.patches import Patch
        from shapely import make_valid
        from shapely.ops import unary_union

        raw_union = unary_union([make_valid(p) for p in raw_polygons])
        cleaned_union = unary_union(cleaned_polygons)
        added = cleaned_union.difference(raw_union)
        removed = raw_union.difference(cleaned_union)
        plot_polygon_geometries(
            axes[2],
            [raw_union.union(cleaned_union)],
            facecolors=["#7d8590"],
            edgecolor="#7d8590",
            alpha=0.25,
            linewidth=0.3,
        )
        legend = []
        for label, geometry, color in (
            ("Added", added, "#2ecc71"),
            ("Removed", removed, "#ff6b6b"),
        ):
            plot_polygon_geometries(
                axes[2],
                [geometry],
                facecolors=[color],
                edgecolor=color,
                linewidth=0.5,
                alpha=1.0,
            )
            legend.append(
                Patch(facecolor=color, label=f"{label}: {geometry.area:,.2f} m²")
            )
        axes[2].set_title("Coverage changes")
        axes[2].legend(handles=legend, loc="lower left")
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
