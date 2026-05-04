from __future__ import annotations

import colorsys
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np


DTCC_COLORS = {
    "background": "#FAFAFA",
    "surface": "#FFFFFF",
    "ink": "#27252A",
    "dark": "#0F0F14",
    "dark_surface": "#1B1B22",
    "dark_panel": "#24242B",
    "muted": "#54545A",
    "muted_light": "#C9C9D1",
    "yellow": "#FADA36",
    "teal": "#78C8BE",
    "orange": "#E35A1D",
    "teal_light": "#C1F1E8",
    "yellow_light": "#F8E36C",
    "grid": "#D7D7DC",
    "grid_dark": "#3A3A42",
}

DTCC_CATEGORY_PALETTE = [
    DTCC_COLORS["teal"],
    DTCC_COLORS["yellow"],
    DTCC_COLORS["orange"],
    DTCC_COLORS["teal_light"],
    DTCC_COLORS["yellow_light"],
    DTCC_COLORS["muted_light"],
]

DTCC_NUMERIC_PALETTE = [
    DTCC_COLORS["teal"],
    DTCC_COLORS["yellow"],
    DTCC_COLORS["orange"],
]

DTCC_THEMES = {
    "light": {
        "figure": DTCC_COLORS["background"],
        "axes": DTCC_COLORS["background"],
        "panel": DTCC_COLORS["surface"],
        "panel_edge": DTCC_COLORS["grid"],
        "text": DTCC_COLORS["ink"],
        "muted": DTCC_COLORS["muted"],
        "grid": DTCC_COLORS["grid"],
    },
    "dark": {
        "figure": DTCC_COLORS["dark"],
        "axes": "#101016",
        "panel": DTCC_COLORS["dark_panel"],
        "panel_edge": DTCC_COLORS["grid_dark"],
        "text": DTCC_COLORS["surface"],
        "muted": DTCC_COLORS["muted_light"],
        "grid": DTCC_COLORS["grid_dark"],
    },
}


def require_matplotlib(feature: str = "plotting"):
    """Import matplotlib lazily and raise a helpful error if unavailable."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(f"Matplotlib is required for {feature}.") from exc
    return plt


def get_theme(theme: str | Mapping[str, str] | None = "dark") -> dict[str, str]:
    """Return a DTCC plotting theme."""
    if theme is None:
        return dict(DTCC_THEMES["dark"])
    if isinstance(theme, Mapping):
        base = dict(DTCC_THEMES["dark"])
        base.update(theme)
        return base
    return dict(DTCC_THEMES.get(theme, DTCC_THEMES["dark"]))


def dtcc_colormap(name: str = "dtcc"):
    """Return the DTCC sequential colormap used for numeric plot values."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name, DTCC_NUMERIC_PALETTE)


def resolve_colormap(cmap=None):
    """Return the requested colormap, defaulting to the DTCC brand colormap."""
    if cmap is None or cmap == "dtcc":
        return dtcc_colormap()
    return cmap


def get_axes(ax=None, *, figsize: tuple[float, float] = (8, 8)):
    """Return an existing axes or create a DTCC-styled matplotlib axes."""
    plt = require_matplotlib()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


def apply_dtcc_style(
    ax,
    *,
    theme: str | Mapping[str, str] | None = "dark",
    equal_aspect: bool = False,
    axis: str = "on",
    grid: bool = False,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    facecolor: str | None = None,
):
    """Apply the DTCC visual style to an axes."""
    theme_values = get_theme(theme)
    ax.figure.patch.set_facecolor(theme_values["figure"])
    ax.set_facecolor(facecolor or theme_values["axes"])

    if title is not None:
        ax.set_title(title, color=theme_values["text"], fontweight=600)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    ax.title.set_color(theme_values["text"])
    ax.xaxis.label.set_color(theme_values["text"])
    ax.yaxis.label.set_color(theme_values["text"])
    ax.tick_params(colors=theme_values["muted"], labelsize=9)

    for spine in ax.spines.values():
        spine.set_color(theme_values["grid"])
        spine.set_linewidth(0.8)

    if grid:
        ax.grid(
            True,
            color=theme_values["grid"],
            linewidth=0.7,
            alpha=0.7,
        )
    else:
        ax.grid(False)

    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")
    if axis == "off":
        ax.set_axis_off()
    return ax


def show_plot(show: bool):
    """Show the current matplotlib figure when requested."""
    if show:
        require_matplotlib().show()


def format_bounds(bounds: Any, *, precision: int = 2) -> list[str]:
    """Format a 2D/3D bounding box as compact metadata lines."""
    if bounds is None:
        return []

    if all(hasattr(bounds, attr) for attr in ("xmin", "ymin", "xmax", "ymax")):
        xmin = bounds.xmin
        ymin = bounds.ymin
        xmax = bounds.xmax
        ymax = bounds.ymax
        zmin = getattr(bounds, "zmin", None)
        zmax = getattr(bounds, "zmax", None)
    else:
        try:
            values = list(bounds)
        except TypeError:
            return []
        if len(values) < 4:
            return []
        xmin, ymin, xmax, ymax = values[:4]
        zmin = values[4] if len(values) > 4 else None
        zmax = values[5] if len(values) > 5 else None

    lines = [
        "Bounds:",
        f"x {_format_number(xmin, precision)} -> {_format_number(xmax, precision)}",
        f"y {_format_number(ymin, precision)} -> {_format_number(ymax, precision)}",
    ]
    if zmin is not None and zmax is not None and (
        float(zmin) != 0.0 or float(zmax) != 0.0
    ):
        lines.append(
            f"z {_format_number(zmin, precision)} -> {_format_number(zmax, precision)}"
        )
    return lines


def add_metadata_box(
    ax,
    metadata: Mapping[str, Any] | Sequence[tuple[str, Any] | str] | None = None,
    *,
    bounds: Any = None,
    theme: str | Mapping[str, str] | None = "dark",
    loc: str = "upper left",
    fontsize: int = 8,
):
    """Add a compact metadata panel to a plot."""
    lines = _metadata_lines(metadata)
    lines.extend(format_bounds(bounds))
    if not lines:
        return None

    theme_values = get_theme(theme)
    x, y, va, ha = _metadata_location(loc)
    return ax.text(
        x,
        y,
        "\n".join(lines),
        transform=ax.transAxes,
        va=va,
        ha=ha,
        fontsize=fontsize,
        color=theme_values["text"],
        zorder=20,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": theme_values["panel"],
            "alpha": 0.9,
            "edgecolor": theme_values["panel_edge"],
            "linewidth": 0.7,
        },
    )


def add_plot_context(
    ax,
    *,
    title: str | None = None,
    metadata: Mapping[str, Any] | Sequence[tuple[str, Any] | str] | None = None,
    bounds: Any = None,
    theme: str | Mapping[str, str] | None = "dark",
    metadata_loc: str = "upper left",
):
    """Add a title and compact metadata panel to an axes."""
    if title:
        theme_values = get_theme(theme)
        ax.set_title(title, color=theme_values["text"], fontweight=600, pad=10)
    return add_metadata_box(
        ax,
        metadata,
        bounds=bounds,
        theme=theme,
        loc=metadata_loc,
    )


def as_numeric_values(values: Sequence[Any]) -> np.ndarray | None:
    """Return values as floats, or None when any non-missing value is categorical."""
    numeric_values = []
    for value in values:
        if _is_missing(value):
            numeric_values.append(np.nan)
            continue
        try:
            numeric_values.append(float(value))
        except (TypeError, ValueError):
            return None
    return np.asarray(numeric_values, dtype=float)


def categorical_labels(
    values: Sequence[Any],
    *,
    missing_label: str = "N/A",
) -> list[str]:
    """Return stable string labels for categorical plot values."""
    return [missing_label if _is_missing(value) else str(value) for value in values]


def category_color_lookup(
    categories: Sequence[str],
    *,
    palette: Sequence[str] | None = None,
    cmap=None,
) -> dict[str, Any]:
    """Return a category-to-color lookup using DTCC colors by default."""
    categories = list(categories)
    if cmap is not None:
        plt = require_matplotlib()
        colormap = plt.get_cmap(cmap, max(len(categories), 1))
        return {category: colormap(i) for i, category in enumerate(categories)}

    colors = list(palette or DTCC_CATEGORY_PALETTE)
    return {
        category: colors[index % len(colors)]
        for index, category in enumerate(categories)
    }


def palette_color(
    index: int,
    *,
    palette: Sequence[str] | None = None,
    lightness_variants: Sequence[float] = (0.0, 0.10, -0.08, 0.16, -0.14),
) -> tuple[float, float, float]:
    """Return a palette color with repeat-cycle lightness variation."""
    colors = list(palette or DTCC_CATEGORY_PALETTE)
    base_color = colors[index % len(colors)]
    delta = lightness_variants[(index // len(colors)) % len(lightness_variants)]
    return adjust_color_lightness(base_color, delta)


def adjust_color_lightness(color: str, delta: float) -> tuple[float, float, float]:
    """Adjust a hex color's lightness in HLS space."""
    red, green, blue = tuple(
        int(color[index : index + 2], 16) / 255.0 for index in (1, 3, 5)
    )
    hue, lightness, saturation = colorsys.rgb_to_hls(red, green, blue)
    lightness = min(0.85, max(0.28, lightness + delta))
    return colorsys.hls_to_rgb(hue, lightness, saturation)


def add_categorical_legend(
    ax,
    color_lookup: dict[str, Any],
    *,
    title: str | None = None,
    linewidth: float = 2.0,
    marker: str = "line",
    theme: str | Mapping[str, str] | None = "dark",
):
    """Add a compact DTCC-styled categorical legend to an axes."""
    theme_values = get_theme(theme)
    if marker == "patch":
        from matplotlib.patches import Patch

        handles = [
            Patch(
                facecolor=color,
                edgecolor=theme_values["text"],
                linewidth=0.6,
                label=category,
            )
            for category, color in color_lookup.items()
        ]
    else:
        from matplotlib.lines import Line2D

        handles = [
            Line2D([0], [0], color=color, lw=linewidth, label=category)
            for category, color in color_lookup.items()
        ]

    legend = ax.legend(handles=handles, title=title, frameon=False)
    if legend.get_title() is not None:
        legend.get_title().set_color(theme_values["text"])
    for text in legend.get_texts():
        text.set_color(theme_values["text"])
    return legend


def plot_line_segments(
    ax,
    segments: Sequence[Any],
    *,
    values: Sequence[Any] | None = None,
    column_label: str | None = None,
    color: str | None = None,
    linewidth: float = 1.0,
    cmap=None,
    legend: bool = True,
    theme: str | Mapping[str, str] | None = "dark",
    **kwargs,
):
    """Plot 2D line segments with shared DTCC color and legend handling."""
    from matplotlib.collections import LineCollection

    if values is None:
        collection = LineCollection(
            segments,
            colors=color or DTCC_COLORS["teal"],
            linewidths=linewidth,
            **kwargs,
        )
        ax.add_collection(collection)
        return collection

    numeric_values = as_numeric_values(values)
    if numeric_values is not None:
        collection = LineCollection(
            segments,
            linewidths=linewidth,
            cmap=resolve_colormap(cmap),
            **kwargs,
        )
        collection.set_array(numeric_values)
        ax.add_collection(collection)
        if legend:
            colorbar = ax.figure.colorbar(collection, ax=ax, label=column_label)
            theme_values = get_theme(theme)
            colorbar.ax.set_facecolor(theme_values["figure"])
            colorbar.ax.yaxis.label.set_color(theme_values["text"])
            colorbar.ax.tick_params(colors=theme_values["muted"], labelsize=9)
            colorbar.outline.set_edgecolor(theme_values["grid"])
        return collection

    labels = categorical_labels(values)
    categories = sorted(set(labels))
    color_lookup = category_color_lookup(categories, cmap=cmap)
    colors = [color_lookup[label] for label in labels]
    collection = LineCollection(
        segments,
        colors=colors,
        linewidths=linewidth,
        **kwargs,
    )
    ax.add_collection(collection)
    if legend:
        add_categorical_legend(
            ax,
            color_lookup,
            title=column_label,
            linewidth=linewidth,
            theme=theme,
        )
    return collection


def plot_geodataframe(
    ax,
    gdf,
    *,
    column: str | None = None,
    edgecolor: str | None = None,
    linewidth: float = 0.6,
    facecolor: str | None = None,
    cmap=None,
    legend: bool = False,
    theme: str | Mapping[str, str] | None = "dark",
    **kwargs,
):
    """Plot a GeoDataFrame with the shared DTCC matplotlib styling."""
    theme_values = get_theme(theme)
    plot_kwargs = {
        "edgecolor": edgecolor or theme_values["text"],
        "linewidth": linewidth,
        **kwargs,
    }
    if column is not None and column in gdf.columns:
        gdf.plot(
            column=column,
            ax=ax,
            legend=legend,
            cmap=resolve_colormap(cmap),
            **plot_kwargs,
        )
    else:
        gdf.plot(
            ax=ax,
            facecolor=facecolor if facecolor is not None else "none",
            **plot_kwargs,
        )
    if legend:
        style_plot_extras(ax, theme=theme)
    return ax


def style_plot_extras(ax, *, theme: str | Mapping[str, str] | None = "dark"):
    """Style legends and colorbar axes created by plotting backends."""
    theme_values = get_theme(theme)
    legend = ax.get_legend()
    if legend is not None:
        legend.get_frame().set_facecolor(theme_values["panel"])
        legend.get_frame().set_edgecolor(theme_values["panel_edge"])
        legend.get_frame().set_alpha(0.9)
        if legend.get_title() is not None:
            legend.get_title().set_color(theme_values["text"])
        for text in legend.get_texts():
            text.set_color(theme_values["text"])

    for extra_ax in ax.figure.axes:
        if extra_ax is ax:
            continue
        extra_ax.set_facecolor(theme_values["figure"])
        extra_ax.tick_params(colors=theme_values["muted"], labelsize=9)
        extra_ax.xaxis.label.set_color(theme_values["text"])
        extra_ax.yaxis.label.set_color(theme_values["text"])
        for spine in extra_ax.spines.values():
            spine.set_color(theme_values["grid"])


def polygon_parts(geometry: Any | None) -> list[Any]:
    """Return non-empty Polygon parts from a Shapely polygonal geometry."""
    from shapely.geometry import GeometryCollection, MultiPolygon, Polygon

    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return [polygon for polygon in geometry.geoms if not polygon.is_empty]
    if isinstance(geometry, GeometryCollection):
        polygons: list[Any] = []
        for part in geometry.geoms:
            polygons.extend(polygon_parts(part))
        return polygons
    return []


def polygon_patch(polygon, *, path_cls=None, patch_cls=None, **kwargs):
    """Create a matplotlib PathPatch for a Shapely polygon, including holes."""
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    from shapely.geometry.polygon import orient

    path_cls = path_cls or Path
    patch_cls = patch_cls or PathPatch
    polygon = orient(polygon, sign=1.0)
    vertices: list[tuple[float, float]] = []
    codes: list[int] = []
    shell_vertices, shell_codes = _ring_path(path_cls, polygon.exterior.coords)
    vertices.extend(shell_vertices)
    codes.extend(shell_codes)
    for hole in polygon.interiors:
        hole_vertices, hole_codes = _ring_path(path_cls, hole.coords)
        vertices.extend(hole_vertices)
        codes.extend(hole_codes)
    path = path_cls(vertices, codes)
    return patch_cls(path, **kwargs)


def plot_polygon_geometries(
    ax,
    geometries: Sequence[Any],
    *,
    edgecolor: str | None = None,
    linewidth: float = 0.8,
    alpha: float = 0.9,
    palette: Sequence[str] | None = None,
    facecolors: Sequence[Any] | None = None,
    **kwargs,
) -> list[Any]:
    """Plot Shapely polygonal geometries as matplotlib patches."""
    plotted_parts = []
    for geometry_index, geometry in enumerate(geometries):
        facecolor = (
            facecolors[geometry_index]
            if facecolors is not None
            else palette_color(geometry_index, palette=palette)
        )
        for polygon in polygon_parts(geometry):
            patch = polygon_patch(
                polygon,
                facecolor=facecolor,
                edgecolor=edgecolor or DTCC_COLORS["ink"],
                linewidth=linewidth,
                alpha=alpha,
                **kwargs,
            )
            ax.add_patch(patch)
            plotted_parts.append(polygon)
    return plotted_parts


def set_axes_extent(
    axes,
    geometries: Sequence[Any],
    *,
    padding_fraction: float = 0.05,
    equal_aspect: bool = True,
    hide_ticks: bool = True,
) -> None:
    """Set matching axes extents around a sequence of Shapely geometries."""
    bounds = geometry_bounds(geometries)
    if bounds is None:
        return

    minx, miny, maxx, maxy = bounds
    span_x = max(maxx - minx, 1.0)
    span_y = max(maxy - miny, 1.0)
    padding = padding_fraction * max(span_x, span_y)

    for ax in axes:
        ax.set_xlim(minx - padding, maxx + padding)
        ax.set_ylim(miny - padding, maxy + padding)
        if equal_aspect:
            ax.set_aspect("equal", adjustable="box")
        if hide_ticks:
            ax.set_xticks([])
            ax.set_yticks([])


def geometry_bounds(geometries: Sequence[Any]) -> tuple[float, float, float, float] | None:
    """Return the 2D bounds of Shapely polygonal geometries."""
    polygons: list[Any] = []
    for geometry in geometries:
        polygons.extend(polygon_parts(geometry))
    if not polygons:
        return None
    return (
        min(polygon.bounds[0] for polygon in polygons),
        min(polygon.bounds[1] for polygon in polygons),
        max(polygon.bounds[2] for polygon in polygons),
        max(polygon.bounds[3] for polygon in polygons),
    )


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    try:
        return bool(np.isscalar(value) and np.isnan(value))
    except (TypeError, ValueError):
        return False


def _metadata_lines(
    metadata: Mapping[str, Any] | Sequence[tuple[str, Any] | str] | None,
) -> list[str]:
    if metadata is None:
        return []
    if isinstance(metadata, Mapping):
        items = metadata.items()
    else:
        items = metadata

    lines: list[str] = []
    for item in items:
        if isinstance(item, str):
            if item:
                lines.append(item)
            continue
        key, value = item
        if value is None or value == "":
            continue
        lines.append(f"{key}: {value}")
    return lines


def _metadata_location(loc: str) -> tuple[float, float, str, str]:
    locations = {
        "upper left": (0.02, 0.98, "top", "left"),
        "upper right": (0.98, 0.98, "top", "right"),
        "lower left": (0.02, 0.02, "bottom", "left"),
        "lower right": (0.98, 0.02, "bottom", "right"),
    }
    return locations.get(loc, locations["upper left"])


def _format_number(value: Any, precision: int) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{number:.{precision}f}"


def _ring_path(path_cls: type[Any], coords) -> tuple[list[tuple[float, float]], list[int]]:
    points = np.asarray(coords, dtype=float)
    if points.ndim != 2 or points.shape[0] < 3:
        return [], []
    xy = points[:, :2]
    if not np.allclose(xy[0], xy[-1]):
        xy = np.vstack([xy, xy[0]])
    vertices = [tuple(xy[0])]
    codes = [path_cls.MOVETO]
    for point in xy[1:-1]:
        vertices.append(tuple(point))
        codes.append(path_cls.LINETO)
    vertices.append(tuple(xy[-1]))
    codes.append(path_cls.CLOSEPOLY)
    return vertices, codes


__all__ = [
    "DTCC_COLORS",
    "DTCC_CATEGORY_PALETTE",
    "DTCC_NUMERIC_PALETTE",
    "DTCC_THEMES",
    "add_categorical_legend",
    "add_metadata_box",
    "add_plot_context",
    "adjust_color_lightness",
    "apply_dtcc_style",
    "as_numeric_values",
    "categorical_labels",
    "category_color_lookup",
    "dtcc_colormap",
    "format_bounds",
    "get_axes",
    "get_theme",
    "geometry_bounds",
    "palette_color",
    "plot_geodataframe",
    "plot_line_segments",
    "plot_polygon_geometries",
    "polygon_parts",
    "polygon_patch",
    "require_matplotlib",
    "resolve_colormap",
    "set_axes_extent",
    "show_plot",
    "style_plot_extras",
]
