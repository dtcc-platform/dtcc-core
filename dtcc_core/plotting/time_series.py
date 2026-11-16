"""Matplotlib previews of temporal values."""

import datetime as dt

import numpy as np

from .style import (
    apply_dtcc_style,
    get_axes,
    get_theme,
    palette_color,
    show_plot,
)

_MAX_CURVES = 20
_MAX_VALUES = 1_000_000


def _index(value, size, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be a nonnegative integer index")
    if not 0 <= value < size:
        raise ValueError(f"{name} index {value} is out of range for size {size}")
    return int(value)


def _elements(elements, size):
    if elements is None:
        return range(size)
    if isinstance(elements, (int, np.integer)):
        elements = [elements]
    try:
        selected = [_index(value, size, "element") for value in elements]
    except TypeError as exc:
        raise ValueError(
            "elements must be an integer or a sequence of indices"
        ) from exc
    if len(set(selected)) != len(selected):
        raise ValueError("elements must not contain duplicate indices")
    return selected


def _components(component, dim):
    if component is None:
        return range(dim)
    if isinstance(component, str) and component == "magnitude":
        return ["magnitude"]
    return [_index(component, dim, "component")]


def _curves(series, elements, components):
    for element in elements:
        for component in components:
            values = series.values[:, element]
            if component == "magnitude":
                # Floating point and hypot avoid integer and squared overflow.
                values = values.astype(float)
                values = (
                    np.abs(values)
                    if series.dim == 1
                    else np.hypot.reduce(values, axis=-1)
                )
            elif series.dim > 1:
                values = values[:, component]
            values = np.where(np.isfinite(values), values, np.nan)
            label = series.name or "Value"
            if series.role:
                label += f" ({series.role})"
            if series.num_elements > 1:
                label += f" / element {element}"
            if component == "magnitude":
                label += " / magnitude"
            elif series.dim > 1:
                label += f" / component {component}"
            yield values, label


def _style_axes(ax, series, component, title, theme):
    from matplotlib.dates import AutoDateLocator, ConciseDateFormatter

    locator = AutoDateLocator(tz=dt.UTC)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(ConciseDateFormatter(locator, tz=dt.UTC))
    ylabel = series.name or "Value"
    if component == "magnitude":
        ylabel += " magnitude"
    if series.unit:
        ylabel += f" [{series.unit}]"
    apply_dtcc_style(
        ax,
        theme=theme,
        grid=True,
        xlabel="Time (UTC)",
        ylabel=ylabel,
        title=title if title is not None else series.name or "Time series",
    )
    ax.xaxis.get_offset_text().set_color(get_theme(theme)["muted"])
    ax.yaxis.get_offset_text().set_color(get_theme(theme)["muted"])
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        legend = ax.legend(handles, labels, frameon=False)
        for text in legend.get_texts():
            text.set_color(get_theme(theme)["text"])


def _plot_time_series(series, *, ax, elements, component, title, theme, show, **kwargs):
    selected = _elements(elements, series.num_elements)
    components = _components(component, series.dim)
    count = len(selected) * len(components)
    if not len(series) or not count:
        raise ValueError("No time series values selected for plotting")
    if count > _MAX_CURVES:
        raise ValueError(
            "Time series plots allow at most 20 curves; select elements/components"
        )
    if len(series) * count > _MAX_VALUES:
        raise ValueError(
            "Time series plots allow at most 1,000,000 values; "
            "select elements/components or use between() / resample()"
        )
    if ax is not None and getattr(ax, "name", None) != "rectilinear":
        raise ValueError("Time series plots require a Matplotlib 2D rectilinear axes")
    curves = list(_curves(series, selected, components))
    if not any(np.isfinite(values).any() for values, _ in curves):
        raise ValueError("No finite time series values selected for plotting")

    ax = get_axes(ax, figsize=(10, 5))
    offset = len(ax.lines)
    for index, (values, label) in enumerate(curves):
        options = dict(kwargs)
        options.setdefault("label", label)
        if "color" not in options and "c" not in options:
            options["color"] = palette_color(offset + index)
        if len(series) == 1:
            options.setdefault("marker", "o")
        ax.plot(series.times, values, **options)
    _style_axes(ax, series, component, title, theme)
    show_plot(show)
    return ax
