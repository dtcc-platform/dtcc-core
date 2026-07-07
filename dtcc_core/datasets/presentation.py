"""Human-readable Dataset v2 presentation helpers."""

from __future__ import annotations

from io import StringIO
import json
import textwrap
from typing import Any

from rich.console import Console
from rich.markup import escape

from dtcc_core.common.dtcc_logging import make_table
from dtcc_core.plotting.style import DTCC_COLORS


MISSING = "Not specified"
PANEL_BACKGROUND = DTCC_COLORS["dark_surface"]
PANEL_TEXT = DTCC_COLORS["surface"]
PANEL_ACCENT = "#FADA36"
PREVIEW_VISUAL_PANEL_BOUNDS = (0.04, 0.08, 0.62, 0.84)
PREVIEW_STORY_PANEL_BOUNDS = (0.695, 0.08, 0.265, 0.84)
PREVIEW_VISUAL_AXES_BOUNDS = (0.055, 0.16, 0.59, 0.69)
PREVIEW_STORY_AXES_BOUNDS = (0.715, 0.11, 0.225, 0.78)
PREVIEW_COLOR_RAMP_AXES_BOUNDS = (0.085, 0.105, 0.31, 0.028)


def format_dataset_context(context, obj: Any | None = None) -> str:
    """Return a readable metadata and presentation summary for a Dataset v2 object."""
    console = Console(
        file=StringIO(),
        record=True,
        width=120,
        color_system=None,
        force_terminal=False,
        soft_wrap=False,
    )

    console.print("Metadata", style="bold")
    console.print(
        make_table(
            [("Field", "left"), ("Value", "left")],
            _metadata_rows(context, obj),
            overflow="fold",
        )
    )

    console.print()
    console.print("Presentation", style="bold")
    console.print(
        make_table(
            [("Field", "left"), ("Value", "left")],
            _presentation_rows(context),
            overflow="fold",
        )
    )

    console.print()
    console.print("Provenance", style="bold")
    console.print(
        make_table(
            [("Field", "left"), ("Value", "left")],
            _provenance_rows(context),
            overflow="fold",
        )
    )

    return console.export_text(styles=False).rstrip()


def plot_product_with_presentation(
    product,
    options,
    *,
    context,
    obj: Any | None = None,
    ax=None,
    show: bool = True,
):
    """Plot a product with a side panel showing Dataset v2 presentation metadata."""
    from dtcc_core.plotting.renderers import plot_product

    if context is None:
        return plot_product(product, options, ax=ax, show=show)

    data_ax, panel_ax = presentation_plot_axes(
        context,
        ax=ax,
        presentation=True,
        figsize=options.figsize,
    )
    plot_product(product, options, ax=data_ax, show=False)
    return finalize_presentation_plot(
        data_ax,
        context=context,
        obj=obj,
        panel_ax=panel_ax,
        presentation=True,
        show=show,
        theme=options.theme,
        title=getattr(options, "title", None),
        subtitle=_preview_subtitle(context, product),
    )


def presentation_plot_axes(
    context,
    *,
    ax=None,
    presentation: bool = True,
    figsize: tuple[float, float] = (8.0, 8.0),
):
    """Return axes for a Dataset v2 presentation-aware plot."""
    from dtcc_core.plotting.style import get_axes, require_matplotlib

    if not presentation or context is None:
        return get_axes(ax, figsize=figsize), None

    if ax is not None:
        return ax, None

    plt = require_matplotlib("dataset presentation plotting")
    width, height = figsize
    fig = plt.figure(
        figsize=(max(width * 1.55, 10.5), max(height, 5.6)),
        facecolor=DTCC_COLORS["dark"],
    )
    add_preview_panel(
        fig,
        PREVIEW_VISUAL_PANEL_BOUNDS,
        facecolor="#101016",
        edgecolor=DTCC_COLORS["grid_dark"],
        shadow=True,
    )
    add_preview_panel(
        fig,
        PREVIEW_STORY_PANEL_BOUNDS,
        facecolor=DTCC_COLORS["dark_surface"],
        edgecolor=DTCC_COLORS["grid_dark"],
        shadow=True,
    )
    data_ax = fig.add_axes(PREVIEW_VISUAL_AXES_BOUNDS, zorder=2)
    data_ax.set_facecolor("#101016")
    panel_ax = fig.add_axes(PREVIEW_STORY_AXES_BOUNDS, zorder=3)
    return data_ax, panel_ax


def finalize_presentation_plot(
    ax,
    *,
    context,
    obj: Any | None = None,
    panel_ax=None,
    presentation: bool = True,
    show: bool = True,
    theme: str = "dark",
    title: str | None = None,
    subtitle: str | None = None,
):
    """Add Dataset v2 presentation chrome to a plot and optionally show it."""
    from dtcc_core.plotting.style import show_plot

    if presentation and context is not None:
        if panel_ax is not None:
            draw_presentation_panel(panel_ax, context, obj=obj)
        else:
            draw_dataset_preview_header(
                ax,
                context,
                theme=theme,
                title=title,
                subtitle=subtitle,
            )
            ax.figure.subplots_adjust(right=0.72)
            ax.text(
                1.03,
                1.0,
                presentation_panel_text(context, obj=obj),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=8.5,
                color=PANEL_TEXT,
                bbox={
                    "boxstyle": "round,pad=0.45",
                    "facecolor": PANEL_BACKGROUND,
                    "edgecolor": "#d0d0d0",
                    "alpha": 0.96,
                },
                clip_on=False,
            )

    show_plot(show)
    return ax


def draw_presentation_panel(ax, context, obj: Any | None = None) -> None:
    """Draw a compact Dataset v2 presentation panel on a Matplotlib axes."""
    ax.set_facecolor(PANEL_BACKGROUND)
    ax.patch.set_visible(True)
    ax.patch.set_alpha(1.0)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(
        left=False,
        bottom=False,
        labelleft=False,
        labelbottom=False,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    _draw_story_panel(ax, context, obj=obj)


def draw_empty_preview_state(
    ax,
    context,
    *,
    object_label: str | None = None,
    message: str | None = None,
) -> None:
    """Draw an explicit empty-domain message inside a presentation visual panel."""
    label = object_label or _object_label(context, plural=True)
    bounds = _requested_domain_label(context)
    ax.text(
        0.5,
        0.55,
        message or f"No {label} in the requested domain",
        ha="center",
        va="center",
        color=DTCC_COLORS["surface"],
        fontsize=13.0,
        fontweight=700,
        transform=ax.transAxes,
    )
    ax.text(
        0.5,
        0.47,
        bounds,
        ha="center",
        va="center",
        color=DTCC_COLORS["muted_light"],
        fontsize=8.8,
        transform=ax.transAxes,
    )
    ax.text(
        0.5,
        0.405,
        "The selected bounds returned an empty layer for this dataset.",
        ha="center",
        va="center",
        color=DTCC_COLORS["muted_light"],
        fontsize=8.0,
        transform=ax.transAxes,
    )


def _draw_story_panel(ax, context, obj: Any | None = None) -> None:
    presentation = context.presentation
    headline = _wrap_text(
        _first_text(
            presentation.headline,
            context.identity.title,
            context.identity.name,
        ),
        width=24,
        max_lines=2,
    )
    ax.text(
        0.0,
        0.985,
        headline,
        ha="left",
        va="top",
        color=DTCC_COLORS["surface"],
        fontsize=13.5,
        fontweight=700,
        linespacing=1.08,
        transform=ax.transAxes,
    )
    y = 0.985 - _text_block_height(headline, line_height=0.058) - 0.035

    if presentation.summary:
        summary = _wrap_text(presentation.summary, width=41, max_lines=4)
        ax.text(
            0.0,
            y,
            summary,
            ha="left",
            va="top",
            color=DTCC_COLORS["muted_light"],
            fontsize=8.0,
            linespacing=1.28,
            transform=ax.transAxes,
        )
        y -= _text_block_height(summary, line_height=0.031) + 0.045
    else:
        y -= 0.025

    y = min(y, 0.705)
    y = draw_preview_chips(ax, _preview_chips(context), y=y) - 0.045

    facts = _preview_facts(context, obj=obj)
    fact_top = 0.18
    narrative_bottom = fact_top + 0.035
    for item in _narrative_items(context)[:2]:
        if y < narrative_bottom + 0.11:
            break
        body = _wrap_text(item[1], width=43, max_lines=7)
        card_height = min(
            0.24,
            max(0.125, 0.067 + _line_count(body) * 0.024),
        )
        if y - card_height < narrative_bottom:
            remaining = y - narrative_bottom
            max_body_lines = max(2, int((remaining - 0.067) / 0.024))
            body = _wrap_text(item[1], width=43, max_lines=max_body_lines)
            card_height = max(0.11, remaining)
        _draw_narrative_card(
            ax,
            heading=item[0],
            body=body,
            top=y,
            height=card_height,
        )
        y -= card_height + 0.02

    ax.text(
        0.0,
        fact_top,
        "Preview facts",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9.0,
        fontweight=700,
        color=DTCC_COLORS["surface"],
    )
    draw_preview_fact_rows(ax, facts, y=fact_top - 0.05)


def add_preview_panel(
    fig,
    bounds: tuple[float, float, float, float],
    *,
    facecolor: str,
    edgecolor: str,
    shadow: bool,
) -> None:
    from matplotlib.patches import FancyBboxPatch

    x, y, width, height = bounds
    if shadow:
        fig.patches.append(
            FancyBboxPatch(
                (x + 0.006, y - 0.008),
                width,
                height,
                boxstyle="round,pad=0.008,rounding_size=0.026",
                transform=fig.transFigure,
                linewidth=0,
                facecolor="#000000",
                alpha=0.24,
                zorder=0,
            )
        )
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.008,rounding_size=0.026",
            transform=fig.transFigure,
            linewidth=0.9,
            edgecolor=edgecolor,
            facecolor=facecolor,
            zorder=1,
        )
    )


def draw_preview_chips(ax, chips: list[str], *, y: float) -> float:
    if len(chips) <= 3:
        x_positions = (0.0, 0.37, 0.72)
        for x, chip in zip(x_positions, chips):
            ax.text(
                x,
                y,
                chip,
                ha="left",
                va="center",
                color=DTCC_COLORS["dark"],
                fontsize=7.5,
                fontweight=700,
                transform=ax.transAxes,
                bbox={
                    "boxstyle": "round,pad=0.26",
                    "facecolor": DTCC_COLORS["teal"],
                    "edgecolor": "none",
                    "alpha": 0.96,
                },
            )
        return y - 0.035

    x = 0.0
    row_y = y
    for chip in chips:
        width = min(0.34, 0.055 + 0.018 * len(chip))
        if x > 0.0 and x + width > 1.0:
            x = 0.0
            row_y -= 0.052
        ax.text(
            x,
            row_y,
            chip,
            ha="left",
            va="center",
            color=DTCC_COLORS["dark"],
            fontsize=7.5,
            fontweight=700,
            transform=ax.transAxes,
            bbox={
                "boxstyle": "round,pad=0.26",
                "facecolor": DTCC_COLORS["teal"],
                "edgecolor": "none",
                "alpha": 0.96,
            },
        )
        x += width + 0.035
    return row_y - 0.035


def draw_preview_fact_strip(ax, facts: list[str], *, y: float) -> float:
    if not facts:
        return y
    if len(facts) <= 3:
        x_positions = (0.0, 0.37, 0.72)
        for x, fact in zip(x_positions, facts):
            ax.text(
                x,
                y,
                fact,
                ha="left",
                va="center",
                color=DTCC_COLORS["surface"],
                fontsize=7.4,
                transform=ax.transAxes,
                bbox={
                    "boxstyle": "round,pad=0.24",
                    "facecolor": "#101016",
                    "edgecolor": DTCC_COLORS["grid_dark"],
                    "linewidth": 0.7,
                },
            )
        return y - 0.035

    x = 0.0
    row_y = y
    for fact in facts[:4]:
        width = min(0.32, 0.07 + 0.012 * len(fact))
        if x > 0.0 and x + width > 1.0:
            x = 0.0
            row_y -= 0.052
        ax.text(
            x,
            row_y,
            fact,
            ha="left",
            va="center",
            color=DTCC_COLORS["surface"],
            fontsize=7.4,
            transform=ax.transAxes,
            bbox={
                "boxstyle": "round,pad=0.24",
                "facecolor": "#101016",
                "edgecolor": DTCC_COLORS["grid_dark"],
                "linewidth": 0.7,
            },
        )
        x += width + 0.035
    return row_y - 0.035


def draw_preview_fact_rows(
    ax,
    facts: list[tuple[str, str]],
    *,
    y: float,
) -> float:
    row_y = y
    for index, (label, value) in enumerate(facts[:4]):
        column = index % 2
        if index and column == 0:
            row_y -= 0.076
        x = 0.0 if column == 0 else 0.52
        ax.text(
            x,
            row_y,
            label,
            ha="left",
            va="top",
            color=DTCC_COLORS["muted_light"],
            fontsize=6.9,
            fontweight=700,
            transform=ax.transAxes,
        )
        ax.text(
            x,
            row_y - 0.029,
            _wrap_text(value, width=20, max_lines=2),
            ha="left",
            va="top",
            color=DTCC_COLORS["surface"],
            fontsize=7.1,
            transform=ax.transAxes,
        )
    return row_y - 0.076


def _draw_narrative_card(
    ax,
    *,
    heading: str,
    body: str,
    top: float,
    height: float,
) -> None:
    from matplotlib.patches import FancyBboxPatch

    ax.add_patch(
        FancyBboxPatch(
            (0.0, top - height),
            1.0,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.035",
            transform=ax.transAxes,
            facecolor=DTCC_COLORS["dark_panel"],
            edgecolor=DTCC_COLORS["grid_dark"],
            linewidth=0.75,
        )
    )
    ax.text(
        0.045,
        top - 0.022,
        heading,
        ha="left",
        va="top",
        color=DTCC_COLORS["yellow"],
        fontsize=8.5,
        fontweight=700,
        transform=ax.transAxes,
    )
    ax.text(
        0.045,
        top - 0.059,
        body,
        ha="left",
        va="top",
        color=DTCC_COLORS["muted_light"],
        fontsize=7.25,
        linespacing=1.13,
        transform=ax.transAxes,
    )


def _preview_chips(context) -> list[str]:
    metadata = context.metadata
    presentation = context.presentation
    view_hints = presentation.view_hints if isinstance(presentation.view_hints, dict) else {}
    values = [
        view_hints.get("table_role") or view_hints.get("preferred_geometry"),
        _plain_value(metadata.crs),
        metadata.data_category,
    ]
    return [
        str(value).replace("_", " ")
        for value in values
        if value and str(value) != MISSING
    ][:3]


def _narrative_items(context) -> list[tuple[str, str]]:
    presentation = context.presentation
    items: list[tuple[str, str]] = []
    for item in presentation.narrative or ():
        if isinstance(item, dict):
            heading = _plain_value(item.get("heading"))
            body = _plain_value(item.get("body"))
        else:
            heading = "Context"
            body = _plain_value(item)
        if heading != MISSING and body != MISSING:
            items.append((heading, body))

    if items:
        return items

    if presentation.summary:
        items.append(("What you are seeing", presentation.summary))
    if presentation.warnings:
        items.append(("Warnings", _plain_value(presentation.warnings[0])))
    if presentation.limitations:
        items.append(("Limitations", _plain_value(presentation.limitations[0])))
    return items


def _preview_facts(context, obj: Any | None = None) -> list[tuple[str, str]]:
    metadata = context.metadata
    facts: list[tuple[str, str]] = []
    facts.append(("Domain", _requested_domain_label(context)))
    facts.append(("Records", _record_count_label(context, obj)))
    crs = _plain_value(metadata.crs)
    if crs != MISSING:
        facts.append(("CRS", _short_list(crs)))
    formats = _plain_value(metadata.formats)
    if formats != MISSING:
        facts.append(("Formats", _short_list(formats)))
    return facts[:4]


def _short_chip(value: Any, *, max_length: int = 14) -> str:
    text = str(value).replace("_", " ").strip()
    if len(text) <= max_length:
        return text
    return text[: max_length - 3].rstrip() + "..."


def _short_list(value: Any, *, max_length: int = 32) -> str:
    text = str(value).replace("_", " ").strip()
    if len(text) <= max_length:
        return text
    return text[: max_length - 3].rstrip() + "..."


def _requested_domain_label(context) -> str:
    bounds = getattr(context.request, "bounds", None)
    if not bounds or len(bounds) < 4:
        return "Request bounds not specified"
    xmin, ymin, xmax, ymax = (float(value) for value in bounds[:4])
    width = max(xmax - xmin, 0.0)
    height = max(ymax - ymin, 0.0)
    return f"{_format_distance(width)} x {_format_distance(height)}"


def _format_distance(value: float) -> str:
    if value >= 1000.0:
        text = f"{value / 1000.0:.1f}".rstrip("0").rstrip(".")
        return f"{text} km"
    text = f"{value:.0f}" if value >= 10.0 else f"{value:.1f}"
    return f"{text} m"


def _record_count_label(context, obj: Any | None = None) -> str:
    count = _object_count(obj)
    label = _object_label(context, plural=count != 1)
    if count is None:
        return label
    if count == 0:
        return f"0 {label} in bounds"
    return f"{count:,} {label}"


def _object_count(obj: Any | None) -> int | None:
    if obj is None:
        return None
    for attr in ("stations", "vehicles"):
        value = getattr(obj, attr, None)
        if callable(value):
            try:
                return len(value())
            except TypeError:
                pass
    if hasattr(obj, "length"):
        try:
            return len(getattr(obj, "length"))
        except TypeError:
            pass
    if hasattr(obj, "points"):
        try:
            return len(getattr(obj, "points"))
        except TypeError:
            pass
    try:
        return len(obj)
    except TypeError:
        return None


def _object_label(context, *, plural: bool) -> str:
    result_kind = str(context.metadata.result_kind or "").replace("_", " ").lower()
    dataset_name = str(context.identity.name or "").replace("_", " ").lower()
    if "sensor" in result_kind or any(
        token in dataset_name
        for token in ("weather", "hydrology", "ocean", "air quality")
    ):
        return "stations" if plural else "station"
    if "vehicle" in result_kind or "transit" in dataset_name:
        return "vehicles" if plural else "vehicle"
    if "road" in result_kind or "road" in dataset_name or "space syntax" in dataset_name:
        return "segments" if plural else "segment"
    if "tree" in result_kind or "tree" in dataset_name:
        return "trees" if plural else "tree"
    if "footprint" in result_kind or "footprint" in dataset_name:
        return "footprints" if plural else "footprint"
    if "grid" in dataset_name:
        return "lines" if plural else "line"
    if "deso" in dataset_name:
        return "areas" if plural else "area"
    if "streamline" in result_kind:
        return "streamlines" if plural else "streamline"
    return "records" if plural else "record"


def _wrap_text(text: Any, *, width: int, max_lines: int | None = None) -> str:
    lines = textwrap.wrap(str(text), width=width)
    if max_lines is not None and len(lines) > max_lines:
        lines = lines[:max_lines]
        last_line = lines[-1].rstrip()
        lines[-1] = f"{last_line[: max(width - 3, 1)].rstrip()}..."
    return "\n".join(lines)


def _line_count(text: str) -> int:
    return max(text.count("\n") + 1, 1)


def _text_block_height(text: str, *, line_height: float) -> float:
    return _line_count(text) * line_height


def draw_preview_header(ax, context, product, options) -> None:
    """Add a polished title treatment above the data preview axes."""
    draw_dataset_preview_header(
        ax,
        context,
        theme=options.theme,
        title=getattr(options, "title", None),
        subtitle=_preview_subtitle(context, product),
    )


def draw_dataset_preview_header(
    ax,
    context,
    *,
    theme: str = "dark",
    title: str | None = None,
    subtitle: str | None = None,
) -> None:
    """Add the Dataset v2 tangible-table preview header above a data axes."""
    from dtcc_core.plotting.style import get_theme

    theme_values = get_theme(theme)
    ax.set_title("")

    header = _first_text(
        title,
        context.presentation.headline,
        context.identity.title,
        context.identity.name,
    )
    ax.text(
        0.0,
        1.12,
        header,
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=19,
        fontweight=700,
        color=theme_values["text"],
        clip_on=False,
    )
    ax.text(
        0.0,
        1.065,
        subtitle or _dataset_preview_subtitle(context),
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9.5,
        color=theme_values["muted"],
        clip_on=False,
    )
    ax.plot(
        [0.0, 0.18],
        [1.035, 1.035],
        transform=ax.transAxes,
        color=PANEL_ACCENT,
        linewidth=3.2,
        solid_capstyle="round",
        clip_on=False,
    )


def presentation_panel_text(
    context,
    obj: Any | None = None,
    *,
    width: int = 38,
) -> str:
    """Return concise plain text for a plot-side Dataset v2 presentation panel."""
    identity = context.identity
    metadata = context.metadata
    presentation = context.presentation

    lines = [
        _first_text(presentation.headline, identity.title, identity.name),
        "",
    ]
    if presentation.summary:
        lines.extend(_wrap(presentation.summary, width=width))
        lines.append("")

    key_points = list(presentation.key_points or ())
    if key_points:
        lines.append("Key points")
        for item in key_points[:3]:
            lines.extend(_bullet_lines(item, width=width))
        lines.append("")

    lines.append("Metadata")
    for label, value in (
        ("Provider", metadata.provider),
        ("License", metadata.license),
        ("CRS", metadata.crs),
        ("Formats", metadata.formats),
    ):
        lines.extend(_key_value_lines(label, _plain_value(value), width=width))

    fields = _field_names(obj)
    if fields:
        lines.extend(_key_value_lines("Fields", ", ".join(fields), width=width))

    return "\n".join(lines).strip()


def _preview_subtitle(context, product) -> str:
    crs = _plain_value(context.metadata.crs)
    if crs == MISSING:
        crs = "local coordinates"

    product_label = "field preview"
    field_name = None
    if hasattr(product, "field_name"):
        product_label = f"{getattr(product, 'axis', 'plane')}-slice preview"
        field_name = product.field_name
    elif hasattr(product, "value_name"):
        product_label = "streamline preview"
        field_name = product.value_name

    parts = [product_label]
    if field_name:
        parts.append(str(field_name))
    parts.append(crs)
    parts.append("synthetic city-flow fixture")
    return " | ".join(parts)


def _dataset_preview_subtitle(context) -> str:
    metadata = context.metadata
    presentation = context.presentation
    view_hints = presentation.view_hints if isinstance(presentation.view_hints, dict) else {}
    role = (
        view_hints.get("table_role")
        or view_hints.get("preferred_geometry")
        or metadata.result_kind
        or "dataset"
    )
    crs = _plain_value(metadata.crs)
    if crs == MISSING:
        crs = "local coordinates"
    category = metadata.data_category or "dataset"
    return " | ".join(str(item) for item in (role, crs, category) if item)


def _metadata_rows(context, obj: Any | None) -> list[tuple[str, str]]:
    metadata = context.metadata
    provenance = context.provenance
    fields = _field_names(obj)
    rows = [
        ("Description", _value(metadata.description)),
        ("Provider", _value(metadata.provider)),
        ("Source", _value(metadata.source)),
        ("License", _value(metadata.license)),
        ("Collection period", _value(metadata.collection_period)),
        ("CRS", _value(metadata.crs)),
        ("LOD", _value(metadata.lod)),
        ("Data types", _value(metadata.data_types)),
        ("Formats", _value(metadata.formats)),
        ("Geographic coverage", _value(metadata.geographic_coverage)),
        ("Update frequency", _value(metadata.update_frequency)),
        ("Data category", _value(metadata.data_category)),
        ("Result kind", _value(metadata.result_kind)),
        ("Python return type", _value(metadata.python_return_type)),
        ("Fields", _value(fields)),
        (
            "Resource links",
            _value(_resource_links(metadata.source, provenance.sources)),
        ),
    ]
    return [(escape(label), escape(value)) for label, value in rows]


def _presentation_rows(context) -> list[tuple[str, str]]:
    presentation = context.presentation
    rows = [
        ("Headline", _value(presentation.headline)),
        ("Summary", _value(presentation.summary)),
        ("Narrative", _value(presentation.narrative)),
        ("Key points", _value(presentation.key_points)),
        ("Legend", _value(presentation.legend)),
        ("Annotations", _value(presentation.annotations)),
        ("View hints", _value(presentation.view_hints)),
        ("Warnings", _value(presentation.warnings)),
        ("Limitations", _value(presentation.limitations)),
    ]
    return [(escape(label), escape(value)) for label, value in rows]


def _provenance_rows(context) -> list[tuple[str, str]]:
    provenance = context.provenance
    rows = [
        ("Sources", _value(provenance.sources)),
        ("Processing steps", _value(provenance.processing_steps)),
        ("Generated by", _value(provenance.generated_by)),
        ("Generated at", _value(provenance.generated_at)),
        ("Derived from", _value(provenance.derived_from)),
    ]
    return [(escape(label), escape(value)) for label, value in rows]


def _field_names(obj: Any | None) -> list[str]:
    if obj is None:
        return []
    names = getattr(obj, "field_names", None)
    if callable(names):
        names = names()
    if names is None:
        return []
    return [str(name) for name in names if str(name)]


def _resource_links(*sources) -> list[str]:
    links = []
    seen = set()
    for source_group in sources:
        for source in source_group or ():
            if isinstance(source, str) and source.startswith(("http://", "https://")):
                if source not in seen:
                    links.append(source)
                    seen.add(source)
            elif isinstance(source, dict):
                for key in ("url", "href", "link", "resource"):
                    value = source.get(key)
                    if isinstance(value, str) and value.startswith(
                        ("http://", "https://")
                    ):
                        if value not in seen:
                            links.append(value)
                            seen.add(value)
    return links


def _value(value: Any) -> str:
    if value is None:
        return MISSING
    if isinstance(value, str):
        return value.strip() or MISSING
    if isinstance(value, (list, tuple, set)):
        items = [_plain_value(item) for item in value if _plain_value(item) != MISSING]
        return "\n".join(items) if items else MISSING
    return _plain_value(value)


def _plain_value(value: Any) -> str:
    if value is None:
        return MISSING
    if isinstance(value, str):
        return value.strip() or MISSING
    if isinstance(value, dict):
        name = value.get("name")
        role = value.get("role")
        url = value.get("url") or value.get("href") or value.get("link")
        if name and role:
            return f"{name} ({role})"
        if name:
            return str(name)
        if url:
            return str(url)
    if isinstance(value, (list, tuple, set)):
        items = [_plain_value(item) for item in value if _plain_value(item) != MISSING]
        return ", ".join(items) if items else MISSING
    try:
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    except TypeError:
        return str(value)


def _first_text(*values: Any) -> str:
    for value in values:
        text = _plain_value(value)
        if text != MISSING:
            return text
    return MISSING


def _wrap(text: str, *, width: int) -> list[str]:
    return textwrap.wrap(text, width=width) or [text]


def _bullet_lines(value: Any, *, width: int) -> list[str]:
    wrapped = _wrap(_plain_value(value), width=max(width - 2, 10))
    return [f"- {wrapped[0]}", *[f"  {line}" for line in wrapped[1:]]]


def _key_value_lines(label: str, value: str, *, width: int) -> list[str]:
    prefix = f"{label}: "
    available = max(width - len(prefix), 10)
    wrapped = _wrap(value, width=available)
    return [
        f"{prefix}{wrapped[0]}",
        *[f"{' ' * len(prefix)}{line}" for line in wrapped[1:]],
    ]
