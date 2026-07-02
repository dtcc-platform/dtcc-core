"""Human-readable Dataset v2 presentation helpers."""

from __future__ import annotations

from io import StringIO
import json
import textwrap
from typing import Any

from rich.console import Console
from rich.markup import escape

from dtcc_core.common.dtcc_logging import make_table


MISSING = "Not specified"
PANEL_BACKGROUND = "#F7F7F7"
PANEL_TEXT = "#111111"
PANEL_ACCENT = "#FADA36"


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
    from dtcc_core.plotting.style import require_matplotlib, show_plot

    if context is None:
        return plot_product(product, options, ax=ax, show=show)

    if ax is None:
        plt = require_matplotlib("dataset presentation plotting")
        width, height = options.figsize
        fig_width = max(width * 1.55, 10.5)
        fig_height = max(height, 5.6)
        fig, (data_ax, panel_ax) = plt.subplots(
            1,
            2,
            figsize=(fig_width, fig_height),
            gridspec_kw={"width_ratios": [2.85, 1.35]},
        )
        fig.subplots_adjust(
            left=0.06,
            right=0.97,
            top=0.82,
            bottom=0.09,
            wspace=0.12,
        )
        plot_product(product, options, ax=data_ax, show=False)
        draw_preview_header(data_ax, context, product, options)
        draw_presentation_panel(panel_ax, context, obj=obj)
        show_plot(show)
        return data_ax

    data_ax = plot_product(product, options, ax=ax, show=False)
    data_ax.figure.subplots_adjust(right=0.72)
    draw_preview_header(data_ax, context, product, options)
    data_ax.text(
        1.03,
        1.0,
        presentation_panel_text(context, obj=obj),
        transform=data_ax.transAxes,
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
    return data_ax


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
    ax.text(
        0.03,
        0.97,
        presentation_panel_text(context, obj=obj),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.0,
        color=PANEL_TEXT,
        wrap=True,
    )


def draw_preview_header(ax, context, product, options) -> None:
    """Add a polished title treatment above the data preview axes."""
    from dtcc_core.plotting.style import get_theme

    theme = get_theme(options.theme)
    ax.set_title("")

    title = _first_text(
        getattr(options, "title", None),
        context.presentation.headline,
        context.identity.title,
        context.identity.name,
    )
    ax.text(
        0.0,
        1.12,
        title,
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=19,
        fontweight=700,
        color=theme["text"],
        clip_on=False,
    )
    ax.text(
        0.0,
        1.065,
        _preview_subtitle(context, product),
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9.5,
        color=theme["muted"],
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
