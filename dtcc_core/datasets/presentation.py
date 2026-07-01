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

    console.print("Tangible Table Metadata", style="bold")
    console.print(
        make_table(
            [("Criterion", "left"), ("Value", "left")],
            _criteria_rows(context, obj),
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
        fig, (data_ax, panel_ax) = plt.subplots(
            1,
            2,
            figsize=(width * 1.42, height),
            gridspec_kw={"width_ratios": [3.0, 1.15]},
            constrained_layout=True,
        )
        plot_product(product, options, ax=data_ax, show=False)
        draw_presentation_panel(panel_ax, context, obj=obj)
        show_plot(show)
        return data_ax

    data_ax = plot_product(product, options, ax=ax, show=False)
    data_ax.figure.subplots_adjust(right=0.72)
    data_ax.text(
        1.03,
        1.0,
        presentation_panel_text(context, obj=obj),
        transform=data_ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        color="#111111",
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "#f7f7f7",
            "edgecolor": "#d0d0d0",
            "alpha": 0.96,
        },
        clip_on=False,
    )
    show_plot(show)
    return data_ax


def draw_presentation_panel(ax, context, obj: Any | None = None) -> None:
    """Draw a compact Dataset v2 presentation panel on a Matplotlib axes."""
    ax.set_axis_off()
    ax.set_facecolor("#f7f7f7")
    ax.text(
        0.03,
        0.97,
        presentation_panel_text(context, obj=obj),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        color="#111111",
        wrap=True,
    )


def presentation_panel_text(context, obj: Any | None = None) -> str:
    """Return concise plain text for a plot-side Dataset v2 presentation panel."""
    identity = context.identity
    metadata = context.metadata
    presentation = context.presentation
    provenance = context.provenance

    lines = [
        _first_text(presentation.headline, identity.title, identity.name),
        "",
    ]
    if presentation.summary:
        lines.extend(_wrap(presentation.summary, width=34))
        lines.append("")

    key_points = list(presentation.key_points or ())
    if key_points:
        lines.append("Key points")
        lines.extend(f"- {line}" for item in key_points for line in _wrap(item, width=32))
        lines.append("")

    lines.extend(
        [
            f"Provider: {_plain_value(metadata.provider)}",
            f"Source: {_plain_value(metadata.source)}",
            f"License: {_plain_value(metadata.license)}",
            f"CRS: {_plain_value(metadata.crs)}",
            f"Type: {_data_type_value(metadata)}",
            f"Formats: {_plain_value(metadata.formats)}",
        ]
    )

    fields = _field_names(obj)
    if fields:
        lines.append(f"Fields: {', '.join(fields)}")

    if provenance.processing_steps:
        lines.append("")
        lines.append("Method")
        for step in provenance.processing_steps[:3]:
            lines.extend(f"- {line}" for line in _wrap(_plain_value(step), width=32))

    return "\n".join(lines).strip()


def _criteria_rows(context, obj: Any | None) -> list[tuple[str, str]]:
    metadata = context.metadata
    provenance = context.provenance
    fields = _field_names(obj)
    rows = [
        ("C-P1 Provider / Source", _join_values(metadata.provider, metadata.source)),
        ("C-P2 Data Collection Year", _value(metadata.collection_period)),
        ("C-P3 License", _value(metadata.license)),
        ("C-P4 CRS", _value(metadata.crs)),
        ("C-P5 LOD", _value(metadata.lod)),
        ("C-P6 Data Type", _data_type_value(metadata)),
        ("C-P7 Available Formats", _value(metadata.formats)),
        ("C-S1 Scale", MISSING),
        ("C-S2 Geographic Coverage", _value(metadata.geographic_coverage)),
        ("C-S3 Description", _value(metadata.description)),
        ("C-S4 Attributes", _value(fields)),
        ("C-S5 Update Frequency", _value(metadata.update_frequency)),
        (
            "C-S6 Link to Original Resource",
            _value(_resource_links(metadata.source, provenance.sources)),
        ),
        ("C-S7 Processing / Methodology", _value(provenance.processing_steps)),
        ("C-S8 Machine-Readable Fields", _machine_fields_value(fields)),
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


def _field_names(obj: Any | None) -> list[str]:
    if obj is None:
        return []
    names = getattr(obj, "field_names", None)
    if callable(names):
        names = names()
    if names is None:
        return []
    return [str(name) for name in names if str(name)]


def _machine_fields_value(fields: list[str]) -> str:
    if not fields:
        return "Manifest v2 is machine-readable; present fields not specified."
    return f"Manifest v2 is machine-readable; present fields: {', '.join(fields)}"


def _data_type_value(metadata) -> str:
    values = [
        *list(metadata.data_types or ()),
        metadata.data_category,
        metadata.result_kind,
        metadata.python_return_type,
    ]
    deduped = []
    for value in values:
        if value and str(value) not in deduped:
            deduped.append(str(value))
    return _value(deduped)


def _resource_links(*sources) -> list[str]:
    links = []
    for source_group in sources:
        for source in source_group or ():
            if isinstance(source, str) and source.startswith(("http://", "https://")):
                links.append(source)
            elif isinstance(source, dict):
                for key in ("url", "href", "link", "resource"):
                    value = source.get(key)
                    if isinstance(value, str) and value.startswith(("http://", "https://")):
                        links.append(value)
    return links


def _join_values(*values) -> str:
    parts = [_value(value) for value in values if _value(value) != MISSING]
    return "\n".join(parts) if parts else MISSING


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
