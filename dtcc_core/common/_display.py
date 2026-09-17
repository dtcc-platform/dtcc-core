"""Shared plain-text inspection formatting, independent of log configuration."""

from io import StringIO
from itertools import islice
import math
from numbers import Number
import reprlib

from rich.console import Console
from rich.markup import escape
from rich.text import Text

from .dtcc_logging import make_table


def compact(value):
    """Describe a value without expanding arrays, mappings or model trees."""
    if hasattr(value, "shape") and hasattr(value, "dtype"):
        if isinstance(value, Number) and hasattr(value, "item"):
            return repr(value.item())
        return f"array(shape={value.shape}, dtype={value.dtype})"
    if isinstance(value, (dict, list, set)):
        return f"{type(value).__name__}(len={len(value)})"
    formatter = reprlib.Repr()
    formatter.maxstring = 100
    formatter.maxother = 100
    formatter.maxlong = 100
    return formatter.repr(value)


def is_literal_number(value):
    """Whether compact() emits a complete finite built-in numeric literal."""
    if isinstance(value, Number) and hasattr(value, "item"):
        value = value.item()
    return (type(value) in (int, float, bool) and len(repr(value)) <= 100
            and (not isinstance(value, float) or math.isfinite(value)))


def format_repr(name, items, *, complete=False):
    body = ", ".join(f"{key}={compact(value)}" for key, value in items)
    result = f"{name}({body})"
    return result if complete else f"<{result}>"


def label(name):
    return {"id": "ID", "srs": "CRS", "crs": "CRS", "dtype": "Data type"}.get(
        name, name.removeprefix("num_").replace("_", " ").capitalize()
    )


def value_text(value):
    if value is None:
        return "None"
    if isinstance(value, str):
        return value
    return compact(value)


def format_info(title, sections, *, max_rows=None):
    """Render (heading, columns, rows) sections; plain prose has columns=None.

    Table cells are literal text, including user attributes containing Rich
    markup. A row limit may be used for model collections; parameter help and
    dataset metadata remain complete. Sections supply sized row sequences.
    """
    console = Console(file=StringIO(), width=120, color_system=None,
                      force_terminal=False, highlight=False)
    if title:
        console.print(Text(title))
    for heading, columns, rows in sections:
        console.print()
        if heading:
            console.print(Text(heading))
        if columns is None:
            console.print(Text(rows))
            continue
        visible = list(islice(rows, max_rows))
        remaining = max(0, len(rows) - len(visible))
        # make_table normally interprets string cells as markup. Escape at this
        # boundary rather than changing logging's existing markup contract.
        console.print(make_table(
            columns,
            ([escape(str(cell)) for cell in row] for row in visible),
            overflow="fold",
        ))
        if remaining:
            console.print(Text(f"... {remaining} more row(s)"))
    return console.file.getvalue().rstrip()
