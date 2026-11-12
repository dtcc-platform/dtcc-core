"""Bounded XY maps of scalar values on native triangular meshes."""

from numbers import Real

import numpy as np

from ..model import Mesh
from .style import (
    apply_dtcc_style,
    get_theme,
    require_matplotlib,
    resolve_colormap,
    style_colorbar,
)


def _mask(value, count, name, default):
    if value is None:
        return np.full(count, default, dtype=bool)
    result = np.asarray(value)
    if result.shape != (count,) or result.dtype != np.bool_:
        raise ValueError(f"{name} must be a boolean array with one entry per face")
    return result


def _limits(values, vmin, vmax):
    for name, value in (("vmin", vmin), ("vmax", vmax)):
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not np.isfinite(value)
        ):
            raise ValueError(f"{name} must be a finite real number")
    lo = float(values.min()) if vmin is None else float(vmin)
    hi = float(values.max()) if vmax is None else float(vmax)
    if lo == hi and (vmin is None or vmax is None):
        delta = max(abs(lo) * 0.01, 1e-12)
        lo = lo - delta if vmin is None else lo
        hi = hi + delta if vmax is None else hi
    if lo >= hi:
        raise ValueError("vmin must be less than vmax")
    return lo, hi


def _mesh_values(mesh, field):
    if not isinstance(mesh, Mesh):
        raise ValueError("projection='xy' requires a standalone triangular Mesh")
    matches = [f for f in mesh.fields if f.name == field]
    if len(matches) != 1:
        raise ValueError("Select exactly one existing face field with field=...")
    f = matches[0]
    values = np.asarray(f.values)
    if f.association != "face" or len(values) != len(mesh.faces):
        raise ValueError("Field count must match its face association")
    if f.dim != 1 or values.shape not in ((len(mesh.faces),), (len(mesh.faces), 1)):
        raise ValueError("XY face maps require a scalar field (dim=1)")
    values = values.reshape(-1)
    if not np.isrealobj(values) or not np.all(np.isfinite(values)):
        raise ValueError("Face values must be finite real numbers")
    return f, values


def _mesh_xy(mesh):
    faces, vertices = np.asarray(mesh.faces), np.asarray(mesh.vertices)
    if (
        faces.ndim != 2
        or faces.shape[1] != 3
        or not len(faces)
        or not np.issubdtype(faces.dtype, np.integer)
        or faces.min() < 0
        or faces.max() >= len(vertices)
    ):
        raise ValueError("XY face maps require nonempty valid triangular connectivity")
    if vertices.ndim != 2 or vertices.shape[1] != 3 or not np.isrealobj(vertices):
        raise ValueError("Mesh vertices must be real XYZ coordinates")
    matrix = mesh.transform._validate_affine(mesh.transform.affine)
    xy = (vertices @ matrix[:3, :3].T + matrix[:3, 3])[:, :2]
    if not np.all(np.isfinite(xy)):
        raise ValueError("Transformed XY coordinates must be finite")
    triangles = xy[faces]
    a, b = triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    if np.any(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0] == 0):
        raise ValueError("Mesh contains degenerate triangles in XY projection")
    return xy, faces


def _perimeter(xy, faces):
    """Outline the XY union, including holes and nonconforming triangle joins."""
    from shapely import polygons, union_all

    geometry = union_all(polygons(xy[faces]))
    pieces = [geometry] if geometry.geom_type == "Polygon" else geometry.geoms
    segments = []
    for piece in pieces:
        for ring in [piece.exterior, *piece.interiors]:
            coordinates = np.asarray(ring.coords)
            segments.extend(np.stack([coordinates[:-1], coordinates[1:]], axis=1))
    return np.asarray(segments)


def _sample(ids, remaining):
    if len(ids) > remaining:
        return ids[np.linspace(0, len(ids) - 1, remaining, dtype=int)]
    return ids


def _scale(ax, span, unit, colors):
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    target = span / 5
    magnitude = 10.0 ** np.floor(np.log10(target))
    size = (
        5 if target / magnitude >= 5 else 2 if target / magnitude >= 2 else 1
    ) * magnitude
    bar = AnchoredSizeBar(
        ax.transData,
        size,
        f"{size:g} {unit}",
        "lower left",
        pad=0.3,
        borderpad=0.5,
        frameon=True,
        color=colors["text"],
        size_vertical=span * 0.003,
        fontproperties={"size": 8},
    )
    bar.patch.set(facecolor=colors["panel"], edgecolor=colors["grid"], alpha=0.9)
    ax.add_artist(bar)


def _decorate(
    ax, mesh, xy, *, theme, title, coordinate_unit, xlabel, ylabel, scale_bar
):
    from matplotlib.ticker import MaxNLocator, ScalarFormatter

    colors = get_theme(theme)
    suffix = f" [{coordinate_unit}]" if coordinate_unit else ""
    apply_dtcc_style(
        ax,
        theme=theme,
        equal_aspect=True,
        title=title,
        xlabel=xlabel if xlabel is not None else "X" + suffix,
        ylabel=ylabel if ylabel is not None else "Y" + suffix,
    )
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = float((hi - lo).max())
    padding = 0.02 * span
    ax.set(
        xlim=(lo[0] - padding, hi[0] + padding), ylim=(lo[1] - padding, hi[1] + padding)
    )
    for axis in (ax.xaxis, ax.yaxis):
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)
        axis.set_major_locator(MaxNLocator(nbins=5))
    if scale_bar:
        _scale(ax, span, coordinate_unit, colors)
    if mesh.transform.srs:
        ax.text(
            0.0,
            1.01,
            mesh.transform.srs,
            transform=ax.transAxes,
            fontsize=8,
            color=colors["muted"],
            ha="left",
            va="bottom",
        )


def _layers(
    ax, xy, faces, values, selected, context, *, max_elements, cmap, limits, theme
):
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.colors import Normalize

    ids = np.flatnonzero(selected | context)
    ground = _sample(np.flatnonzero(selected), max_elements)
    background = _sample(np.flatnonzero(context), max_elements - len(ground))
    remaining = max_elements - len(ground) - len(background)
    artist = PolyCollection(
        xy[faces[ground]],
        array=values[ground],
        cmap=resolve_colormap(cmap),
        norm=Normalize(*limits),
        edgecolors="none",
        antialiaseds=False,
    )
    artist.set_gid("face_field")
    ax.add_collection(artist)
    if context.any():
        colors = get_theme(theme)
        neutral = PolyCollection(
            xy[faces[background]],
            facecolors=colors["muted"],
            edgecolors="none",
            antialiaseds=False,
        )
        neutral.set_gid("face_context")
        ax.add_collection(neutral)
        edges = _perimeter(xy, faces[context]) if remaining else np.empty((0, 2, 2))
        edge_ids = _sample(np.arange(len(edges)), remaining)
        outline = LineCollection(edges[edge_ids], colors=colors["text"], linewidths=0.7)
        outline.set_gid("context_outline")
        ax.add_collection(outline)
    else:
        edges = []
    if len(ids) + len(edges) > max_elements or (context.any() and not remaining):
        ax.text(
            0.01,
            0.01,
            "Large geometry sampled; context outline may be incomplete",
            transform=ax.transAxes,
            fontsize=8,
            color=get_theme(theme)["muted"],
        )
    return artist


def _plot_face_map(
    mesh,
    *,
    ax,
    field,
    face_mask,
    context_mask,
    vmin,
    vmax,
    cmap,
    title,
    coordinate_unit,
    xlabel,
    ylabel,
    scale_bar,
    max_elements,
    theme,
    show,
):
    f, values = _mesh_values(mesh, field)
    xy, faces = _mesh_xy(mesh)
    context = _mask(context_mask, len(faces), "context_mask", False)
    selected = (
        _mask(face_mask, len(faces), "face_mask", True)
        if face_mask is not None
        else ~context
    )
    if np.any(selected & context):
        raise ValueError("face_mask and context_mask must be disjoint")
    if not selected.any():
        raise ValueError("Select at least one field face")
    limits = _limits(values[selected], vmin, vmax)
    if coordinate_unit is not None and (
        not isinstance(coordinate_unit, str) or not coordinate_unit.strip()
    ):
        raise ValueError("coordinate_unit must be a nonempty string")
    if scale_bar and not coordinate_unit:
        raise ValueError("scale_bar requires an explicit coordinate_unit")
    if ax is not None and getattr(ax, "name", None) != "rectilinear":
        raise ValueError("XY face maps require a Matplotlib 2D axes")
    plt = require_matplotlib("XY face maps")
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7), layout="constrained")
    artist = _layers(
        ax,
        xy,
        faces,
        values,
        selected,
        context,
        max_elements=max_elements,
        cmap=cmap,
        limits=limits,
        theme=theme,
    )
    label = f.name + (f" [{f.unit}]" if f.unit else "")
    style_colorbar(
        ax.figure.colorbar(artist, ax=ax, pad=0.03, fraction=0.045, label=label),
        theme=theme,
    )
    _decorate(
        ax,
        mesh,
        xy[faces[selected | context]].reshape(-1, 2),
        theme=theme,
        title=title if title is not None else f.name,
        coordinate_unit=coordinate_unit,
        xlabel=xlabel,
        ylabel=ylabel,
        scale_bar=scale_bar,
    )
    ax.set_gid("face_field_map")
    if np.all(values[selected] == 0):
        ax.text(
            0.5,
            0.96,
            "All selected values are zero",
            transform=ax.transAxes,
            ha="center",
            va="top",
            color=get_theme(theme)["text"],
        )
    if show:
        plt.show()
    return ax
