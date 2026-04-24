"""Matplotlib helpers for inspecting footprint conditioning results."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient


def _load_plot_modules():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required to show footprint cleaning plots."
        ) from exc

    return plt, PathPatch, Path


def _iter_polygon_parts(geometry: BaseGeometry | None) -> list[Polygon]:
    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return [polygon for polygon in geometry.geoms if not polygon.is_empty]
    if isinstance(geometry, GeometryCollection):
        polygons: list[Polygon] = []
        for part in geometry.geoms:
            polygons.extend(_iter_polygon_parts(part))
        return polygons
    return []


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


def _polygon_patch(polygon: Polygon, *, path_cls: type[Any], patch_cls: type[Any], **kwargs):
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


def _plot_polygons(ax, polygons: Sequence[BaseGeometry], *, facecolor: str, edgecolor: str, linewidth: float, alpha: float):
    plotted_parts: list[Polygon] = []
    for geometry in polygons:
        for polygon in _iter_polygon_parts(geometry):
            patch = _polygon_patch(
                polygon,
                path_cls=ax._dtcc_path_cls,
                patch_cls=ax._dtcc_patch_cls,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=linewidth,
                alpha=alpha,
            )
            ax.add_patch(patch)
            plotted_parts.append(polygon)
    return plotted_parts


def _set_axes_extent(axes, polygons: Sequence[Polygon]) -> None:
    if not polygons:
        return
    minx = min(polygon.bounds[0] for polygon in polygons)
    miny = min(polygon.bounds[1] for polygon in polygons)
    maxx = max(polygon.bounds[2] for polygon in polygons)
    maxy = max(polygon.bounds[3] for polygon in polygons)
    span_x = max(maxx - minx, 1.0)
    span_y = max(maxy - miny, 1.0)
    padding = 0.05 * max(span_x, span_y)
    xmin = minx - padding
    xmax = maxx + padding
    ymin = miny - padding
    ymax = maxy + padding
    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])


def plot_footprint_cleaning_comparison(
    raw_polygons: Sequence[BaseGeometry],
    cleaned_polygons: Sequence[BaseGeometry],
    *,
    title: str = "Footprint cleaning",
    raw_title: str = "Input footprints",
    cleaned_title: str = "Conditioned footprints",
    show: bool = True,
    block: bool = True,
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

    plt, path_patch_cls, path_cls = _load_plot_modules()
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    for ax in axes:
        ax._dtcc_patch_cls = path_patch_cls  # type: ignore[attr-defined]
        ax._dtcc_path_cls = path_cls  # type: ignore[attr-defined]

    raw_parts = _plot_polygons(
        axes[0],
        raw_polygons,
        facecolor="#cbd5e1",
        edgecolor="#0f172a",
        linewidth=0.7,
        alpha=0.9,
    )
    cleaned_parts = _plot_polygons(
        axes[1],
        cleaned_polygons,
        facecolor="#8ecae6",
        edgecolor="#023047",
        linewidth=0.7,
        alpha=0.95,
    )

    axes[0].set_title(raw_title)
    axes[1].set_title(cleaned_title)

    all_parts = [*raw_parts, *cleaned_parts]
    _set_axes_extent(axes, all_parts)

    axes[0].text(
        0.02,
        0.02,
        f"Polygons: {len(raw_parts)}",
        transform=axes[0].transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "alpha": 0.9,
            "edgecolor": "#334155",
            "linewidth": 0.6,
        },
    )
    axes[1].text(
        0.02,
        0.02,
        f"Polygons: {len(cleaned_parts)}",
        transform=axes[1].transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "alpha": 0.9,
            "edgecolor": "#334155",
            "linewidth": 0.6,
        },
    )

    fig.suptitle(title)
    if show:
        plt.show(block=block)
    return fig, axes
