from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..geometry.bounds import Bounds
from .raster import Raster

POPULATION_GRID_CREDIT = "Source: Statistics Sweden (SCB), 1 km population grid"
POPULATION_GRID_UNIT = "Persons per 1 km cell"
ZERO_COLOR = "#54545A"
# Grids with at most this many cells get their counts written in the cells.
MAX_LABELLED_CELLS = 64


@dataclass(repr=False)
class PopulationGrid(Raster):
    """Registered population per 1 km cell, as a :class:`Raster` with plot metadata.

    ``data`` holds persons per cell for one variable. Cells without registered
    residents are 0; cells outside Sweden are ``nodata`` (NaN).

    Attributes
    ----------
    year : int or None
        Reference year; the population registered on 31 December.
    variable : str
        Population variable, such as ``"total"``, ``"women"`` or ``"age_20_24"``.
    requested_bounds : Bounds or None
        Bounds that were requested. The whole 1 km cells covering them extend
        beyond these bounds.
    area_boundaries : shapely geometry or None
        Statistical area (DeSO) boundaries within the grid, drawn by
        :meth:`plot` for geographic context.
    """

    year: int | None = field(default=None, kw_only=True)
    variable: str = field(default="total", kw_only=True)
    requested_bounds: Bounds | None = field(default=None, kw_only=True)
    area_boundaries: Any = field(default=None, kw_only=True)

    def _info_sections(self):
        sections = super()._info_sections()
        rows = sections[0][2]
        rows.extend([("Variable", self.variable), ("Year", self.year)])
        data = self.data if self.data.ndim == 2 else np.empty((0, 0))
        rows.extend(
            [
                ("Cells", data.size),
                ("Cells with residents", int(np.sum(data > 0))),
                ("Cells outside Sweden (nodata)", int(np.sum(np.isnan(data)))),
                ("Persons in cells", int(np.nansum(data))),
            ]
        )
        return sections

    def to_raster(self, nodata: float | None = None) -> Raster:
        """Return the data as a plain :class:`Raster`, for example for Protobuf export.

        Parameters
        ----------
        nodata : float, optional
            Value that replaces the NaN cells outside Sweden and becomes the
            raster's ``nodata``, for formats that need a numeric nodata value
            such as GeoTIFF. NaN is kept by default.
        """
        data = self.data.copy()
        if nodata is None:
            nodata = self.nodata
        else:
            data[np.isnan(data)] = nodata
        return Raster(data=data, georef=self.georef, nodata=nodata, crs=self.crs)

    def plot(
        self,
        ax=None,
        *,
        cmap: str = "viridis",
        labels: bool = True,
        title: str | None = None,
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
    ):
        """Plot the 1 km cells as a 2D map.

        Inhabited cells are coloured on a sequential scale of persons per cell,
        cells with no residents are grey, and cells outside Sweden are hatched.
        The requested bounds are dashed, since whole cells extend beyond them.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw into. A new figure is created when omitted.
        cmap : str, default "viridis"
            Sequential Matplotlib colormap for inhabited cells.
        labels : bool, default True
            Write the count in each cell for small grids (at most 64 cells).
        title, theme, presentation, show
            As for the other DTCC dataset plots.

        Returns
        -------
        matplotlib.axes.Axes
            Axes containing the map.
        """
        from dtcc_core.datasets.presentation import (
            finalize_presentation_plot,
            presentation_plot_axes,
        )
        from dtcc_core.plotting.style import apply_dtcc_style, get_theme

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        own_figure = ax is None
        ax, panel_ax = presentation_plot_axes(context, ax=ax, presentation=presentation)
        colors = get_theme(theme)

        data = np.asarray(self.data, dtype=float)
        if data.ndim != 2 or data.size == 0:
            raise ValueError("PopulationGrid.plot requires a 2D grid of cells")
        image = self._draw_cells(ax, data, cmap, colors)
        self._draw_context(ax, data.shape, colors)
        if labels and data.size <= MAX_LABELLED_CELLS:
            self._draw_counts(ax, data, colors)
        self._draw_key(
            ax,
            image,
            data,
            theme,
            colors,
            below=-0.02 if presentation_enabled else -0.1,
        )

        bounds = self.bounds
        pad = 0.03 * max(bounds.xmax - bounds.xmin, bounds.ymax - bounds.ymin)
        ax.set_xlim(bounds.xmin - pad, bounds.xmax + pad)
        ax.set_ylim(bounds.ymin - pad, bounds.ymax + pad)
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=True,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x (m, EPSG:3006)",
            ylabel=None if presentation_enabled else "y (m, EPSG:3006)",
        )
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(
            title or self._default_title(), color=colors["text"], fontweight=600, pad=10
        )
        self._draw_notes(
            ax, data, colors, below=-0.08 if presentation_enabled else -0.17
        )
        if own_figure and not presentation_enabled:
            ax.figure.tight_layout()
        return finalize_presentation_plot(
            ax,
            context=context,
            obj=self,
            panel_ax=panel_ax,
            presentation=presentation,
            show=show,
            theme=theme,
            title=title,
        )

    def _draw_cells(self, ax, data, cmap, colors):
        """Draw inhabited cells on ``cmap``, zero cells grey and nodata cells hatched."""
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Rectangle

        bounds = self.bounds
        extent = (bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax)
        cell = abs(self.georef.a)
        outside = np.isnan(data)
        inhabited = data > 0
        vmax = float(np.max(data[inhabited])) if inhabited.any() else 1.0
        # Masked cells are transparent, leaving zero and nodata cells to their own layers.
        image = ax.imshow(
            np.ma.masked_where(~inhabited, data),
            cmap=cmap,
            vmin=0.0,
            vmax=vmax,
            extent=extent,
            origin="upper",
            interpolation="nearest",
            zorder=1,
        )
        zero = np.where(~inhabited & ~outside, 1.0, np.nan)
        ax.imshow(
            np.ma.masked_invalid(zero),
            cmap=ListedColormap([ZERO_COLOR]),
            vmin=0,
            vmax=1,
            extent=extent,
            origin="upper",
            interpolation="nearest",
            zorder=1,
        )
        rows, cols = np.nonzero(outside)
        for row, col in zip(rows, cols, strict=True):
            ax.add_patch(
                Rectangle(
                    (bounds.xmin + col * cell, bounds.ymax - (row + 1) * cell),
                    cell,
                    cell,
                    facecolor="none",
                    edgecolor=colors["muted"],
                    hatch="///",
                    linewidth=0.0,
                    zorder=1,
                )
            )
        return image

    def _draw_context(self, ax, shape, colors):
        """Draw the 1 km cell lines, DeSO boundaries and the requested bounds."""
        bounds = self.bounds
        height, width = shape
        for x in np.linspace(bounds.xmin, bounds.xmax, width + 1):
            ax.axvline(x, color=colors["axes"], linewidth=0.5, zorder=2)
        for y in np.linspace(bounds.ymin, bounds.ymax, height + 1):
            ax.axhline(y, color=colors["axes"], linewidth=0.5, zorder=2)
        boundaries = self.area_boundaries
        if boundaries is not None and not boundaries.is_empty:
            for line in getattr(boundaries, "geoms", [boundaries]):
                if hasattr(line, "xy"):
                    xs, ys = line.xy
                    ax.plot(
                        xs,
                        ys,
                        color=colors["text"],
                        linewidth=0.4,
                        alpha=0.35,
                        zorder=3,
                    )
        if self.requested_bounds is not None:
            r = self.requested_bounds
            ax.plot(
                [r.xmin, r.xmax, r.xmax, r.xmin, r.xmin],
                [r.ymin, r.ymin, r.ymax, r.ymax, r.ymin],
                color=colors["text"],
                linewidth=1.2,
                linestyle="--",
                zorder=4,
            )

    def _draw_counts(self, ax, data, colors):
        """Write each cell's count at its centre, except for cells outside Sweden."""
        import matplotlib.patheffects as patheffects

        bounds = self.bounds
        cell = abs(self.georef.a)
        outline = [patheffects.withStroke(linewidth=2.2, foreground=colors["axes"])]
        for (row, col), value in np.ndenumerate(data):
            if np.isnan(value):
                continue
            ax.text(
                bounds.xmin + (col + 0.5) * cell,
                bounds.ymax - (row + 0.5) * cell,
                f"{value:,.0f}",
                ha="center",
                va="center",
                fontsize=8,
                color=colors["text"],
                path_effects=outline,
                zorder=5,
            )

    @staticmethod
    def _draw_key(ax, image, data, theme, colors, below):
        """Add the persons-per-cell colour scale and the zero, nodata and bounds legend."""
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        from dtcc_core.plotting.style import style_colorbar

        if np.any(data > 0):
            colorbar = ax.figure.colorbar(
                image, cax=ax.inset_axes([1.03, 0.1, 0.035, 0.8])
            )
            colorbar.set_label(POPULATION_GRID_UNIT)
            style_colorbar(colorbar, theme=theme)
        handles = [
            Patch(
                facecolor=ZERO_COLOR,
                edgecolor=colors["text"],
                linewidth=0.6,
                label="No registered residents (0)",
            ),
            Patch(
                facecolor="none",
                edgecolor=colors["muted"],
                hatch="///",
                label="Outside Sweden (no data)",
            ),
            Line2D(
                [0], [0], color=colors["text"], linestyle="--", label="Requested bounds"
            ),
        ]
        legend = ax.legend(
            handles=handles,
            loc="upper center",
            frameon=False,
            ncol=3,
            fontsize=8,
            handlelength=1.6,
            bbox_to_anchor=(0.5, below),
        )
        for text in legend.get_texts():
            text.set_color(colors["text"])

    def _draw_notes(self, ax, data, colors, below):
        """Add the whole-cell note, the SCB credit and any empty-grid message."""
        ax.text(
            0.0,
            below,
            "Counts are for whole 1 km cells, which extend beyond the requested bounds.\n"
            + POPULATION_GRID_CREDIT,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color=colors["muted"],
        )
        message = self._empty_message(data > 0, np.isnan(data))
        if message:
            ax.text(
                0.5,
                0.5,
                message,
                transform=ax.transAxes,
                ha="center",
                va="center",
                color=colors["text"],
                fontsize=11,
                zorder=6,
                bbox={
                    "boxstyle": "round,pad=0.4",
                    "facecolor": colors["panel"],
                    "edgecolor": colors["panel_edge"],
                },
            )

    def _default_title(self) -> str:
        variable = self.variable or "total"
        if variable.startswith("age_"):
            parts = variable[4:].split("_")
            label = (
                "Aged 100+"
                if parts == ["100", "plus"]
                else f"Aged {parts[0]}–{parts[1]}"
            )
        else:
            label = {"total": "Total population", "men": "Men", "women": "Women"}.get(
                variable, variable
            )
        return label if self.year is None else f"{label}, 31 December {self.year}"

    @staticmethod
    def _empty_message(inhabited, outside) -> str | None:
        if outside.all():
            return "All cells are outside Sweden (no data)"
        if not inhabited.any():
            return "No registered residents in these cells"
        return None
