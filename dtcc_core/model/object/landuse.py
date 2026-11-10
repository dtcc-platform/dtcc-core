from collections import Counter
from dataclasses import dataclass, field
from enum import Enum, auto

import numpy as np

from ..geometry import Surface
from .object import GeometryType, Object


class LanduseClasses(Enum):
    """
    Enumeration of land use classification types for geographic or urban modeling applications.

    This enum defines various categories of land cover or land utilization, useful in GIS,
    remote sensing, simulation, or 3D urban environments. Each class represents a general
    type of surface usage, ranging from natural to heavily developed areas.

    Attributes
    ----------
    WATER
        Bodies of water such as lakes, rivers, and oceans.
    GRASS
        Grass-covered areas including parks, fields, and lawns.
    FOREST
        Forested or densely vegetated regions.
    FARMLAND
        Agricultural land used for crops or pasture.
    LIGHT_URBAN
        Low-density urban development, such as suburbs or residential zones.
    URBAN
        Medium-density urban areas with mixed land use.
    HEAVY_URBAN
        High-density urban cores, often with large buildings or infrastructure.
    INDUSTRIAL
        Zones designated for industrial activity, including factories and
        warehouses.
    MILITARY
        Areas reserved for military use or restricted access.
    ROAD
        Land primarily used for roads and transportation corridors.
    RAIL
        Railway infrastructure and corridors.
    UNKNOWN
        Land use is unknown or unclassified (explicitly set to 9999).
    """

    WATER = auto()
    GRASS = auto()
    FOREST = auto()
    FARMLAND = auto()
    LIGHT_URBAN = auto()
    URBAN = auto()
    HEAVY_URBAN = auto()
    INDUSTRIAL = auto()
    MILITARY = auto()
    ROAD = auto()
    RAIL = auto()
    UNKNOWN = 9999


# Fixed colours so a class keeps its colour across plots. Urban Atlas classes are
# grouped by family: reds for housing, purples for industry and transport hubs,
# greys for roads and rail, greens for vegetation, yellows for agriculture and
# blues for water.
URBAN_ATLAS_CLASS_COLORS = {
    "11100": "#7A0A0A",
    "11210": "#B3261E",
    "11220": "#D9534F",
    "11230": "#EE8C83",
    "11240": "#F7BDB5",
    "11300": "#C9A39B",
    "12100": "#9C4FC4",
    "12210": "#5E5E66",
    "12220": "#9A9AA2",
    "12230": "#3F3A36",
    "12300": "#7E57C2",
    "12400": "#B39DDB",
    "13100": "#8D6E63",
    "13300": "#BCAAA4",
    "13400": "#E6D3A3",
    "14100": "#7CC46B",
    "14200": "#B5E08C",
    "21000": "#F5E663",
    "22000": "#E8B84B",
    "23000": "#D8E07A",
    "24000": "#F0D58C",
    "25000": "#C9A227",
    "31000": "#2E7D32",
    "32000": "#A5C882",
    "33000": "#D7D2C4",
    "40000": "#8FB8D9",
    "50000": "#3B82C4",
}
# Short Urban Atlas class names for legends; ``ua_class`` keeps the full names.
URBAN_ATLAS_SHORT_NAMES = {
    "11100": "Continuous urban fabric (>80% sealed)",
    "11210": "Dense urban fabric (50-80%)",
    "11220": "Medium density urban fabric (30-50%)",
    "11230": "Low density urban fabric (10-30%)",
    "11240": "Very low density urban fabric (<10%)",
    "11300": "Isolated structures",
    "12100": "Industrial, commercial, public units",
    "12210": "Fast transit roads",
    "12220": "Other roads",
    "12230": "Railways",
    "12300": "Port areas",
    "12400": "Airports",
    "13100": "Mineral extraction and dump sites",
    "13300": "Construction sites",
    "13400": "Land without current use",
    "14100": "Green urban areas",
    "14200": "Sports and leisure facilities",
    "21000": "Arable land",
    "22000": "Permanent crops",
    "23000": "Pastures",
    "24000": "Complex and mixed cultivation",
    "25000": "Orchards",
    "31000": "Forests",
    "32000": "Herbaceous vegetation",
    "33000": "Open spaces with little vegetation",
    "40000": "Wetlands",
    "50000": "Water",
}
LANDUSE_CLASS_COLORS = {
    "WATER": "#3B82C4",
    "GRASS": "#7CC46B",
    "FOREST": "#2E7D32",
    "FARMLAND": "#F5E663",
    "LIGHT_URBAN": "#EE8C83",
    "URBAN": "#D9534F",
    "HEAVY_URBAN": "#7A0A0A",
    "INDUSTRIAL": "#9C4FC4",
    "MILITARY": "#6D4C41",
    "ROAD": "#9A9AA2",
    "RAIL": "#3F3A36",
    "UNKNOWN": "#54545A",
}
UNKNOWN_LANDUSE_COLOR = "#54545A"
# Legends list at most this many classes; the rest are grouped as "Other".
MAX_LEGEND_CLASSES = 14


@dataclass(repr=False)
class Landuse(Object):
    """
    Represents a land use object with associated land use classifications and geometric surfaces.

    This class models areas of land with specified types of usage (e.g., forest, urban, water),
    along with their associated geometric representation as multi-surfaces.

    Attributes
    ----------
    landuses : List[LanduseClasses]
        A list of land use classes describing how the land is used.
    """

    landuses: list[LanduseClasses] = field(default_factory=list)

    def _summary_items(self):
        return super()._summary_items() + [("num_landuses", len(self.landuses))]

    def _info_sections(self):
        from collections import Counter

        sections = super()._info_sections()
        if self.landuses:
            counts = Counter(value.name for value in self.landuses)
            sections.append(
                ("Land use classes", ("Class", "Count"), list(counts.items()))
            )
        return sections

    @property
    def surfaces(self) -> list[Surface]:
        """
        Access the list of surfaces representing land use polygons.

        Returns
        -------
        list[Surface]
            Surfaces stored under the MultiSurface geometry, or an empty list.
        """
        geom = self.get_geometry(GeometryType.MULTISURFACE)
        if geom is None:
            return []
        return geom.surfaces

    def plot(
        self,
        ax=None,
        *,
        legend: bool = True,
        title: str | None = None,
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
    ):
        """Plot the land use polygons as a 2D map coloured by class.

        Urban Atlas results (with a ``ua_code`` attribute) are coloured by
        their 27 Urban Atlas classes, others by :class:`LanduseClasses`. Each
        class keeps a fixed colour, and the legend lists the classes by share
        of the mapped area.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw into. A new figure is created when omitted.
        legend : bool, default True
            Whether to show the class legend with area shares.
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
        from dtcc_core.plotting.style import (
            apply_dtcc_style,
            get_theme,
            plot_polygon_geometries,
        )

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        own_figure = ax is None
        ax, panel_ax = presentation_plot_axes(context, ax=ax, presentation=presentation)
        colors = get_theme(theme)

        polygons = [surface.to_polygon() for surface in self.surfaces]
        labels, facecolors = self._plot_classes()
        if polygons:
            plot_polygon_geometries(
                ax,
                polygons,
                facecolors=facecolors,
                edgecolor=colors["axes"],
                linewidth=0.3,
                alpha=1.0,
            )
            if legend:
                self._draw_legend(
                    ax, polygons, labels, facecolors, colors, panel=presentation_enabled
                )
        else:
            ax.text(
                0.5,
                0.5,
                self._empty_message(),
                transform=ax.transAxes,
                ha="center",
                va="center",
                color=colors["text"],
                fontsize=11,
                wrap=True,
            )

        self._set_extent(ax, polygons, colors)
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=True,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x (m, EPSG:3006)",
            ylabel=None if presentation_enabled else "y (m, EPSG:3006)",
        )
        ax.ticklabel_format(useOffset=False, style="plain")
        _draw_scale_bar(ax, colors)
        ax.set_title(
            title or self._default_title(), color=colors["text"], fontweight=600, pad=10
        )
        attribution = self.attributes.get("attribution")
        if attribution:
            ax.text(
                0.995,
                0.005,
                attribution,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=5.5,
                color=colors["text"],
                zorder=6,
                bbox={
                    "boxstyle": "square,pad=0.25",
                    "facecolor": colors["axes"],
                    "edgecolor": "none",
                    "alpha": 0.75,
                },
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

    def _plot_classes(self):
        """Return a legend label and colour per surface."""
        codes = self.attributes.get("ua_code")
        if isinstance(codes, list) and len(codes) == len(self.surfaces):
            names = self.attributes.get("ua_class") or [None] * len(codes)
            labels = [
                URBAN_ATLAS_SHORT_NAMES.get(str(code)) or name or str(code)
                for code, name in zip(codes, names, strict=True)
            ]
            facecolors = [
                URBAN_ATLAS_CLASS_COLORS.get(str(code), UNKNOWN_LANDUSE_COLOR)
                for code in codes
            ]
            return labels, facecolors
        names = [value.name for value in self.landuses]
        names += ["UNKNOWN"] * (len(self.surfaces) - len(names))
        labels = [name.replace("_", " ").capitalize() for name in names]
        return labels, [
            LANDUSE_CLASS_COLORS.get(name, UNKNOWN_LANDUSE_COLOR) for name in names
        ]

    @staticmethod
    def _draw_legend(ax, polygons, labels, facecolors, colors, panel):
        from matplotlib.patches import Patch

        area = Counter()
        color_of = {}
        for polygon, label, color in zip(polygons, labels, facecolors, strict=True):
            area[label] += polygon.area
            color_of[label] = color
        total = sum(area.values()) or 1.0
        ranked = area.most_common()
        shown, rest = ranked[:MAX_LEGEND_CLASSES], ranked[MAX_LEGEND_CLASSES:]
        handles = [
            Patch(
                facecolor=color_of[label],
                edgecolor=colors["text"],
                linewidth=0.4,
                label=f"{label} ({_share(value, total)})",
            )
            for label, value in shown
        ]
        if rest:
            other = sum(value for _, value in rest)
            handles.append(
                Patch(
                    facecolor="none",
                    edgecolor=colors["text"],
                    linewidth=0.4,
                    label=f"{len(rest)} other classes ({_share(other, total)})",
                )
            )
        ncol = 3 if panel else 2
        if panel:
            # The panel axes have a fixed box; make room below the map for the legend rows.
            rows = -(-len(handles) // ncol)
            box = ax.get_position()
            shift = min(0.025 * rows, 0.3 * box.height)
            ax.set_position([box.x0, box.y0 + shift, box.width, box.height - shift])
        legend = ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.02 if panel else -0.1),
            frameon=False,
            ncol=ncol,
            fontsize=6.5 if panel else 7,
            handlelength=1.2,
            columnspacing=1.0,
        )
        for text in legend.get_texts():
            text.set_color(colors["text"])

    def _set_extent(self, ax, polygons, colors):
        requested = self.attributes.get("requested_bounds")
        if requested is not None:
            xmin, ymin, xmax, ymax = requested
            ax.plot(
                [xmin, xmax, xmax, xmin, xmin],
                [ymin, ymin, ymax, ymax, ymin],
                color=colors["muted"],
                linewidth=0.8,
                linestyle="--",
                zorder=3,
            )
        elif polygons:
            xmin = min(p.bounds[0] for p in polygons)
            ymin = min(p.bounds[1] for p in polygons)
            xmax = max(p.bounds[2] for p in polygons)
            ymax = max(p.bounds[3] for p in polygons)
        else:
            return
        pad = 0.03 * max(xmax - xmin, ymax - ymin, 1.0)
        ax.set_xlim(xmin - pad, xmax + pad)
        ax.set_ylim(ymin - pad, ymax + pad)

    def _default_title(self) -> str:
        source = self.attributes.get("source") or "Land use"
        places = sorted(
            {name for name in self.attributes.get("fua_name") or [] if name}
        )
        return f"{source}: {', '.join(places)}" if places else source

    def _empty_message(self) -> str:
        if self.attributes.get("partial_result"):
            return "No land use: the download failed (see upstream_errors)"
        if "ua_code" in self.attributes:
            return "No Urban Atlas land use here.\nUrban Atlas covers Functional Urban Areas only."
        return "No land use polygons"


def _share(value: float, total: float) -> str:
    percent = 100 * value / total
    return "<1%" if 0 < percent < 1 else f"{percent:.0f}%"


def _draw_scale_bar(ax, colors):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    target = (xmax - xmin) / 5
    if target <= 0:
        return
    magnitude = 10 ** np.floor(np.log10(target))
    length = max(step * magnitude for step in (1, 2, 5) if step * magnitude <= target)
    x0 = xmin + 0.04 * (xmax - xmin)
    y0 = ymin + 0.05 * (ymax - ymin)
    ax.plot([x0, x0 + length], [y0, y0], color=colors["text"], linewidth=2, zorder=5)
    label = f"{length / 1000:g} km" if length >= 1000 else f"{length:g} m"
    ax.text(
        x0 + length / 2,
        y0 + 0.012 * (ymax - ymin),
        label,
        ha="center",
        va="bottom",
        fontsize=7.5,
        color=colors["text"],
        zorder=5,
    )
