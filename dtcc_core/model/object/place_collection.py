from collections import Counter
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ...common import warning
from ...plotting.style import (
    add_categorical_legend,
    add_plot_context,
    apply_dtcc_style,
    get_theme,
    plot_line_segments,
)
from ..geometry import Point
from .object import Object
from .roadnetwork import RoadNetwork

# Fixed colour per place category, so a category keeps its colour across plots.
PLACE_CATEGORY_COLORS = {
    "food_and_drink": "#E35A1D",
    "shopping": "#FADA36",
    "education": "#4E9FE6",
    "health": "#E5484D",
    "culture": "#A27FE0",
    "accommodation": "#78C8BE",
    "tourism": "#4CB050",
    "leisure": "#B5DB6E",
    "public_service": "#B5835A",
    "finance": "#E377C2",
    "office": "#7C8CA3",
    "craft": "#A8A832",
    "other": "#C9C9D1",
}
UNKNOWN_CATEGORY_COLOR = "#54545A"
DEFAULT_PLACE_ATTRIBUTION = "© OpenStreetMap contributors"


@dataclass(repr=False)
class PlaceCollection(Object):
    """Represents a collection of places of interest, such as shops and cafés.

    Each place is a child ``Object`` with a point geometry for its location
    and metadata such as name, category and source tags stored in attributes.

    Attributes
    ----------
    context_roads : RoadNetwork or None
        Roads for the same area, drawn underneath the places by :meth:`plot`.
        They give the map geographic context and are not places themselves.
    """

    context_roads: RoadNetwork | None = field(default=None, kw_only=True)

    def _info_sections(self):
        sections = super()._info_sections()
        places = self.places()
        counts = Counter(
            place.attributes.get("category", "unknown") for place in places
        )
        if counts:
            sections.append(
                (
                    "Places by category",
                    ("Category", "Count"),
                    [(category, count) for category, count in counts.most_common()],
                )
            )
        rows = []
        for place in places[:3]:
            point = self._point_geometry(place)
            rows.append(
                (
                    place.attributes.get("name", place.id),
                    place.attributes.get("category", "unknown"),
                    "N/A" if point is None else f"({point.x:.2f}, {point.y:.2f})",
                )
            )
        if rows:
            sections.append(
                ("Sample places (first 3)", ("Name", "Category", "Location"), rows)
            )
        return sections

    def _summary_items(self):
        return [("num_places", sum(len(group) for group in self.children.values()))]

    def add_place(self, place: Object) -> None:
        """Add a place as a child object."""
        if not isinstance(place, Object):
            raise TypeError("Place must be an Object instance")
        self.add_child(place)

    def places(self) -> list[Object]:
        """Return all places in the collection."""
        result = []
        for child_list in self.children.values():
            result.extend(child_list)
        return result

    def to_geodataframe(self, crs: str = "EPSG:3006"):
        """Return the places as a GeoDataFrame with one point row per place.

        Parameters
        ----------
        crs : str, default "EPSG:3006"
            Coordinate reference system of the place coordinates.

        Returns
        -------
        geopandas.GeoDataFrame
            One row per place with its attributes as columns.
        """
        import geopandas as gpd
        from shapely.geometry import Point as ShapelyPoint

        rows, geometries = [], []
        for place in self.places():
            point = self._point_geometry(place)
            if point is None:
                continue
            rows.append(dict(place.attributes))
            geometries.append(ShapelyPoint(point.x, point.y))
        return gpd.GeoDataFrame(rows, geometry=geometries, crs=crs)

    def plot(
        self,
        ax=None,
        *,
        roads=None,
        labels: int = 12,
        size: float = 26.0,
        legend: bool = True,
        interactive: bool = True,
        equal_aspect: bool = True,
        title: str | None = None,
        metadata: bool | dict[str, Any] = True,
        metadata_loc: str = "upper left",
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
    ):
        """Plot the places as a 2D map coloured by category.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw into. A new figure is created when omitted.
        roads : RoadNetwork or False, optional
            Roads drawn underneath the places for geographic context. Defaults
            to :attr:`context_roads`; pass False to draw no roads.
        labels : int, default 12
            Maximum number of place names written on the map. Labels are
            spread out so they do not overlap; 0 disables them.
        size : float, default 26.0
            Marker size.
        legend : bool, default True
            Whether to show the category legend with place counts.
        interactive : bool, default True
            Show a place's name, category and OSM tag when hovering over it
            in an interactive Matplotlib window.
        equal_aspect : bool, default True
            Whether to use the same scale on both axes.
        title, metadata, metadata_loc, theme, presentation, show
            As for the other DTCC dataset plots.

        Returns
        -------
        matplotlib.axes.Axes
            Axes containing the map.
        """
        from dtcc_core.datasets.presentation import (
            draw_empty_preview_state,
            finalize_presentation_plot,
            presentation_plot_axes,
        )

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        own_figure = ax is None
        ax, panel_ax = presentation_plot_axes(context, ax=ax, presentation=presentation)
        theme_values = get_theme(theme)

        places = [(place, self._point_geometry(place)) for place in self.places()]
        places = [(place, point) for place, point in places if point is not None]
        requested = self.attributes.get("requested_bounds")

        if roads is None:
            roads = self.context_roads
        if roads is not None and roads is not False:
            segments = roads._plot_segments()
            if segments:
                plot_line_segments(
                    ax,
                    segments,
                    color=theme_values["grid"],
                    linewidth=0.7,
                    theme=theme,
                    zorder=1,
                )
        if requested is not None:
            xmin, ymin, xmax, ymax = requested
            ax.plot(
                [xmin, xmax, xmax, xmin, xmin],
                [ymin, ymin, ymax, ymax, ymin],
                color=theme_values["muted"],
                linewidth=0.8,
                linestyle="--",
                zorder=1,
            )

        if places:
            xy = np.array([(point.x, point.y) for _, point in places])
            categories = [
                place.attributes.get("category") or "unknown" for place, _ in places
            ]
            colors = [
                PLACE_CATEGORY_COLORS.get(c, UNKNOWN_CATEGORY_COLOR) for c in categories
            ]
            scatter = ax.scatter(
                xy[:, 0],
                xy[:, 1],
                s=size,
                c=colors,
                edgecolors=theme_values["axes"],
                linewidths=0.6,
                zorder=3,
            )
            if legend:
                counts = Counter(categories)
                add_categorical_legend(
                    ax,
                    {
                        f"{category.replace('_', ' ')} ({count})": PLACE_CATEGORY_COLORS.get(
                            category, UNKNOWN_CATEGORY_COLOR
                        )
                        for category, count in counts.most_common()
                    },
                    title=None,
                    marker="patch",
                    theme=theme,
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.02 if presentation_enabled else -0.1),
                    ncol=min(4, len(counts)),
                    fontsize=8,
                    handlelength=1.2,
                    columnspacing=1.2,
                )
            names = [place.attributes.get("name") or "" for place, _ in places]
            self._draw_labels(ax, xy, names, labels, requested, theme_values)
            if interactive:
                self._connect_hover(ax, scatter, places, theme_values)
        else:
            warning("PlaceCollection has no places to plot.")
            if presentation_enabled:
                draw_empty_preview_state(ax, context)
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No places in the requested bounds",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    color=theme_values["text"],
                    fontsize=12,
                )

        self._set_extent(ax, places, requested)
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=equal_aspect,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x (m, EPSG:3006)",
            ylabel=None if presentation_enabled else "y (m, EPSG:3006)",
            grid=not presentation_enabled,
        )
        ax.ticklabel_format(useOffset=False, style="plain")
        self._draw_scale_bar(ax, theme_values)
        ax.text(
            0.99,
            0.01,
            self.attributes.get("attribution", DEFAULT_PLACE_ATTRIBUTION),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7,
            color=theme_values["muted"],
            zorder=5,
        )
        add_plot_context(
            ax,
            title=title if presentation_enabled else title or "DTCC Places of Interest",
            metadata=None
            if presentation_enabled
            else self._plot_metadata(len(places), metadata),
            bounds=None,
            theme=theme,
            metadata_loc=metadata_loc,
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

    def _plot_metadata(self, place_count: int, metadata):
        if metadata is False:
            return None
        values = {"Places": place_count}
        if isinstance(metadata, dict):
            values.update(metadata)
        return values

    @staticmethod
    def _set_extent(ax, places, requested):
        if requested is not None:
            xmin, ymin, xmax, ymax = requested
        elif places:
            xs = [point.x for _, point in places]
            ys = [point.y for _, point in places]
            xmin, ymin, xmax, ymax = min(xs), min(ys), max(xs), max(ys)
        else:
            return
        padding = 0.04 * max(xmax - xmin, ymax - ymin, 1.0)
        ax.set_xlim(xmin - padding, xmax + padding)
        ax.set_ylim(ymin - padding, ymax + padding)

    @staticmethod
    def _draw_labels(ax, xy, names, max_labels, requested, theme_values):
        if max_labels <= 0:
            return
        import matplotlib.patheffects as patheffects

        if requested is not None:
            span = max(requested[2] - requested[0], requested[3] - requested[1])
        else:
            span = float(np.ptp(xy, axis=0).max()) if len(xy) > 1 else 1.0
        spacing = 0.12 * max(span, 1.0)
        chosen = []
        for index, name in enumerate(names):
            if not name:
                continue
            if all(np.hypot(*(xy[index] - xy[other])) >= spacing for other in chosen):
                chosen.append(index)
                if len(chosen) == max_labels:
                    break
        outline = [
            patheffects.withStroke(linewidth=2.5, foreground=theme_values["axes"])
        ]
        for index in chosen:
            ax.annotate(
                names[index],
                xy[index],
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7.5,
                color=theme_values["text"],
                path_effects=outline,
                zorder=4,
            )

    @staticmethod
    def _connect_hover(ax, scatter, places, theme_values):
        tooltip = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(10, 10),
            textcoords="offset points",
            fontsize=8,
            color=theme_values["text"],
            zorder=6,
            visible=False,
            bbox={
                "boxstyle": "round,pad=0.3",
                "facecolor": theme_values["panel"],
                "edgecolor": theme_values["panel_edge"],
            },
        )

        def on_move(event):
            if event.inaxes is not ax:
                return
            hit, details = scatter.contains(event)
            if hit:
                place, point = places[details["ind"][0]]
                attributes = place.attributes
                tooltip.xy = (point.x, point.y)
                tooltip.set_text(
                    f"{attributes.get('name', '')}\n{attributes.get('category', '')}"
                    f" · {attributes.get('osm_key', '')}={attributes.get('osm_value', '')}"
                )
                tooltip.set_visible(True)
            elif tooltip.get_visible():
                tooltip.set_visible(False)
            else:
                return
            ax.figure.canvas.draw_idle()

        ax.figure.canvas.mpl_connect("motion_notify_event", on_move)

    @staticmethod
    def _draw_scale_bar(ax, theme_values):
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        target = (xmax - xmin) / 5
        if target <= 0:
            return
        magnitude = 10 ** np.floor(np.log10(target))
        length = max(
            step * magnitude for step in (1, 2, 5) if step * magnitude <= target
        )
        x0 = xmin + 0.04 * (xmax - xmin)
        y0 = ymin + 0.05 * (ymax - ymin)
        ax.plot(
            [x0, x0 + length],
            [y0, y0],
            color=theme_values["text"],
            linewidth=2,
            zorder=5,
        )
        label = f"{length / 1000:g} km" if length >= 1000 else f"{length:g} m"
        ax.text(
            x0 + length / 2,
            y0 + 0.012 * (ymax - ymin),
            label,
            ha="center",
            va="bottom",
            fontsize=7.5,
            color=theme_values["text"],
            zorder=5,
        )

    @staticmethod
    def _point_geometry(place: Object):
        return next((g for g in place.get_geometries() if isinstance(g, Point)), None)
