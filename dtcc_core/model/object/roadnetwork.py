from collections import Counter
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any

import numpy as np

from ...common import warning
from ...plotting.style import (
    DTCC_COLORS,
    add_categorical_legend,
    add_plot_context,
    apply_dtcc_style,
    plot_line_segments,
)
from ..geometry import Bounds, LineString, MultiLineString
from .object import GeometryType, Object


class RoadType(Enum):
    """Enumeration representing different road types."""

    MOTORWAY = auto()
    PRIMARY = auto()
    SECONDARY = auto()
    TERTIARY = auto()
    RESIDENTIAL = auto()
    SERVICE = auto()
    TRACK = auto()
    PEDESTRIAN = auto()
    CYCLEWAY = auto()
    FOOTWAY = auto()
    BRIDLEWAY = auto()
    PATH = auto()


@dataclass(repr=False)
class RoadNetwork(Object):
    """
    Represents a road network as a graph of vertices and edges with associated lengths.

    Attributes
    ----------
    vertices : np.ndarray
        An array of vertex coordinates in the road network (shape: [n_vertices, dim]).

    edges : np.ndarray
        An array of edge indices representing connections between vertices
        (shape: [n_edges, 2], each row is [start_idx, end_idx]).

    length : np.ndarray
        An array of lengths corresponding to each edge in the network
        (shape: [n_edges]).
    """

    vertices: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    edges: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )  # each edge (start_idx, end_idx)
    length: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))

    def _summary_items(self):
        return [
            ("num_vertices", len(self.vertices)),
            ("num_edges", len(self.edges)),
            ("num_segments", len(self.length)),
        ]

    @property
    def linestrings(self) -> list[LineString]:
        """
        Access individual LineString geometries if present.

        Returns
        -------
        list[LineString]
            Line strings stored under the MultiLineString geometry, or an empty list.
        """
        geom = self.get_geometry(GeometryType.MULTILINESTRING)
        if geom is None:
            return []
        return geom.linestrings

    @property
    def multilinestrings(self) -> MultiLineString:
        """
        Access the MultiLineString geometry representation.

        Returns
        -------
        MultiLineString or None
            Stored multilinestring geometry if available.
        """
        geom = self.get_geometry(GeometryType.MULTILINESTRING)
        return geom

    @property
    def bounds(self):
        """
        Compute bounds from geometry or vertex coordinates.

        Returns
        -------
        Bounds
            Bounding box derived from the MultiLineString geometry when present,
            otherwise computed from vertex extrema.
        """
        if self._bounds is not None:
            return self._bounds
        geom = self.get_geometry(GeometryType.MULTILINESTRING)
        if geom is None:
            if len(self.vertices) == 0:
                return Bounds()
            xmin = np.min(self.vertices[:, 0])
            ymin = np.min(self.vertices[:, 1])
            xmax = np.max(self.vertices[:, 0])
            ymax = np.max(self.vertices[:, 1])
            return Bounds(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax)
        return geom.bounds

    def to_arrays(self, include_attributes=True, include_geometry=False) -> dict:
        """
        Convert the road network to low-level NumPy-friendly arrays.

        Parameters
        ----------
        include_attributes : bool, default True
            Include edge-aligned attributes converted to NumPy arrays.
        include_geometry : bool, default False
            Include flattened line geometry as ``line_vertices`` and
            ``line_offsets``. The offsets array has length n_lines + 1.

        Returns
        -------
        dict
            Dictionary containing at least ``vertices``, ``edges`` and
            ``lengths``.
        """
        arrays = {
            "vertices": np.asarray(self.vertices),
            "edges": np.asarray(self.edges),
            "lengths": np.asarray(self.length),
        }

        if include_attributes:
            arrays["attributes"] = {
                key: np.asarray(value) for key, value in self.attributes.items()
            }

        if include_geometry:
            line_vertices = []
            line_offsets = [0]
            dim = (
                self.vertices.shape[1] if getattr(self.vertices, "ndim", 0) == 2 else 2
            )
            for line in self.linestrings:
                vertices = np.asarray(line.vertices)
                if vertices.size == 0:
                    vertices = np.empty((0, dim))
                line_vertices.append(vertices)
                line_offsets.append(line_offsets[-1] + len(vertices))
            if line_vertices:
                arrays["line_vertices"] = np.vstack(line_vertices)
            else:
                arrays["line_vertices"] = np.empty((0, dim))
            arrays["line_offsets"] = np.asarray(line_offsets, dtype=np.int64)

        return arrays

    def _info_sections(self):
        sections = super()._info_sections()
        sections[0][2].append(("Line geometries", len(self.linestrings)))
        geometry = self.get_geometry(GeometryType.MULTILINESTRING)
        if not self.transform.srs and geometry is not None:
            sections[0][2].append(
                ("Geometry CRS", geometry.transform.srs or "Not specified")
            )
        if len(self.length):
            lengths = np.asarray(self.length, dtype=float)
            sections.append(
                (
                    "Length statistics",
                    ("Statistic", "Value"),
                    [
                        ("Count", len(lengths)),
                        ("Total", f"{lengths.sum():.2f}"),
                        ("Min", f"{lengths.min():.2f}"),
                        ("Max", f"{lengths.max():.2f}"),
                        ("Mean", f"{lengths.mean():.2f}"),
                    ],
                )
            )
        if "highway" in self.attributes:
            counts = Counter(
                value for value in self.attributes["highway"] if value not in (None, "")
            )
            if counts:
                sections.append(
                    ("Highway classes", ("Class", "Count"), counts.most_common())
                )
        if "oneway" in self.attributes:
            count = sum(
                value is True or str(value).strip().lower() in {"true", "yes", "1"}
                for value in self.attributes["oneway"]
            )
            sections[0][2].append(("One-way segments", count))
        return sections

    def plot(
        self,
        ax=None,
        column=None,
        color=DTCC_COLORS["teal"],
        linewidth=1.0,
        cmap=None,
        legend=True,
        equal_aspect=True,
        title: str | None = None,
        metadata: bool | dict[str, Any] = True,
        metadata_loc: str = "upper left",
        theme: str = "dark",
        presentation: bool = True,
        show=True,
        **kwargs,
    ):
        """
        Plot the road network using matplotlib.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw into. A new figure and axes are created when
            omitted.
        column : str, optional
            Edge-aligned attribute to color by. ``"length"`` and ``"lengths"``
            refer to road segment lengths.
        color : str, default DTCC teal
            Line color when no column is given.
        linewidth : float, default 1.0
            Width of road lines.
        cmap : str or matplotlib colormap, optional
            Matplotlib colormap for numeric or categorical columns. Defaults to
            the DTCC brand colormap.
        legend : bool, default True
            Whether to show a colorbar or categorical legend when coloring by a
            column.
        equal_aspect : bool, default True
            Whether to set equal axis scaling.
        show : bool, default True
            Whether to call ``matplotlib.pyplot.show()`` before returning.
        **kwargs
            Additional keyword arguments passed to ``LineCollection``.

        Returns
        -------
        matplotlib.axes.Axes
            Axes containing the road plot.
        """
        from dtcc_core.datasets.presentation import (
            draw_empty_preview_state,
            finalize_presentation_plot,
            presentation_plot_axes,
        )

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        segments = self._plot_segments()
        ax, panel_ax = presentation_plot_axes(
            context,
            ax=ax,
            presentation=presentation,
        )

        if len(segments) == 0:
            warning("RoadNetwork has no road segments to plot.")
            if presentation_enabled:
                draw_empty_preview_state(ax, context)
            apply_dtcc_style(
                ax,
                theme=theme,
                equal_aspect=equal_aspect,
                axis="off" if presentation_enabled else "on",
                xlabel=None if presentation_enabled else "x",
                ylabel=None if presentation_enabled else "y",
                grid=False if presentation_enabled else True,
            )
            add_plot_context(
                ax,
                title=None if presentation_enabled else title or "DTCC Road Network",
                metadata=(
                    None
                    if presentation_enabled
                    else self._plot_metadata(column, 0, metadata)
                ),
                bounds=None if presentation_enabled else self.bounds,
                theme=theme,
                metadata_loc=metadata_loc,
            )
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

        if column is None:
            plot_line_segments(
                ax,
                segments,
                color=color,
                linewidth=linewidth,
                theme=theme,
                **kwargs,
            )
            if legend and not presentation_enabled:
                add_categorical_legend(
                    ax,
                    {"Road segments": color},
                    title="Layers",
                    linewidth=linewidth,
                    theme=theme,
                )
        else:
            values = self._plot_column_values(column, len(segments))
            plot_line_segments(
                ax,
                segments,
                values=values,
                column_label=column,
                linewidth=linewidth,
                cmap=cmap,
                legend=legend and not presentation_enabled,
                theme=theme,
                **kwargs,
            )

        ax.autoscale()
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=equal_aspect,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x",
            ylabel=None if presentation_enabled else "y",
            grid=False if presentation_enabled else True,
        )
        add_plot_context(
            ax,
            title=None if presentation_enabled else title or "DTCC Road Network",
            metadata=(
                None
                if presentation_enabled
                else self._plot_metadata(column, len(segments), metadata)
            ),
            bounds=None if presentation_enabled else self.bounds,
            theme=theme,
            metadata_loc=metadata_loc,
        )
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

    def _plot_metadata(self, column, segment_count: int, metadata):
        if metadata is False:
            return None

        plot_metadata: dict[str, Any] = {
            "Segments": segment_count,
            "Vertices": len(self.vertices),
            "CRS": self._plot_crs(),
            "Column": column,
        }
        if len(self.length) > 0:
            plot_metadata["Total length"] = f"{float(np.sum(self.length)):.2f}"
        highway_values = self.attributes.get("highway")
        if highway_values is not None:
            highway_classes = {
                value for value in highway_values if value not in (None, "")
            }
            if highway_classes:
                plot_metadata["Highway classes"] = len(highway_classes)

        if isinstance(metadata, dict):
            plot_metadata.update(metadata)
        return plot_metadata

    def _plot_crs(self):
        crs = self.transform.srs
        geom = self.get_geometry(GeometryType.MULTILINESTRING)
        if not crs and geom is not None:
            crs = geom.transform.srs
        return crs

    def _plot_segments(self):
        segments = []
        if len(self.linestrings) > 0:
            for line in self.linestrings:
                vertices = np.asarray(line.vertices)
                if vertices.ndim == 2 and len(vertices) > 1:
                    segments.append(vertices[:, :2])
            return segments

        vertices = np.asarray(self.vertices)
        edges = np.asarray(self.edges, dtype=np.int64)
        if vertices.ndim != 2 or edges.size == 0:
            return segments

        edges = edges.reshape((-1, 2))
        for start, end in edges:
            if start < len(vertices) and end < len(vertices):
                segments.append(vertices[[start, end], :2])
        return segments

    def _plot_column_values(self, column, segment_count):
        if column in ("length", "lengths"):
            values = self.length
        else:
            if column not in self.attributes:
                available = ["length"] + sorted(self.attributes.keys())
                raise KeyError(
                    f"RoadNetwork has no attribute column '{column}'. "
                    f"Available columns: {available}"
                )
            values = self.attributes[column]

        if len(values) != segment_count:
            raise ValueError(
                f"Column '{column}' has {len(values)} value(s), "
                f"but the plot has {segment_count} segment(s)."
            )
        return values

    def to_shapely(self):
        """
        Convert the road network geometry to a Shapely MultiLineString.

        Returns
        -------
        shapely.geometry.MultiLineString or None
            Shapely representation when geometry exists, else ``None``.
        """
        multilinestring = self.get_geometry(GeometryType.MULTILINESTRING)
        if multilinestring is None:
            return None
        return multilinestring.to_shapely()
