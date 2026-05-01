from dataclasses import dataclass, field
import builtins
from collections import Counter
from typing import Union, List, Tuple
from enum import Enum, auto
from .object import Object, GeometryType
from ..geometry import LineString, MultiLineString
from ..geometry import Bounds
from .. import dtcc_pb2 as proto

import numpy as np


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


@dataclass
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

    def __str__(self):
        """Return a compact human-readable summary of the road network."""
        return (
            f"DTCC RoadNetwork with {len(self.vertices)} vertices, "
            f"{len(self.edges)} edge(s) and {len(self.length)} road segment(s)"
        )

    @property
    def linestrings(self) -> List[LineString]:
        """
        Access individual LineString geometries if present.

        Returns
        -------
        list[LineString]
            Line strings stored under the MultiLineString geometry, or an empty list.
        """
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
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
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
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
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
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
                self.vertices.shape[1]
                if getattr(self.vertices, "ndim", 0) == 2
                else 2
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

    def info(self, print: bool = True) -> str | None:
        """Print or return a human-readable multi-line summary."""
        lines = []
        lines.append("=" * 70)
        lines.append("DTCC RoadNetwork")
        lines.append("=" * 70)
        lines.append(f"Vertices: {len(self.vertices)}")
        lines.append(f"Edges: {len(self.edges)}")
        lines.append(f"Road segments: {len(self.length)}")
        lines.append(f"Line geometries: {len(self.linestrings)}")

        bounds = self.bounds
        if bounds is not None:
            lines.append(f"Bounds: {bounds}")

        crs = self.transform.srs
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
        if not crs and geom is not None:
            crs = geom.transform.srs
        if crs:
            lines.append(f"CRS: {crs}")

        if len(self.length) > 0:
            lengths = np.asarray(self.length, dtype=float)
            lines.append("")
            lines.append("Length Statistics:")
            lines.append(f"  Count: {len(lengths)}")
            lines.append(f"  Total: {np.sum(lengths):.2f}")
            lines.append(f"  Min: {np.min(lengths):.2f}")
            lines.append(f"  Max: {np.max(lengths):.2f}")
            lines.append(f"  Mean: {np.mean(lengths):.2f}")

        if self.attributes:
            lines.append("")
            lines.append("Attributes:")
            for key in sorted(self.attributes.keys()):
                value = self.attributes[key]
                try:
                    count = len(value)
                except TypeError:
                    count = 1
                lines.append(f"  {key}: {count} value(s)")

        if "highway" in self.attributes:
            highway_counts = Counter(
                value
                for value in self.attributes["highway"]
                if value not in (None, "")
            )
            if highway_counts:
                lines.append("")
                lines.append("Highway Classes:")
                for highway, count in highway_counts.most_common():
                    lines.append(f"  {highway}: {count}")

        if "oneway" in self.attributes:
            oneway_count = sum(
                value is True or str(value).strip().lower() in {"true", "yes", "1"}
                for value in self.attributes["oneway"]
            )
            lines.append("")
            lines.append(f"One-way segments: {oneway_count}")

        lines.append("=" * 70)
        summary = "\n".join(lines)
        if print:
            builtins.print(summary)
            return None
        return summary

    def plot(
        self,
        ax=None,
        column=None,
        color="black",
        linewidth=1.0,
        cmap="viridis",
        legend=True,
        equal_aspect=True,
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
        color : str, default "black"
            Line color when no column is given.
        linewidth : float, default 1.0
            Width of road lines.
        cmap : str, default "viridis"
            Matplotlib colormap for numeric or categorical columns.
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
        try:
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.lines import Line2D
        except ImportError as exc:
            raise ImportError(
                "matplotlib is required to plot RoadNetwork objects."
            ) from exc

        segments = self._plot_segments()
        if ax is None:
            _, ax = plt.subplots()

        if len(segments) == 0:
            warning("RoadNetwork has no road segments to plot.")
            return ax

        if column is None:
            collection = LineCollection(
                segments, colors=color, linewidths=linewidth, **kwargs
            )
            ax.add_collection(collection)
        else:
            values = self._plot_column_values(column, len(segments))
            numeric_values = self._as_numeric_plot_values(values)
            if numeric_values is not None:
                collection = LineCollection(
                    segments, linewidths=linewidth, cmap=cmap, **kwargs
                )
                collection.set_array(numeric_values)
                ax.add_collection(collection)
                if legend:
                    ax.figure.colorbar(collection, ax=ax, label=column)
            else:
                labels = [str(value) if value not in (None, "") else "N/A" for value in values]
                categories = sorted(set(labels))
                colormap = plt.get_cmap(cmap, max(len(categories), 1))
                color_lookup = {
                    category: colormap(i) for i, category in enumerate(categories)
                }
                colors = [color_lookup[label] for label in labels]
                collection = LineCollection(
                    segments, colors=colors, linewidths=linewidth, **kwargs
                )
                ax.add_collection(collection)
                if legend:
                    handles = [
                        Line2D([0], [0], color=color_lookup[category], lw=linewidth)
                        for category in categories
                    ]
                    ax.legend(handles, categories, title=column)

        ax.autoscale()
        if equal_aspect:
            ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        if show:
            plt.show()
        return ax

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

    @staticmethod
    def _as_numeric_plot_values(values):
        numeric_values = []
        for value in values:
            if value is None or value == "":
                numeric_values.append(np.nan)
                continue
            try:
                numeric_values.append(float(value))
            except (TypeError, ValueError):
                return None
        return np.asarray(numeric_values, dtype=float)

    def to_shapely(self):
        """
        Convert the road network geometry to a Shapely MultiLineString.

        Returns
        -------
        shapely.geometry.MultiLineString or None
            Shapely representation when geometry exists, else ``None``.
        """
        multilinestring = self.geometry.get(GeometryType.MULTILINESTRING)
        if multilinestring is None:
            return None
        return multilinestring.to_shapely()

    def to_proto(self):
        """
        Convert the road network to a protobuf Object message.

        Returns
        -------
        proto.Object
            Serialized road network including vertices, edges, and lengths.
        """
        pb = Object.to_proto(self)
        _pb = proto.RoadNetwork()
        dim = self.vertices.shape[1] if self.vertices.ndim == 2 else 0
        _pb.vertices.extend(self.vertices.flatten())
        _pb.dim = dim
        _pb.edges.extend(self.edges.flatten())
        _pb.lengths.extend(self.length.flatten())
        pb.road_network.CopyFrom(_pb)

        return pb

    def from_proto(self, pb):
        """
        Populate the road network from a protobuf Object message.

        Parameters
        ----------
        pb : proto.Object or bytes
            Protobuf message or serialized bytes containing a road network.
        """
        if isinstance(pb, bytes):
            pb = proto.Object.FromString(pb)
        Object.from_proto(self, pb)
        _pb = pb.road_network
        dim = _pb.dim
        if dim == 0:
            self.vertices = np.empty((0, 0))
        else:
            self.vertices = np.array(_pb.vertices).reshape(-1, dim)
        self.edges = np.array(_pb.edges).reshape(-1, 2)
        self.length = np.array(_pb.lengths)
