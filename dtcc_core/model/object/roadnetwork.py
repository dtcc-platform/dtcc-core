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
