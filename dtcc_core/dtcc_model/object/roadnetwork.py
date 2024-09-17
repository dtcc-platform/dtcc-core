from dataclasses import dataclass, field
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
    vertices: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    edges: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )  # each edge (start_idx, end_idx)
    length: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))

    @property
    def linestrings(self) -> List[LineString]:
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
        if geom is None:
            return []
        return geom.linestrings

    @property
    def multilinestrings(self) -> MultiLineString:
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
        return geom

    @property
    def bounds(self):
        geom = self.geometry.get(GeometryType.MULTILINESTRING)
        if geom is None:
            xmin = np.min(self.vertices[:, 0])
            ymin = np.min(self.vertices[:, 1])
            xmax = np.max(self.vertices[:, 0])
            ymax = np.max(self.vertices[:, 1])
            return Bounds(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax)
        return geom.bounds

    def to_shapely(self):
        multilinestring = self.geometry.get(GeometryType.MULTILINESTRING)
        if multilinestring is None:
            return None
        return multilinestring.to_shapely()

    def to_proto(self):
        pb = Object.to_proto(self)
        _pb = proto.RoadNetwork()
        dim = self.vertices.shape[1]
        _pb.vertices.extend(self.vertices.flatten())
        _pb.dim = dim
        _pb.edges.extend(self.edges.flatten())
        _pb.lengths.extend(self.length.flatten())
        pb.road_network.CopyFrom(_pb)

        return pb

    def from_proto(self, pb):
        if isinstance(pb, bytes):
            pb = proto.Object.FromString(pb)
        Object.from_proto(self, pb)
        _pb = pb.road_network
        dim = _pb.dim
        self.vertices = np.array(_pb.vertices).reshape(-1, dim)
        self.edges = np.array(_pb.edges).reshape(-1, 2)
        self.length = np.array(_pb.lengths)
