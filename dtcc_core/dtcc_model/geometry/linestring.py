# Copyright(C) 2024 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union
from .geometry import Geometry, Bounds
import numpy as np

from shapely.geometry import (
    LineString as ShapelyLineString,
    MultiLineString as ShapelyMultiLineString,
)

from .. import dtcc_pb2 as proto


@dataclass
class LineString(Geometry):
    vertices: np.ndarray = field(default_factory=lambda: np.empty(0))

    def calculate_bounds(self):
        """Calculate the bounding box of the line string."""
        if len(self.vertices) == 0:
            self._bounds = Bounds()
            return self._bounds
        xmin = np.min(self.vertices[:, 0])
        ymin = np.min(self.vertices[:, 1])
        xmax = np.max(self.vertices[:, 0])
        ymax = np.max(self.vertices[:, 1])
        if self.vertices.shape[1] == 3:
            zmin = np.min(self.vertices[:, 2])
            zmax = np.max(self.vertices[:, 2])
        else:
            zmin = 0
            zmax = 0

        self._bounds = Bounds(
            xmin=xmin,
            ymin=ymin,
            zmin=zmin,
            xmax=xmax,
            ymax=ymax,
            zmax=zmax,
        )
        return self._bounds

    @property
    def length(self):
        """Calculate the length of the line string."""
        return np.sum(np.linalg.norm(np.diff(self.vertices, axis=0), axis=1))

    def to_shapely(self):
        """Convert the LineString to a Shapely LineString."""
        return ShapelyLineString(self.vertices)

    def from_shapely(self, shape):
        """Initialize the LineString from a Shapely LineString."""
        self.vertices = np.array(shape.coords)
        return self

    def to_proto(self) -> proto.Geometry:
        """Return a protobuf representation of the LineString.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the MultiSurface as a Geometry.
        """

        # Handle Geometry fields
        pb = Geometry.to_proto(self)
        _pb = proto.LineString()
        _pb.vertices.extend(self.vertices.flatten())
        _pb.dim = self.vertices.shape[1]
        pb.line_string.CopyFrom(_pb)
        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes], only_linestring_field=False):
        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)

        # Handle Geometry fields
        if not only_linestring_field:
            Geometry.from_proto(self, pb)

        _pb = pb if only_linestring_field else pb.line_string
        dim = _pb.dim
        self.vertices = np.array(_pb.vertices).reshape(-1, dim)


@dataclass
class MultiLineString(Geometry):
    linestrings: list[LineString] = field(default_factory=lambda: [])

    def calculate_bounds(self):
        """Calculate the bounding box of the multi line string."""
        if len(self.linestrings) == 0:
            return Bounds()
        bounds = [line.calculate_bounds() for line in self.linestrings]
        self._bounds = Bounds(
            xmin=np.min([b.xmin for b in bounds]),
            ymin=np.min([b.ymin for b in bounds]),
            zmin=np.min([b.zmin for b in bounds]),
            xmax=np.max([b.xmax for b in bounds]),
            ymax=np.max([b.ymax for b in bounds]),
            zmax=np.max([b.zmax for b in bounds]),
        )
        return self._bounds

    @property
    def length(self):
        """Calculate the length of the multi line string."""
        return sum([line.length for line in self.linestrings])

    def to_shapely(self):
        """Convert the MultiLineString to a Shapely MultiLineString."""
        return ShapelyMultiLineString([line.to_shapely() for line in self.linestrings])

    def from_shapely(self, shape):
        """Initialize the MultiLineString from a Shapely MultiLineString."""
        self.linestrings = [LineString().from_shapely(l) for l in shape]
        return self

    def to_proto(self):
        """Return a protobuf representation of the MultiSurface.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the MultiSurface as a Geometry.
        """

        # Handle Geometry fields
        pb = Geometry.to_proto(self)
        _pb = proto.MultiLineString()
        _pb.line_strings.extend(
            [line.to_proto().line_string for line in self.linestrings]
        )
        pb.multi_line_string.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes]):
        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)

        # Handle Geometry fields
        Geometry.from_proto(self, pb)
        # Handle specific fields
        _pb = pb.multi_line_string
        for line_string in _pb.line_strings:
            _linestring = LineString()
            _linestring.from_proto(line_string, only_linestring_field=True)
            self.linestrings.append(_linestring)
