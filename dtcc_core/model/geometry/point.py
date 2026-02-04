# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass
from typing import Union
import numpy as np

from .geometry import Geometry
from .bounds import Bounds
from .. import dtcc_pb2 as proto


@dataclass
class Point(Geometry):
    """Represents a single point in 3D space.

    Attributes
    ----------
    x : float
        The x-coordinate of the point.
    y : float
        The y-coordinate of the point.
    z : float
        The z-coordinate of the point (default: 0.0).
    """

    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def __str__(self):
        """Return a string representation of the Point.

        Returns
        -------
        str
            A string representation of the Point.
        """
        return f"DTCC Point at ({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"

    def __repr__(self):
        return f"Point(x={self.x}, y={self.y}, z={self.z})"

    def calculate_bounds(self):
        """Calculate the bounds of the point and update the bounds attribute.

        For a single point, the min and max bounds are the same.
        """
        self._bounds = Bounds(
            xmin=self.x,
            xmax=self.x,
            ymin=self.y,
            ymax=self.y,
            zmin=self.z,
            zmax=self.z,
        )

    def offset(self, dx: float = 0.0, dy: float = 0.0, dz: float = 0.0):
        """Offset the point by the given amounts.

        Parameters
        ----------
        dx : float
            Offset in x-direction.
        dy : float
            Offset in y-direction.
        dz : float
            Offset in z-direction.
        """
        self.x += dx
        self.y += dy
        self.z += dz
        self.calculate_bounds()

    def to_proto(self) -> proto.Geometry:
        """Return a protobuf representation of the Point.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the Point as a Geometry.
        """
        # Handle Geometry fields (bounds, transform, fields)
        pb = Geometry.to_proto(self)

        # Handle specific fields for Point
        _pb = proto.Point()
        _pb.x = self.x
        _pb.y = self.y
        _pb.z = self.z
        pb.point.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes]):
        """Initialize Point from a protobuf representation.

        Parameters
        ----------
        pb : Union[proto.Geometry, bytes]
            The protobuf message or its serialized bytes representation.
        """
        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)

        # Handle Geometry fields (bounds, transform, fields)
        Geometry.from_proto(self, pb)

        # Handle specific fields for Point
        _pb = pb.point
        self.x = _pb.x
        self.y = _pb.y
        self.z = _pb.z
