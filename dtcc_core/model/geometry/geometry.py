# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union
from abc import abstractmethod

from .. import dtcc_pb2 as proto
from ..model import Model
from ..values import Field
from .bounds import Bounds
from .transform import Transform


@dataclass
class Geometry(Model):
    """Base class for all geometry classes.

    Geometry classes represent geometric objects such as point clouds,
    surfaces, polygons, and meshes. They are used to represent the geometry of
    city objects.

    All geometries are stored in a local coordinate system, which may be
    different for each geometry. The transform attribute is used to transform
    the geometry from the local coordinate system to a global coordinate system.

    A geometry may have a list of fields which are values (scalars or vectors)
    defined on the entities of the geometry (vertices, edges, faces, etc.).

    Attributes
    ----------
    bounds : Bounds
        Bounding box of the geometry in the local coordinate system.
    transform : Transform
        Affine transform to a global coordinate system.
    fields: list[Field]
    """

    _bounds: Bounds = field(default_factory=Bounds)
    transform: Transform = field(default_factory=Transform)
    fields: list[Field] = field(default_factory=list)

    @abstractmethod
    def calculate_bounds(self):
        pass

    @property
    def bounds(self) -> Bounds:
        if self._bounds is None or self._bounds.area == 0:
            self.calculate_bounds()
        return self._bounds

    @bounds.setter
    def bounds(self, bounds: Bounds):
        self._bounds = bounds

    def add_field(self, field: Field):
        """Add a field to the geometry.

        Parameters
        ----------
        field : Field
            The field to add to the geometry.
        """
        self.fields.append(field)

    def to_proto(self) -> proto.Geometry:
        """Return a protobuf representation of the Geometry.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the Geometry.
        """
        pb = proto.Geometry()
        pb.bounds.CopyFrom(self.bounds.to_proto())
        pb.transform.CopyFrom(self.transform.to_proto())
        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes]):
        """Initialize Geometry from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Geometry, bytes]
            The protobuf message or its serialized bytes representation.
        """
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)
        self.bounds.from_proto(pb.bounds)
        self.transform.from_proto(pb.transform)
