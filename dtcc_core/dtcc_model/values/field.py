# Copyright(C) 2024 Anders Logg
# Licensed under the MIT License


import numpy as np
from typing import Union
from dataclasses import dataclass, field

from ..model import Model
from .. import dtcc_pb2 as proto


@dataclass
class Field(Model):
    """Represents a field (scalar or vector) defined on a geometry.

    A field has a value and a unit of measurement and can represent e.g.
    physical properties such as a temperature or a velocity vector field. The
    field takes a value on each element of a geometry and the value may be
    either scalar or vector-valued.


    Attributes
    ----------
    name: str
        Name of the field.
    unit: str
        Unit of measurement of the field.
    description: str
        Description of the field.
    values: np.ndarray
        An array of values (scalar or vector-valued) of the field. The
        dimension is n x d, where n is the number of elements in the geometry
        and d is the dimension of the field.
    dim: int
        The dimension of the field.
    """

    name: str = ""
    unit: str = ""
    description: str = ""
    values: np.ndarray = field(default_factory=lambda: np.empty(0))
    dim: int = 0

    def to_proto(self) -> proto.Field:
        """Return a protobuf representation of the Field.

        Returns
        -------
        proto.Field
            A protobuf representation of the Field.
        """

        pb = proto.Field()
        pb.name = self.name
        pb.unit = self.unit
        pb.description = self.description
        pb.values.extend(self.values.flatten())
        pb.dim = self.dim

        return pb

    def from_proto(self, pb: Union[proto.Field, bytes]):
        """Initialize Field from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Field, bytes]
            The protobuf message or its serialized bytes representation.
        """

        if isinstance(pb, bytes):
            pb = proto.Field.FromString(pb)
        self.name = pb.name
        self.unit = pb.unit
        self.description = pb.description
        self.values = np.array(pb.values).reshape((-1, pb.dim))
        self.dim = pb.dim
