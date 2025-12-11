# Copyright(C) 2025 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
import numpy as np
from typing import Union

from .object import Object
from .. import dtcc_pb2 as proto


@dataclass
class Tree(Object):
    """
    Represents a single tree with minimal geometric attributes.

    Attributes
    ----------
    position : numpy.ndarray
        XYZ coordinate of the tree, expected as a 3-element array in the
        project's spatial reference system.
    height : float
        Height of the tree crown apex above ground, in meters.
    crown_radius : float
        Plan-view radius of the tree crown, in meters.
    """
    position: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    height: float = 0.0
    crown_radius: float = 0.0

    def to_proto(self) -> proto.Object:
        """
        Convert the tree to a protobuf Object message.

        Returns
        -------
        proto.Object
            Serialized tree representation.
        """
        pass

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """
        Populate the tree from a protobuf Object message.

        Parameters
        ----------
        pb : proto.Object or bytes
            Protobuf message or serialized bytes containing a tree.
        """
        pass
