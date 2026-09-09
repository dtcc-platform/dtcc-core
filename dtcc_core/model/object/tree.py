# Copyright(C) 2025 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
import math
import numpy as np
from typing import Union
from numbers import Real

from .object import Object
from .. import dtcc_pb2 as proto


def _tree_values(position, height, crown_radius):
    """Validate the scalar/tree-position contract before writing or reading it."""
    try:
        position = np.asarray(position)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Tree.position must contain three numeric XYZ coordinates or be empty"
        ) from exc
    if position.shape not in ((3,), (0,), (0, 3)):
        raise ValueError("Tree.position must contain three XYZ coordinates or be empty")
    if position.dtype.kind not in "fiu" or not np.all(np.isfinite(position)):
        raise ValueError("Tree.position must contain finite numeric XYZ coordinates")
    float_limit = float(np.finfo(np.float32).max)
    if np.any(np.abs(position) > float_limit):
        raise ValueError("Tree.position exceeds the protobuf float32 range")
    for name, value in (("height", height), ("crown_radius", crown_radius)):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, Real)
            or value < 0
            or value > float_limit
            or not math.isfinite(value)
        ):
            raise ValueError(f"Tree.{name} must be a finite, nonnegative float32 number")
    if position.size == 0:
        # An unlocated tree is the intentional default of Tree().
        position = np.empty((0, 3))
    return position, float(height), float(crown_radius)


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
        position, height, crown_radius = _tree_values(
            self.position, self.height, self.crown_radius
        )
        pb = Object.to_proto(self)
        tree_pb = proto.Tree(height=height, crown_radius=crown_radius)
        tree_pb.position.extend(position.flatten())
        pb.tree.CopyFrom(tree_pb)
        return pb

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """
        Populate the tree from a protobuf Object message.

        Parameters
        ----------
        pb : proto.Object or bytes
            Protobuf message or serialized bytes containing a tree.
        """
        if isinstance(pb, bytes):
            pb = proto.Object.FromString(pb)
        if pb.WhichOneof("type") != "tree":
            raise ValueError("Tree.from_proto requires an Object message of type tree")
        position, height, crown_radius = _tree_values(
            np.asarray(pb.tree.position), pb.tree.height, pb.tree.crown_radius
        )
        Object.from_proto(self, pb)
        self.position = position
        self.height = height
        self.crown_radius = crown_radius
