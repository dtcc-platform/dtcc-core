# Copyright(C) 2025 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
import numpy as np
from typing import Union

from .object import Object
from .. import dtcc_pb2 as proto


@dataclass
class Tree(Object):
    position: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    height: float = 0.0
    crown_radius: float = 0.0

    def to_proto(self) -> proto.Object:
        pass

    def from_proto(self, pb: Union[proto.Object, bytes]):
        pass
