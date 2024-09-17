from dataclasses import dataclass, field
from typing import Union, List, Tuple
from enum import Enum, auto
from .object import Object, GeometryType
from dtcc_model.geometry import Surface, MultiSurface
from dtcc_model.geometry import Bounds
from dtcc_model import dtcc_pb2 as proto

import numpy as np


class LanduseClasses(Enum):
    WATER = auto()
    GRASS = auto()
    FOREST = auto()
    FARMLAND = auto()
    LIGHT_URBAN = auto()
    URBAN = auto()
    HEAVY_URBAN = auto()
    INDUSTRIAL = auto()
    MILITARY = auto()
    ROAD = auto()
    RAIL = auto()
    UNKNOWN = 9999


@dataclass
class Landuse(Object):

    landuses: List[LanduseClasses] = field(default_factory=list)

    @property
    def surfaces(self) -> List[Surface]:
        geom = self.geometry.get(GeometryType.MULTISURFACE)
        if geom is None:
            return []
        return geom.surfaces

    def to_proto(self) -> proto.Object:
        pass

    def from_proto(self, pb: Union[proto.Object, bytes]):
        pass
