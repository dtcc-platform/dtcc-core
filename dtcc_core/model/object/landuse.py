from dataclasses import dataclass, field
from typing import Union, List, Tuple
from enum import Enum, auto
from .object import Object, GeometryType
from ..geometry import Surface, MultiSurface
from ..geometry import Bounds
from .. import dtcc_pb2 as proto

import numpy as np


class LanduseClasses(Enum):
    """
    Enumeration of land use classification types for geographic or urban modeling applications.

    This enum defines various categories of land cover or land utilization, useful in GIS, 
    remote sensing, simulation, or 3D urban environments. Each class represents a general 
    type of surface usage, ranging from natural to heavily developed areas.

    Attributes:
        WATER: Bodies of water such as lakes, rivers, and oceans.
        GRASS: Grass-covered areas including parks, fields, and lawns.
        FOREST: Forested or densely vegetated regions.
        FARMLAND: Agricultural land used for crops or pasture.
        LIGHT_URBAN: Low-density urban development, such as suburbs or residential zones.
        URBAN: Medium-density urban areas with mixed land use.
        HEAVY_URBAN: High-density urban cores, often with large buildings or infrastructure.
        INDUSTRIAL: Zones designated for industrial activity, including factories and warehouses.
        MILITARY: Areas reserved for military use or restricted access.
        ROAD: Land primarily used for roads and transportation corridors.
        RAIL: Railway infrastructure and corridors.
        UNKNOWN: Land use is unknown or unclassified (explicitly set to 9999).
    """
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
    """
    Represents a land use object with associated land use classifications and geometric surfaces.

    This class models areas of land with specified types of usage (e.g., forest, urban, water),
    along with their associated geometric representation as multi-surfaces.

    Attributes:
        landuses (List[LanduseClasses]): A list of land use classes describing how the land is used.
    """
    landuses: List[LanduseClasses] = field(default_factory=list)

    @property
    def surfaces(self) -> List[Surface]:
        """
        Access the list of surfaces representing land use polygons.

        Returns
        -------
        list[Surface]
            Surfaces stored under the MultiSurface geometry, or an empty list.
        """
        geom = self.geometry.get(GeometryType.MULTISURFACE)
        if geom is None:
            return []
        return geom.surfaces

    def to_proto(self) -> proto.Object:
        """
        Convert the Landuse object to a protobuf representation.

        Returns
        -------
        proto.Object
            Protobuf message encoding the land use data.
        """
        pass

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """
        Populate the Landuse object from a protobuf message.

        Parameters
        ----------
        pb : proto.Object or bytes
            Protobuf message or serialized bytes containing land use data.
        """
        pass
