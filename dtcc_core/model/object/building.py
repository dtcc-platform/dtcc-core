# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union
from shapely.geometry import Polygon
from shapely.ops import unary_union

from .object import Object, GeometryType
from ..geometry import Bounds, Surface
from .. import dtcc_pb2 as proto
from ..logging import warning


@dataclass
class Building(Object):
    """Represents a building in a city."""

    @property
    def building_parts(self):
        """Return list of building parts in building."""
        return self.children[BuildingPart] if BuildingPart in self.children else []

    @property
    def height(self):
        """
        Get the height of the building.
        
        Returns the height from the building's attributes, or calculates it
        from the bounds if no height attribute is set.
        
        Returns
        -------
        float
            The height of the building in meters.
        """
        height = self.attributes.get("height", None)
        if height is None:
            height = self.bounds.zmax - self.bounds.zmin
        return height

    def footprint(self, geom_type: GeometryType = None) -> Surface | None:
        """
        Extract a footprint surface from the building geometry.

        Parameters
        ----------
        geom_type : GeometryType, optional
            Specific geometry type to use. If omitted, the highest available
            LOD geometry is used.

        Returns
        -------
        Surface or None
            The footprint as a surface, or ``None`` when the building has no
            usable LOD geometry.
        """
        lod_levels = [
            GeometryType.LOD0,
            GeometryType.LOD1,
            GeometryType.LOD2,
            GeometryType.LOD3,
        ]

        geom = self.flatten_geometry(geom_type) if geom_type is not None else None
        if geom is None:
            for lod in lod_levels:
                geom = self.flatten_geometry(lod)
                if geom is not None:
                    break

        if geom is None:
            warning(f"Building {self.id} has no LOD geometry.")
            return None

        footprint = geom.to_polygon()
        if footprint is None or footprint.is_empty:
            warning(f"Building {self.id} has no footprint polygon.")
            return None
        if footprint.geom_type == "MultiPolygon":
            merged = unary_union(footprint.geoms)
            if merged.geom_type == "MultiPolygon":
                footprint = max(merged.geoms, key=lambda p: p.area)
            else:
                footprint = merged
        if footprint.geom_type != "Polygon":
            warning(f"Building {self.id} footprint is not a polygon.")
            return None

        surface = Surface()
        surface.from_polygon(footprint, geom.bounds.zmax)
        return surface

    def to_proto(self) -> proto.Object:
        """Return a protobuf representation of the Building.

        Returns
        -------
        proto.Object
            A protobuf representation of the Building as an Object.
        """

        # Handle Object fields
        pb = Object.to_proto(self)

        # Handle specific fields (currently none)
        _pb = proto.Building()
        pb.building.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """Initialize Building from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Object, bytes]
            The protobuf message or its serialized bytes representation.
        """

        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Object.FromString(pb)

        # Handle Object fields
        Object.from_proto(self, pb)

        # Handle specific fields (currently none)
        pass


class BuildingPart(Object):
    """Represents a building part object with protobuf serialization support.
    
        A specialized Object subclass that provides conversion methods for 
        protobuf serialization and deserialization of building part data.
    """
    def to_proto(self) -> proto.Object:
        """Return a protobuf representation of the BuildingPart.

        Returns
        -------
        proto.Object
            A protobuf representation of the BuildingPart as an Object.
        """

        # Handle Object fields
        pb = Object.to_proto(self)

        # Handle specific fields (currently none)
        _pb = proto.BuildingPart()
        pb.building_part.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """Initialize BuildingPart from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Object, bytes]
            The protobuf message or its serialized bytes representation.
        """

        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Object.FromString(pb)

        # Handle Object fields
        Object.from_proto(self, pb)

        # Handle specific fields (currently none)
        pass
