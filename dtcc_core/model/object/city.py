# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union
from collections import defaultdict

from .object import Object, GeometryType
from .building import Building
from .terrain import Terrain
from ..values.raster import Raster
from .. import geometry
from ..geometry import Bounds

from .. import dtcc_pb2 as proto
from ..logging import info, warning, error, debug


@dataclass
class City(Object):
    """Represents a city, the top-level container class for city models."""

    @property
    def buildings(self):
        """Return list of buildings in city."""
        return self.children[Building] if Building in self.children else []

    @property
    def terrain(self):
        """Return terrain in city."""
        if Terrain in self.children:
            return self.children[Terrain][0]
        else:
            return Terrain()

    def has_terrain(self) -> bool:
        return Terrain in self.children

    @property
    def num_buildings(self):
        """Return number of buildings in city."""
        return len(self.buildings)

    def get_building_attribute(self, attribute):
        """Return list of values for a specific building attribute."""
        return [b.attributes.get(attribute, None) for b in self.buildings]

    def set_building_attribute(self, attribute, values):
        """Set specific building attribute for all buildings."""
        if len(values) != self.num_buildings:
            raise ValueError(
                f"Number of values ({len(values)}) does not match number of buildings ({self.num_buildings})"
            )
        for b, v in zip(self.buildings, values):
            b.attributes[attribute] = v

    def get_building_attributes(self):

        city_buildings = self.buildings
        if len(city_buildings) == 0:
            return {}

        building_attributes = defaultdict(list)
        # assuming all buildings have the same attributes
        # TODO: handle buildings with different attributes
        attribute_keys = city_buildings[0].attributes.keys()
        for b in city_buildings:
            for key in attribute_keys:
                building_attributes[key].append(b.attributes.get(key, None))
        return dict(building_attributes)

    def add_terrain(self, terrain):
        """Add terrain to city."""
        if isinstance(terrain, Terrain):
            self.add_child(terrain)
        else:
            terrain_object = Terrain()
            if isinstance(terrain, Raster):
                terrain_object.add_geometry(terrain, GeometryType.RASTER)
            elif isinstance(terrain, geometry.Mesh):
                terrain_object.add_geometry(terrain, GeometryType.MESH)
            else:
                raise ValueError(f"Invalid terrain type {type(terrain)}.")
            self.add_child(terrain_object)

    def remove_terrain(self):
        if Terrain in self.children:
            self.children[Terrain] = []

    def remove_buildings(self):
        if Building in self.children:
            self.children[Building] = []

    def add_building(self, building: Building):
        """Add building to city."""
        self.add_child(building)

    def add_buildings(
        self, buildings: list[Building], remove_outside_terrain: bool = False
    ):
        """Add building to city.
        args:
            buildings: list[Building]
                List of buildings to add to the city.
            remove_outside_terrain: bool
                If True, remove buildings that are outside the terrain.
        """
        if remove_outside_terrain:
            initial_count = len(buildings)
            terrain = self.terrain
            if terrain is None:
                warning("City has no terrain. Cannot remove buildings outside terrain.")
            else:
                terrain_bounds = terrain.bounds
                buildings = [
                    b for b in buildings if terrain_bounds.contains_bounds(b.bounds)
                ]
                info(
                    f"Removed {initial_count - len(buildings)} buildings outside terrain."
                )
        self.add_children(buildings)

    def to_proto(self) -> proto.Object:
        """Return a protobuf representation of the City.

        Returns
        -------
        proto.Object
            A protobuf representation of the City as an Object.
        """

        # Handle Object fields
        pb = Object.to_proto(self)

        # Handle specific fields (currently none)
        _pb = proto.City()
        pb.city.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """Initialize City from a protobuf representation.

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


@dataclass
class CityObject(Object):
    """Represents a generic object in a city."""

    def to_proto(self) -> proto.Object:
        """Return a protobuf representation of the CityObject.

        Returns
        -------
        proto.Object
            A protobuf representation of the CityObject as an Object.
        """

        # Handle Object fields
        pb = Object.to_proto(self)

        # Handle specific fields (currently none)
        _pb = proto.CityObject()
        pb.city.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """Initialize CityObject from a protobuf representation.

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
