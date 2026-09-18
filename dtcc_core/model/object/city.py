# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Literal, Union

from .. import geometry
from ..geometry import Bounds
from ..logging import debug, error, info, warning
from ..mixins.city import (
    CityBuilderMixin,
    CityDownloadMixin,
    CityLoaderMixin,
    CityModifyingMixin,
    CitySaveMixin,
)
from ..values.raster import Raster
from .building import Building
from .dataset_collections import BuildingCollection, FootprintCollection
from .object import GeometryType, Object
from .terrain import Terrain
from .tree import Tree


@dataclass(repr=False)
class City(
    CityLoaderMixin,
    CitySaveMixin,
    CityDownloadMixin,
    CityBuilderMixin,
    CityModifyingMixin,
    Object,
):
    """Represents a city, the top-level container class for city models."""

    @property
    def buildings(self) -> list[Building]:
        """Return list of buildings in city."""
        return self.children[Building] if Building in self.children else []

    @property
    def terrain(self):
        """Return terrain in city."""
        if Terrain in self.children:
            return self.children[Terrain][0]
        else:
            return Terrain()

    @property
    def trees(self):
        """Return trees in city."""
        return self.children[Tree] if Tree in self.children else []

    def has_terrain(self) -> bool:
        """
        Check whether the city has a terrain child.

        Returns
        -------
        bool
            ``True`` if a terrain object is present, otherwise ``False``.
        """
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
        """
        Collect attributes from all buildings in the city.

        Returns
        -------
        dict[str, list]
            Mapping of all attribute names to building-aligned lists. A missing
            attribute is represented by None, and keys follow first-seen order.
        """

        city_buildings = self.buildings
        if len(city_buildings) == 0:
            return {}

        building_attributes = defaultdict(list)
        attribute_keys = dict.fromkeys(
            key for building in city_buildings for key in building.attributes
        )
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
        """Remove any terrain objects from the city."""
        if Terrain in self.children:
            self.children[Terrain] = []

    def remove_buildings(self):
        """Remove all buildings from the city."""
        if Building in self.children:
            self.children[Building] = []

    def add_building(self, building: Building):
        """Add building to city."""
        self.add_child(building)

    def building_collection(self) -> BuildingCollection:
        """Return city buildings as a semantic collection."""
        return BuildingCollection(self.buildings)

    def building_footprints(
        self,
        geom_type: GeometryType | None = None,
        *,
        z: Literal["geometry", "ground"] | float = "geometry",
    ) -> FootprintCollection:
        """Return city building footprints as a semantic collection.

        By default, this extracts canonical LOD0 footprints only. Passing
        ``geom_type`` requests an advanced/derived extraction from that explicit
        geometry type.
        """
        return self.building_collection().footprints(geom_type, z=z)

    def replace_buildings(self, buildings: list[Building]):
        """Replace all buildings in city with new list of buildings."""
        self.remove_buildings()
        self.add_buildings(buildings)

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

    def add_trees(self, trees: list[Tree]):
        """Add trees to city
        args:
            trees: list[Tree]
                List of tree objects to add to the city.
        """
        self.add_children(trees)


@dataclass(repr=False)
class CityObject(Object):
    """Represents a generic object in a city."""
