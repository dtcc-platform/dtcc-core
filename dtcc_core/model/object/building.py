# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass
from numbers import Real
from typing import Literal

from shapely.ops import unary_union

from ..geometry import Surface
from ..logging import warning
from .object import GeometryType, Object


@dataclass(repr=False)
class Building(Object):
    """Represents a building in a city."""

    @property
    def building_parts(self):
        """Return list of building parts in building."""
        return self.children[BuildingPart] if BuildingPart in self.children else []

    @property
    def measured_height(self) -> float | None:
        """Stored measured height in metres, or None when absent.

        This is direct access to ``attributes["measured_height"]``. Validation
        occurs at save/load; neither estimates nor geometry are substituted.
        """
        return self.attributes.get("measured_height")

    @measured_height.setter
    def measured_height(self, value: float | None):
        if value is None:
            self.attributes.pop("measured_height", None)
        else:
            self.attributes["measured_height"] = value

    height = measured_height

    @property
    def estimated_height(self) -> float | None:
        """Stored modelling height in metres, or None; never a measurement fallback."""
        return self.attributes.get("estimated_height")

    @estimated_height.setter
    def estimated_height(self, value: float | None):
        if value is None:
            self.attributes.pop("estimated_height", None)
        else:
            self.attributes["estimated_height"] = value

    def footprint(
        self,
        geom_type: GeometryType | None = None,
        *,
        z: Literal["geometry", "ground"] | float = "geometry",
    ) -> Surface | None:
        """
        Extract a footprint surface from the building geometry.

        Parameters
        ----------
        geom_type : GeometryType, optional
            Geometry type to use. If omitted, only ``GeometryType.LOD0`` is
            considered because LOD0 is the canonical building footprint
            geometry. Passing another geometry type is an advanced/derived
            extraction path and only that explicit type is considered.
        z : {"geometry", "ground"} or float, default "geometry"
            Height for the returned footprint surface. ``"geometry"`` uses the
            source geometry ``zmax`` and preserves previous behavior.
            ``"ground"`` uses the source geometry ``zmin``. A numeric value
            places the footprint at that exact z height.

        Returns
        -------
        Surface or None
            The footprint as a surface, or ``None`` when the building has no
            usable LOD geometry.
        """
        selected_geom_type = GeometryType.LOD0 if geom_type is None else geom_type
        geom = self.flatten_geometry(selected_geom_type)
        if geom is not None and not _is_footprint_compatible(geom):
            warning(f"Building {self.id} geometry cannot produce a footprint.")
            return None

        if geom is None:
            warning(f"Building {self.id} has no {selected_geom_type.name} geometry.")
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
        surface.from_polygon(footprint, _resolve_footprint_z(geom, z))
        return surface


class BuildingPart(Object):
    """Represents a building part object with protobuf serialization support.

    A specialized Object subclass that provides conversion methods for
    protobuf serialization and deserialization of building part data.
    """

    building_parts = Building.building_parts
    measured_height = Building.measured_height
    height = Building.height
    estimated_height = Building.estimated_height


def _resolve_footprint_z(
    geom,
    z: Literal["geometry", "ground"] | float,
) -> float:
    if z == "geometry":
        return float(geom.bounds.zmax)
    if z == "ground":
        return float(geom.bounds.zmin)
    if isinstance(z, Real) and not isinstance(z, bool):
        return float(z)
    raise ValueError("z must be 'geometry', 'ground', or a numeric height.")


def _is_footprint_compatible(geom) -> bool:
    return geom is not None and hasattr(geom, "to_polygon")
