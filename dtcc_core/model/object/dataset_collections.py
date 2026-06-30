"""Semantic Dataset v2 model collections for city-domain returns."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from ..geometry import Bounds, Surface
from ..model import Model
from .building import Building
from .object import GeometryType
from .tree import Tree


@dataclass
class FootprintCollection(Model):
    """Collection of building footprint surfaces."""

    footprints: list[Surface] = field(default_factory=list)
    source_ids: list[str | None] = field(default_factory=list)
    source_indices: list[int] = field(default_factory=list)

    @classmethod
    def from_buildings(
        cls,
        buildings: list[Building],
        geom_type: GeometryType | None = None,
        *,
        z: Literal["geometry", "ground"] | float = "geometry",
    ) -> "FootprintCollection":
        """Build a footprint collection from buildings with available geometry."""
        footprints = []
        source_ids: list[str | None] = []
        source_indices: list[int] = []
        for source_index, building in enumerate(buildings):
            footprint = building.footprint(geom_type, z=z)
            if footprint is not None:
                footprints.append(footprint)
                source_ids.append(getattr(building, "id", None))
                source_indices.append(source_index)
        return cls(
            footprints=footprints,
            source_ids=source_ids,
            source_indices=source_indices,
        )

    def __len__(self) -> int:
        return len(self.footprints)

    def __iter__(self) -> Iterator[Surface]:
        return iter(self.footprints)

    def __getitem__(self, index):
        return self.footprints[index]

    def to_list(self) -> list[Surface]:
        """Return footprints as a plain list."""
        return list(self.footprints)

    @property
    def bounds(self) -> Bounds:
        """Return bounds spanning all footprints."""
        return _bounds_for_models(self.footprints)

    def to_arrays(self) -> list[np.ndarray]:
        """Return footprint vertex arrays."""
        return [footprint.vertices.copy() for footprint in self.footprints]

    def to_shapely(self):
        """Return footprints as Shapely polygons."""
        return [footprint.to_polygon(simplify=0.0) for footprint in self.footprints]

    def to_geojson(self, crs: str | None = None) -> dict[str, Any]:
        """Return footprints as a GeoJSON FeatureCollection."""
        features: list[dict[str, Any]] = []
        for index, polygon in enumerate(self.to_shapely()):
            if polygon is None or polygon.is_empty:
                continue
            features.append(
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": _polygon_coordinates(polygon),
                    },
                    "properties": self._feature_properties(index),
                }
            )

        payload: dict[str, Any] = {
            "type": "FeatureCollection",
            "features": features,
        }
        if crs is not None:
            payload["crs"] = {
                "type": "name",
                "properties": {"name": crs},
            }
        return payload

    def to_proto(self):
        raise NotImplementedError(
            "FootprintCollection protobuf serialization is not implemented."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "FootprintCollection protobuf deserialization is not implemented."
        )

    def _feature_properties(self, index: int) -> dict[str, Any]:
        properties: dict[str, Any] = {"index": index}
        if index < len(self.source_indices):
            properties["source_index"] = self.source_indices[index]
        if index < len(self.source_ids) and self.source_ids[index] is not None:
            properties["source_id"] = self.source_ids[index]
        return properties


@dataclass
class BuildingCollection(Model):
    """Collection of DTCC buildings."""

    buildings: list[Building] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.buildings)

    def __iter__(self) -> Iterator[Building]:
        return iter(self.buildings)

    def __getitem__(self, index):
        return self.buildings[index]

    def to_list(self) -> list[Building]:
        """Return buildings as a plain list."""
        return list(self.buildings)

    @property
    def bounds(self) -> Bounds:
        """Return bounds spanning all buildings."""
        return _bounds_for_models(self.buildings)

    def footprints(
        self,
        geom_type: GeometryType | None = None,
        *,
        z: Literal["geometry", "ground"] | float = "geometry",
    ) -> FootprintCollection:
        """Return building footprints as a semantic collection."""
        return FootprintCollection.from_buildings(self.buildings, geom_type, z=z)

    def to_proto(self):
        raise NotImplementedError(
            "BuildingCollection protobuf serialization is not implemented."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "BuildingCollection protobuf deserialization is not implemented."
        )


@dataclass
class TreeCollection(Model):
    """Collection of DTCC tree objects."""

    trees: list[Tree] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.trees)

    def __iter__(self) -> Iterator[Tree]:
        return iter(self.trees)

    def __getitem__(self, index):
        return self.trees[index]

    def to_list(self) -> list[Tree]:
        """Return trees as a plain list."""
        return list(self.trees)

    def to_arrays(self) -> np.ndarray:
        """Return tree positions as an ``N x 3`` array when available."""
        positions = []
        for tree in self.trees:
            position = np.asarray(tree.position, dtype=float)
            if position.size < 3:
                continue
            positions.append(position.reshape(-1, 3)[0])
        if not positions:
            return np.empty((0, 3))
        return np.asarray(positions, dtype=float)

    def to_proto(self):
        raise NotImplementedError(
            "TreeCollection protobuf serialization is not implemented."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "TreeCollection protobuf deserialization is not implemented."
        )


@dataclass
class CalibrationGrid(Model):
    """Semantic model for a synthetic calibration grid."""

    bounds: Bounds = field(default_factory=Bounds)
    divisions: int = 0
    crs: str | None = None
    features: list[dict[str, Any]] = field(default_factory=list)
    name: str = "calibration_grid"
    metadata_payload: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_geojson(cls, geojson: dict[str, Any]) -> "CalibrationGrid":
        """Build a calibration grid model from its GeoJSON representation."""
        metadata = dict(geojson.get("metadata") or {})
        bounds_values = metadata.get("bounds") or [0.0, 0.0, 0.0, 0.0]
        bounds = Bounds(
            xmin=float(bounds_values[0]),
            ymin=float(bounds_values[1]),
            xmax=float(bounds_values[2]),
            ymax=float(bounds_values[3]),
        )
        crs = metadata.get("crs")
        if crs is None:
            crs = (
                geojson.get("crs", {})
                .get("properties", {})
                .get("name")
            )
        return cls(
            bounds=bounds,
            divisions=int(metadata.get("divisions") or 0),
            crs=crs,
            features=list(geojson.get("features") or []),
            name=str(geojson.get("name") or "calibration_grid"),
            metadata_payload=metadata,
        )

    def __len__(self) -> int:
        return len(self.features)

    def __iter__(self) -> Iterator[dict[str, Any]]:
        return iter(self.features)

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            return self.features[key]
        return self.to_geojson()[key]

    def __contains__(self, key) -> bool:
        return key in self.to_geojson()

    def get(self, key, default=None):
        """Return GeoJSON member by key with mapping-style fallback."""
        return self.to_geojson().get(key, default)

    def keys(self):
        """Return GeoJSON member keys."""
        return self.to_geojson().keys()

    def items(self):
        """Return GeoJSON member items."""
        return self.to_geojson().items()

    def values(self):
        """Return GeoJSON member values."""
        return self.to_geojson().values()

    def to_python(self) -> dict[str, Any]:
        """Return the GeoJSON representation."""
        return self.to_geojson()

    def to_geojson(self) -> dict[str, Any]:
        """Return the calibration grid as a GeoJSON FeatureCollection."""
        metadata = (
            dict(self.metadata_payload)
            if self.metadata_payload
            else self._default_metadata()
        )
        payload: dict[str, Any] = {
            "type": "FeatureCollection",
            "name": self.name,
            "features": list(self.features),
            "metadata": metadata,
        }
        if self.crs is not None:
            payload["crs"] = {
                "type": "name",
                "properties": {"name": self.crs},
            }
            payload["metadata"].setdefault("crs", self.crs)
        return payload

    def _default_metadata(self) -> dict[str, Any]:
        spacing = [0.0, 0.0]
        if self.divisions:
            spacing = [
                self.bounds.width / self.divisions,
                self.bounds.height / self.divisions,
            ]
        metadata = {
            "dataset": "calibration_grid",
            "divisions": self.divisions,
            "line_count": len(self.features),
            "spacing": spacing,
            "bounds": [
                self.bounds.xmin,
                self.bounds.ymin,
                self.bounds.xmax,
                self.bounds.ymax,
            ],
        }
        if self.crs is not None:
            metadata["crs"] = self.crs
        return metadata

    def to_proto(self):
        raise NotImplementedError(
            "CalibrationGrid protobuf serialization is not implemented."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "CalibrationGrid protobuf deserialization is not implemented."
        )


def _bounds_for_models(models) -> Bounds:
    bounds = None
    for model in models:
        model_bounds = getattr(model, "bounds", None)
        if model_bounds is None:
            continue
        if bounds is None:
            bounds = model_bounds.copy()
        else:
            bounds.union(model_bounds)
    return bounds or Bounds()


def _polygon_coordinates(polygon) -> list[list[list[float]]]:
    coordinates = [
        [[float(x), float(y)] for x, y in polygon.exterior.coords],
    ]
    for interior in polygon.interiors:
        coordinates.append([[float(x), float(y)] for x, y in interior.coords])
    return coordinates
