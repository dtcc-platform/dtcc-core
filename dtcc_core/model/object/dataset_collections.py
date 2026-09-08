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
from ...plotting.style import (
    add_categorical_legend,
    add_plot_context,
    apply_dtcc_style,
    plot_line_segments,
    plot_polygon_geometries,
    resolve_colormap,
    style_colorbar,
)


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
        """Build a footprint collection from buildings with available geometry.

        By default, this extracts canonical LOD0 building footprints only.
        Passing ``geom_type`` requests an advanced/derived extraction from that
        explicit geometry type.
        """
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

    def __str__(self) -> str:
        return f"DTCC FootprintCollection with {len(self)} footprint(s)"

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

    def plot(
        self,
        ax=None,
        edgecolor: str | None = None,
        linewidth: float = 0.8,
        alpha: float = 0.45,
        title: str | None = None,
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
        **kwargs,
    ):
        """Plot building footprint polygons with matplotlib."""
        from dtcc_core.datasets.presentation import (
            draw_empty_preview_state,
            finalize_presentation_plot,
            presentation_plot_axes,
        )

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        ax, panel_ax = presentation_plot_axes(
            context,
            ax=ax,
            presentation=presentation,
        )
        geometries = self.to_shapely()
        if presentation_enabled and len(geometries) == 0:
            draw_empty_preview_state(ax, context)
        plot_polygon_geometries(
            ax,
            geometries,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=alpha,
            **kwargs,
        )
        ax.autoscale()
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=True,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x",
            ylabel=None if presentation_enabled else "y",
            grid=False if presentation_enabled else True,
        )
        add_plot_context(
            ax,
            title=None if presentation_enabled else title or "DTCC Building Footprints",
            metadata=None if presentation_enabled else {"Footprints": len(self)},
            bounds=None if presentation_enabled else self.bounds,
            theme=theme,
        )
        return finalize_presentation_plot(
            ax,
            context=context,
            obj=self,
            panel_ax=panel_ax,
            presentation=presentation,
            show=show,
            theme=theme,
            title=title,
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

    def __str__(self) -> str:
        return f"DTCC BuildingCollection with {len(self)} building(s)"

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
        """Return building footprints as a semantic collection.

        By default, this extracts canonical LOD0 footprints only. Passing
        ``geom_type`` requests an advanced/derived extraction from that explicit
        geometry type.
        """
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

    def __str__(self) -> str:
        return f"DTCC TreeCollection with {len(self)} tree(s)"

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

    @property
    def bounds(self) -> Bounds:
        """Return bounds spanning all tree positions."""
        return _bounds_for_points(self.to_arrays())

    def to_proto(self):
        raise NotImplementedError(
            "TreeCollection protobuf serialization is not implemented."
        )

    def from_proto(self, pb):
        raise NotImplementedError(
            "TreeCollection protobuf deserialization is not implemented."
        )

    def plot(
        self,
        ax=None,
        column: str = "height",
        size: float = 42.0,
        cmap=None,
        title: str | None = None,
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
        **kwargs,
    ):
        """Plot tree positions with matplotlib."""
        from dtcc_core.datasets.presentation import (
            draw_empty_preview_state,
            finalize_presentation_plot,
            presentation_plot_axes,
        )

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        ax, panel_ax = presentation_plot_axes(
            context,
            ax=ax,
            presentation=presentation,
        )
        points = self.to_arrays()
        values = np.asarray([float(getattr(tree, column, np.nan)) for tree in self])

        if len(points) == 0:
            if presentation_enabled:
                draw_empty_preview_state(ax, context)
            apply_dtcc_style(
                ax,
                theme=theme,
                equal_aspect=True,
                axis="off" if presentation_enabled else "on",
                xlabel=None if presentation_enabled else "x",
                ylabel=None if presentation_enabled else "y",
            )
            add_plot_context(
                ax,
                title=None if presentation_enabled else title or "DTCC Trees",
                metadata=None if presentation_enabled else {"Trees": 0, "Column": column},
                bounds=None if presentation_enabled else self.bounds,
                theme=theme,
            )
            return finalize_presentation_plot(
                ax,
                context=context,
                obj=self,
                panel_ax=panel_ax,
                presentation=presentation,
                show=show,
                theme=theme,
                title=title,
            )

        if np.all(np.isnan(values)):
            ax.scatter(points[:, 0], points[:, 1], s=size, **kwargs)
        else:
            scatter = ax.scatter(
                points[:, 0],
                points[:, 1],
                s=size,
                c=values,
                cmap=resolve_colormap(cmap),
                **kwargs,
            )
            if not presentation_enabled:
                colorbar = ax.figure.colorbar(scatter, ax=ax, shrink=0.78)
                colorbar.set_label(column)
                style_colorbar(colorbar, theme=theme)

        ax.autoscale()
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=True,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x",
            ylabel=None if presentation_enabled else "y",
            grid=False if presentation_enabled else True,
        )
        add_plot_context(
            ax,
            title=None if presentation_enabled else title or "DTCC Trees",
            metadata=None if presentation_enabled else {"Trees": len(self), "Column": column},
            bounds=None if presentation_enabled else self.bounds,
            theme=theme,
        )
        return finalize_presentation_plot(
            ax,
            context=context,
            obj=self,
            panel_ax=panel_ax,
            presentation=presentation,
            show=show,
            theme=theme,
            title=title,
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

    def __str__(self) -> str:
        return (
            f"DTCC CalibrationGrid with {len(self.features)} line feature(s), "
            f"{self.divisions} division(s), {self.crs or 'unknown CRS'}"
        )

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

    def plot(
        self,
        ax=None,
        color: str | None = None,
        linewidth: float = 0.8,
        title: str | None = None,
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
        **kwargs,
    ):
        """Plot calibration grid line features with matplotlib."""
        from dtcc_core.datasets.presentation import (
            draw_empty_preview_state,
            finalize_presentation_plot,
            presentation_plot_axes,
        )

        context = self.dataset_context
        presentation_enabled = presentation and context is not None
        ax, panel_ax = presentation_plot_axes(
            context,
            ax=ax,
            presentation=presentation,
        )
        segments = []
        for feature in self.features:
            geometry = feature.get("geometry") or {}
            if geometry.get("type") != "LineString":
                continue
            coordinates = geometry.get("coordinates") or []
            if len(coordinates) >= 2:
                segments.append([(float(x), float(y)) for x, y, *_ in coordinates])

        effective_color = color or "#4cc9f0"
        if presentation_enabled and len(segments) == 0:
            draw_empty_preview_state(ax, context)
        plot_line_segments(
            ax,
            segments,
            color=effective_color,
            linewidth=linewidth,
            theme=theme,
            **kwargs,
        )
        if segments and not presentation_enabled:
            add_categorical_legend(
                ax,
                {"Calibration grid": effective_color},
                title="Layers",
                linewidth=linewidth,
                theme=theme,
            )
        ax.autoscale()
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=True,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x",
            ylabel=None if presentation_enabled else "y",
            grid=False if presentation_enabled else True,
        )
        add_plot_context(
            ax,
            title=None if presentation_enabled else title or "DTCC Calibration Grid",
            metadata=(
                None
                if presentation_enabled
                else {"Lines": len(self.features), "Divisions": self.divisions}
            ),
            bounds=None if presentation_enabled else self.bounds,
            theme=theme,
        )
        return finalize_presentation_plot(
            ax,
            context=context,
            obj=self,
            panel_ax=panel_ax,
            presentation=presentation,
            show=show,
            theme=theme,
            title=title,
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


def _bounds_for_points(points: np.ndarray) -> Bounds:
    if len(points) == 0:
        return Bounds()
    return Bounds(
        xmin=float(np.min(points[:, 0])),
        ymin=float(np.min(points[:, 1])),
        xmax=float(np.max(points[:, 0])),
        ymax=float(np.max(points[:, 1])),
    )


def _polygon_coordinates(polygon) -> list[list[list[float]]]:
    coordinates = [
        [[float(x), float(y)] for x, y in polygon.exterior.coords],
    ]
    for interior in polygon.interiors:
        coordinates.append([[float(x), float(y)] for x, y in interior.coords])
    return coordinates
