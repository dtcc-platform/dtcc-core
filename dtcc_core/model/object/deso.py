from dataclasses import dataclass
import builtins
from collections import Counter
from typing import Any

import numpy as np

from .object import Object, GeometryType
from ..geometry import MultiSurface
from ..values import Field


@dataclass
class DeSO(Object):
    """Collection of Swedish DeSO demographic statistical areas."""

    @property
    def areas(self) -> list[Object]:
        """Return DeSO area objects."""
        return self.get_children(Object)

    @property
    def codes(self) -> list[str]:
        """Return DeSO area codes."""
        return [area.attributes.get("desokod", area.id) for area in self.areas]

    @property
    def fields(self) -> dict[str, Field]:
        """Return area-aligned fields aggregated from DeSO area geometries."""
        fields: dict[str, Field] = {}
        metadata: dict[str, Field] = {}
        for area in self.areas:
            geometry = area.geometry.get(GeometryType.LOD0)
            if geometry is None:
                continue
            for field in geometry.fields:
                if field.name not in fields:
                    metadata[field.name] = field
                    fields[field.name] = []

        for name, rows in fields.items():
            dim = metadata[name].dim
            for area in self.areas:
                value = self._area_field_value(area, name, dim)
                rows.append(value)

        return {
            name: Field(
                name=name,
                unit=metadata[name].unit,
                description=metadata[name].description,
                values=np.asarray(values, dtype=float).reshape(
                    (-1, metadata[name].dim)
                ),
                dim=metadata[name].dim,
            )
            for name, values in fields.items()
        }

    @property
    def field_names(self) -> list[str]:
        """Return names of area-aligned fields attached to the DeSO areas."""
        return list(self.fields.keys())

    def __len__(self):
        return len(self.areas)

    def __str__(self):
        year = self.attributes.get("year")
        year_text = f" {year}" if year else ""
        return f"DTCC DeSO{year_text} with {len(self)} area(s)"

    def info(self, print: bool = True) -> str | None:
        """Print or return a human-readable multi-line summary."""
        lines = []
        lines.append("=" * 70)
        lines.append("DTCC DeSO")
        lines.append("=" * 70)
        lines.append(f"Areas: {len(self)}")

        year = self.attributes.get("year")
        if year is not None:
            lines.append(f"Year: {year}")

        source = self.attributes.get("source")
        if source:
            lines.append(f"Source: {source}")

        crs = self.transform.srs
        if crs:
            lines.append(f"CRS: {crs}")

        if self.bounds is not None:
            lines.append(f"Bounds: {self.bounds}")

        if len(self) > 0:
            fields = self.fields
            if fields:
                lines.append("")
                lines.append("Fields:")
                for field in fields.values():
                    lines.append(f"  {field.name} ({field.unit})")

            area_types = Counter(
                code[4]
                for code in self.codes
                if isinstance(code, str) and len(code) > 4
            )
            if area_types:
                lines.append("")
                lines.append("Area Types:")
                for area_type, count in sorted(area_types.items()):
                    lines.append(f"  {area_type}: {count}")

            attribute_keys = sorted(
                {key for area in self.areas for key in area.attributes.keys()}
            )
            if attribute_keys:
                lines.append("")
                lines.append("Attributes:")
                for key in attribute_keys:
                    lines.append(f"  {key}")

        lines.append("=" * 70)
        summary = "\n".join(lines)
        if print:
            builtins.print(summary)
            return None
        return summary

    def to_dataframe(self):
        """Convert DeSO areas to a GeoPandas GeoDataFrame."""
        try:
            import geopandas as gpd
            from shapely.geometry import MultiPolygon
        except ImportError as exc:
            raise ImportError("GeoPandas is required for DeSO.to_dataframe().") from exc

        rows: list[dict[str, Any]] = []
        geometries = []
        fields = self.fields
        for index, area in enumerate(self.areas):
            geometry = area.geometry.get(GeometryType.LOD0)
            if not isinstance(geometry, MultiSurface):
                continue

            polygons = [
                surface.to_polygon(simplify=0.0)
                for surface in geometry.surfaces
                if len(surface.vertices) >= 3
            ]
            polygons = [polygon for polygon in polygons if not polygon.is_empty]
            if len(polygons) == 0:
                continue

            row = dict(area.attributes)
            row["id"] = area.id
            for name, field in fields.items():
                row.setdefault(name, self._to_python_value(field.values[index]))
            rows.append(row)
            geometry = polygons[0] if len(polygons) == 1 else MultiPolygon(polygons)
            geometries.append(geometry)

        return gpd.GeoDataFrame(
            rows,
            geometry=geometries,
            crs=self.transform.srs or None,
        )

    def to_df(self):
        """Alias for :meth:`to_dataframe`."""
        return self.to_dataframe()

    def to_arrays(self, include_attributes: bool = True) -> dict:
        """Return low-level arrays with area centroids and optional attributes."""
        centroids = []
        codes = []
        for area in self.areas:
            geometry = area.geometry.get(GeometryType.LOD0)
            if not isinstance(geometry, MultiSurface) or len(geometry.surfaces) == 0:
                continue
            centroid = np.asarray(geometry.centroid(), dtype=float)
            centroids.append(centroid)
            codes.append(area.attributes.get("desokod", area.id))

        arrays = {
            "centroids": np.asarray(centroids, dtype=float),
            "codes": np.asarray(codes, dtype=object),
        }

        if include_attributes:
            attribute_keys = sorted(
                {key for area in self.areas for key in area.attributes.keys()}
            )
            attributes: dict[str, list[Any]] = {key: [] for key in attribute_keys}
            for area in self.areas:
                for key in attribute_keys:
                    attributes[key].append(area.attributes.get(key))
            arrays["attributes"] = {
                key: np.asarray(value, dtype=object)
                for key, value in attributes.items()
            }

        fields = self.fields
        if fields:
            arrays["fields"] = {
                name: field.values.copy()
                for name, field in fields.items()
            }

        return arrays

    def attach_field(self, field: Field, attribute_name: str | None = None):
        """Attach an area-aligned field to the DeSO areas.

        The field must have one row per DeSO area. Each row is stored on the
        corresponding area's LOD0 geometry and mirrored to area attributes so
        tabular export and viewer picking can expose the values directly.
        """
        if not isinstance(field, Field):
            raise TypeError("field must be a dtcc_core.model.Field.")
        if len(field.values) != len(self):
            raise ValueError(
                f"Field '{field.name}' has {len(field.values)} values, "
                f"but DeSO has {len(self)} areas."
            )

        values = np.asarray(field.values, dtype=float).reshape((-1, field.dim))
        name = attribute_name or field.name
        for area, value in zip(self.areas, values):
            geometry = area.geometry.get(GeometryType.LOD0)
            if geometry is None:
                continue

            area_field = Field(
                name=field.name,
                unit=field.unit,
                description=field.description,
                values=np.asarray(value, dtype=float).reshape((1, field.dim)),
                dim=field.dim,
            )
            geometry.fields = [f for f in geometry.fields if f.name != field.name]
            geometry.fields.append(area_field)
            area.attributes[name] = self._to_python_value(value)

        return self

    def get_field(self, name: str) -> Field | None:
        """Return an area-aligned field by name, if present."""
        return self.fields.get(name)

    def plot(
        self,
        ax=None,
        column: str | None = "desokod",
        edgecolor: str = "black",
        linewidth: float = 0.6,
        facecolor: str = "none",
        legend: bool = False,
        show: bool = True,
        **kwargs,
    ):
        """Plot DeSO polygons using GeoPandas/Matplotlib."""
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("Matplotlib is required for DeSO.plot().") from exc

        gdf = self.to_dataframe()
        if ax is None:
            _, ax = plt.subplots()

        plot_kwargs = dict(edgecolor=edgecolor, linewidth=linewidth, **kwargs)
        if column is not None and column in gdf.columns:
            gdf.plot(column=column, ax=ax, legend=legend, **plot_kwargs)
        else:
            gdf.plot(ax=ax, facecolor=facecolor, **plot_kwargs)

        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()

        if show:
            plt.show()
        return ax

    @staticmethod
    def _area_field_value(area: Object, name: str, dim: int):
        geometry = area.geometry.get(GeometryType.LOD0)
        if geometry is not None:
            for field in geometry.fields:
                if field.name == name and len(field.values) > 0:
                    value = np.asarray(field.values, dtype=float).reshape((-1, dim))[0]
                    return value if dim != 1 else value[0]
        if dim == 1:
            return np.nan
        return np.full(dim, np.nan)

    @staticmethod
    def _to_python_value(value):
        value = np.asarray(value)
        if value.size == 1:
            scalar = value.reshape(-1)[0]
            if np.isnan(scalar):
                return None
            return scalar.item() if hasattr(scalar, "item") else scalar
        return [
            None
            if np.isnan(item)
            else (item.item() if hasattr(item, "item") else item)
            for item in value.reshape(-1)
        ]


__all__ = ["DeSO"]
