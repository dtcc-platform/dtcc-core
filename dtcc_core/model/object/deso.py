from dataclasses import dataclass
import builtins
from collections import Counter
from typing import Any

import numpy as np

from .object import Object, GeometryType
from ..geometry import MultiSurface


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
        for area in self.areas:
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

        return arrays

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


__all__ = ["DeSO"]
