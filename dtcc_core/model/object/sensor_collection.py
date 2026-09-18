# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass
from typing import Any, List

import numpy as np

from ...plotting.style import (
    add_plot_context,
    apply_dtcc_style,
    resolve_colormap,
    style_colorbar,
)
from ..geometry import Point
from ..values import Field
from .object import Object


@dataclass(repr=False)
class SensorCollection(Object):
    """Represents a collection of sensor stations/devices.

    Each sensor station is represented as a child Object with:
    - Point geometry for location
    - Field(s) attached to the geometry with measurement values
    - Metadata stored in attributes (station ID, timestamp, etc.)

    This class provides convenience methods for working with sensor data
    in the DTCC object model.
    """

    def _info_sections(self):
        from ...common._display import value_text

        sections = super()._info_sections()
        stations = self.stations()
        phenomenon = self.attributes.get("phenomenon")
        if phenomenon is None:
            phenomenon = next(iter(self.attributes.get("parameter_fields") or {}), None)
        if phenomenon is None:
            phenomenon = next(
                (
                    f.name
                    for station in stations
                    for geometry in station.get_geometries()
                    for f in geometry.fields
                ),
                None,
            )
        if phenomenon is not None:
            _, values = self.to_arrays(phenomenon)
            if values.dtype.kind in "biuf":
                valid = values[np.isfinite(values)]
                if valid.size:
                    sections.append(
                        (
                            f"Measurement statistics: {phenomenon}",
                            ("Statistic", "Value"),
                            [
                                ("Count", valid.size),
                                ("Min", f"{valid.min():.2f}"),
                                ("Max", f"{valid.max():.2f}"),
                                ("Mean", f"{valid.mean():.2f}"),
                                ("Median", f"{np.median(valid):.2f}"),
                            ],
                        )
                    )
        rows = []
        for station in stations[:3]:
            attrs = station.attributes
            point = next(
                (g for g in station.get_geometries() if isinstance(g, Point)), None
            )
            value, unit = attrs.get("value"), attrs.get("unit", "")
            if value is None and point is not None and point.fields:
                field = next(
                    (f for f in point.fields if f.name == phenomenon), point.fields[0]
                )
                value = field.values[0] if len(field.values) else None
                unit = field.unit
            if isinstance(value, np.ndarray):
                value = np.array2string(value, threshold=6, edgeitems=2)

            rows.append(
                (
                    attrs.get("station_name", station.id),
                    "N/A" if point is None else f"({point.x:.4f}, {point.y:.4f})",
                    value_text(value),
                    unit,
                    attrs.get("timestamp", "N/A"),
                )
            )
        if rows:
            sections.append(
                (
                    "Sample stations (first 3)",
                    ("Station", "Location", "Value", "Unit", "Timestamp"),
                    rows,
                )
            )
        return sections

    def _summary_items(self):
        return [("num_stations", sum(len(group) for group in self.children.values()))]

    def add_station(self, station: Object) -> None:
        """Add a sensor station as a child object.

        Parameters
        ----------
        station : Object
            The sensor station object to add. Should have a Point geometry
            and Field(s) with measurement data.
        """
        if not isinstance(station, Object):
            raise TypeError("Station must be an Object instance")
        self.add_child(station)

    def stations(self) -> list[Object]:
        """Get all sensor stations (child objects).

        Returns
        -------
        List[Object]
            List of all sensor station objects.
        """
        # self.children is a dict[type, list], flatten to single list
        result = []
        for child_list in self.children.values():
            result.extend(child_list)
        return result

    def to_arrays(self, field_name: str = None):
        """Convert sensor data to numpy arrays.

        Extracts locations and values for all stations. If field_name is specified,
        only that field is extracted. Otherwise, the first field is used.

        Parameters
        ----------
        field_name : str, optional
            Name of the field to extract. If None, uses the first field.

        Returns
        -------
        points : np.ndarray
            Nx3 array of station coordinates
        values : np.ndarray
            N array of measurement values
        """
        points = []
        values = []

        for station in self.stations():
            # Find Point geometry
            point_geom = None
            for geom in station.get_geometries():
                if isinstance(geom, Point):
                    point_geom = geom
                    break

            if point_geom is None:
                continue

            # Find field
            field = None
            if field_name:
                # Search for field by name
                for f in point_geom.fields:
                    if f.name == field_name:
                        field = f
                        break
            else:
                # Use first field
                if point_geom.fields:
                    field = point_geom.fields[0]

            if field is None or len(field.values) == 0:
                continue

            # Add to arrays
            points.append([point_geom.x, point_geom.y, point_geom.z])
            values.append(field.values[0] if field.dim == 1 else field.values)

        if not points:
            return np.empty((0, 3)), np.empty(0)

        return np.array(points), np.array(values)

    def plot(
        self,
        field_name: str | None = None,
        ax=None,
        size: float = 44.0,
        cmap=None,
        equal_aspect: bool = True,
        title: str | None = None,
        metadata: bool | dict[str, Any] = True,
        metadata_loc: str = "upper left",
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
        **kwargs,
    ):
        """Plot sensor station points with matplotlib."""
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
        points, values = self.to_arrays(field_name)
        label = field_name or self._first_field_name() or "value"

        if len(points) == 0:
            if presentation_enabled:
                draw_empty_preview_state(ax, context)
            apply_dtcc_style(
                ax,
                theme=theme,
                equal_aspect=equal_aspect,
                axis="off" if presentation_enabled else "on",
                xlabel=None if presentation_enabled else "x",
                ylabel=None if presentation_enabled else "y",
                grid=False if presentation_enabled else True,
            )
            add_plot_context(
                ax,
                title=None if presentation_enabled else title or "DTCC Sensors",
                metadata=(
                    None
                    if presentation_enabled
                    else self._plot_metadata(label, 0, metadata)
                ),
                bounds=None if presentation_enabled else self.bounds,
                theme=theme,
                metadata_loc=metadata_loc,
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

        numeric_values = _numeric_values(values)
        if numeric_values is None:
            scatter = ax.scatter(points[:, 0], points[:, 1], s=size, **kwargs)
        else:
            scatter = ax.scatter(
                points[:, 0],
                points[:, 1],
                s=size,
                c=numeric_values,
                cmap=resolve_colormap(cmap),
                **kwargs,
            )
            if not presentation_enabled:
                colorbar = ax.figure.colorbar(scatter, ax=ax, shrink=0.78)
                colorbar.set_label(label)
                style_colorbar(colorbar, theme=theme)

        ax.autoscale()
        apply_dtcc_style(
            ax,
            theme=theme,
            equal_aspect=equal_aspect,
            axis="off" if presentation_enabled else "on",
            xlabel=None if presentation_enabled else "x",
            ylabel=None if presentation_enabled else "y",
            grid=False if presentation_enabled else True,
        )
        add_plot_context(
            ax,
            title=None if presentation_enabled else title or "DTCC Sensors",
            metadata=(
                None
                if presentation_enabled
                else self._plot_metadata(label, len(points), metadata)
            ),
            bounds=None if presentation_enabled else self.bounds,
            theme=theme,
            metadata_loc=metadata_loc,
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

    def _first_field_name(self) -> str | None:
        for station in self.stations():
            for geom in station.get_geometries():
                fields = getattr(geom, "fields", ())
                if fields:
                    return fields[0].name
        return None

    def _plot_metadata(self, field_name: str, station_count: int, metadata):
        if metadata is False:
            return None
        values = {"Stations": station_count, "Field": field_name}
        if isinstance(metadata, dict):
            values.update(metadata)
        return values


def _numeric_values(values):
    try:
        numeric = np.asarray(values, dtype=float)
    except (TypeError, ValueError):
        return None
    if numeric.ndim != 1:
        return None
    if len(numeric) == 0 or np.all(np.isnan(numeric)):
        return None
    return numeric
