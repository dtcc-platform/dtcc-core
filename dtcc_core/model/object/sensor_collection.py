# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass
from typing import Any, List
import numpy as np

from .object import Object
from ..geometry import Point
from ..values import Field
from ...plotting.style import (
    add_plot_context,
    apply_dtcc_style,
    resolve_colormap,
    style_colorbar,
)


@dataclass
class SensorCollection(Object):
    """Represents a collection of sensor stations/devices.

    Each sensor station is represented as a child Object with:
    - Point geometry for location
    - Field(s) attached to the geometry with measurement values
    - Metadata stored in attributes (station ID, timestamp, etc.)

    This class provides convenience methods for working with sensor data
    in the DTCC object model.
    """

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

    def stations(self) -> List[Object]:
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
            for geom_type, geom in station.geometry.items():
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
            for geom in station.geometry.values():
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

    def __str__(self):
        """Return a pretty-printed representation of the SensorCollection."""
        lines = []
        lines.append("=" * 70)
        lines.append("DTCC SensorCollection")
        lines.append("=" * 70)

        # Get basic info
        stations = self.stations()
        n_stations = len(stations)

        lines.append(f"Number of stations: {n_stations}")

        # Get bounds if available
        if self.bounds:
            lines.append(f"Bounds: {self.bounds}")

        # Get attributes if available
        if self.attributes:
            # Show key attributes
            attrs_to_show = [
                ("source", "Source"),
                ("phenomenon", "Phenomenon"),
                ("crs", "CRS"),
                ("retrieval_time", "Retrieved"),
            ]

            lines.append("")
            lines.append("Dataset Information:")
            for key, label in attrs_to_show:
                if key in self.attributes:
                    lines.append(f"  {label}: {self.attributes[key]}")

        # Show statistics if we have numeric values
        if n_stations > 0:
            try:
                import numpy as np

                # Try to get values for statistics
                # Discover the field name: check 'phenomenon' attribute first,
                # then try 'parameter_fields' keys, then fall back to the first
                # field on the first station's geometry.
                phenomenon = self.attributes.get("phenomenon")
                if phenomenon is None:
                    pf = self.attributes.get("parameter_fields")
                    if pf and isinstance(pf, dict):
                        phenomenon = next(iter(pf))
                if phenomenon is None:
                    for st in stations:
                        for _k, geom in st.geometry.items():
                            if hasattr(geom, "fields") and geom.fields:
                                phenomenon = geom.fields[0].name
                                break
                        if phenomenon:
                            break
                if phenomenon is None:
                    phenomenon = "value"
                points, values = self.to_arrays(phenomenon)

                if len(values) > 0 and not np.all(np.isnan(values)):
                    lines.append("")
                    lines.append("Measurement Statistics:")
                    valid_values = values[~np.isnan(values)]
                    if len(valid_values) > 0:
                        lines.append(f"  Count: {len(valid_values)}")
                        lines.append(f"  Min: {np.min(valid_values):.2f}")
                        lines.append(f"  Max: {np.max(valid_values):.2f}")
                        lines.append(f"  Mean: {np.mean(valid_values):.2f}")
                        lines.append(f"  Median: {np.median(valid_values):.2f}")
            except Exception:
                pass  # Skip statistics if there's any issue

        # Get sample station info
        if n_stations > 0:
            lines.append("")
            lines.append("Sample Stations:")
            for i, station in enumerate(stations[:3]):  # Show first 3
                attrs = station.attributes

                # Get location
                loc_str = "N/A"
                for geom_type, geom in station.geometry.items():
                    if hasattr(geom, "x") and hasattr(geom, "y"):
                        loc_str = f"({geom.x:.4f}, {geom.y:.4f})"
                        break

                # Get value – try attrs['value'] first, then read from
                # the first Field on the Point geometry.
                value_str = attrs.get("value", None)
                unit = attrs.get("unit", "")
                if value_str is None:
                    for _gk, geom in station.geometry.items():
                        if hasattr(geom, "fields") and geom.fields:
                            f = geom.fields[0]
                            v = f.values[0] if len(f.values) > 0 else None
                            if v is not None and not np.isnan(v):
                                value_str = v
                                unit = getattr(f, "unit", unit)
                            break
                if value_str is None:
                    value_str = "N/A"
                if isinstance(value_str, (int, float)):
                    value_str = f"{value_str:.2f}"
                timestamp = attrs.get("timestamp", "N/A")

                station_name = attrs.get("station_name", f"Station {i+1}")
                lines.append(f"  {i+1}. {station_name}")
                lines.append(f"     Location: {loc_str}")
                lines.append(f"     Value: {value_str} {unit}")
                if timestamp and timestamp != "N/A":
                    # Shorten timestamp if it's ISO format
                    if "T" in timestamp:
                        timestamp = timestamp.split("T")[0]
                    lines.append(f"     Timestamp: {timestamp}")

            if n_stations > 3:
                lines.append(f"  ... and {n_stations - 3} more stations")

        lines.append("=" * 70)

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()


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
