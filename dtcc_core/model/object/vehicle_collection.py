# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List

import numpy as np

from .object import Object
from ..geometry import Point
from ...plotting.style import (
    DTCC_CATEGORY_PALETTE,
    DTCC_COLORS,
    add_categorical_legend,
    add_plot_context,
    apply_dtcc_style,
    resolve_colormap,
    style_colorbar,
)


@dataclass
class VehicleCollection(Object):
    """Represents a live snapshot of public transport vehicles.

    Each vehicle is represented as a child ``Object`` with point geometry for
    the latest known location and metadata such as route, line, mode, provider,
    speed, bearing, and timestamp stored in attributes.
    """

    def add_vehicle(self, vehicle: Object) -> None:
        """Add a vehicle as a child object."""
        if not isinstance(vehicle, Object):
            raise TypeError("Vehicle must be an Object instance")
        self.add_child(vehicle)

    def vehicles(self) -> List[Object]:
        """Return all vehicles in the collection."""
        result = []
        for child_list in self.children.values():
            result.extend(child_list)
        return result

    def to_arrays(self, field_name: str | None = None):
        """Convert vehicle locations and optional values to numpy arrays.

        If ``field_name`` is provided, values are read first from point fields
        and then from vehicle attributes. If omitted, the first point field is
        used when available, otherwise the ``mode`` attribute is used.
        """
        points = []
        values = []
        have_values = False

        for vehicle in self.vehicles():
            point = self._point_geometry(vehicle)
            if point is None:
                continue

            points.append([point.x, point.y, point.z])
            value = self._value_for_vehicle(vehicle, point, field_name)
            if value is not None:
                have_values = True
            values.append(value)

        if not points:
            return np.empty((0, 3)), np.empty(0)

        value_array = np.asarray(values, dtype=object)
        if have_values and all(_is_number(value) or value is None for value in values):
            value_array = np.asarray(
                [np.nan if value is None else float(value) for value in values],
                dtype=float,
            )
        return np.asarray(points, dtype=float), value_array

    def plot(
        self,
        ax=None,
        column: str | None = "mode",
        color=DTCC_COLORS["teal"],
        size: float = 42.0,
        cmap=None,
        legend: bool = True,
        equal_aspect: bool = True,
        show_direction: bool = True,
        title: str | None = None,
        metadata: bool | dict[str, Any] = True,
        metadata_loc: str = "upper left",
        theme: str = "dark",
        presentation: bool = True,
        show: bool = True,
        **kwargs,
    ):
        """Plot vehicle positions using matplotlib.

        Parameters
        ----------
        column : str, optional
            Attribute or point field used for coloring. Defaults to ``"mode"``.
        show_direction : bool, default True
            Draw small heading arrows when vehicle bearings are available.
        """
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
        points, values = self.to_arrays(column)

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
                title=None if presentation_enabled else title or "DTCC Transit Vehicles",
                metadata=(
                    None
                    if presentation_enabled
                    else self._plot_metadata(column, 0, metadata)
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

        if column is None or len(values) == 0 or all(value is None for value in values):
            ax.scatter(points[:, 0], points[:, 1], s=size, c=color, **kwargs)
        elif _is_numeric_array(values):
            mappable = ax.scatter(
                points[:, 0],
                points[:, 1],
                s=size,
                c=values.astype(float),
                cmap=resolve_colormap(cmap),
                **kwargs,
            )
            if legend and not presentation_enabled:
                cbar = ax.figure.colorbar(mappable, ax=ax, shrink=0.78)
                cbar.set_label(column)
                style_colorbar(cbar, theme=theme)
        else:
            categories = [str(value) if value is not None else "unknown" for value in values]
            colors = _categorical_colors(categories)
            ax.scatter(points[:, 0], points[:, 1], s=size, c=colors, **kwargs)
            if legend and not presentation_enabled:
                legend_items = {
                    category: _color_for_category(category)
                    for category in sorted(set(categories))
                }
                add_categorical_legend(
                    ax,
                    legend_items,
                    title=column or "Vehicles",
                    marker="patch",
                    theme=theme,
                )

        if show_direction:
            self._plot_bearings(ax, points)

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
            title=None if presentation_enabled else title or "DTCC Transit Vehicles",
            metadata=(
                None
                if presentation_enabled
                else self._plot_metadata(column, len(points), metadata)
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

    def _point_geometry(self, vehicle: Object):
        for geom in vehicle.geometry.values():
            if isinstance(geom, Point):
                return geom
            if hasattr(geom, "x") and hasattr(geom, "y"):
                return geom
        return None

    def _value_for_vehicle(self, vehicle: Object, point, field_name: str | None):
        if field_name is not None:
            for field in getattr(point, "fields", []):
                if field.name == field_name and len(field.values) > 0:
                    return np.ravel(field.values)[0]
            return vehicle.attributes.get(field_name)

        for field in getattr(point, "fields", []):
            if len(field.values) > 0:
                return np.ravel(field.values)[0]
        return vehicle.attributes.get("mode")

    def _plot_bearings(self, ax, points):
        bearings = []
        for vehicle in self.vehicles():
            value = vehicle.attributes.get("bearing")
            bearings.append(np.nan if value is None else float(value))

        if len(bearings) != len(points) or np.all(np.isnan(bearings)):
            return

        bearings_rad = np.deg2rad(np.asarray(bearings, dtype=float))
        valid = ~np.isnan(bearings_rad)
        if not np.any(valid):
            return

        span = max(np.ptp(points[:, 0]), np.ptp(points[:, 1]), 1.0)
        scale = span / 28.0
        u = np.sin(bearings_rad[valid]) * scale
        v = np.cos(bearings_rad[valid]) * scale
        ax.quiver(
            points[valid, 0],
            points[valid, 1],
            u,
            v,
            angles="xy",
            scale_units="xy",
            scale=1,
            width=0.003,
            color=DTCC_COLORS["yellow"],
            alpha=0.9,
            zorder=10,
        )

    def _plot_metadata(self, column, vehicle_count: int, metadata):
        if metadata is False:
            return None

        plot_metadata: dict[str, Any] = {
            "Vehicles": vehicle_count,
            "Column": column,
        }
        for key, label in (
            ("provider", "Provider"),
            ("source", "Source"),
            ("retrieval_time", "Retrieved"),
            ("crs", "CRS"),
        ):
            if key in self.attributes:
                plot_metadata[label] = self.attributes[key]

        modes = sorted(
            {
                str(vehicle.attributes.get("mode"))
                for vehicle in self.vehicles()
                if vehicle.attributes.get("mode")
            }
        )
        if modes:
            plot_metadata["Modes"] = ", ".join(modes)

        if isinstance(metadata, dict):
            plot_metadata.update(metadata)
        return plot_metadata

    def __str__(self):
        lines = []
        vehicles = self.vehicles()
        lines.append("=" * 70)
        lines.append("DTCC VehicleCollection")
        lines.append("=" * 70)
        lines.append(f"Number of vehicles: {len(vehicles)}")
        if self.bounds:
            lines.append(f"Bounds: {self.bounds}")

        if self.attributes:
            lines.append("")
            lines.append("Dataset Information:")
            for key in ("source", "provider", "crs", "retrieval_time"):
                if key in self.attributes:
                    lines.append(f"  {key}: {self.attributes[key]}")

            if "partial_result" in self.attributes:
                lines.append("")
                lines.append("Live Data Status:")
                lines.append(f"  partial_result: {self.attributes['partial_result']}")
                lines.append(
                    f"  upstream_error_count: "
                    f"{self.attributes.get('upstream_error_count', 0)}"
                )
                errors = self.attributes.get("upstream_errors") or []
                for i, error in enumerate(errors[:2]):
                    failure = error.get("failure_class", "upstream")
                    message = error.get("message", "")
                    lines.append(f"  error {i + 1} ({failure}): {message}")
                if len(errors) > 2:
                    lines.append(f"  ... and {len(errors) - 2} more upstream errors")

            help_lines = self.attributes.get("configuration_help") or []
            if help_lines:
                lines.append("")
                lines.append("Configuration Help:")
                for item in help_lines:
                    lines.append(f"  - {item}")

            if not vehicles and not self.attributes.get("partial_result", False):
                lines.append("")
                lines.append("No vehicles matched the selected bounds and modes.")

        if vehicles:
            lines.append("")
            lines.append("Sample Vehicles:")
            for i, vehicle in enumerate(vehicles[:3]):
                attrs = vehicle.attributes
                point = self._point_geometry(vehicle)
                loc = "N/A"
                if point is not None:
                    loc = f"({point.x:.2f}, {point.y:.2f})"
                label = attrs.get("line") or attrs.get("route_id") or attrs.get("vehicle_id")
                lines.append(f"  {i + 1}. {label or 'Vehicle'}")
                lines.append(f"     Mode: {attrs.get('mode', 'unknown')}")
                lines.append(f"     Location: {loc}")
                if attrs.get("timestamp"):
                    lines.append(f"     Timestamp: {attrs['timestamp']}")
            if len(vehicles) > 3:
                lines.append(f"  ... and {len(vehicles) - 3} more vehicles")

        lines.append("=" * 70)
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()


def _is_number(value) -> bool:
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def _is_numeric_array(values: np.ndarray) -> bool:
    if values.dtype.kind in "fiu":
        return True
    return all(_is_number(value) or value is None for value in values)


def _color_for_category(category: str) -> str:
    index = sum(ord(char) for char in category) % len(DTCC_CATEGORY_PALETTE)
    return DTCC_CATEGORY_PALETTE[index]


def _categorical_colors(categories: list[str]) -> list[str]:
    return [_color_for_category(category) for category in categories]
