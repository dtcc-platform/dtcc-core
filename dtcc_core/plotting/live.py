from __future__ import annotations

from collections import deque
from time import monotonic
from typing import Callable

import numpy as np

from .style import DTCC_COLORS, apply_dtcc_style, get_axes, require_matplotlib


class LiveVehiclePlot:
    """Animated Matplotlib plot for live vehicle collections.

    ``VehicleCollection.plot()`` is best for one snapshot. This helper owns the
    Matplotlib artists needed for repeated updates: current vehicle markers,
    per-vehicle trail lines, and optional static road-network context.
    """

    def __init__(
        self,
        *,
        road_network=None,
        bounds=None,
        ax=None,
        title: str = "Live Vehicles",
        window_title: str | None = None,
        trail_length: int = 80,
        stale_after_s: float = 120.0,
        vehicle_color: str = DTCC_COLORS["yellow"],
        trail_color: str = DTCC_COLORS["teal_light"],
        road_color: str = "#4f5b62",
        road_linewidth: float = 0.8,
        vehicle_size: float = 58.0,
        theme: str = "dark",
    ):
        self.title = title
        self.trail_length = int(trail_length)
        self.stale_after_s = float(stale_after_s)
        self.trails: dict[str, deque[tuple[float, float]]] = {}
        self.trail_lines = {}
        self.last_seen: dict[str, float] = {}

        if road_network is not None and ax is None:
            ax = road_network.plot(
                color=road_color,
                linewidth=road_linewidth,
                legend=False,
                metadata=False,
                title=title,
                show=False,
            )
        else:
            ax = get_axes(ax)
            apply_dtcc_style(
                ax,
                theme=theme,
                equal_aspect=True,
                xlabel="x",
                ylabel="y",
                grid=True,
                title=title,
            )

        self.ax = ax
        self.fig = ax.figure
        if window_title and hasattr(self.fig.canvas.manager, "set_window_title"):
            self.fig.canvas.manager.set_window_title(window_title)

        self.current = ax.scatter(
            [],
            [],
            s=vehicle_size,
            color=vehicle_color,
            edgecolors=DTCC_COLORS["dark"],
            linewidths=0.8,
            zorder=5,
            label="Current vehicle position",
        )
        self.trail_proxy = ax.plot(
            [],
            [],
            color=trail_color,
            linewidth=1.4,
            alpha=0.65,
            label="Vehicle trail",
        )[0]
        self.trail_proxy.set_data([], [])
        self.trail_color = trail_color
        self.ax.legend(loc="upper right")
        self._set_bounds(bounds or getattr(road_network, "bounds", None))

    def update(self, vehicles, *, now: float | None = None):
        """Update current markers and trails from a VehicleCollection."""
        now = monotonic() if now is None else float(now)
        positions = vehicle_positions(vehicles)

        for vehicle_id, xy in positions.items():
            self.trails.setdefault(
                vehicle_id, deque(maxlen=self.trail_length)
            ).append(xy)
            self.last_seen[vehicle_id] = now

        self._prune_stale(now)
        self._update_current_markers(positions)
        self._update_trails()

        retrieval_time = getattr(vehicles, "attributes", {}).get("retrieval_time", "")
        self.ax.set_title(
            f"{self.title}: {len(positions)} current, "
            f"{len(self.trails)} tracked\n{retrieval_time}"
        )
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return positions

    def run(
        self,
        fetch_vehicles: Callable[[], object],
        *,
        interval_s: float = 1.0,
    ):
        """Poll ``fetch_vehicles`` until the Matplotlib window is closed."""
        plt = require_matplotlib("live vehicle animation")
        plt.show(block=False)
        try:
            while plt.fignum_exists(self.fig.number):
                self.update(fetch_vehicles())
                plt.pause(interval_s)
        except KeyboardInterrupt:
            return self
        return self

    def _set_bounds(self, bounds):
        if bounds is None:
            return
        if all(hasattr(bounds, attr) for attr in ("xmin", "ymin", "xmax", "ymax")):
            self.ax.set_xlim(bounds.xmin, bounds.xmax)
            self.ax.set_ylim(bounds.ymin, bounds.ymax)
            return
        try:
            xmin, ymin, xmax, ymax = list(bounds)[:4]
        except (TypeError, ValueError):
            return
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)

    def _prune_stale(self, now: float):
        for vehicle_id in list(self.trails):
            if now - self.last_seen.get(vehicle_id, now) <= self.stale_after_s:
                continue
            self.trails.pop(vehicle_id, None)
            self.last_seen.pop(vehicle_id, None)
            line = self.trail_lines.pop(vehicle_id, None)
            if line is not None:
                line.remove()

    def _update_current_markers(self, positions: dict[str, tuple[float, float]]):
        offsets = np.asarray(list(positions.values()), dtype=float)
        if len(offsets) == 0:
            offsets = np.empty((0, 2), dtype=float)
        self.current.set_offsets(offsets)

    def _update_trails(self):
        for vehicle_id, trail in self.trails.items():
            line = self.trail_lines.get(vehicle_id)
            if line is None:
                line = self.ax.plot(
                    [],
                    [],
                    color=self.trail_color,
                    linewidth=1.4,
                    alpha=0.65,
                    zorder=4,
                )[0]
                self.trail_lines[vehicle_id] = line
            xy = np.asarray(trail, dtype=float)
            line.set_data(xy[:, 0], xy[:, 1])


def vehicle_positions(vehicles) -> dict[str, tuple[float, float]]:
    """Extract keyed 2D positions from a VehicleCollection-like object."""
    positions = {}
    for vehicle in vehicles.vehicles():
        point = next(
            (
                geom
                for geom in vehicle.geometry.values()
                if hasattr(geom, "x") and hasattr(geom, "y")
            ),
            None,
        )
        if point is None:
            continue
        vehicle_id = vehicle.attributes.get("vehicle_id") or vehicle.id
        positions[str(vehicle_id)] = (float(point.x), float(point.y))
    return positions


def animate_vehicle_collection(fetch_vehicles, **kwargs) -> LiveVehiclePlot:
    """Create a LiveVehiclePlot and run it immediately."""
    interval_s = kwargs.pop("interval_s", 1.0)
    plotter = LiveVehiclePlot(**kwargs)
    plotter.run(fetch_vehicles, interval_s=interval_s)
    return plotter


__all__ = ["LiveVehiclePlot", "animate_vehicle_collection", "vehicle_positions"]
