"""Animate live buses over Gothenburg streets around Chalmers Johanneberg."""

import dtcc_core as dtcc

# Chalmers Johanneberg centre in EPSG:3006
center_x = 319970.0
center_y = 6397910.0

# Domain size
L = 2000.0

# Define bounds
bounds = dtcc.Bounds(
    center_x - 0.5 * L, center_y - 0.5 * L, center_x + 0.5 * L, center_y + 0.5 * L
)

# Fetch the road network once
roads = dtcc.datasets.roads(bounds=bounds)

# Animate live buses indefinitely
print("Close the matplotlib window or press Ctrl+C to stop.")
dtcc.plotting.LiveVehiclePlot(
    road_network=roads,
    bounds=bounds,
    title="Live buses",
    window_title="DTCC Live Buses",
    trail_length=240,
    stale_after_s=180.0,
).run(
    lambda: dtcc.datasets.buses(bounds=bounds, strict_live=True),
    interval_s=1.0,
)
