# This demo animates live local buses on top of the street network in Gothenburg.

import dtcc_core as dtcc

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Domain size
L = 2000.0

# Define bounds
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

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
    lambda: dtcc.datasets.buses(bounds=bounds),
    interval_s=1.0,
)
