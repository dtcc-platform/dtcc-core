"""Compare hourly SMHI observations and ERA5 in Stockholm (requires internet)."""

import matplotlib.pyplot as plt

import dtcc_core as dtcc

request = dict(
    location=(674572, 6580743),  # SWEREF99 TM; both APIs receive WGS84.
    start="2020-01-01",
    end="2020-01-04",  # Exclusive end; all times are UTC.
    parameters=["air_temperature", "wind_speed", "wind_direction"],
    strict_live=True,
)
print("Fetching SMHI station observations (station discovery can take a few minutes).")
smhi = dtcc.datasets.weather_history(source="smhi", max_distance_m=50_000, **request)
print("Fetching Open-Meteo ERA5 reanalysis.")
era5 = dtcc.datasets.weather_history(source="open_meteo", **request)
station = smhi.attributes["selection"]
station_id = station["station_id"]
label = f"SMHI station {station_id} (observed)"
print(f"{label}: {station['distance_m'] / 1000:.1f} km from the requested location")
fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True, layout="constrained")
for ax, parameter in zip(axes, request["parameters"], strict=True):
    # Direction is circular: points avoid artificial lines across north (0/360).
    style = dict(marker=".", linestyle="none") if parameter == "wind_direction" else {}
    smhi[(f"smhi:{station_id}", parameter)].plot(
        ax=ax, show=False, label=label, color="tab:blue", **style
    )
    era5[("open_meteo:location", parameter)].plot(
        ax=ax, show=False, label="Open-Meteo ERA5", color="tab:orange", **style
    )
    ax.set_title(parameter.replace("_", " ").capitalize())
axes[-1].set_ylim(0, 360)
axes[-1].set_yticks([0, 90, 180, 270, 360])
fig.suptitle(
    f"Historical weather: {request['start']} to {request['end']} (end excluded)\n"
    f"SMHI station distance: {station['distance_m'] / 1000:.1f} km",
    color=axes[0].title.get_color(),
)
plt.show()
