"""Fetch one day of hourly temperature in Stockholm (requires internet)."""

import dtcc_core as dtcc

weather = dtcc.datasets.weather_history(
    source="open_meteo",
    location=(674572, 6580743),
    crs="EPSG:3006",
    start="2020-01-01",
    end="2020-01-02",  # Exclusive end; times are UTC.
    parameters=["air_temperature"],
    strict_live=True,
)

print(weather.to_dataframe()[["time", "value", "unit"]].to_string(index=False))
