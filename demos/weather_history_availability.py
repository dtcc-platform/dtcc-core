"""Discover advertised weather periods near Stockholm; no readings are downloaded."""

import dtcc_core as dtcc

availability = dtcc.datasets.weather_history.availability(
    source="smhi",
    location=(674572, 6580743),
    crs="EPSG:3006",
    max_distance_m=5_000,
    strict_live=True,
)

for parameter, details in availability["parameters"].items():
    print(parameter, details["unit"], details["status"])
    for station in details["candidates"]:
        print(" ", station["station_name"], round(station["distance_m"]), "m")
        for period in station["periods"]:
            print("   ", period["start"], period["end_exclusive"], "(advertised)")

print("Observations verified:", availability["observations_verified"])
