#!/usr/bin/env python3
"""
Air Quality Snapshot Example
=============================

This example demonstrates how to fetch and work with air quality
sensor data using the DTCC air quality dataset.

The dataset fetches snapshot measurements from datavardluft.smhi.se
and returns a SensorCollection with station locations and values.
"""

import sys
from pathlib import Path

from dtcc_core import datasets

# Define bounds for Stockholm area (SWEREF99 TM / EPSG:3006)
stockholm_bounds = (674000, 6580000, 676000, 6582000)

print("Fetching NO2 measurements from Stockholm area...")
print(f"Bounds: {stockholm_bounds}")
print()

try:
    # Fetch air quality data (default CRS is EPSG:3006)
    sensors = datasets.air_quality(bounds=stockholm_bounds, phenomenon="NO2")

    print(f"✓ Fetched sensor collection")
    print(f"  Source: {sensors.attributes.get('source')}")
    print(f"  Phenomenon: {sensors.attributes.get('phenomenon')}")
    print(f"  Stations found: {sensors.attributes.get('total_stations_found')}")
    print(f"  Stations used: {sensors.attributes.get('stations_used')}")
    print()

    # Get data as arrays
    if len(sensors.stations()) > 0:
        points, values = sensors.to_arrays("NO2")

        print(f"Measurement statistics:")
        print(f"  Count: {len(values)}")
        print(f"  Min: {values.min():.2f}")
        print(f"  Max: {values.max():.2f}")
        print(f"  Mean: {values.mean():.2f}")
        print()

        # Show first few stations
        print("Sample stations:")
        for i, station in enumerate(sensors.stations()[:3]):
            attrs = station.attributes
            point = station.geometry.get("location")
            print(f"  {i+1}. {attrs.get('station_name', 'Unknown')}")
            print(f"     Location: ({point.x:.4f}, {point.y:.4f})")
            print(f"     Value: {attrs.get('value', 'N/A')} {attrs.get('unit', '')}")
            print(f"     Timestamp: {attrs.get('timestamp', 'N/A')}")
    else:
        print("No stations found in the specified area.")

    print()
    print("✓ Example completed successfully")

except Exception as e:
    print(f"Error: {e}")
    import traceback

    traceback.print_exc()
    sys.exit(1)
