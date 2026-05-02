#!/usr/bin/env python3
"""Smoke-test prototype LOD2 reconstruction using live DTCC downloads."""

import os
from pathlib import Path

# The downloader reads these when dtcc_core.io is imported, so set defaults first.
os.environ.setdefault("DTCC_LIDAR_URL", "http://13.60.69.202:8001")
os.environ.setdefault("DTCC_GPKG_URL", "http://13.60.69.202:8001")

from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from dtcc_core.io import download_footprints, download_pointcloud
from dtcc_core.model import Bounds, City


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    output_path = repo_root / "sandbox" / "lod2_live.city.json"
    bounds = Bounds(xmin=386075.0, ymin=6174447.0, xmax=386575.0, ymax=6174947.0)

    print(f"DTCC_LIDAR_URL: {os.environ['DTCC_LIDAR_URL']}")
    print(f"DTCC_GPKG_URL: {os.environ['DTCC_GPKG_URL']}")
    print(f"bounds: {bounds.tuple}")

    city = City()
    city.bounds = bounds
    city.add_buildings(download_footprints(bounds))
    city.add_pointcloud(download_pointcloud(bounds))
    city.build_lod2_buildings(calculate_heights=True)

    lod2_buildings = [building for building in city.buildings if building.lod2 is not None]
    lod1_buildings = [building for building in city.buildings if building.lod1 is not None]
    watertight_count = sum(is_watertight(building.lod2) for building in lod2_buildings)

    print(f"buildings: {len(city.buildings)}")
    print(f"lod2: {len(lod2_buildings)}")
    print(f"lod1 fallback: {len(lod1_buildings)}")
    print(f"watertight lod2: {watertight_count}")

    if not lod2_buildings:
        raise SystemExit("Live LOD2 smoke test failed: no LOD2 buildings were created")
    if watertight_count != len(lod2_buildings):
        raise SystemExit("Live LOD2 smoke test failed: at least one LOD2 shell is not watertight")

    city.save_cityjson(output_path)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
