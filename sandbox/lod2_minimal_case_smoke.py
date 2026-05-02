#!/usr/bin/env python3
"""Smoke-test prototype LOD2 reconstruction on bundled MinimalCase data."""

from pathlib import Path

from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from dtcc_core.model import City


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    case_dir = repo_root / "tests" / "data" / "MinimalCase"
    output_path = repo_root / "sandbox" / "lod2_minimal.city.json"

    city = City()
    city.load_footprints(str(case_dir / "PropertyMap.shp"))
    city.load_pointcloud(str(case_dir / "pointcloud.las"))
    city.build_lod2_buildings(calculate_heights=True)

    lod2_buildings = [building for building in city.buildings if building.lod2 is not None]
    watertight_count = sum(is_watertight(building.lod2) for building in lod2_buildings)

    print(f"buildings: {len(city.buildings)}")
    print(f"lod2: {len(lod2_buildings)}")
    print(f"watertight lod2: {watertight_count}")

    if not lod2_buildings:
        raise SystemExit("LOD2 smoke test failed: no LOD2 buildings were created")
    if watertight_count != len(lod2_buildings):
        raise SystemExit("LOD2 smoke test failed: at least one LOD2 shell is not watertight")

    city.save_cityjson(output_path)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
