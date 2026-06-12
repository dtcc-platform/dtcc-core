# This demo builds LOD2 buildings for three areas in Kungsbacka.

from pathlib import Path

import dtcc_core as dtcc
from dtcc_core.io import download_footprints, download_pointcloud
from dtcc_core.model import City


# Bounds are in EPSG:3006 / SWEREF 99 TM
AREAS = (
    {
        "name": "verkstadsgatan_3",
        "bounds": dtcc.Bounds(
            xmin=324861,
            ymin=6374485,
            xmax=325668,
            ymax=6375155,
        ),
    },
    {
        "name": "turkosvagen_4",
        "bounds": dtcc.Bounds(
            xmin=323963,
            ymin=6374269,
            xmax=324426,
            ymax=6374687,
        ),
    },
    {
        "name": "alvsakersvagen_500",
        "bounds": dtcc.Bounds(
            xmin=330476,
            ymin=6381158,
            xmax=331668,
            ymax=6382045,
        ),
    },
)


output_dir = Path("output")
output_dir.mkdir(exist_ok=True)


for area in AREAS:
    # Build LOD2 buildings (watertight LOD2 where supported, LOD1 fallback)
    city = City()
    city.bounds = area["bounds"]
    city.add_buildings(download_footprints(area["bounds"]))
    city.add_pointcloud(download_pointcloud(area["bounds"]))
    city.build_lod2_buildings(calculate_heights=True, log_rejections=True)

    lod2_count = sum(1 for building in city.buildings if building.lod2 is not None)
    print(
        f"{area['name']}: buildings={len(city.buildings)} "
        f"lod2={lod2_count} lod1_fallback={len(city.buildings) - lod2_count}"
    )

    city.save_cityjson(output_dir / f"lod2_kungsbacka_{area['name']}.city.json")
