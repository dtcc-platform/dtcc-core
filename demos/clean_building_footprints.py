"""Clean building footprints and inspect the result."""

import dtcc_core as dtcc

city = dtcc.model.City()
city.bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
city.download_footprints()

# Clean each building independently, keeping buildings separate.
isolated_city = dtcc.builder.clean_footprints(
    city, min_building_detail=0.5, isolated=True, show=True
)

# Clean the original city jointly, merging nearby footprints.
merged_city = dtcc.builder.clean_footprints(
    city, min_building_detail=0.5, isolated=False, show=True
)
