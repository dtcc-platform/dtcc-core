# This demo illustrates how to build a city mesh from a point cloud and footprints.

from dtcc_core.model import City, Bounds
import dtcc_viewer
from dtcc_core.builder.building.modify import (
    merge_building_footprints,
    fix_building_footprint_clearance,
    simplify_building_footprints,
    clean_building_footprints,
)
import sys

city = City()
h = 2000.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h, 0, 200)
city.bounds = bounds
#
# # Download pointcloud and building footprints
# city.download_pointcloud(bounds=bounds, filter_on_z_bounds=True)
city.download_footprints(bounds=bounds)

city.save_building_footprints("footprints_original.gpkg")

merged_footprints = merge_building_footprints(city.buildings, max_distance=1.0)
city.replace_buildings(merged_footprints)
city.save_building_footprints("footprints_merged_sb3.gpkg")

fixed_footprints = clean_building_footprints(
    merged_footprints, clearance=2.0, smallest_hole_area=4.0
)

city.replace_buildings(fixed_footprints)
city.save_building_footprints("footprints_merged_fixed4.gpkg")

simplified_footprints = simplify_building_footprints(fixed_footprints, tolerance=1)
city.replace_buildings(simplified_footprints)

city.save_building_footprints("footprints_merged_fixed_simplifed.gpkg")
sys.exit()

# city.building_heights_from_pointcloud()
# surface_mesh = city.build_surface_mesh()
#
# # View mesh
# surface_mesh.view()
