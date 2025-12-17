# This demo illustrates how to build a city mesh from a point cloud and footprints.

from dtcc_core.model import City, Bounds, GeometryType
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
city.download_pointcloud(bounds=bounds, filter_on_z_bounds=True)
city.download_footprints(bounds=bounds)
city.building_heights_from_pointcloud()


# city.save_building_footprints("footprints_original.gpkg")
lod = GeometryType.LOD0
merge_tolerance = 1.0
min_building_area = 15
min_building_detail = 1.0


# merged_buildings = merge_building_footprints(
#     city.buildings, lod, max_distance=merge_tolerance, min_area=min_building_area
# )
#
# smallest_hole = max(min_building_detail, min_building_detail**2)
# cleaned_footprints = clean_building_footprints(
#     merged_buildings,
#     clearance=min_building_detail,
#     smallest_hole_area=smallest_hole,
# )

# city.replace_buildings(cleaned_footprints)
# city.save_building_footprints("footprints_merged_cleaned_1.0.gpkg")

# merged_footprints = merge_building_footprints(city.buildings, max_distance=1.0)
# city.replace_buildings(merged_footprints)
# city.save_building_footprints("footprints_merged_sb3.gpkg")
#
# fixed_footprints = clean_building_footprints(
#     merged_footprints, clearance=2.0, smallest_hole_area=4.0
# )
#
# city.replace_buildings(fixed_footprints)
# city.save_building_footprints("footprints_merged_fixed4.gpkg")
#
# simplified_footprints = simplify_building_footprints(fixed_footprints, tolerance=1)
# city.replace_buildings(simplified_footprints)
#
# city.save_building_footprints("footprints_merged_fixed_simplifed.gpkg")


surface_mesh = city.build_surface_mesh(min_building_detail=1.0, smoothing=0)
surface_mesh.save("city_surface_mesh7.vtk")

# # View mesh
surface_mesh.view()
