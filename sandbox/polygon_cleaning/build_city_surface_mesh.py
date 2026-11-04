# This demo illustrates how to build a city mesh from a point cloud and footprints.

from dtcc_core.model import City, Bounds
import dtcc_viewer

city = City()
h = 2000.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h, 0, 200)
city.bounds = bounds
#
# # Download pointcloud and building footprints
city.download_pointcloud(bounds=bounds, filter_on_z_bounds=True)
city.download_footprints(bounds=bounds)
city.building_heights_from_pointcloud()

min_building_area = 15
min_building_detail = 1.0

# Clean the footprints once; the mesh builder reuses the recorded cleaning.
city = city.clean_footprints(
    min_building_detail=min_building_detail, min_area=min_building_area
)
# city.save_building_footprints("footprints_cleaned_1.0.gpkg")

surface_mesh = city.build_surface_mesh(
    min_building_detail=min_building_detail, smoothing=0
)
surface_mesh.save("city_surface_mesh7.vtk")

# # View mesh
surface_mesh.view()
