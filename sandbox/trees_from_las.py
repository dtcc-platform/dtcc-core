from dtcc_core import io, builder
from dtcc_core.builder import find_tree_tops, tree_crown_polygons
from dtcc_core.model import Tree, Raster, Bounds, City
import dtcc_viewer

import fiona
import shapely


h = 2000.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h)
#
# # Download point cloud and footprints
pc = io.download_pointcloud(bounds=bounds)
buildings = io.download_footprints(bounds=bounds)
pc = pc.remove_global_outliers(3.0)
#
# pc.save("test_pointcloud.las")
# city = City()
# city.add_buildings(buildings)
#
tree_raster: Raster = builder.tree_raster_from_pointcloud(
    pc,
    buildings=buildings,
    tree_type="urban",
    cell_size=0.5,
    smallest_cluster=4,
    fill_hole_size=1.5,
)
tree_raster.save("tree_demo3_n5.tif")
# city = City()
# city.add_buildings(buildings)
# city.save_building_footprints("tree_demo_footprints2.gpkg")
# print("tree raster built")
# coords, height, radius = find_tree_tops(tree_raster, tree_type="urban")
# tree_polygons, tree_mask = tree_crown_polygons(tree_raster, coords, "urban")
#
# schema = {
#     "geometry": "Polygon",
#     "properties": {
#         "id": "int",
#     },
# }
# with fiona.open("tree_polygons.gpkg", "w", driver="GPKG", schema=schema) as dst:
#     for idx, tree_polygon in enumerate(tree_polygons):
#         dst.write(
#             {
#                 "geometry": shapely.geometry.mapping(tree_polygon),
#                 "properties": {
#                     "id": idx,
#                 },
#             }
#         )

# trees = builder.trees_from_pointcloud(pc, radius_from_segmentation=True)
# print(f"Found {len(trees)} trees")
#
# io.trees.save_trees(trees, "tree_r_seg_area2_rf.gpkg", as_circles=True)


# labels = segment_tree_crowns(tree_raster.data, coords, 2.5)

pass
#
# maxima_footprint = ski.morphology.disk(3, dtype=bool)
# local_max = ski.morphology.local_maxima(
#     tree_raster.data, footprint=np.ones((3, 3))
# ).astype(np.float64)
# local_max *= 100
# local_max_raster = tree_raster.copy(no_data=True)
# local_max_raster.data = local_max
# tree_raster.data = tree_raster.data + local_max
#
# tree_pc = tree_raster.to_pointcloud(nodata=0)

# tree_raster.view()

pass
