from dtcc_core import io, builder
from dtcc_core.model import Tree, Raster, Bounds, City
import dtcc_viewer

import scipy.ndimage
import numpy as np
import skimage as ski

h = 2000.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h)

# Download point cloud and footprints
pc = io.download_pointcloud(bounds=bounds)
buildings = io.download_footprints(bounds=bounds)
pc = pc.remove_global_outliers(3.0)

pc.save("test_pointcloud.las")
city = City()
city.add_buildings(buildings)

tree_raster: Raster = builder.tree_raster_from_pointcloud(
    pc,
    buildings=buildings,
    cell_size=0.5,
    shortest_tree=2.0,
    smallest_cluster=4,
    fill_hole_size=1.5,
    sigma=0.5,
)

pc_veg = pc.get_vegetation()
# pc_veg.save("test_vegetation.las")

# tree_raster.save("test_tree_raster.tif")
# city.save_building_footprints("footprints_original.gpkg")


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

tree_raster.view()

pass
