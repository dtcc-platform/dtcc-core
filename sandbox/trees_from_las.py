from dtcc_core import io, builder
from dtcc_core.model import Tree, Raster
import dtcc_viewer

import scipy.ndimage
import numpy as np
import skimage as ski

pc = io.load_pointcloud("../data/helsingborg-residential-2022/pointcloud.las")
buildings = io.load_footprints("../data/helsingborg-residential-2022/footprints.shp")

pc = pc.remove_global_outliers(3.0)

tree_raster: Raster = builder.tree_raster_from_pointcloud(
    pc,
    cell_size=0.5,
    shortest_tree=2.0,
    smallest_cluster=100,
    fill_hole_size=100,
    sigma=1.0,
)

#
# maxima_footprint = ski.morphology.disk(3, dtype=bool)
# local_max = ski.morphology.local_maxima(
#     tree_raster.data, footprint=np.ones((3, 3))
# ).astype(np.float64)
# local_max *= 100
#
# local_max_raster = tree_raster.copy(no_data=True)
# local_max_raster.data = local_max
# tree_raster.data = tree_raster.data + local_max
#
# tree_pc = tree_raster.to_pointcloud(nodata=0)

tree_raster.view()

pass
