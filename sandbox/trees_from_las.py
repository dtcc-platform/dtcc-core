from dtcc_core import io, builder
from dtcc_core.model import Tree, Raster
import dtcc_viewer

pc = io.load_pointcloud("../data/helsingborg-residential-2022/pointcloud.las")
buildings = io.load_footprints("../data/helsingborg-residential-2022/footprints.shp")

pc = pc.remove_global_outliers(3.0)


terrain_raster = builder.build_terrain_raster(pc, cell_size=0.5, ground_only=True)

trees = pc.get_vegetation()

building_polygon_footprints = [b.get_footprint().to_polygon() for b in buildings]

tree_raster: Raster = trees.rasterize(
    cell_size=0.5,
    bounds=terrain_raster.bounds,
    fill_holes=False,
    window_size=1,
)

shortest_tree = 2


tree_raster = tree_raster.erode_small_lines(neighborhood_size=None, nodata=0)
tree_raster = tree_raster.remove_small_masks(min_size=100, nodata=0)
tree_raster = tree_raster.fill_small_holes(hole_size=100, nodata=0)

tree_raster.data -= terrain_raster.data
tree_raster.data[tree_raster.data < shortest_tree] = 0


tree_raster.view()


pass
