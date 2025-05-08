from dtcc_core import io, builder
from dtcc_core.model import Tree
import dtcc_viewer

pc = io.load_pointcloud("../data/helsingborg-residential-2022/pointcloud.las")
buildings = io.load_footprints("../data/helsingborg-residential-2022/footprints.shp")

pc = pc.remove_global_outliers(3.0)

terrain_raster = builder.build_terrain_raster(pc, cell_size=0.5, ground_only=True)

trees = pc.get_vegetation()

building_polygon_footprints = [b.get_footprint().to_polygon() for b in buildings]

tree_raster = trees.rasterize(
    cell_size=0.5,
    bounds=terrain_raster.bounds,
    fill_holes=False,
    window_size=1,
)

tree_raster.view()

pass
