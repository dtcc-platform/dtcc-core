from dtcc_core.model import City, Bounds

h = 2000.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h)

city = City()
city.bounds = bounds

city.download_footprints()
city.download_pointcloud()

city.build_trees_from_pointcloud(tree_type="urban")

city.save_trees("trees_demo4.gpkg", save_as_circles=False)
