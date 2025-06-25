import dtcc_core
from dtcc_core.model import City, Bounds
import dtcc_viewer

city = City()
h = 500.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h)
city.bounds = bounds

city.download_footprints()
city.download_pointcloud()
city.build_lod1_buildings()
city.view()
