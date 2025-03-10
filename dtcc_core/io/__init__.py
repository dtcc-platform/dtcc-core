from . import pointcloud
from . import meshes

# from . import landuse
# from . import roadnetwork
from . import raster
from . import info
from . import footprints
from . import city
from . import roadnetwork
from . import landuse
from . import pointcloud_directory
from .cityjson import cityjson

load_pointcloud = pointcloud.load
save_pointcloud = pointcloud.save

load_pointcloud_directory = pointcloud_directory.load

load_raster = raster.load
save_raster = raster.save

load_mesh = meshes.load_mesh
load_volume_mesh = meshes.load_volume_mesh
load_mesh_as_city = meshes.load_mesh_as_city

save_mesh = meshes.save
save_volume_mesh = meshes.save
list_mesh_io = meshes.list_io
print_mesh_io = meshes.print_io

load_footprints = footprints.load
save_footprints = footprints.save

load_landuse = landuse.load

load_roadnetwork = roadnetwork.load
# save_roadnetwork = roadnetwork.save

load_cityjson = cityjson.load
load_city = city.load
save_city = city.save


from ..model import City, PointCloud, Raster, Mesh, VolumeMesh, City, RoadNetwork

City.add_methods(load_city, "load")
City.add_methods(save_city, "save")


# City.add_methods(save_city, "save")
PointCloud.add_methods(save_pointcloud, "save")
Raster.add_methods(save_raster, "save")
Mesh.add_methods(save_mesh, "save")
VolumeMesh.add_methods(save_mesh, "save")

RoadNetwork.add_methods(roadnetwork.to_dataframe, "to_df")

__all__ = [
    "load_mesh",
    "save_mesh",
    "load_volume_mesh",
    "save_volume_mesh",
    "load_pointcloud",
    "save_pointcloud",
    "load_raster",
    "save_raster",
    "load_city",
    "load_footprints",
    "save_footprints",
    "load_cityjson",
    "load_roadnetwork",
]
