import dtcc_core
import dtcc_data
import dtcc_viewer
from dtcc_core.model import Bounds, City, Object, GeometryType


from dtcc_core.builder.meshing import tile_surface_mesh
from time import time

# Define bounds (a residential area in Helsingborg)
city = City()
h = 1000.0
bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h, 0, 200)
city.bounds = bounds

# Download pointcloud and building footprints
city.download_pointcloud(bounds=bounds, filter_on_z_bounds=True)
city.download_footprints(bounds=bounds)

city.building_heights_from_pointcloud()

surface_mesh = city.build_surface_mesh()

# move to (0,0,0)
surface_mesh.offset_to_origin()
surface_mesh.offset([0, 0, 1])  # lift 1m above ground since we extrude to 0 by default

start = time()
tiles = surface_mesh.tile(tile_size=100, progress=True)
print(f"Tiling took {time()-start:.2f} seconds for {len(tiles)} tiles")

print_tiles = [t.create_printable_solid() for t in tiles]

# hack to view a lits of meshes
obj = Object()
obj.geometry[GeometryType.MESH] = print_tiles
obj.view()
