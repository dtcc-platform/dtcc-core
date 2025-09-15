import dtcc_core
import dtcc_data
import dtcc_viewer
from dtcc_core.model import Bounds, City


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

surface_mesh = city.build_city_surface_mesh()


start = time()
tiles = tile_surface_mesh(surface_mesh, tile_size=250)
print(f"Tiling took {time()-start:.2f} seconds for {len(tiles)} tiles")

view_tiles = [t for idx, t in enumerate(tiles) if idx % 3 == 0]

# dtcc_core.builder.meshing.merge_meshes(view_tiles).view()

test_tile = tiles[14]
capped_tile = test_tile.create_printable_solid(extrusion_depth=20)

capped_tile.view()
