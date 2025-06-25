from dtcc_core import builder, io, model
import dtcc_viewer
from pathlib import Path

# Define bounds (a residential area in Helsingborg)
h = 2000.0
bounds = model.Bounds(319891, 6399790, 319891 + h, 6399790 + h)

data_dir = Path(__file__).parent / ".." / "data"
data_dir = data_dir / "helsingborg-residential-2022"

# Download pointcloud and building footprints
pointcloud = io.pointcloud.load(data_dir / "pointcloud.las")
buildings = io.load_footprints(data_dir / "footprints.shp")

# Remove global outliers
pointcloud = pointcloud.remove_global_outliers(3.0)

# Build terrain raster and mesh
raster = builder.build_terrain_raster(pointcloud, cell_size=5, ground_only=True)

# Extract roof points and compute building heights
buildings = builder.extract_roof_points(buildings, pointcloud)
buildings = builder.compute_building_heights(buildings, raster, overwrite=True)

# Build LOD1 buildings
buildings = builder.build_lod1_buildings(
    buildings, default_ground_height=raster.min, always_use_default_ground=True
)

b1 = buildings[0].lod1
b1.mesh(weld=True, snap=0.005)

building_meshes = [b.lod1.mesh(weld=True, snap=0.005) for b in buildings]
merged_mesh = builder.meshing.merge_meshes(building_meshes, weld=False, snap=0)

merged_mesh.save("merged_norm_fix.stl")
merged_mesh.view()
