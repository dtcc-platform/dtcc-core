import dtcc_core

# Define bounds (a residential area in Helsingborg)
h = 2000.0
bounds = dtcc_core.model.Bounds(319891, 6399790, 319891 + h, 6399790 + h)

# Download point cloud
pointcloud = dtcc_core.io.download_pointcloud(bounds=bounds)

# Remove global outliers
pointcloud = pointcloud.remove_global_outliers(3.0)

# Build terrain raster
for error in (0, 0.5, 2, 5):
    mesh = dtcc_core.builder.adaptive_terrain_mesh(pointcloud, max_error=error)
    mesh.save(f"adaptive_mesh_{error}.obj")
