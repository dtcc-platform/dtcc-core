"""Rasterize ground points into a terrain elevation model."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
pointcloud = dtcc.io.data.download_pointcloud(bounds=bounds)

raster = dtcc.builder.build_terrain_raster(pointcloud, cell_size=2.0, ground_only=True)
raster.plot()
raster.view()
