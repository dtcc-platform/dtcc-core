"""Triangulate terrain from ground points using the builder API."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)
pointcloud = dtcc.io.data.download_pointcloud(bounds=bounds)
pointcloud = pointcloud.remove_global_outliers(3.0)

mesh = dtcc.builder.build_terrain_surface_mesh(pointcloud, max_mesh_size=10.0)
mesh.plot()
mesh.view()
