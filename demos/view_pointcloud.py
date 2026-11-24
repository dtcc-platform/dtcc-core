"""Inspect LiDAR points in a plot and in the viewer, coloured by elevation."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)
pointcloud = dtcc.io.data.download_pointcloud(bounds=bounds)
pointcloud.plot()
pointcloud.view(data=pointcloud.points[:, 2])
