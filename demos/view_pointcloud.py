"""Inspect LiDAR points in a plot and in the viewer, coloured by elevation."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
pointcloud = dtcc.io.data.download_pointcloud(bounds=bounds)
pointcloud.plot()
pointcloud.view(data=pointcloud.points[:, 2])
