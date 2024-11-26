from ...model import Mesh, PointCloud, Raster
import numpy as np

from ..pointcloud.convert import rasterize


def mesh_to_raster(mesh: Mesh, cell_size: float) -> Raster:
    unique_points = np.unique(mesh.faces.flatten())
    unique_vertices = mesh.vertices[unique_points]
    pc = PointCloud(points=unique_vertices)
    pc.calculate_bounds()
    raster = rasterize(pc, cell_size, radius=cell_size, ground_only=False)
    return raster

