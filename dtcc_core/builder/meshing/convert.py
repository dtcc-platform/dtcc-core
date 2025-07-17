from ...model import Mesh, PointCloud, Raster
import numpy as np

from ..pointcloud.convert import rasterize


def mesh_to_raster(mesh: Mesh, cell_size: float) -> Raster:
    """
    Convert a mesh to a raster representation.
    
    This function converts a triangular mesh into a raster by extracting unique
    vertices and rasterizing them into a grid at the specified cell size.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh to convert to raster.
    cell_size : float
        The size of each raster cell in meters.
        
    Returns
    -------
    Raster
        A raster representation of the mesh.
    """
    unique_points = np.unique(mesh.faces.flatten())
    unique_vertices = mesh.vertices[unique_points]
    pc = PointCloud(points=unique_vertices)
    pc.calculate_bounds()
    raster = rasterize(pc, cell_size, radius=cell_size, ground_only=False)
    return raster

