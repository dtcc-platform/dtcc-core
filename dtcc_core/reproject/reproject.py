import numpy as np
import pyproj
from pyproj import Transformer
from dtcc_core.model import PointCloud, Building
from dtcc_core.common import warning



def reproject_pointcloud(pointcloud:PointCloud, src_crs: str, target_crs:str, override_geometry_crs:bool=False) -> PointCloud:
    if src_crs == target_crs:
        return pointcloud

    if target_crs is None or len(target_crs) == 0:
        raise ValueError("Target CRS not provided")

    if not pointcloud.transform.is_identity:
        pts = pointcloud.transform(pointcloud.points)
    else:
        pts = pointcloud.points

    if src_crs is None or len(src_crs) == 0:
        src_crs = pointcloud.transform.srs
        if src_crs is None or len(src_crs) == 0:
            src_crs = "EPSG:3006"
            warning("Source CRS not provided, defaulting to EPSG:3006")

    else:
        if not override_geometry_crs:
            if pointcloud.transform.srs:
                src_crs = pointcloud.transform.srs
            else:
                warning("Pointcloud CRS not set, defaulting to EPSG:3006")
                src_crs = "EPSG:3006"


    transformer = Transformer.from_crs(src_crs, target_crs, always_xy=True)

    x_trans, y_trans = transformer.transform(pts[:,0], pts[:,1])
    pts = np.column_stack((x_trans, y_trans, pts[:,2]))
    pc = PointCloud(points=pts, classification=pointcloud.classification,
                      intensity=pointcloud.intensity,
                      return_number=pointcloud.return_number,
                      num_returns=pointcloud.num_returns,
                      transform=pointcloud.transform)
    pc.calculate_bounds()
    pc.transform.srs = target_crs
    return pc








