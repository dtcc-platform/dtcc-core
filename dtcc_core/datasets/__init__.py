# __init__.py
from .pointcloud import PointCloudDataset
from .buildings import BuildingDataset
from .terrain import TerrainDataset

pointcloud = PointCloudDataset()
buildings = BuildingDataset()
terrain = TerrainDataset()
