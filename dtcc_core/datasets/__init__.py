# __init__.py
from .pointcloud import PointCloudDataset
from .buildings import BuildingDataset
from .terrain import TerrainDataset

pointcloud = PointCloudDataset()
buildings = BuildingDataset()
terrain = TerrainDataset()

__datasets_registry = [pointcloud, buildings, terrain]


def list():
    return {dataset.name: dataset for dataset in __datasets_registry}


__all__ = [
    "pointcloud",
    "buildings",
    "terrain",
    "list",
]
