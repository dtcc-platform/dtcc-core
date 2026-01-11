import logging
from .dataset import DatasetDescriptor, DatasetBaseArgs

# Import dataset classes to trigger auto-registration
from .pointcloud import PointCloudDataset
from .buildings import BuildingDataset
from .terrain import TerrainDataset

# Import registry infrastructure from separate module to avoid circular imports
from .registry import (
    _register_dataset_class,
    register,
    register_class,
    unregister,
    list_datasets as list,
    get_dataset,
)

# Override module attributes with registered dataset instances
# This is necessary because importing from submodules makes the submodules
# themselves available as attributes, which would shadow the registered instances
pointcloud = get_dataset("pointcloud")
buildings = get_dataset("Buildings LoD1")
terrain = get_dataset("terrain")


def __getattr__(name):
    """
    Dynamic attribute lookup for registered datasets.

    This enables access like `datasets.pointcloud`.

    Args:
        name: Attribute name to look up

    Returns:
        Dataset instance if found

    Raises:
        AttributeError: If dataset not found
    """
    dataset = get_dataset(name)
    if dataset is not None:
        return dataset

    raise AttributeError(f"module 'dtcc_core.datasets' has no attribute '{name}'")


__all__ = [
    "DatasetDescriptor",
    "DatasetBaseArgs",
    "register",
    "register_class",
    "unregister",
    "list",
    "pointcloud",
    "buildings",
    "terrain",
]
