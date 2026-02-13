import logging
from .dataset import DatasetDescriptor, DatasetBaseArgs

# Import dataset classes to trigger auto-registration
from .pointcloud import PointCloudDataset
from .buildings import BuildingDataset
from .terrain_surface_mesh import TerrainSurfaceMeshDataset
from .city_surface_mesh import CitySurfaceMeshDataset
from .city_flat_mesh import CityFlatMeshDataset
from .city_volume_mesh import CityVolumeMeshDataset
from .air_quality import AirQualityDataset
from .trees import TreesDataset
from .weather import WeatherDataset
from .hydrology import HydrologyDataset
from .ocean import OceanDataset

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
point_cloud = get_dataset("point_cloud")
buildings = get_dataset("buildings")
terrain_surface_mesh = get_dataset("terrain_surface_mesh")
city_surface_mesh = get_dataset("city_surface_mesh")
city_flat_mesh = get_dataset("city_flat_mesh")
city_volume_mesh = get_dataset("city_volume_mesh")
air_quality = get_dataset("air_quality")
trees = get_dataset("trees")
weather = get_dataset("weather")
hydrology = get_dataset("hydrology")
ocean = get_dataset("ocean")


def info():
    """
    Print information about all available datasets.

    This function iterates over all registered datasets and prints
    a nicely formatted summary of each one.

    Example:
        >>> import dtcc_core.datasets as datasets
        >>> datasets.info()
    """
    datasets_dict = list()

    if not datasets_dict:
        print("=" * 70)
        print("DTCC Datasets")
        print("=" * 70)
        print("No datasets are currently registered.")
        print("=" * 70)
        return

    print()
    print("=" * 70)
    print(f"DTCC Datasets ({len(datasets_dict)} available)")
    print("=" * 70)
    print()
    print("Use datasets.<name>() to access a dataset.")
    print("Use print(datasets.<name>) to see dataset parameters.")
    print()

    # Print a summary table
    print("Available datasets:")
    print("-" * 70)
    for name, dataset in datasets_dict.items():
        desc = (
            dataset.description[:50] + "..."
            if len(dataset.description) > 50
            else dataset.description
        )
        print(f"  • {name:20s} - {desc}")

    print()
    print("=" * 70)
    print()
    print("For detailed information on a specific dataset, use:")
    print("  print(datasets.<name>)")
    print()


def __getattr__(name):
    """
    Dynamic attribute lookup for registered datasets.

    This enables access like `datasets.pointcloud`.

    Args:
        name (str): The attribute name being accessed.

    Returns:
        The registered dataset instance, if found.

    Raises:
        AttributeError: If the attribute is not a registered dataset.
    """
    try:
        return get_dataset(name)
    except KeyError:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
