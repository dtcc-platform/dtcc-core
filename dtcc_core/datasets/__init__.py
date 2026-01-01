"""
DTCC Datasets - Dynamic Dataset Registration System
===================================================

This module provides a dynamic registration system for datasets in dtcc_core.
Datasets automatically register when they are defined, making them discoverable
and accessible at runtime.

Basic Usage
----------

Access built-in datasets::

    from dtcc_core import datasets

    # Use built-in datasets
    pointcloud = datasets.pointcloud(bounds=[minx, miny, maxx, maxy])
    buildings = datasets.buildings(bounds=[minx, miny, maxx, maxy])
    terrain = datasets.terrain(bounds=[minx, miny, maxx, maxy])

    # List all available datasets
    available = datasets.list()
    print(available.keys())  # ['pointcloud', 'buildings', 'terrain', ...]

Creating Custom Datasets
-----------------------

Define a custom dataset by inheriting from DatasetDescriptor.
The dataset will automatically register when the class is defined::

    from dtcc_core.datasets import DatasetDescriptor, DatasetBaseArgs
    from pydantic import Field

    class MyDatasetArgs(DatasetBaseArgs):
        custom_param: str = Field(..., description="Custom parameter")

    class MyDataset(DatasetDescriptor):
        name = "my_dataset"
        description = "My custom dataset"
        ArgsModel = MyDatasetArgs

        def build(self, args):
            # Implement data fetching/processing
            bounds = self.parse_bounds(args.bounds)
            # ... fetch data using bounds and custom_param
            return result

    # No registration needed - it's automatic!
    # Access it like any other dataset
    from dtcc_core import datasets
    result = datasets.my_dataset(
        bounds=[minx, miny, maxx, maxy],
        custom_param="value"
    )

Abstract Base Classes
-------------------

To create abstract base classes that should NOT register, use register=False::

    class AbstractSpatialDataset(DatasetDescriptor, register=False):
        \"\"\"Abstract base class - won't be registered.\"\"\"

        def common_spatial_logic(self):
            # Shared logic for spatial datasets
            pass

    class ConcreteDataset(AbstractSpatialDataset):
        name = "concrete"  # This WILL register
        ArgsModel = MyArgs

        def build(self, args):
            self.common_spatial_logic()
            return result

Explicit Registration API
-------------------------

For advanced use cases, explicit registration is also available::

    from dtcc_core.datasets import register, register_class, unregister

    # Register an instance
    my_dataset = MyDataset()
    register("my_custom_name", my_dataset)

    # Register a class (will instantiate it)
    register_class("another_name", AnotherDataset)

    # Unregister (useful for testing)
    unregister("my_custom_name")

External Plugins
--------------

External packages can provide datasets by simply defining them::

    # In your external package: my_dtcc_plugin/datasets.py
    from dtcc_core.datasets import DatasetDescriptor, DatasetBaseArgs
    from pydantic import Field

    class PluginArgs(DatasetBaseArgs):
        plugin_param: str = Field(..., description="Plugin parameter")

    class PluginDataset(DatasetDescriptor):
        name = "plugin_dataset"
        description = "Dataset from external plugin"
        ArgsModel = PluginArgs

        def build(self, args):
            # Plugin implementation
            return result

    # In user code:
    from dtcc_core import datasets
    import my_dtcc_plugin.datasets  # Import triggers registration

    result = datasets.plugin_dataset(bounds=[...], plugin_param="value")

Dataset Discovery
---------------

All registered datasets can be discovered programmatically::

    from dtcc_core import datasets

    # Get all datasets
    all_datasets = datasets.list()

    # Iterate through datasets
    for name, dataset in all_datasets.items():
        print(f"{name}: {dataset.description}")

        # Get parameter schema
        schema = dataset.show_options()
        print(schema)

Notes
-----
- Datasets are instantiated eagerly when classes are imported
- The `name` attribute determines the dataset's registered name
- Datasets without a name attribute won't register
- Later registrations with the same name replace earlier ones (with a warning)
- All datasets support `.show_options()` to get Pydantic JSON schema
- Datasets are callable: `dataset(bounds=[...], param=value)`

See Also
--------
- DatasetDescriptor : Base class for all datasets
- DatasetBaseArgs : Base class for dataset arguments with bounds validation
"""

import logging
from .dataset import DatasetDescriptor, DatasetBaseArgs

# Import dataset classes to trigger auto-registration
from .pointcloud import PointCloudDataset
from .buildings import BuildingDataset
from .terrain import TerrainDataset

# Registry infrastructure
__dataset_classes = {}  # Maps name -> class (for class-level access)
__datasets_by_name = {}  # Maps name -> instance (cached instances)
__datasets_registry = []  # List of instances (for compatibility)


def _register_dataset_class(name: str, cls: type):
    """
    Internal function called by __init_subclass__ hook to register dataset classes.

    Args:
        name: Name to register the dataset under
        cls: Dataset class to register
    """
    if name in __dataset_classes:
        logging.warning(
            f"Dataset '{name}' already registered (class: {__dataset_classes[name].__name__}), "
            f"replacing with {cls.__name__}"
        )

    # Store class reference
    __dataset_classes[name] = cls

    # Eagerly instantiate and cache
    instance = cls()
    __datasets_by_name[name] = instance

    # Add to registry list if not already present
    if instance not in __datasets_registry:
        __datasets_registry.append(instance)


def register(name: str, instance: DatasetDescriptor):
    """
    Explicitly register a dataset instance.

    Args:
        name: Name to register the dataset under
        instance: Instance of DatasetDescriptor to register

    Raises:
        TypeError: If instance is not a DatasetDescriptor

    Example:
        >>> my_dataset = MyDataset()
        >>> register("my_dataset", my_dataset)
    """
    if not isinstance(instance, DatasetDescriptor):
        raise TypeError(
            f"Expected DatasetDescriptor instance, got {type(instance).__name__}"
        )

    if name in __datasets_by_name:
        logging.warning(f"Dataset '{name}' already registered, replacing it.")

    __datasets_by_name[name] = instance

    # Add to registry list if not already present
    if instance not in __datasets_registry:
        __datasets_registry.append(instance)


def register_class(name: str, cls: type, **init_kwargs):
    """
    Register a dataset class (instantiates it).

    Args:
        name: Name to register the dataset under
        cls: Subclass of DatasetDescriptor to register
        **init_kwargs: Keyword arguments to pass to class constructor

    Raises:
        TypeError: If cls is not a DatasetDescriptor subclass

    Example:
        >>> register_class("my_dataset", MyDataset, some_param="value")
    """
    if not issubclass(cls, DatasetDescriptor):
        raise TypeError(
            f"Expected DatasetDescriptor subclass, got {cls.__name__}"
        )

    instance = cls(**init_kwargs)
    register(name, instance)


def unregister(name: str):
    """
    Remove a dataset from the registry.

    Args:
        name: Name of the dataset to unregister

    Example:
        >>> unregister("my_dataset")
    """
    if name in __datasets_by_name:
        instance = __datasets_by_name.pop(name)
        if instance in __datasets_registry:
            __datasets_registry.remove(instance)

    if name in __dataset_classes:
        __dataset_classes.pop(name)


def list():
    """
    Return all registered datasets.

    Returns:
        dict: Dictionary mapping dataset names to dataset instances

    Example:
        >>> available_datasets = list()
        >>> print(available_datasets.keys())
        dict_keys(['pointcloud', 'buildings', 'terrain'])
    """
    return {name: dataset for name, dataset in __datasets_by_name.items()}


def __getattr__(name):
    """
    Dynamic attribute lookup for registered datasets.

    This enables backward-compatible access like `datasets.pointcloud`.

    Args:
        name: Attribute name to look up

    Returns:
        Dataset instance if found

    Raises:
        AttributeError: If dataset not found
    """
    if name in __datasets_by_name:
        return __datasets_by_name[name]

    raise AttributeError(f"module 'dtcc_core.datasets' has no attribute '{name}'")


# Legacy function for backward compatibility (now deprecated)
def register_dataset(dataset: DatasetDescriptor):
    """
    Register a dataset (deprecated - use register() instead).

    Args:
        dataset: Dataset instance to register
    """
    logging.warning(
        "register_dataset() is deprecated, use register() instead"
    )
    if hasattr(dataset, 'name'):
        register(dataset.name, dataset)
    else:
        raise ValueError("Dataset must have a 'name' attribute")


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
