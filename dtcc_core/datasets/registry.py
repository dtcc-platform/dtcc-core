"""
Dataset registry infrastructure.

This module contains the registry and registration functions for datasets.
It's separate from __init__.py to avoid circular import issues during
module initialization.
"""

from ..common import info, warning

# Registry storage
_dataset_classes = {}  # Maps name -> class (for class-level access)
_datasets_by_name = {}  # Maps name -> instance (cached instances)
_datasets_registry = []  # List of instances (for compatibility)


def _register_dataset_class(name: str, cls: type):
    """
    Internal function called by __init_subclass__ hook to register dataset classes.

    Args:
        name: Name to register the dataset under
        cls: Dataset class to register
    """
    if name in _dataset_classes:
        warning(
            f"Dataset '{name}' already registered (class: {_dataset_classes[name].__name__}), "
            f"replacing with {cls.__name__}"
        )

    # Store class reference
    _dataset_classes[name] = cls

    # Eagerly instantiate and cache
    instance = cls()
    _datasets_by_name[name] = instance
    info(f"Registered dataset: '{name}' ({cls.__name__})")

    # Add to registry list if not already present
    if instance not in _datasets_registry:
        _datasets_registry.append(instance)


def register(name: str, instance):
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
    from .dataset import DatasetDescriptor

    if not isinstance(instance, DatasetDescriptor):
        raise TypeError(
            f"Expected DatasetDescriptor instance, got {type(instance).__name__}"
        )

    if name in _datasets_by_name:
        warning(f"Dataset '{name}' already registered, replacing it.")

    _datasets_by_name[name] = instance
    info(f"Registered dataset instance: '{name}' ({type(instance).__name__})")

    # Add to registry list if not already present
    if instance not in _datasets_registry:
        _datasets_registry.append(instance)


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
    from .dataset import DatasetDescriptor

    if not issubclass(cls, DatasetDescriptor):
        raise TypeError(f"Expected DatasetDescriptor subclass, got {cls.__name__}")

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
    if name in _datasets_by_name:
        instance = _datasets_by_name.pop(name)
        if instance in _datasets_registry:
            _datasets_registry.remove(instance)

    if name in _dataset_classes:
        _dataset_classes.pop(name)


def list_datasets():
    """
    Return all registered datasets.

    Returns:
        dict: Dictionary mapping dataset names to dataset instances

    Example:
        >>> available_datasets = list_datasets()
        >>> print(available_datasets.keys())
        dict_keys(['pointcloud', 'buildings', 'terrain'])
    """
    return {name: dataset for name, dataset in _datasets_by_name.items()}


def get_dataset(name: str):
    """
    Get a dataset by name.

    Args:
        name: Name of the dataset

    Returns:
        Dataset instance if found

    Raises:
        KeyError: If dataset name is not registered

    Example:
        >>> dataset = get_dataset('pointcloud')
        >>> result = dataset(bounds=[...])
    """
    if name not in _datasets_by_name:
        raise KeyError(f"Dataset '{name}' not found in registry")
    return _datasets_by_name[name]
