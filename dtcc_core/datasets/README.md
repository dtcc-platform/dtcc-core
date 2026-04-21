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

Graceful Upstream Handling
--------------------------

Live/network-backed datasets should follow this contract::

    from dtcc_core import datasets

    sensors = datasets.weather(
        bounds=[minx, miny, maxx, maxy],
        strict_live=False,  # default
    )

    if sensors.attributes.get("partial_result"):
        print("Result is incomplete due to upstream failures")
        print(sensors.attributes["upstream_errors"])

Default mode (`strict_live=False`) degrades gracefully on upstream failures and
returns an empty or partial `SensorCollection` with health metadata in
`.attributes`, including:

- `partial_result`
- `upstream_error_count`
- `upstream_errors`
- `stations_skipped_upstream`

Datasets that fetch multiple parameters may also report:

- `requested_parameters`
- `fetched_parameters`

Use `strict_live=True` when you want upstream failures to raise a typed
`DatasetUpstreamError` instead of returning a degraded result. User input
errors such as malformed bounds or unknown parameters still fail fast with a
clear exception.

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
