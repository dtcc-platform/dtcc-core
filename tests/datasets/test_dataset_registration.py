import pytest
import logging

try:
    from pydantic import BaseModel, Field, ValidationError
except ImportError:
    # Pydantic should be installed as dtcc_core dependency
    pytest.skip("pydantic not available", allow_module_level=True)

import dtcc_core.datasets as datasets
from dtcc_core.datasets import (
    DatasetDescriptor,
    DatasetBaseArgs,
    register,
    register_class,
    unregister,
    list as list_datasets,
)


# Test Fixtures - Custom Dataset Classes for Testing


class BaseTestArgs(DatasetBaseArgs):
    """Test arguments for custom datasets."""

    test_param: str = Field(..., description="Test parameter")


class BaseTestDataset(DatasetDescriptor):
    """Test dataset that auto-registers."""

    name = "test_dataset"
    description = "A test dataset"
    ArgsModel = BaseTestArgs

    def build(self, args):
        return f"test_result_{args.test_param}"


class AnotherTestArgs(DatasetBaseArgs):
    """Arguments for another test dataset."""

    value: int = Field(default=42, description="Test value")


class AnotherTestDataset(DatasetDescriptor):
    """Another test dataset."""

    name = "another_test"
    description = "Another test dataset"
    ArgsModel = AnotherTestArgs

    def build(self, args):
        return args.value * 2


class AbstractTestDataset(DatasetDescriptor, register=False):
    """Abstract test dataset that should NOT register."""

    name = "abstract_test"
    description = "Should not be registered"
    ArgsModel = BaseTestArgs

    def build(self, args):
        return "abstract_result"


class ConcreteTestDataset(AbstractTestDataset):
    """Concrete dataset inheriting from abstract - should register."""

    name = "concrete_test"
    description = "Concrete test dataset"

    def build(self, args):
        return f"concrete_{args.test_param}"


class NoNameDataset(DatasetDescriptor):
    """Dataset without a name - should NOT register."""

    name = ""
    ArgsModel = BaseTestArgs

    def build(self, args):
        return "no_name_result"


# Auto-Registration Tests


def test_auto_registration_on_class_definition():
    """Test that datasets auto-register when class is defined."""
    # TestDataset should have auto-registered during import
    available = list_datasets()
    assert "test_dataset" in available
    assert isinstance(available["test_dataset"], BaseTestDataset)


def test_multiple_datasets_register():
    """Test that multiple datasets can be registered."""
    available = list_datasets()
    assert "test_dataset" in available
    assert "another_test" in available


def test_abstract_class_does_not_register():
    """Test that datasets with register=False don't auto-register."""
    available = list_datasets()
    assert "abstract_test" not in available


def test_concrete_subclass_registers():
    """Test that concrete subclasses of abstract classes register."""
    available = list_datasets()
    assert "concrete_test" in available
    assert isinstance(available["concrete_test"], ConcreteTestDataset)


def test_dataset_without_name_does_not_register():
    """Test that datasets without a name attribute don't register."""
    available = list_datasets()
    # NoNameDataset has empty name, should not register
    assert "" not in available


def test_dataset_base_args_reject_short_bounds():
    """DatasetBaseArgs should fail fast on malformed bounds length."""
    with pytest.raises(ValueError, match="Bounds must be 4 or 6 floats"):
        BaseTestArgs(bounds=(0.0, 1.0, 2.0), test_param="x")


def test_dataset_base_args_reject_inverted_bounds():
    """DatasetBaseArgs should fail fast on inverted 2D bounds."""
    with pytest.raises(ValueError, match="Invalid bounds: xmin < xmax, ymin < ymax"):
        BaseTestArgs(bounds=(2.0, 1.0, 1.0, 3.0), test_param="x")


def test_dataset_base_args_reject_unknown_arguments():
    """DatasetBaseArgs should fail fast on misspelled or unsupported options."""
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        BaseTestArgs(
            bounds=(0.0, 1.0, 2.0, 3.0),
            test_param="x",
            unknown_option=True,
        )


def test_descriptor_extracts_format_metadata_from_schema():
    """Descriptor metadata should expose formats without Atlas parsing schemas."""
    ds = datasets.point_cloud

    assert ds.list_supported_formats() == ["copc", "las", "laz"]
    assert ds.format_metadata()[0] == {
        "format": "copc",
        "extension": "copc",
        "media_type": "application/octet-stream",
        "data_kind": "point_cloud",
        "multi_file": False,
    }


def test_descriptor_describe_returns_dataset_contract():
    """describe() should return the stable dataset contract for clients."""
    ds = datasets.city_volume_mesh
    metadata = ds.describe()

    assert metadata["name"] == "city_volume_mesh"
    assert metadata["title"] == "City Volume Mesh"
    assert metadata["data_category"] == "derived"
    assert metadata["result_kind"] == "mesh"
    assert metadata["python_return_type"] == "dtcc_core.model.VolumeMesh"
    assert metadata["supported_formats"] == ["xdmf", "vtu"]
    assert metadata["multi_file_formats"] == ["xdmf"]
    assert metadata["serialization"]["python_object_when_format_omitted"] is True
    assert metadata["serialization"]["bytes_when_format_is_set"] is True
    assert metadata["serialization"]["format_parameter"] is True

    xdmf = next(item for item in metadata["formats"] if item["format"] == "xdmf")
    assert xdmf["data_kind"] == "mesh"
    assert xdmf["multi_file"] is True


def test_dataset_info_groups_rows_by_category():
    """Dataset info summary rows should be grouped in the public category order."""
    grouped = datasets._group_dataset_summary_rows(
        {
            "smoke": datasets.smoke,
            "city": datasets.city,
            "point_cloud": datasets.point_cloud,
        }
    )

    assert list(grouped) == ["raw", "derived", "simulation"]
    assert grouped["raw"][0][0] == "point_cloud"
    assert grouped["derived"][0][0] == "city"
    assert grouped["simulation"][0][0] == "smoke"


def test_dataset_info_prints_one_table_per_category(monkeypatch):
    """datasets.info() should delegate grouped output to the common table logger."""
    messages = []
    tables = []

    monkeypatch.setattr(datasets, "log_info", messages.append)

    def capture_table(log_fn, message, columns, rows, **kwargs):
        tables.append((message, columns, rows, kwargs))

    monkeypatch.setattr(datasets, "log_table", capture_table)

    datasets.info()

    table_titles = [table[0] for table in tables]
    assert "Raw Datasets" in table_titles[0]
    assert any(title.startswith("Derived Datasets") for title in table_titles)
    assert any(title.startswith("Simulation Datasets") for title in table_titles)
    assert all(len(table[2]) > 0 for table in tables)
    assert any("DTCC Datasets" in message for message in messages)


# Explicit API Tests


def test_register_instance():
    """Test explicit registration of a dataset instance."""
    # Create a custom dataset instance
    custom_dataset = BaseTestDataset()

    # Register it with a custom name
    register("custom_test", custom_dataset)

    # Verify it's registered
    available = list_datasets()
    assert "custom_test" in available
    assert available["custom_test"] is custom_dataset

    # Cleanup
    unregister("custom_test")


def test_register_class():
    """Test registration by class."""
    # Register the class
    register_class("class_test", AnotherTestDataset)

    # Verify it's registered and instantiated
    available = list_datasets()
    assert "class_test" in available
    assert isinstance(available["class_test"], AnotherTestDataset)

    # Cleanup
    unregister("class_test")


def test_register_invalid_instance():
    """Test that registering invalid instance raises TypeError."""
    with pytest.raises(TypeError):
        register("invalid", "not a dataset")


def test_register_class_invalid():
    """Test that registering invalid class raises TypeError."""
    with pytest.raises(TypeError):
        register_class("invalid", str)


def test_unregister():
    """Test unregistering a dataset."""
    # Register a dataset
    register_class("temp_test", BaseTestDataset)
    assert "temp_test" in list_datasets()

    # Unregister it
    unregister("temp_test")
    assert "temp_test" not in list_datasets()


def test_unregister_nonexistent():
    """Test that unregistering a non-existent dataset doesn't raise error."""
    # Should not raise
    unregister("nonexistent_dataset")


# Discovery Tests


def test_list_returns_all_datasets():
    """Test that list() returns all registered datasets."""
    available = list_datasets()

    # Should include built-in datasets
    assert "point_cloud" in available
    assert "buildings" in available
    assert "roads" in available
    assert "terrain_surface_mesh" in available

    # Should include test datasets
    assert "test_dataset" in available
    assert "another_test" in available


def test_list_returns_dict():
    """Test that list() returns a dictionary."""
    available = list_datasets()
    assert isinstance(available, dict)


def test_list_values_are_instances():
    """Test that list() values are dataset instances."""
    available = list_datasets()
    for name, dataset in available.items():
        assert isinstance(dataset, DatasetDescriptor)


# Backward Compatibility Tests


def test_module_attribute_access_pointcloud():
    """Test access to point_cloud dataset."""
    # Should be able to access via module attribute
    pc = datasets.point_cloud
    assert isinstance(pc, DatasetDescriptor)
    assert pc.name == "point_cloud"


def test_module_attribute_access_buildings():
    """Test access to buildings dataset."""
    buildings = datasets.buildings
    assert isinstance(buildings, DatasetDescriptor)
    assert buildings.name == "buildings"


def test_module_attribute_access_terrain():
    """Test backward-compatible access to terrain surface mesh dataset."""
    terrain_surface_mesh = datasets.terrain_surface_mesh
    assert isinstance(terrain_surface_mesh, DatasetDescriptor)
    assert terrain_surface_mesh.name == "terrain_surface_mesh"


def test_module_attribute_access_roads():
    """Test access to roads dataset."""
    roads = datasets.roads
    assert isinstance(roads, DatasetDescriptor)
    assert roads.name == "roads"


def test_module_attribute_access_custom():
    """Test that custom datasets are accessible via module attribute."""
    # TestDataset should be accessible
    test_ds = datasets.test_dataset
    assert isinstance(test_ds, BaseTestDataset)


def test_module_attribute_nonexistent_raises():
    """Test that accessing non-existent dataset raises AttributeError."""
    with pytest.raises(AttributeError):
        _ = datasets.nonexistent_dataset


def test_callable_interface():
    """Test that datasets remain callable."""
    # Create a test dataset call (with mock data to avoid actual API calls)
    test_ds = datasets.test_dataset

    # Should be callable (will fail validation since we don't provide bounds)
    assert callable(test_ds)


def test_show_options():
    """Test that datasets have show_options() method."""
    test_ds = datasets.test_dataset
    schema = test_ds.show_options()

    # Should return a schema dict
    assert isinstance(schema, dict)
    assert "properties" in schema


# Edge Cases


def test_inheritance_hierarchy():
    """Test that inheritance hierarchies work correctly."""
    # Concrete class should register
    assert "concrete_test" in list_datasets()

    # Abstract parent should not
    assert "abstract_test" not in list_datasets()

    # Concrete should work
    concrete = datasets.concrete_test
    assert isinstance(concrete, ConcreteTestDataset)


def test_dataset_name_attribute_used():
    """Test that dataset name comes from the name attribute."""
    test_ds = datasets.test_dataset
    assert test_ds.name == "test_dataset"

    another_ds = datasets.another_test
    assert another_ds.name == "another_test"


def test_multiple_instances_same_class():
    """Test registering multiple instances of the same class."""
    # Create two instances
    instance1 = BaseTestDataset()
    instance2 = BaseTestDataset()

    # Register with different names
    register("instance1", instance1)
    register("instance2", instance2)

    # Both should be registered
    available = list_datasets()
    assert "instance1" in available
    assert "instance2" in available

    # They should be different instances
    assert available["instance1"] is instance1
    assert available["instance2"] is instance2

    # Cleanup
    unregister("instance1")
    unregister("instance2")


# External Plugin Simulation


def test_external_plugin_pattern():
    """Test that external plugins can register datasets."""

    # Simulate an external plugin defining a dataset
    class ExternalPluginArgs(DatasetBaseArgs):
        plugin_param: str = Field(..., description="Plugin parameter")

    class ExternalPluginDataset(DatasetDescriptor):
        name = "external_plugin"
        description = "External plugin dataset"
        ArgsModel = ExternalPluginArgs

        def build(self, args):
            return f"external_{args.plugin_param}"

    # Just defining the class should register it (via __init_subclass__)
    available = list_datasets()
    assert "external_plugin" in available

    # Should be accessible
    plugin_ds = datasets.external_plugin
    assert isinstance(plugin_ds, ExternalPluginDataset)

    # Cleanup
    unregister("external_plugin")


def test_runtime_dataset_definition():
    """Test defining and using a dataset at runtime."""

    # Define a dataset class at runtime
    class RuntimeArgs(DatasetBaseArgs):
        runtime_value: int = Field(default=100, description="Runtime value")

    class RuntimeDataset(DatasetDescriptor):
        name = "runtime_dataset"
        description = "Runtime defined dataset"
        ArgsModel = RuntimeArgs

        def build(self, args):
            return args.runtime_value * 3

    # Should be registered
    assert "runtime_dataset" in list_datasets()

    # Should be accessible and usable
    runtime_ds = datasets.runtime_dataset
    assert isinstance(runtime_ds, RuntimeDataset)

    # Cleanup
    unregister("runtime_dataset")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
