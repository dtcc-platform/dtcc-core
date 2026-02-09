import pytest
import logging

try:
    from pydantic import BaseModel, Field
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


def test_duplicate_registration_warns(caplog):
    """Test that duplicate registration logs a warning."""
    # Register first time
    register_class("dup_test", BaseTestDataset)

    # Register again with same name
    with caplog.at_level(logging.WARNING):
        register_class("dup_test", AnotherTestDataset)

    # Should have logged a warning
    assert "already registered" in caplog.text.lower()

    # Latest registration should win
    available = list_datasets()
    assert isinstance(available["dup_test"], AnotherTestDataset)

    # Cleanup
    unregister("dup_test")


# Discovery Tests


def test_list_returns_all_datasets():
    """Test that list() returns all registered datasets."""
    available = list_datasets()

    # Should include built-in datasets
    assert "point_cloud" in available
    assert "buildings" in available
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
