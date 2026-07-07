"""Test SensorCollection protobuf roundtrip."""

import pytest
import numpy as np

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg", force=True)

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context
from dtcc_core.model.object import SensorCollection, Object
from dtcc_core.model.geometry import Point
from dtcc_core.model.values import Field


def test_sensor_collection_creation():
    """Test basic SensorCollection creation."""
    sc = SensorCollection()
    assert len(sc.stations()) == 0


def test_sensor_collection_add_station():
    """Test adding stations to collection."""
    sc = SensorCollection()
    
    # Create station
    station = Object()
    station.attributes = {"id": "station1", "name": "Test Station"}
    
    # Add geometry
    point = Point(x=10.0, y=20.0, z=0.0)
    station.geometry["location"] = point
    
    # Add station
    sc.add_station(station)
    
    assert len(sc.stations()) == 1
    assert sc.stations()[0] == station


def test_sensor_collection_with_field():
    """Test SensorCollection with field values."""
    sc = SensorCollection()
    
    # Create station with field
    station = Object()
    station.attributes = {"id": "s1", "value": 42.5}
    
    point = Point(x=100.0, y=200.0, z=0.0)
    
    field = Field()
    field.name = "NO2"
    field.unit = "µg/m³"
    field.dim = 1
    field.values = np.array([42.5], dtype=np.float32)
    
    point.fields = [field]
    station.geometry["location"] = point
    
    sc.add_station(station)
    
    # Test to_arrays
    points, values = sc.to_arrays("NO2")
    
    assert points.shape == (1, 3)
    assert values.shape == (1,)
    assert values[0] == pytest.approx(42.5)


def test_sensor_collection_plot_returns_axes():
    sc = SensorCollection()
    station = Object()
    point = Point(x=100.0, y=200.0, z=0.0)
    field = Field()
    field.name = "NO2"
    field.unit = "ug/m3"
    field.dim = 1
    field.values = np.array([42.5], dtype=np.float32)
    point.fields = [field]
    station.geometry["location"] = point
    sc.add_station(station)

    ax = sc.plot("NO2", show=False)

    assert ax.get_title()


def test_sensor_collection_plot_uses_presentation_panel_by_default():
    sc = SensorCollection()
    station = Object()
    point = Point(x=100.0, y=200.0, z=0.0)
    field = Field()
    field.name = "air_temperature"
    field.unit = "celsius"
    field.dim = 1
    field.values = np.array([12.5], dtype=np.float32)
    point.fields = [field]
    station.geometry["location"] = point
    sc.add_station(station)
    context = datasets.weather.create_context(
        datasets.weather.validate(
            {"bounds": (0.0, 0.0, 1.0, 1.0), "parameters": ["temperature"]}
        )
    )
    attach_dataset_context(sc, context)

    ax = sc.plot("air_temperature", show=False)
    simple_ax = sc.plot("air_temperature", show=False, presentation=False)

    assert len(ax.figure.axes) >= 2
    assert ax.get_title() == ""
    assert ax.figure.get_figwidth() > simple_ax.figure.get_figwidth()
    assert simple_ax.get_title()


def test_sensor_collection_protobuf_roundtrip():
    """Test SensorCollection protobuf roundtrip with multiple stations."""
    # Create collection
    sc1 = SensorCollection()
    sc1.attributes = {"source": "test", "phenomenon": "PM10"}
    
    # Add two stations
    for i in range(2):
        station = Object()
        station.attributes = {
            "station_id": f"s{i}",
            "value": float(i * 10),
        }
        
        point = Point(x=float(i), y=float(i * 2), z=0.0)
        
        field = Field()
        field.name = "PM10"
        field.unit = "µg/m³"
        field.dim = 1
        field.values = np.array([float(i * 10)], dtype=np.float32)
        
        point.fields = [field]
        station.geometry["location"] = point
        
        sc1.add_station(station)
    
    # Serialize to protobuf
    pb = sc1.to_proto()
    
    # Deserialize
    sc2 = SensorCollection()
    sc2.from_proto(pb)
    
    # Verify structure
    assert len(sc2.stations()) == 2
    assert sc2.attributes.get("phenomenon") == "PM10"
    
    # Verify stations
    for i, station in enumerate(sc2.stations()):
        assert station.attributes["station_id"] == f"s{i}"
        
        # Check geometry
        assert "location" in station.geometry
        point = station.geometry["location"]
        assert point.x == float(i)
        assert point.y == float(i * 2)
        
        # Check field
        assert len(point.fields) == 1
        field = point.fields[0]
        assert field.name == "PM10"
        assert np.ravel(field.values)[0] == pytest.approx(float(i * 10))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
