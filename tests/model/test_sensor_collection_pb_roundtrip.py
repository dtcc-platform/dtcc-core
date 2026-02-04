"""Test SensorCollection protobuf roundtrip."""

import pytest
import numpy as np

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


@pytest.mark.skip(reason="SensorCollection protobuf serialization not yet implemented")
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
        assert field.values[0] == pytest.approx(float(i * 10))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
