"""Test VehicleCollection protobuf roundtrip."""

import numpy as np
import pytest

from dtcc_core.model.object import Object, VehicleCollection
from dtcc_core.model.geometry import Point
from dtcc_core.model.values import Field


def _vehicle(vehicle_id="bus-1", x=10.0, y=20.0, mode="bus"):
    vehicle = Object()
    vehicle.attributes = {
        "vehicle_id": vehicle_id,
        "mode": mode,
        "provider": "test",
        "speed": 8.5,
    }
    point = Point(x=x, y=y, z=0.0)
    field = Field()
    field.name = "speed"
    field.unit = "m/s"
    field.dim = 1
    field.values = np.array([8.5], dtype=np.float32)
    point.fields = [field]
    vehicle.geometry["location"] = point
    return vehicle


def test_vehicle_collection_creation():
    vehicles = VehicleCollection()
    assert vehicles.vehicles() == []


def test_vehicle_collection_add_vehicle_and_arrays():
    vehicles = VehicleCollection()
    vehicles.add_vehicle(_vehicle())

    points, speeds = vehicles.to_arrays("speed")

    assert points.shape == (1, 3)
    assert speeds.shape == (1,)
    assert speeds[0] == pytest.approx(8.5)


def test_vehicle_collection_attribute_arrays():
    vehicles = VehicleCollection()
    vehicles.add_vehicle(_vehicle(mode="tram"))

    points, modes = vehicles.to_arrays("mode")

    assert points.shape == (1, 3)
    assert modes.tolist() == ["tram"]


def test_vehicle_collection_protobuf_roundtrip():
    first = VehicleCollection()
    first.attributes = {"source": "test", "provider": "test"}
    first.add_vehicle(_vehicle("bus-1", 10.0, 20.0))
    first.add_vehicle(_vehicle("bus-2", 11.0, 21.0))

    payload = first.to_proto()

    second = VehicleCollection()
    second.from_proto(payload)

    assert second.attributes["provider"] == "test"
    assert len(second.vehicles()) == 2
    for vehicle in second.vehicles():
        assert vehicle.attributes["mode"] == "bus"
        assert "location" in vehicle.geometry
        point = vehicle.geometry["location"]
        assert len(point.fields) == 1
        assert point.fields[0].name == "speed"
