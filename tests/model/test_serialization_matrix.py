"""Coverage of the legacy protobuf contract, not a claim of lossless exchange.

The property and migration limits are recorded in docs/design/model-contract.md.
Every exported Model subclass must be classified here when the public API grows.
"""

import inspect
import json

import numpy as np
import pytest
from google.protobuf.json_format import Parse

from dtcc_core import model
from dtcc_core.model.model import Model


OBJECTS = (
    "Object", "Building", "BuildingPart", "City", "CityObject", "Terrain",
    "Tree", "RoadNetwork",
)
ROOT_ONLY_OBJECTS = ("SensorCollection", "VehicleCollection", "DeSO")
GEOMETRIES = (
    "Grid", "VolumeGrid", "Mesh", "VolumeMesh", "MultiSurface", "Surface",
    "Point", "PointCloud", "LineString", "MultiLineString",
)
VALUES = ("Bounds", "Transform", "Field", "Raster")
UNSUPPORTED = (
    "Landuse", "BuildingCollection", "CalibrationGrid", "FootprintCollection",
    "TreeCollection", "FieldSlice", "StreamlineCollection", "DatasetCollection",
    "DatasetValue",
)
SUPPORTED = OBJECTS + ROOT_ONLY_OBJECTS + GEOMETRIES + VALUES


def test_every_public_model_has_a_serialization_classification():
    public_models = {
        name for name in model.__all__
        if inspect.isclass(getattr(model, name))
        and issubclass(getattr(model, name), Model)
    }
    assert public_models == set(SUPPORTED + UNSUPPORTED + ("Geometry",))
    assert inspect.isabstract(model.Geometry)


@pytest.mark.parametrize("name", SUPPORTED)
def test_minimal_model_round_trips_bytes_and_protobuf_json(name):
    source = getattr(model, name)()
    message = source.to_proto()
    restored = type(source)()
    restored.from_proto(message.SerializeToString())
    # Compare the defined legacy wire state; dtype/shape and metadata limits
    # cannot be inferred from equality of protobuf messages alone.
    assert json.loads(restored.to_json()) == json.loads(source.to_json())
    json_message = Parse(source.to_json(), type(message)())
    restored.from_proto(json_message.SerializeToString())
    assert json.loads(restored.to_json()) == json.loads(source.to_json())


@pytest.mark.parametrize("name", UNSUPPORTED)
def test_unsupported_model_serialization_fails_explicitly(name):
    value = getattr(model, name)()
    for operation in (value.to_proto, value.to_json, lambda: value.from_proto(b"")):
        with pytest.raises(NotImplementedError, match="[Pp]rotobuf"):
            operation()


@pytest.mark.parametrize("name", OBJECTS)
def test_nested_objects_preserve_type_attributes_and_children(name):
    root = model.City(id="city")
    child = getattr(model, name)(id="feature")
    child.attributes = {"name": "A", "optional": None, "tags": ["urban", 2, True]}
    part = model.BuildingPart(id="part")
    part.add_geometry(model.Point(x=2, y=3, z=4), "location")
    child.add_child(part)
    root.add_child(child)
    payload = root.to_proto().SerializeToString()
    restored = model.City()
    # Reuse is replacement, including nested content.
    restored.from_proto(payload)
    restored.from_proto(payload)
    assert list(restored.children) == [type(child)]
    assert len(restored.children[type(child)]) == 1
    result = restored.children[type(child)][0]
    assert type(result) is type(child)
    assert result.id == child.id
    assert result.attributes == child.attributes
    assert len(result.children[model.BuildingPart]) == 1
    assert result.children[model.BuildingPart][0].geometry["location"].z == 4


@pytest.mark.parametrize("name", ROOT_ONLY_OBJECTS)
def test_collections_cannot_silently_lose_type_when_nested(name):
    root = model.City()
    root.add_child(getattr(model, name)())
    with pytest.raises(NotImplementedError, match=name):
        root.to_proto()


def _representative_geometry(name):
    vertices = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.]])
    constructors = {
        "Grid": lambda: model.Grid(width=2, height=3),
        "VolumeGrid": lambda: model.VolumeGrid(width=2, height=3, depth=4),
        "Mesh": lambda: model.Mesh(vertices=vertices[:3], faces=np.array([[0, 1, 2]])),
        "VolumeMesh": lambda: model.VolumeMesh(vertices=vertices, cells=np.array([[0, 1, 2, 3]])),
        "Surface": lambda: model.Surface(vertices=vertices[:3]),
        "MultiSurface": lambda: model.MultiSurface(surfaces=[model.Surface(vertices=vertices[:3])]),
        "LineString": lambda: model.LineString(vertices=vertices[:2]),
        "MultiLineString": lambda: model.MultiLineString(linestrings=[model.LineString(vertices=vertices[:2])]),
        "PointCloud": lambda: model.PointCloud(points=vertices),
        "Point": lambda: model.Point(x=2, y=3, z=4),
    }
    geometry = constructors[name]()
    geometry.transform.srs = "EPSG:3006"
    geometry.transform.set_translation(10, 20, 30)
    geometry.add_field(model.Field(name="temperature", unit="K", values=np.array([273.5])))
    return geometry


@pytest.mark.parametrize("name", GEOMETRIES)
def test_nested_geometry_preserves_shape_payload_and_base_metadata(name):
    geometry = _representative_geometry(name)
    root = model.Object()
    root.add_geometry(geometry, "representation")
    restored = model.Object()
    payload = root.to_proto().SerializeToString()
    restored.from_proto(payload)
    restored.from_proto(payload)
    result = restored.geometry["representation"]
    assert type(result) is type(geometry)
    assert result.transform.srs == "EPSG:3006"
    np.testing.assert_array_equal(result.transform.affine, geometry.transform.affine)
    assert len(result.fields) == 1
    assert result.fields[0].name == "temperature"
    assert result.fields[0].unit == "K"
    np.testing.assert_array_equal(result.fields[0].values, [[273.5]])
    assert json.loads(result.to_json()) == json.loads(geometry.to_json())


@pytest.mark.parametrize("name", ("MultiSurface", "MultiLineString"))
def test_multi_geometry_decode_replaces_shape_children(name):
    source = _representative_geometry(name)
    payload = source.to_proto().SerializeToString()
    source.from_proto(payload)
    source.from_proto(payload)
    children = source.surfaces if name == "MultiSurface" else source.linestrings
    assert len(children) == 1


def test_nested_geometry_subclass_cannot_silently_lose_type():
    class CustomPoint(model.Point):
        pass

    root = model.Object()
    root.add_geometry(CustomPoint(), "location")
    with pytest.raises(NotImplementedError, match="CustomPoint.*restore Point"):
        root.to_proto()
