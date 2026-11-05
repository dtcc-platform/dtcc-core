"""Regression coverage for the public Object and Tree protobuf contract."""

import math

import numpy as np
import pytest

from dtcc_core.model import (
    Bounds,
    Building,
    City,
    DeSO,
    Field,
    GeometryType,
    Grid,
    LineString,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    Raster,
    SensorCollection,
    Tree,
    VehicleCollection,
)
from dtcc_core.model import dtcc_pb2 as proto


@pytest.mark.parametrize("as_bytes", [False, True])
def test_tree_nested_protobuf_roundtrip(as_bytes):
    tree = Tree(
        id="oak",
        attributes={"species": "oak", "survey": {"year": 2026}},
        position=np.array([1.25, 2.5, 3.75]),
        height=8.5,
        crown_radius=2.25,
    )
    tree.add_geometry(Point(x=1.25, y=2.5, z=3.75), "location")
    tree.add_child(Object(id="survey"))
    city = City()
    city.add_trees([tree])
    payload = city.to_proto()
    restored = City()
    restored.from_proto(payload.SerializeToString() if as_bytes else payload)

    assert len(restored.trees) == 1
    result = restored.trees[0]
    assert type(result) is Tree
    assert result.id == tree.id
    assert result.attributes == tree.attributes
    np.testing.assert_array_equal(result.position, tree.position)
    assert result.height == tree.height
    assert result.crown_radius == tree.crown_radius
    assert result.get_children(Object)[0].id == "survey"
    assert result.get_geometry("location").x == 1.25


def test_empty_tree_roundtrip():
    restored = Tree()
    restored.from_proto(Tree().to_proto())
    assert restored.position.shape == (0, 3)
    assert restored.height == 0.0
    assert restored.crown_radius == 0.0


@pytest.mark.parametrize(
    "kwargs, field",
    [
        ({"position": [1, 2]}, "position"),
        ({"position": [[1, 2, 3]]}, "position"),
        ({"position": [[1], [2, 3]]}, "position"),
        ({"position": [1, 2, math.nan]}, "position"),
        ({"position": [1, 2, 1e40]}, "position"),
        ({"position": ["1", "2", "3"]}, "position"),
        ({"height": -1}, "height"),
        ({"height": math.inf}, "height"),
        ({"height": 10**400}, "height"),
        ({"crown_radius": True}, "crown_radius"),
        ({"crown_radius": math.nan}, "crown_radius"),
    ],
)
def test_tree_rejects_invalid_values(kwargs, field):
    with pytest.raises(ValueError, match=f"Tree.{field}"):
        Tree(**kwargs).to_proto()


def test_tree_rejects_invalid_payload_without_replacing_state():
    tree = Tree(id="original")
    malformed = Tree(id="replacement", position=np.array([1.0, 2.0, 3.0])).to_proto()
    malformed.object.tree.position.shape[:] = [1, 3]
    with pytest.raises(ValueError, match="Tree.position"):
        tree.from_proto(malformed)
    assert tree.id == "original"
    with pytest.raises(ValueError, match="Expected Tree"):
        tree.from_proto(Object().to_proto())


@pytest.mark.parametrize("collection_type", [SensorCollection, VehicleCollection, DeSO])
def test_top_level_collection_still_roundtrips_when_type_is_known(collection_type):
    collection = collection_type(attributes={"source": "test"})
    collection.add_child(Object(id="entry"))
    restored = collection_type()
    restored.from_proto(collection.to_proto())
    assert restored.attributes == {"source": "test"}
    assert restored.get_children(Object)[0].id == "entry"


@pytest.mark.parametrize("base_type", [Building, Tree])
def test_nested_inherited_concrete_serializer_cannot_erase_subclass(base_type):
    class CustomObject(base_type):
        pass

    parent = Object()
    parent.add_child(CustomObject())
    with pytest.raises(
        NotImplementedError,
        match="CustomObject",
    ):
        parent.to_proto()


@pytest.mark.parametrize("as_bytes", [False, True])
def test_object_from_proto_replaces_children_geometry_attributes_and_bounds(as_bytes):
    source = Object(id="source", attributes={"active": True})
    source.add_child(Building(id="building"))
    source.add_geometry(Point(x=3, y=4, z=5), "location")
    payload = source.to_proto()
    if as_bytes:
        payload = payload.SerializeToString()

    restored = Object(attributes={"stale": True})
    restored.add_child(Tree())
    restored.add_geometry(Point(), "old")
    restored.from_proto(payload)
    restored.from_proto(payload)
    assert restored.id == "source"
    assert restored.attributes == {"active": True}
    assert list(restored.children) == [Building]
    assert len(restored.get_children(Building)) == 1
    assert list(restored.geometry) == ["location"]

    restored.bounds = Bounds(xmin=1, ymin=2, xmax=3, ymax=4)
    restored.from_proto(Object(id="empty").to_proto())
    assert restored.children == {}
    assert restored.geometry == {}
    assert restored.attributes == {}
    assert restored._bounds is None


@pytest.mark.parametrize(
    "name, expected",
    [
        ("MultiSurface", GeometryType.MULTISURFACE),
        ("multi_surface", GeometryType.MULTISURFACE),
        ("LineString", GeometryType.LINESTRING),
        ("multi_line_string", GeometryType.MULTILINESTRING),
        ("pointcloud", GeometryType.POINT_CLOUD),
        ("volume_mesh", GeometryType.VOLUME_MESH),
    ],
)
def test_geometry_names_normalize_to_existing_enum_members(name, expected):
    assert GeometryType.from_str(name) is expected
    assert GeometryType.from_class_name(name) is expected


@pytest.mark.parametrize("geometry_type", [MultiSurface, LineString, MultiLineString])
def test_geometry_class_name_can_infer_existing_representation(geometry_type):
    obj = Object()
    geometry = geometry_type()
    obj.add_geometry(geometry)
    assert obj.get_geometry(GeometryType.from_class(geometry_type)) is geometry


@pytest.mark.parametrize("geometry", [Point(), Grid()])
def test_geometry_without_builtin_representation_requires_explicit_role(geometry):
    with pytest.raises(ValueError, match="explicit geometry_type string role"):
        Object().add_geometry(geometry)


@pytest.mark.parametrize("key", ["location", "centerlines", "grid", "volume_grid"])
def test_custom_geometry_roles_are_preserved(key):
    obj = Object()
    obj.add_geometry(Point(x=1, y=2, z=3), key)
    restored = Object()
    restored.from_proto(obj.to_proto())
    assert list(restored.geometry) == [key]
    restored.remove_geometry(key)
    assert restored.geometry == {}


@pytest.mark.parametrize("geometry", [None, 12, object()])
def test_invalid_geometry_values_fail_before_mutating_populated_object(geometry):
    obj = Object(attributes={"source": "test"})
    point = Point(x=1, y=2, z=3)
    obj.add_geometry(point, "location")
    with pytest.raises(TypeError, match="[Gg]eometry"):
        obj.add_geometry(geometry, "invalid")
    assert list(obj.geometry) == ["location"]
    assert obj.get_geometry("location") is point
    assert obj.attributes == {"source": "test"}

    obj.geometry["invalid"] = geometry
    with pytest.raises(ValueError, match="GeometryRepresentation"):
        obj.to_proto()


@pytest.mark.parametrize(
    "geometry",
    [Raster(data=np.array([[1.0, 2.0]])), Bounds(xmin=1, ymin=2, xmax=3, ymax=4)],
)
def test_raster_and_bounds_representations_roundtrip(geometry):
    obj = Object()
    obj.add_geometry(geometry)
    assert obj.get_geometry(GeometryType.from_class(type(geometry))) is geometry
    restored = Object()
    restored.from_proto(obj.to_proto())
    result = restored.get_geometries()[0]
    assert type(result) is type(geometry)
    if isinstance(geometry, Raster):
        np.testing.assert_array_equal(result.data, geometry.data)
    else:
        assert result == geometry


@pytest.mark.parametrize("payload", [None, proto.Point(), proto.Geometry()])
def test_geometry_serialization_requires_concrete_geometry_message(payload):
    class InvalidGeometry(Point):
        def to_proto(self):
            return payload

    obj = Object()
    obj.add_geometry(InvalidGeometry(), "location")
    with pytest.raises(NotImplementedError, match="InvalidGeometry"):
        obj.to_proto()


@pytest.mark.parametrize("key", ["mesh", "MESH", "GeometryType.MESH"])
def test_add_remove_and_field_use_same_geometry_key_normalization(key):
    obj = Object()
    geometry = Point()
    obj.add_geometry(geometry, key)
    assert obj.get_geometry(GeometryType.MESH) is geometry
    field = Field(name="test")
    obj.add_field(field, key)
    assert geometry.fields == [field]
    obj.remove_geometry(key)
    assert obj.geometry == {}


@pytest.mark.parametrize(
    "key", [1, None, (), "", "  ", "GeometryType.missing", "GeometryTypeMesh"]
)
def test_invalid_geometry_keys_fail(key):
    obj = Object()
    # Direct dictionary mutation also receives validation at the wire boundary.
    obj.geometry[key] = Point()
    with pytest.raises((TypeError, ValueError), match="[Gg]eometry|Representation"):
        obj.to_proto()
    with pytest.raises((TypeError, ValueError), match="[Gg]eometry"):
        obj.remove_geometry(key)


def test_defined_geometries_handles_enums_and_custom_roles():
    obj = Object()
    keys = [GeometryType.MESH, "location", GeometryType.LOD0]
    for key in keys:
        obj.add_geometry(Point(), key)
    assert obj.defined_geometries() == sorted(map(str, keys))


def test_missing_geometry_discriminator_has_actionable_error():
    payload = Object().to_proto()
    payload.object.representations.add(id="location").geometry.SetInParent()
    with pytest.raises(NotImplementedError, match="geometry kind"):
        Object().from_proto(payload)


def test_json_attributes_roundtrip():
    attributes = {
        "nested": [None, True, 4, 1.5, "oak", {"labels": ["a", "b"]}],
        "height": np.float64(8.5),
    }
    restored = Object()
    restored.from_proto(Object(attributes=attributes).to_proto())
    assert restored.attributes == attributes


@pytest.mark.parametrize(
    "value",
    [(1, 2), {1, 2}, np.array([1, 2]), np.int64(3), math.inf, math.nan],
)
def test_invalid_attribute_value_reports_its_path(value):
    obj = Object(attributes={"survey": [{"value": value}]})
    with pytest.raises(
        (TypeError, ValueError), match=r"attributes\['survey'\]\[0\]\['value'\]"
    ):
        obj.to_proto()


def test_attribute_keys_and_cycles_fail_explicitly():
    with pytest.raises(TypeError, match="non-string dictionary key"):
        Object(attributes={"nested": {1: "value"}}).to_proto()
    attributes = {}
    attributes["self"] = attributes
    with pytest.raises(ValueError, match="circular reference"):
        Object(attributes=attributes).to_proto()


def test_city_collects_all_building_attributes_with_aligned_missing_values():
    city = City()
    city.add_buildings(
        [
            Building(attributes={"name": "A", "height": 10}),
            Building(attributes={"height": 20, "year": 2026}),
            Building(),
        ]
    )
    assert city.get_building_attributes() == {
        "name": ["A", None, None],
        "height": [10, 20, None],
        "year": [None, 2026, None],
    }
    assert list(city.get_building_attributes()) == ["name", "height", "year"]
    assert City().get_building_attributes() == {}
