"""Coverage of the sole DTCC Protobuf contract.

The property and migration limits are recorded in docs/design/model-contract.md.
Every exported Model subclass must be classified here when the public API grows.
"""

import inspect
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pytest
from google.protobuf.json_format import Parse

from dtcc_core import model
from dtcc_core.model.model import Model
from dtcc_core.model.values import raster_tiles

OBJECTS = (
    "Object",
    "Building",
    "BuildingPart",
    "City",
    "CityObject",
    "Terrain",
    "Tree",
    "RoadNetwork",
    "Landuse",
    "SensorCollection",
    "VehicleCollection",
    "DeSO",
)
GEOMETRIES = (
    "Grid",
    "VolumeGrid",
    "Mesh",
    "VolumeMesh",
    "MultiSurface",
    "Surface",
    "Point",
    "PointCloud",
    "LineString",
    "MultiLineString",
    "Solid",
)
VALUES = ("Bounds", "Transform", "Field", "Raster")
# Root-only Dataset wrappers exchanged through a transient Object.
ROOT_WRAPPERS = (
    "BuildingCollection",
    "CalibrationGrid",
    "FieldSlice",
    "FootprintCollection",
    "RasterTileCollection",
    "StreamlineCollection",
    "TreeCollection",
)
# Root wrappers that dtcc_core.model does not export; built from their own module.
UNEXPORTED_WRAPPERS = ("RasterTileCollection",)
UNSUPPORTED = (
    "DatasetCollection",
    "DatasetValue",
    "TimeSeries",
)
SUPPORTED = OBJECTS + GEOMETRIES + VALUES


def test_every_public_model_has_a_serialization_classification():
    public_models = {
        name
        for name in model.__all__
        if inspect.isclass(getattr(model, name))
        and issubclass(getattr(model, name), Model)
    }
    assert public_models == set(
        SUPPORTED + ROOT_WRAPPERS + UNSUPPORTED + ("Geometry",)
    ) - set(UNEXPORTED_WRAPPERS)
    assert inspect.isabstract(model.Geometry)


@pytest.mark.parametrize("name", SUPPORTED)
def test_minimal_model_round_trips_bytes_and_protobuf_json(name):
    source = getattr(model, name)()
    if name == "Field":
        source.association = "sample"
    if name == "Solid":
        source = model.Solid(
            surfaces=[
                model.Surface(
                    vertices=np.array(
                        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
                    )
                )
            ],
            shells=[np.array([0])],
        )
    message = source.to_proto()
    restored = type(source)()
    restored.from_proto(message.SerializeToString())
    # Compare the defined legacy wire state; dtype/shape and metadata limits
    # cannot be inferred from equality of protobuf messages alone.
    assert json.loads(restored.to_json()) == json.loads(source.to_json())
    json_message = Parse(source.to_json(), type(message)())
    restored.from_proto(json_message.SerializeToString())
    assert json.loads(restored.to_json()) == json.loads(source.to_json())


def _populated_wrapper(name):
    square = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0]])
    if name == "TreeCollection":
        moved = model.Tree(id="moved", position=np.array([1.0, 2.0, 0.0]), height=6.0)
        moved.transform.set_translation(10, 20, 0)
        return model.TreeCollection(trees=[model.Tree(id="tree"), moved])
    if name == "BuildingCollection":
        building = model.Building(id="building", attributes={"height": 12.5})
        building.add_geometry(model.Surface(vertices=square), "lod0")
        return model.BuildingCollection(buildings=[building, model.Building(id="b2")])
    if name == "FootprintCollection":
        return model.FootprintCollection(
            footprints=[model.Surface(vertices=square)] * 2,
            source_ids=["building", None],
            source_indices=[0, 3],
        )
    if name == "FieldSlice":
        source = model.FieldSlice(
            points=square,
            slice_position=0.5,
            resolution=3,
            domain_bounds=model.Bounds(xmin=0, ymin=0, xmax=4, ymax=4),
            metadata_payload={"source": "matrix"},
        )
        source.transform.set_translation(10, 20, 30)
        source.add_field(
            model.Field(
                name="speed", unit="m/s", values=np.ones(3), association="sample"
            )
        )
        return source
    if name == "StreamlineCollection":
        line = model.LineString(vertices=square)
        line.add_field(
            model.Field(
                name="speed", unit="m/s", values=np.ones(3), association="vertex"
            )
        )
        return model.StreamlineCollection(
            lines=[line], requested_line_count=1, streamline_steps=2
        )
    if name == "RasterTileCollection":
        tile = raster_tiles.RasterTile(
            path=Path("cache/tile.tif"),
            id="tile",
            collection="orto",
            datetime=datetime(2025, 5, 31, tzinfo=UTC),
            extent=model.Bounds(xmin=0, ymin=0, xmax=4, ymax=4),
            crs="EPSG:3006",
            spektraltyp="rgb",
            resolution=0.5,
            size_bytes=None,
        )
        return raster_tiles.RasterTileCollection(tiles=[tile, tile])
    return model.CalibrationGrid(
        bounds=model.Bounds(xmin=0, ymin=0, xmax=4, ymax=4),
        divisions=2,
        crs="EPSG:3006",
        features=[{"type": "Feature", "properties": {"index": 0}}],
        metadata_payload={"divisions": 2},
    )


@pytest.mark.parametrize("name", ROOT_WRAPPERS)
def test_populated_root_wrapper_round_trips_bytes_and_protobuf_json(name):
    source = _populated_wrapper(name)
    message = source.to_proto()
    restored = type(source)()
    restored.from_proto(message.SerializeToString())
    assert type(restored) is type(source)
    assert json.loads(restored.to_json()) == json.loads(source.to_json())
    json_message = Parse(source.to_json(), type(message)())
    restored.from_proto(json_message.SerializeToString())
    assert json.loads(restored.to_json()) == json.loads(source.to_json())
    assert restored.to_proto().SerializeToString(
        deterministic=True
    ) == message.SerializeToString(deterministic=True)


@pytest.mark.parametrize("name", UNSUPPORTED)
def test_unsupported_model_serialization_fails_explicitly(name):
    value = getattr(model, name)()
    for operation in (value.to_proto, value.to_json, lambda: value.from_proto(b"")):
        with pytest.raises(NotImplementedError, match=name):
            operation()


@pytest.mark.parametrize("name", OBJECTS)
def test_nested_objects_preserve_type_attributes_and_children(name):
    root = model.City(id="city")
    child = getattr(model, name)(id="feature")
    child.attributes = {"name": "A", "optional": None, "tags": ["urban", 2, True]}
    part = model.Object(id="part")
    part.add_geometry(model.Point(x=2, y=3, z=4), "location")
    child.add_child(part)
    root.add_child(child)
    payload = root.to_proto(validate_schema=False).SerializeToString()
    restored = model.City()
    # Reuse is replacement, including nested content.
    restored.from_proto(payload, validate_schema=False)
    restored.from_proto(payload, validate_schema=False)
    assert list(restored.children) == [type(child)]
    assert len(restored.children[type(child)]) == 1
    result = restored.children[type(child)][0]
    assert type(result) is type(child)
    assert result.id == child.id
    assert result.attributes == child.attributes
    assert len(result.children[model.Object]) == 1
    assert result.children[model.Object][0].get_geometry("location").z == 4


def _representative_geometry(name):
    vertices = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    )
    constructors = {
        "Grid": lambda: model.Grid(width=2, height=3),
        "VolumeGrid": lambda: model.VolumeGrid(width=2, height=3, depth=4),
        "Mesh": lambda: model.Mesh(vertices=vertices[:3], faces=np.array([[0, 1, 2]])),
        "VolumeMesh": lambda: model.VolumeMesh(
            vertices=vertices, cells=np.array([[0, 1, 2, 3]])
        ),
        "Surface": lambda: model.Surface(vertices=vertices[:3]),
        "MultiSurface": lambda: model.MultiSurface(
            surfaces=[model.Surface(vertices=vertices[:3])]
        ),
        "LineString": lambda: model.LineString(vertices=vertices[:2]),
        "MultiLineString": lambda: model.MultiLineString(
            linestrings=[model.LineString(vertices=vertices[:2])]
        ),
        "PointCloud": lambda: model.PointCloud(points=vertices),
        "Point": lambda: model.Point(x=2, y=3, z=4),
        "Solid": lambda: model.Solid(
            surfaces=[model.Surface(vertices=vertices[:3])], shells=[np.array([0])]
        ),
    }
    geometry = constructors[name]()
    geometry.transform.srs = "EPSG:3006"
    geometry.transform.set_translation(10, 20, 30)
    geometry.add_field(
        model.Field(
            name="temperature",
            unit="K",
            values=np.array([273.5]),
            association="geometry",
        )
    )
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
    result = restored.get_geometry("representation")
    assert type(result) is type(geometry)
    assert result.transform.srs == "EPSG:3006"
    np.testing.assert_array_equal(result.transform.affine, geometry.transform.affine)
    assert len(result.fields) == 1
    assert result.fields[0].name == "temperature"
    assert result.fields[0].unit == "K"
    np.testing.assert_array_equal(result.fields[0].values, [273.5])
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
    with pytest.raises(NotImplementedError, match="CustomPoint"):
        root.to_proto()
