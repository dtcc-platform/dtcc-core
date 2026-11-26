"""Root-only Dataset wrappers through every canonical exchange boundary."""

from importlib.resources import files

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (
    DatasetContext,
    DatasetIdentity,
    DatasetMetadata,
    DatasetPresentation,
    DatasetProvenance,
    DatasetRequest,
)
from dtcc_core.model import (
    Bounds,
    Building,
    BuildingCollection,
    BuildingPart,
    CalibrationGrid,
    City,
    Field,
    FootprintCollection,
    GeometryType,
    Object,
    Surface,
    Tree,
    TreeCollection,
    exchange,
)
from dtcc_core.model import dtcc_pb2 as wire
from dtcc_core.model.profiles import SemanticProfile

WRAPPERS = (TreeCollection, BuildingCollection, FootprintCollection, CalibrationGrid)
KINDS = {
    TreeCollection: wire.Object.TREE_COLLECTION,
    BuildingCollection: wire.Object.BUILDING_COLLECTION,
    FootprintCollection: wire.Object.FOOTPRINT_COLLECTION,
    CalibrationGrid: wire.Object.CALIBRATION_GRID,
}
SCHEMA = [True, False]


def square(x=0.0, size=10.0):
    return np.array(
        [[x, 0.0, 0.0], [x + size, 0.0, 0.0], [x + size, size, 0.0], [x, size, 0.0]]
    )


def populated(cls):
    """Two members, one carrying its own non-identity transform."""
    if cls is TreeCollection:
        first = Tree(
            id="tree-1",
            position=np.array([1.0, 2.0, 0.0]),
            height=8.5,
            crown_radius=2.0,
            attributes={"species": "oak"},
        )
        second = Tree(id="tree-2", position=np.array([5.0, 6.0, 0.0]), height=4.0)
        second.transform.set_translation(100.0, 200.0, 0.0)
        return TreeCollection(trees=[first, second])
    if cls is BuildingCollection:
        first = Building(id="building-1", attributes={"height": 12.5, "levels": 4})
        first.add_geometry(Surface(vertices=square()), GeometryType.LOD0)
        second = Building(id="building-2")
        second.add_geometry(Surface(vertices=square(20.0)), GeometryType.LOD0)
        second.transform.srs = "EPSG:3006"
        second.transform.set_translation(319000.0, 6399000.0, 0.0)
        return BuildingCollection(buildings=[first, second])
    if cls is FootprintCollection:
        holed = Surface(vertices=square(), holes=[square(2.0, 2.0)])
        moved = Surface(vertices=square(20.0))
        moved.transform.set_translation(5.0, 5.0, 0.0)
        moved.add_field(
            Field(name="area", unit="m2", values=np.array([100.0]), association="face")
        )
        return FootprintCollection(
            footprints=[holed, moved],
            source_ids=["building-1", None],
            source_indices=[0, 2],
        )
    return CalibrationGrid(
        bounds=Bounds(xmin=0.0, ymin=0.0, xmax=100.0, ymax=50.0),
        divisions=2,
        crs="EPSG:3006",
        features=[
            {
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": [[0, 0], [0, 50]]},
                "properties": {"axis": "x", "index": 0},
            },
            {
                "type": "Feature",
                "geometry": {
                    "type": "LineString",
                    "coordinates": [[0.5, 0], [0.5, 50]],
                },
                "properties": {"axis": "x", "index": 1},
            },
        ],
        name="calibration_grid",
        metadata_payload={"dataset": "custom", "spacing": [50.0, 25.0], "note": None},
    )


def context():
    return DatasetContext(
        identity=DatasetIdentity(name="example", title="Example"),
        metadata=DatasetMetadata(crs=["EPSG:3006"]),
        provenance=DatasetProvenance(sources=["example"]),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name="example", bounds=[0.0, 0.0, 1.0, 1.0]),
        health={"status": "ok"},
        warnings=["synthetic"],
    )


def assert_geometry_equal(a, b):
    assert type(a) is type(b)
    assert a.transform.srs == b.transform.srs
    np.testing.assert_array_equal(a.transform.affine, b.transform.affine)
    np.testing.assert_array_equal(a.vertices, b.vertices)
    assert len(a.holes) == len(b.holes)
    for x, y in zip(a.holes, b.holes, strict=True):
        np.testing.assert_array_equal(x, y)
    assert [(f.name, f.unit, f.association) for f in a.fields] == [
        (f.name, f.unit, f.association) for f in b.fields
    ]
    for x, y in zip(a.fields, b.fields, strict=True):
        np.testing.assert_array_equal(x.values, y.values)


def assert_object_equal(a, b):
    assert type(a) is type(b)
    assert a.id == b.id
    assert a.attributes == b.attributes
    assert a.transform.srs == b.transform.srs
    np.testing.assert_array_equal(a.transform.affine, b.transform.affine)
    assert a.geometry.keys() == b.geometry.keys()
    for key in a.geometry:
        x, y = a.geometry[key], b.geometry[key]
        assert (x.lod, x.role) == (y.lod, y.role)
        assert_geometry_equal(x.geometry, y.geometry)
    if type(a) is Tree:
        np.testing.assert_array_equal(a.position, b.position)
        assert (a.height, a.crown_radius) == (b.height, b.crown_radius)


def assert_wrapper_equal(source, restored):
    assert type(restored) is type(source)
    if type(source) is TreeCollection:
        assert len(restored.trees) == len(source.trees)
        for a, b in zip(source.trees, restored.trees, strict=True):
            assert_object_equal(a, b)
    elif type(source) is BuildingCollection:
        assert len(restored.buildings) == len(source.buildings)
        for a, b in zip(source.buildings, restored.buildings, strict=True):
            assert_object_equal(a, b)
    elif type(source) is FootprintCollection:
        assert restored.source_ids == source.source_ids
        assert restored.source_indices == source.source_indices
        assert len(restored.footprints) == len(source.footprints)
        for a, b in zip(source.footprints, restored.footprints, strict=True):
            assert_geometry_equal(a, b)
    else:
        for name in ("xmin", "ymin", "zmin", "xmax", "ymax", "zmax"):
            assert getattr(restored.bounds, name) == getattr(source.bounds, name)
        assert restored.divisions == source.divisions
        assert restored.crs == source.crs
        assert restored.features == source.features
        assert restored.name == source.name
        assert restored.metadata_payload == source.metadata_payload


@pytest.fixture(params=WRAPPERS, ids=lambda cls: cls.__name__)
def wrapper(request):
    return populated(request.param)


@pytest.mark.parametrize("validate_schema", SCHEMA)
def test_populated_wrapper_round_trips_bytes_and_proto(wrapper, validate_schema):
    payload = exchange.dumps(wrapper, validate_schema=validate_schema)
    pb = wire.ModelFile.FromString(payload)
    assert pb.version == 7
    assert pb.object.kind == KINDS[type(wrapper)]
    assert pb.object.id == "collection"

    restored = exchange.loads(payload, validate_schema=validate_schema)
    assert_wrapper_equal(wrapper, restored)
    assert exchange.dumps(restored, validate_schema=validate_schema) == payload

    message = wrapper.to_proto(validate_schema=validate_schema)
    target = type(wrapper)()
    target.from_proto(message, validate_schema=validate_schema)
    assert_wrapper_equal(wrapper, target)
    target.from_proto(message.SerializeToString(), validate_schema=validate_schema)
    assert_wrapper_equal(wrapper, target)
    assert exchange.dumps(target) == payload


@pytest.mark.parametrize("validate_schema", SCHEMA)
def test_populated_wrapper_round_trips_dtcc_file(wrapper, validate_schema, tmp_path):
    path = tmp_path / "wrapper.dtcc"
    io.save_model(wrapper, path, validate_schema=validate_schema)
    restored = io.load_model(
        path, expected_type=type(wrapper), validate_schema=validate_schema
    )
    assert_wrapper_equal(wrapper, restored)
    assert path.read_bytes() == exchange.dumps(restored)


@pytest.mark.parametrize("validate_schema", SCHEMA)
@pytest.mark.parametrize("name", ["package", "package.dtccpkg"])
def test_populated_wrapper_round_trips_canonical_package(
    wrapper, name, validate_schema, tmp_path
):
    wrapper.dataset_context = context()
    package = wrapper.export(
        tmp_path / name, canonical=True, validate_schema=validate_schema
    )
    assert package.artifacts[0].model_type == type(wrapper).__name__
    restored = load_model_package(package.path, validate_schema=validate_schema)
    assert_wrapper_equal(wrapper, restored)
    assert restored.dataset_context == wrapper.dataset_context
    assert exchange.dumps(restored) == exchange.dumps(wrapper)


@pytest.mark.parametrize("cls", WRAPPERS, ids=lambda cls: cls.__name__)
@pytest.mark.parametrize("validate_schema", SCHEMA)
def test_empty_wrapper_round_trips(cls, validate_schema, tmp_path):
    source = cls()
    payload = exchange.dumps(source, validate_schema=validate_schema)
    restored = exchange.loads(payload, validate_schema=validate_schema)
    assert_wrapper_equal(source, restored)
    assert exchange.dumps(restored, validate_schema=validate_schema) == payload
    source.dataset_context = context()
    package = source.export(
        tmp_path / "package.dtccpkg", canonical=True, validate_schema=validate_schema
    )
    restored = load_model_package(package.path, validate_schema=validate_schema)
    assert_wrapper_equal(source, restored)
    assert restored.dataset_context == source.dataset_context


def test_bundled_schema_validates_wrapper_as_transient_object(wrapper):
    profile = SemanticProfile(files("dtcc_core").joinpath("schemas", "dtcc.yaml"))
    assert profile.validate(wrapper).valid
    invalid = TreeCollection(trees=[Building(id="building")])
    with pytest.raises(ValueError, match="children must be Tree"):
        profile.validate(invalid)


# Expected rejection reason per failure, at encode and at the bytes boundary.
REASONS = {
    "building_in_tree_collection": "TreeCollection children must be Tree objects",
    "wrapper_under_city": "kind 13|TreeCollection is unsupported",
    "member_id_collection": "member id must not be 'collection'",
    "member_relation_collection": "member relations must not target 'collection'",
    "nested_relation_collection": "member relations must not target 'collection'",
    "missing_source_indices": "exactly the attributes",
    "non_list_source_indices": "source_indices must be a list",
    "negative_divisions": "int divisions >= 0",
    "negative_stored_divisions": "int divisions >= 0",
    "short_source_ids": "source_ids must be empty or have one entry per footprint",
    "long_source_indices": "source_indices must be empty or have one entry per",
    "negative_source_index": "source_indices must be a list of int >= 0",
    "missing_bounds": "exactly the representations",
    "root_transform": "transform must be identity",
    "root_srs": "transform must be identity",
    "root_relation": "no relations",
    "root_semantic_type": "no semantic type",
    "root_profile": "no semantic type or profile",
    "root_id": "transient Object 'collection'",
    "root_tree_state": "state must match its concrete kind",
    "extra_attribute": "exactly the attributes",
    "representation_lod": "exactly the representations",
    "representation_role": "exactly the representations",
    "extra_representation": "exactly the representations",
    "carrier_field": "carrier MultiSurface",
    "carrier_transform": "carrier MultiSurface",
    "carrier_srs": "carrier MultiSurface",
    "carrier_regions": "carrier MultiSurface",
}


FOOTPRINT_LISTS = {
    "short_source_ids": {"source_ids": ["a"], "source_indices": [0, 1]},
    "long_source_indices": {"source_ids": ["a", None], "source_indices": [0, 1, 2]},
    "negative_source_index": {"source_ids": ["a", None], "source_indices": [0, -1]},
}


def native_failure(failure):
    """Invalid wrapper state constructible through the native classes."""
    if failure == "building_in_tree_collection":
        return TreeCollection(trees=[Tree(id="tree"), Building(id="building")])
    if failure == "member_id_collection":
        return TreeCollection(trees=[Tree(id="collection")])
    if failure == "member_relation_collection":
        return TreeCollection(
            trees=[Tree(id="tree", relations={"member_of": ["collection"]})]
        )
    if failure == "non_list_source_indices":
        return FootprintCollection(
            footprints=[Surface(vertices=square())],
            source_ids=["a"],
            source_indices=(0,),
        )
    if failure == "negative_divisions":
        return CalibrationGrid(divisions=-1)
    if failure == "negative_stored_divisions":
        return CalibrationGrid(divisions=2, metadata_payload={"divisions": -3})
    if failure in FOOTPRINT_LISTS:
        return FootprintCollection(
            footprints=[Surface(vertices=square()), Surface(vertices=square(20.0))],
            **FOOTPRINT_LISTS[failure],
        )
    city = City(id="city")
    city.children[TreeCollection] = [populated(TreeCollection)]
    return city


@pytest.mark.parametrize("validate_schema", SCHEMA)
@pytest.mark.parametrize(
    "failure",
    [
        "building_in_tree_collection",
        "member_id_collection",
        "member_relation_collection",
        "non_list_source_indices",
        "negative_divisions",
        "negative_stored_divisions",
        *FOOTPRINT_LISTS,
        "wrapper_under_city",
    ],
)
def test_native_rejection_precedes_any_file_write(failure, validate_schema, tmp_path):
    model = native_failure(failure)
    path = tmp_path / "existing.dtcc"
    path.write_bytes(b"previous")
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        exchange.dumps(model, validate_schema=validate_schema)
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        io.save_model(model, path, validate_schema=validate_schema)
    assert path.read_bytes() == b"previous"
    model.dataset_context = context()
    target = tmp_path / "package.dtccpkg"
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        model.export(target, canonical=True, validate_schema=validate_schema)
    assert not target.exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["existing.dtcc"]


def carrier_field(root):
    field = Field(name="n", values=np.array([1.0, 2.0]), association="face")
    root.representations[0].geometry.fields.add().CopyFrom(
        exchange._encode_field(field)
    )


def set_affine(transform):
    transform.affine[3] = 1.0


def nested_relation(root):
    part = BuildingPart(id="part", relations={"member_of": ["collection"]})
    exchange._encode_object(part, root.children[1].children.add())


def carrier_region(root):
    root.representations[0].geometry.regions.add(
        semantic_type="https://example.org/Roof",
        indices=exchange._encode_array(np.array([0])),
    )


# Payload edits on an otherwise valid wrapper; reachable only through bytes.
MUTATIONS = {
    "member_id_collection": (
        TreeCollection,
        lambda root: setattr(root.children[1], "id", "collection"),
    ),
    "nested_relation_collection": (BuildingCollection, nested_relation),
    "missing_source_indices": (
        FootprintCollection,
        lambda root: root.attributes.pop("source_indices"),
    ),
    "non_list_source_indices": (
        FootprintCollection,
        lambda root: root.attributes["source_indices"].CopyFrom(wire.Value(text="0")),
    ),
    "negative_divisions": (
        CalibrationGrid,
        lambda root: root.attributes["divisions"].CopyFrom(wire.Value(integer="-1")),
    ),
    "negative_stored_divisions": (
        CalibrationGrid,
        lambda root: (
            root.attributes["metadata"]
            .object.values["divisions"]
            .CopyFrom(wire.Value(integer="-3"))
        ),
    ),
    "short_source_ids": (
        FootprintCollection,
        lambda root: root.attributes["source_ids"].list.values.pop(),
    ),
    "long_source_indices": (
        FootprintCollection,
        lambda root: root.attributes["source_indices"].list.values.add(integer="3"),
    ),
    "negative_source_index": (
        FootprintCollection,
        lambda root: (
            root.attributes["source_indices"]
            .list.values[1]
            .CopyFrom(wire.Value(integer="-1"))
        ),
    ),
    "missing_bounds": (
        CalibrationGrid,
        lambda root: root.ClearField("representations"),
    ),
    "root_transform": (TreeCollection, lambda root: set_affine(root.transform)),
    "root_srs": (
        TreeCollection,
        lambda root: setattr(root.transform, "srs", "EPSG:3006"),
    ),
    "root_relation": (
        TreeCollection,
        lambda root: root.relations["observes"].ids.append("tree-1"),
    ),
    "root_semantic_type": (
        TreeCollection,
        lambda root: setattr(root, "semantic_type", "https://example.org/Collection"),
    ),
    "root_profile": (
        TreeCollection,
        lambda root: (
            setattr(root, "profile_id", "https://example.org/profile"),
            setattr(root, "profile_version", "1"),
        ),
    ),
    "root_id": (TreeCollection, lambda root: setattr(root, "id", "other")),
    "root_tree_state": (TreeCollection, lambda root: setattr(root.tree, "height", 1.0)),
    "extra_attribute": (
        TreeCollection,
        lambda root: root.attributes["extra"].CopyFrom(wire.Value(text="x")),
    ),
    "representation_lod": (
        CalibrationGrid,
        lambda root: setattr(root.representations[0], "lod", "0"),
    ),
    "representation_role": (
        CalibrationGrid,
        lambda root: setattr(root.representations[0], "role", "x"),
    ),
    "extra_representation": (
        CalibrationGrid,
        lambda root: root.representations.add(id="extra").bounds.CopyFrom(
            root.representations[0].bounds
        ),
    ),
    "carrier_field": (FootprintCollection, carrier_field),
    "carrier_transform": (
        FootprintCollection,
        lambda root: set_affine(root.representations[0].geometry.transform),
    ),
    "carrier_srs": (
        FootprintCollection,
        lambda root: setattr(
            root.representations[0].geometry.transform, "srs", "EPSG:3006"
        ),
    ),
    "carrier_regions": (FootprintCollection, carrier_region),
}


def payload_failure(failure):
    """Invalid wrapper payloads reachable through bytes."""
    if failure == "building_in_tree_collection":
        pb = wire.ModelFile.FromString(exchange.dumps(populated(BuildingCollection)))
        pb.object.kind = wire.Object.TREE_COLLECTION
        return pb
    if failure == "wrapper_under_city":
        city = City(id="city")
        city.add_child(Object(id="child"))
        pb = wire.ModelFile.FromString(exchange.dumps(city))
        pb.object.children[0].kind = wire.Object.TREE_COLLECTION
        return pb
    cls, mutate = MUTATIONS[failure]
    pb = wire.ModelFile.FromString(exchange.dumps(populated(cls)))
    mutate(pb.object)
    return pb


@pytest.mark.parametrize("validate_schema", SCHEMA)
@pytest.mark.parametrize(
    "failure",
    [
        "building_in_tree_collection",
        "wrapper_under_city",
        *MUTATIONS,
    ],
)
def test_bytes_boundary_rejects_invalid_wrapper_payload(
    failure, validate_schema, tmp_path
):
    pb = payload_failure(failure)
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        exchange.loads(pb.SerializeToString(), validate_schema=validate_schema)
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        exchange.loads(pb, validate_schema=validate_schema)
    path = tmp_path / "invalid.dtcc"
    path.write_bytes(pb.SerializeToString())
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        io.load_model(path, validate_schema=validate_schema)


@pytest.mark.parametrize("cls", WRAPPERS, ids=lambda cls: cls.__name__)
def test_nested_dispatch_rejects_wrapper_kinds_below_the_root(cls):
    """Nested-dispatch behaviour: wrapper kinds decode only at the root.

    This exercises the current decoder's nested dispatch in _decode_object; it is
    not a run of a decoder that predates the wrapper kinds.
    """
    root = Object(id="root")
    member = Object(id="member")
    member.add_child(Object(id="nested"))
    root.add_child(member)
    pb = wire.ModelFile.FromString(exchange.dumps(root))
    pb.object.children[0].children[0].kind = KINDS[cls]
    message = f"Unsupported canonical object kind {KINDS[cls]}"
    with pytest.raises(NotImplementedError, match=message):
        exchange._decode_object(pb.object.children[0].children[0], {}, depth=2)
    for validate_schema in SCHEMA:
        with pytest.raises(NotImplementedError, match=message):
            exchange.loads(pb.SerializeToString(), validate_schema=validate_schema)
