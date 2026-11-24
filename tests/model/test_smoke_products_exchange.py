"""FieldSlice and StreamlineCollection through every canonical exchange boundary."""

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
    City,
    Field,
    FieldSlice,
    LineString,
    Object,
    Point,
    StreamlineCollection,
    exchange,
)
from dtcc_core.model import dtcc_pb2 as wire
from dtcc_core.model.profiles import SemanticProfile

PRODUCTS = (FieldSlice, StreamlineCollection)
KINDS = {
    FieldSlice: wire.Object.FIELD_SLICE,
    StreamlineCollection: wire.Object.STREAMLINE_COLLECTION,
}
SCHEMA = [True, False]
SLICE_ATTRIBUTES = (
    "slice_axis",
    "slice_position",
    "resolution",
    "axes",
    "name",
    "time",
    "period",
    "crs",
    "include_z",
    "domain_bounds",
    "metadata_payload",
)
STREAMLINE_ATTRIBUTES = (
    "seed_axis",
    "seed_position",
    "requested_line_count",
    "streamline_steps",
    "streamline_step_size",
    "time",
    "period",
    "crs",
    "include_z",
    "axes",
    "name",
    "domain_bounds",
    "metadata_payload",
)
DOMAIN = Bounds(xmin=0.0, ymin=0.0, zmin=0.0, xmax=10.0, ymax=20.0, zmax=40.0)


def populated(cls):
    """Populated product with Fields, a domain and a metadata payload."""
    if cls is FieldSlice:
        points = np.array(
            [[0.0, 0.0, 5.0], [10.0, 0.0, 5.0], [0.0, 20.0, 5.0], [10.0, 20.0, 5.0]]
        )
        source = FieldSlice(
            points=points,
            classification=np.array([1, 2, 2, 1], dtype=np.uint8),
            intensity=np.array([10, 20, 30, 40], dtype=np.uint16),
            return_number=np.array([1, 1, 2, 1], dtype=np.uint8),
            num_returns=np.array([1, 2, 2, 1], dtype=np.uint8),
            slice_axis="z",
            slice_position=0.125,
            resolution=2,
            axes=(0, 1),
            name="slice",
            time=1.5,
            period=8.0,
            crs="EPSG:3006",
            include_z=False,
            domain_bounds=DOMAIN,
            metadata_payload={"source": "synthetic", "levels": [1, 2], "note": None},
        )
        source.transform.srs = "EPSG:3006"
        source.transform.set_rotation(
            np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        )
        source.transform.set_translation(319000.25, 6399000.5, 2.0)
        for name, unit, dim in (("velocity", "m/s", 3), ("speed", "m/s", 1)):
            source.add_field(
                Field(
                    name=name,
                    unit=unit,
                    values=np.arange(4.0 * dim).reshape((4, dim)),
                    dim=dim,
                    association="sample",
                )
            )
        source.add_field(
            Field(
                name="pressure",
                unit="Pa",
                description="gauge",
                values=np.array([[1.0], [-2.0], [3.0], [-4.0]]),
                association="sample",
            )
        )
        return source
    lines = []
    for index in range(3):
        vertices = np.array([[0.0, index, 1.0], [1.0, index, 2.0], [2.0, index, 3.0]])
        line = LineString(vertices=vertices)
        line.add_field(
            Field(
                name="speed",
                unit="m/s",
                values=np.array([[1.0], [2.0], [index + 0.5]]),
                association="vertex",
            )
        )
        lines.append(line)
    lines[1].transform.set_translation(100.0, 200.0, 0.0)
    return StreamlineCollection(
        lines=lines,
        seed_axis="x",
        seed_position=0.25,
        requested_line_count=4,
        streamline_steps=3,
        streamline_step_size=0.5,
        time=2.0,
        period=6.0,
        crs="EPSG:3006",
        include_z=True,
        axes=(1, 2),
        name="lines",
        domain_bounds=DOMAIN,
        metadata_payload={"source": "synthetic", "seeds": 3},
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


def assert_fields_equal(a, b):
    assert [(f.name, f.unit, f.description, f.dim, f.association) for f in a] == [
        (f.name, f.unit, f.description, f.dim, f.association) for f in b
    ]
    for x, y in zip(a, b, strict=True):
        np.testing.assert_array_equal(x.values, y.values)


def assert_transform_equal(a, b):
    assert a.srs == b.srs
    np.testing.assert_array_equal(a.affine, b.affine)


def assert_product_equal(source, restored):
    assert type(restored) is type(source)
    names = SLICE_ATTRIBUTES if type(source) is FieldSlice else STREAMLINE_ATTRIBUTES
    for name in names:
        assert getattr(restored, name) == getattr(source, name), name
    assert type(restored.axes) is tuple
    assert restored.default_artifact_format == source.default_artifact_format
    if type(source) is FieldSlice:
        for name in (
            "points",
            "classification",
            "intensity",
            "return_number",
            "num_returns",
        ):
            np.testing.assert_array_equal(
                getattr(restored, name), getattr(source, name)
            )
        assert_transform_equal(source.transform, restored.transform)
        assert_fields_equal(source.fields, restored.fields)
        return
    assert len(restored.lines) == len(source.lines)
    for a, b in zip(source.lines, restored.lines, strict=True):
        assert type(b) is LineString
        np.testing.assert_array_equal(a.vertices, b.vertices)
        assert_transform_equal(a.transform, b.transform)
        assert_fields_equal(a.fields, b.fields)


@pytest.fixture(params=PRODUCTS, ids=lambda cls: cls.__name__)
def product(request):
    return populated(request.param)


@pytest.mark.parametrize("validate_schema", SCHEMA)
def test_populated_product_round_trips_bytes_and_proto(product, validate_schema):
    payload = exchange.dumps(product, validate_schema=validate_schema)
    pb = wire.ModelFile.FromString(payload)
    assert pb.version == 7
    assert pb.object.kind == KINDS[type(product)]
    assert pb.object.id == "collection"

    restored = exchange.loads(payload, validate_schema=validate_schema)
    assert_product_equal(product, restored)
    assert exchange.dumps(restored, validate_schema=validate_schema) == payload

    message = product.to_proto(validate_schema=validate_schema)
    target = type(product)()
    target.from_proto(message, validate_schema=validate_schema)
    assert_product_equal(product, target)
    target.from_proto(message.SerializeToString(), validate_schema=validate_schema)
    assert_product_equal(product, target)
    assert exchange.dumps(target) == payload


def test_slice_samples_keep_the_slice_transform_exactly():
    source = populated(FieldSlice)
    pb = wire.ModelFile.FromString(exchange.dumps(source))
    samples = {r.id: r for r in pb.object.representations}["samples"].geometry
    assert samples.WhichOneof("shape") == "point_cloud"
    assert len(samples.fields) == 3
    assert samples.transform.srs == "EPSG:3006"
    restored = exchange.loads(pb.SerializeToString())
    assert_transform_equal(source.transform, restored.transform)
    assert not restored.transform.is_identity


@pytest.mark.parametrize("validate_schema", SCHEMA)
def test_populated_product_round_trips_dtcc_file(product, validate_schema, tmp_path):
    path = tmp_path / "product.dtcc"
    io.save_model(product, path, validate_schema=validate_schema)
    restored = io.load_model(
        path, expected_type=type(product), validate_schema=validate_schema
    )
    assert_product_equal(product, restored)
    assert path.read_bytes() == exchange.dumps(restored)


@pytest.mark.parametrize("validate_schema", SCHEMA)
@pytest.mark.parametrize("name", ["package", "package.dtccpkg"])
def test_populated_product_round_trips_canonical_package(
    product, name, validate_schema, tmp_path
):
    product.dataset_context = context()
    package = product.export(
        tmp_path / name, canonical=True, validate_schema=validate_schema
    )
    assert package.artifacts[0].model_type == type(product).__name__
    restored = load_model_package(package.path, validate_schema=validate_schema)
    assert_product_equal(product, restored)
    assert restored.dataset_context == product.dataset_context
    assert exchange.dumps(restored) == exchange.dumps(product)


@pytest.mark.parametrize("cls", PRODUCTS, ids=lambda cls: cls.__name__)
@pytest.mark.parametrize("validate_schema", SCHEMA)
def test_empty_product_round_trips(cls, validate_schema, tmp_path):
    source = cls()
    payload = exchange.dumps(source, validate_schema=validate_schema)
    pb = wire.ModelFile.FromString(payload)
    assert "domain" not in {r.id for r in pb.object.representations}
    restored = exchange.loads(payload, validate_schema=validate_schema)
    assert_product_equal(source, restored)
    assert restored.domain_bounds is None
    assert exchange.dumps(restored, validate_schema=validate_schema) == payload
    source.dataset_context = context()
    package = source.export(
        tmp_path / "package.dtccpkg", canonical=True, validate_schema=validate_schema
    )
    restored = load_model_package(package.path, validate_schema=validate_schema)
    assert_product_equal(source, restored)
    assert restored.dataset_context == source.dataset_context


def test_bundled_schema_validates_product_as_transient_object(product):
    profile = SemanticProfile(files("dtcc_core").joinpath("schemas", "dtcc.yaml"))
    assert profile.validate(product).valid
    invalid = StreamlineCollection(lines=[Point()])
    with pytest.raises(ValueError, match="must contain LineString"):
        profile.validate(invalid)


# Expected rejection reason per failure, at encode and at the bytes boundary.
REASONS = {
    "slice_field_without_association": "requires an explicit",
    "line_field_without_association": "requires an explicit",
    "non_line_child": "must contain LineString objects",
    "slice_under_city": "FieldSlice is unsupported|kind 17",
    "streamlines_under_city": "StreamlineCollection is unsupported|kind 18",
    "float_resolution": "int resolution",
    "negative_resolution": "int resolution >= 0",
    "negative_requested_line_count": "requested_line_count >= 0",
    "negative_streamline_steps": "streamline_steps >= 0",
    "negative_streamline_step_size": "streamline_step_size >= 0",
    "three_axes": "two int axes",
    "axes_out_of_range": "two int axes in 0..2",
    "samples_not_point_cloud": "exactly the representations",
    "domain_not_bounds": "exactly the representations",
    "samples_lod": "exactly the representations",
    "extra_representation": "exactly the representations",
    "missing_samples": "exactly the representations",
    "missing_resolution": "exactly the attributes",
    "missing_requested_line_count": "exactly the attributes",
    "extra_attribute": "exactly the attributes",
    "text_time": "numeric",
    "root_transform": "transform must be identity",
    "root_id": "transient Object 'collection'",
    "carrier_field": "carrier MultiLineString",
    "carrier_transform": "carrier MultiLineString",
    "carrier_srs": "carrier MultiLineString",
}


NEGATIVE_STREAMLINE_COUNTS = {
    "negative_requested_line_count": {"requested_line_count": -1},
    "negative_streamline_steps": {"streamline_steps": -1},
    "negative_streamline_step_size": {"streamline_step_size": -0.5},
}


def native_failure(failure):
    """Invalid product state constructible through the native classes."""
    if failure == "slice_field_without_association":
        source = populated(FieldSlice)
        source.fields[1].association = None
        return source
    if failure == "line_field_without_association":
        source = populated(StreamlineCollection)
        source.lines[2].fields[0].association = None
        return source
    if failure == "non_line_child":
        source = populated(StreamlineCollection)
        source.lines.append(Point())
        return source
    if failure == "float_resolution":
        return FieldSlice(resolution=2.5)
    if failure == "negative_resolution":
        return FieldSlice(resolution=-1)
    if failure in NEGATIVE_STREAMLINE_COUNTS:
        return StreamlineCollection(**NEGATIVE_STREAMLINE_COUNTS[failure])
    if failure == "three_axes":
        return StreamlineCollection(axes=(0, 1, 2))
    if failure == "axes_out_of_range":
        return FieldSlice(axes=(5, -7))
    city = City(id="city")
    cls = FieldSlice if failure == "slice_under_city" else StreamlineCollection
    city.children[cls] = [populated(cls)]
    return city


@pytest.mark.parametrize("validate_schema", SCHEMA)
@pytest.mark.parametrize(
    "failure",
    [
        "slice_field_without_association",
        "line_field_without_association",
        "non_line_child",
        "float_resolution",
        "negative_resolution",
        *NEGATIVE_STREAMLINE_COUNTS,
        "three_axes",
        "axes_out_of_range",
        "slice_under_city",
        "streamlines_under_city",
    ],
)
def test_native_rejection_precedes_any_file_write(failure, validate_schema, tmp_path):
    model = native_failure(failure)
    with pytest.raises((ValueError, NotImplementedError), match=REASONS[failure]):
        exchange.validate(model)
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


def representation(root, key):
    return next(r for r in root.representations if r.id == key)


def replace_geometry(root, key, value):
    target = representation(root, key)
    target.ClearField("geometry")
    if type(value) is Bounds:
        target.bounds.CopyFrom(exchange._encode_bounds(value))
    else:
        exchange._encode_geometry(value, target.geometry)


def drop_representation(root, key):
    kept = [r for r in root.representations if r.id != key]
    root.ClearField("representations")
    root.representations.extend(kept)


def carrier_field(root):
    field = Field(name="n", values=np.arange(3.0), association="sample")
    representation(root, "lines").geometry.fields.add().CopyFrom(
        exchange._encode_field(field)
    )


def set_affine(transform):
    transform.affine[3] = 1.0


def non_line_child(root):
    lines = representation(root, "lines").geometry.multi_line_string.linestrings
    exchange._encode_geometry(Point(x=1.0), lines.add())


# Payload edits on an otherwise valid product; reachable only through bytes.
MUTATIONS = {
    "samples_not_point_cloud": (
        FieldSlice,
        lambda root: replace_geometry(
            root, "samples", LineString(vertices=np.zeros((2, 3)))
        ),
    ),
    "domain_not_bounds": (
        StreamlineCollection,
        lambda root: replace_geometry(root, "domain", Point()),
    ),
    "samples_lod": (
        FieldSlice,
        lambda root: setattr(representation(root, "samples"), "lod", "0"),
    ),
    "extra_representation": (
        StreamlineCollection,
        lambda root: root.representations.add(id="extra").bounds.CopyFrom(
            representation(root, "domain").bounds
        ),
    ),
    "missing_samples": (FieldSlice, lambda root: drop_representation(root, "samples")),
    "missing_resolution": (
        FieldSlice,
        lambda root: root.attributes.pop("resolution"),
    ),
    "missing_requested_line_count": (
        StreamlineCollection,
        lambda root: root.attributes.pop("requested_line_count"),
    ),
    "extra_attribute": (
        FieldSlice,
        lambda root: root.attributes["extra"].CopyFrom(wire.Value(text="x")),
    ),
    "float_resolution": (
        FieldSlice,
        lambda root: root.attributes["resolution"].CopyFrom(wire.Value(number=2.5)),
    ),
    "negative_resolution": (
        FieldSlice,
        lambda root: root.attributes["resolution"].CopyFrom(wire.Value(integer="-1")),
    ),
    "negative_requested_line_count": (
        StreamlineCollection,
        lambda root: root.attributes["requested_line_count"].CopyFrom(
            wire.Value(integer="-1")
        ),
    ),
    "negative_streamline_steps": (
        StreamlineCollection,
        lambda root: root.attributes["streamline_steps"].CopyFrom(
            wire.Value(integer="-1")
        ),
    ),
    "negative_streamline_step_size": (
        StreamlineCollection,
        lambda root: root.attributes["streamline_step_size"].CopyFrom(
            wire.Value(number=-0.5)
        ),
    ),
    "three_axes": (
        StreamlineCollection,
        lambda root: root.attributes["axes"].list.values.add(integer="2"),
    ),
    "axes_out_of_range": (
        FieldSlice,
        lambda root: (
            root.attributes["axes"].list.values[1].CopyFrom(wire.Value(integer="3"))
        ),
    ),
    "text_time": (
        StreamlineCollection,
        lambda root: root.attributes["time"].CopyFrom(wire.Value(text="1")),
    ),
    "root_transform": (FieldSlice, lambda root: set_affine(root.transform)),
    "root_id": (StreamlineCollection, lambda root: setattr(root, "id", "other")),
    "non_line_child": (StreamlineCollection, non_line_child),
    "carrier_field": (StreamlineCollection, carrier_field),
    "carrier_transform": (
        StreamlineCollection,
        lambda root: set_affine(representation(root, "lines").geometry.transform),
    ),
    "carrier_srs": (
        StreamlineCollection,
        lambda root: setattr(
            representation(root, "lines").geometry.transform, "srs", "EPSG:3006"
        ),
    ),
}


def payload_failure(failure):
    """Invalid product payloads reachable through bytes."""
    if failure in ("slice_under_city", "streamlines_under_city"):
        city = City(id="city")
        city.add_child(Object(id="child"))
        pb = wire.ModelFile.FromString(exchange.dumps(city))
        cls = FieldSlice if failure == "slice_under_city" else StreamlineCollection
        pb.object.children[0].kind = KINDS[cls]
        return pb
    cls, mutate = MUTATIONS[failure]
    pb = wire.ModelFile.FromString(exchange.dumps(populated(cls)))
    mutate(pb.object)
    return pb


@pytest.mark.parametrize("validate_schema", SCHEMA)
@pytest.mark.parametrize(
    "failure", ["slice_under_city", "streamlines_under_city", *MUTATIONS]
)
def test_bytes_boundary_rejects_invalid_product_payload(
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


@pytest.mark.parametrize("cls", PRODUCTS, ids=lambda cls: cls.__name__)
def test_nested_dispatch_rejects_product_kinds_below_the_root(cls):
    """Product kinds decode only at the root of the current decoder."""
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
