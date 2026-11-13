"""Canonical TimeSeries wire shapes, ownership and transactional I/O."""

import numpy as np
import pytest

from dtcc_core.io import load_model, save_model
from dtcc_core.model import (
    Field,
    Mesh,
    MultiSurface,
    Object,
    Surface,
    TimeSeries,
    exchange,
)


def series(dim=1, steps=3, elements=2, dtype="float64", association="vertex"):
    shape = (steps, elements) if dim == 1 else (steps, elements, dim)
    values = np.zeros(shape, dtype=dtype)
    if values.size:
        if values.dtype.kind == "u":
            values.flat[0] = np.iinfo(dtype).max
        elif values.dtype.kind == "i":
            values.flat[0] = np.iinfo(dtype).max
        elif values.dtype.kind == "f":
            values.flat[0] = np.nan
    times = np.array([-(2**63) + 1, -1, 2**63 - 1][:steps], dtype="int64")
    return TimeSeries(
        "temperature",
        "K",
        "測定",
        times.view("datetime64[ns]"),
        values,
        dim,
        association=association,
        role="observed",
    )


def assert_series(actual, expected):
    for name in ("name", "unit", "description", "dim", "association", "role"):
        assert getattr(actual, name) == getattr(expected, name)
    np.testing.assert_array_equal(actual.times, expected.times)
    np.testing.assert_array_equal(actual.values, expected.values)
    assert actual.values.dtype == expected.values.dtype.newbyteorder("<")
    assert actual.values.shape == expected.values.shape
    assert not actual.times.flags.writeable
    assert actual.values.flags.writeable


@pytest.mark.parametrize(
    "dtype",
    ["bool", "int8", "uint8", ">i2", ">u2", ">i4", ">u4", ">i8", ">u8", ">f4", ">f8"],
)
@pytest.mark.parametrize("dim", [1, 3])
@pytest.mark.parametrize("steps", [0, 3])
def test_standalone_round_trip(tmp_path, dtype, dim, steps):
    source = series(dim, steps, dtype=dtype)
    pb = source.to_proto()
    assert pb.version == 7
    assert pb.WhichOneof("root") == "time_series"
    assert pb.time_series.times.dtype == "<i8"
    assert list(pb.time_series.times.shape) == [steps]
    assert_series(exchange.loads(pb.SerializeToString()), source)
    target = TimeSeries()
    target.from_proto(pb)
    assert_series(target, source)
    path = tmp_path / "series.dtcc"
    save_model(source, path)
    assert_series(load_model(path, expected_type=TimeSeries), source)
    assert target.schema_id == pb.schema_id
    assert target.schema_version == pb.schema_version


def test_geometry_and_object_attachments():
    source = Mesh(vertices=np.zeros((2, 3)), time_series=[series(), series(dim=3)])
    source.fields = [Field(association="vertex", values=np.ones(2))]
    restored = exchange.loads(exchange.dumps(source))
    for a, b in zip(restored.time_series, source.time_series, strict=True):
        assert_series(a, b)
    np.testing.assert_array_equal(restored.fields[0].values, source.fields[0].values)
    obj = Object(id="owner")
    obj.add_geometry(source, id="mesh")
    restored = exchange.loads(exchange.dumps(obj))
    assert_series(
        restored.geometry["mesh"].geometry.time_series[0], source.time_series[0]
    )
    surface = Surface(vertices=np.zeros((3, 3)), time_series=[series(elements=3)])
    restored = exchange.loads(exchange.dumps(MultiSurface(surfaces=[surface])))
    assert_series(restored.surfaces[0].time_series[0], surface.time_series[0])


@pytest.mark.parametrize("steps", [0, 3])
@pytest.mark.parametrize(
    "association,elements", [("vertex", 3), ("face", 2), ("geometry", 2), ("edge", 2)]
)
def test_owner_count_and_association_rejected(steps, association, elements):
    source = Mesh(
        vertices=np.zeros((2, 3)),
        time_series=[series(steps=steps, elements=elements, association=association)],
    )
    with pytest.raises(ValueError, match="association/count"):
        exchange.dumps(source, validate_schema=False)


def test_empty_attached_series_keeps_element_count():
    source = Mesh(vertices=np.zeros((2, 3)), time_series=[series(dim=3, steps=0)])
    restored = exchange.loads(exchange.dumps(source))
    assert_series(restored.time_series[0], source.time_series[0])


@pytest.mark.parametrize(
    "damage",
    [
        "dtype",
        "rank",
        "nat",
        "duplicate",
        "descending",
        "length",
        "values_rank",
        "dim",
        "association",
        "missing_times",
        "missing_values",
        "bytes",
    ],
)
def test_invalid_wire_is_rejected_atomically(tmp_path, damage):
    pb = series().to_proto()
    ts = pb.time_series
    if damage == "dtype":
        ts.times.dtype = "<u8"
    elif damage == "rank":
        ts.times.shape[:] = [3, 1]
    elif damage in ("nat", "duplicate", "descending"):
        ticks = {
            "nat": [-(2**63), 0, 1],
            "duplicate": [0, 0, 1],
            "descending": [1, 0, -1],
        }[damage]
        ts.times.data = np.array(ticks, dtype="<i8").tobytes()
    elif damage == "length":
        ts.times.shape[:] = [2]
        ts.times.data = ts.times.data[:16]
    elif damage == "values_rank":
        ts.values.shape[:] = [6]
    elif damage == "dim":
        ts.dim = 0
    elif damage == "association":
        ts.association = ""
    elif damage == "missing_times":
        ts.ClearField("times")
    elif damage == "missing_values":
        ts.ClearField("values")
    elif damage == "bytes":
        ts.values.data = b""
    target = series()
    previous = target.copy()
    with pytest.raises(ValueError):
        target.from_proto(pb, validate_schema=False)
    assert_series(target, previous)


def test_append_capacity_and_failed_save(tmp_path):
    source = TimeSeries(times=["2020-01-01"], values=[[1, 2]], association="vertex")
    source.append("2020-01-02", [3, 4])
    assert len(source._times) > len(source)
    pb = source.to_proto()
    assert len(pb.time_series.times.data) == 16
    restored = exchange.loads(pb.SerializeToString())
    restored.append("2020-01-03", [5, 6])
    assert len(restored) == 3
    path = tmp_path / "series.dtcc"
    save_model(source, path)
    before = path.read_bytes()
    source.association = None
    with pytest.raises(ValueError, match="explicit"):
        save_model(source, path)
    assert path.read_bytes() == before


def test_field_slice_preserves_attachments():
    from dtcc_core.model import FieldSlice

    source = FieldSlice(points=np.zeros((2, 3)), time_series=[series()])
    restored = exchange.loads(exchange.dumps(source))
    assert_series(restored.time_series[0], source.time_series[0])


def test_semantic_rules_reach_standalone_and_nested_series(tmp_path):
    from importlib.resources import files

    import yaml

    from dtcc_core.model.profiles import SemanticProfile

    schema = yaml.safe_load(
        files("dtcc_core").joinpath("schemas/dtcc.yaml").read_text()
    )
    schema["classes"]["TimeSeries"]["attributes"]["role"]["pattern"] = "^predicted$"
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.safe_dump(schema))
    profile = SemanticProfile(path)
    source = series(elements=3)
    assert profile.validate(source).issues[0].path == "time_series.role"
    surface = Surface(vertices=np.zeros((3, 3)), time_series=[source])
    report = profile.validate(MultiSurface(surfaces=[surface]))
    assert report.issues[0].path == "geometry.surfaces[0].time_series[0].role"
    source.role = "predicted"
    assert profile.validate(source).valid


def test_invalid_native_attachments():
    mesh = Mesh(vertices=np.zeros((2, 3)), time_series=[series()])
    mesh.time_series[0].dataset_context = object()
    with pytest.raises(NotImplementedError, match="Nested Dataset Context"):
        exchange.dumps(mesh)
    mesh.time_series = [Field()]
    with pytest.raises(NotImplementedError, match="TimeSeries subtype"):
        exchange.dumps(mesh)
    mesh.time_series = ()
    with pytest.raises(ValueError, match="must be a list"):
        exchange.dumps(mesh)


def test_strict_cityjson_rejects_time_series(tmp_path):
    from dtcc_core.model import Building, City

    city = City(id="city")
    child = Building(id="building")
    child.add_geometry(
        MultiSurface(
            surfaces=[
                Surface(
                    vertices=np.array(
                        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
                    )
                )
            ],
            time_series=[series(elements=1, association="geometry")],
        ),
        id="cityjson-0",
        lod="2",
    )
    city.add_child(child)
    with pytest.raises(NotImplementedError, match="Time series"):
        city.save(tmp_path / "city.json", strict=True)


@pytest.mark.parametrize(
    "corrupt, message",
    [
        (
            lambda ts: ts._times.__setitem__(slice(0, 2), ts._times[1::-1].copy()),
            r"^TimeSeries\.times must be strictly increasing\.$",
        ),
        (
            lambda ts: ts.__dict__.__setitem__("_values", np.zeros((3, 2, 2))),
            r"^TimeSeries\.values must have shape",
        ),
    ],
    ids=["order", "shape"],
)
def test_save_checks_the_time_series_invariants(corrupt, message):
    source = series()
    corrupt(source)
    with pytest.raises(ValueError, match=message):
        exchange.dumps(source, validate_schema=False)
