"""Model == contract: exact values, identical shape and dtype, NaN equals NaN."""

import copy
import dataclasses
import inspect

import numpy as np
import pytest

import dtcc_core.model as model
from dtcc_core.model import (
    Bounds,
    Building,
    DatasetValue,
    Field,
    GeometryRepresentation,
    Grid,
    Mesh,
    PointCloud,
    Raster,
    SemanticRegion,
    Transform,
    VolumeGrid,
    exchange,
)
from dtcc_core.model._equality import model_equal
from dtcc_core.model.model import Model


def _exported_dataclasses():
    return sorted(
        (
            cls
            for name in model.__all__
            if inspect.isclass(cls := getattr(model, name))
            and dataclasses.is_dataclass(cls)
        ),
        key=lambda cls: cls.__name__,
    )


def _fresh(cls):
    if cls is SemanticRegion:
        return SemanticRegion("urn:test:roof", indices=np.array([0, 2]))
    if cls is GeometryRepresentation:
        return GeometryRepresentation(Mesh(vertices=np.zeros((3, 3))), lod="lod1")
    return cls()


def test_every_exported_model_class_uses_the_shared_comparator():
    classes = _exported_dataclasses()
    assert {Mesh, Grid, Transform, Field, SemanticRegion} <= set(classes)
    for cls in classes:
        assert cls.__hash__ is None, cls
        if cls is DatasetValue:
            assert cls.__eq__ is not model_equal
        else:
            # Defined explicitly, not by the dataclass code generator, so
            # behaviour does not depend on the Python version.
            assert cls.__eq__ is model_equal, cls


@pytest.mark.parametrize(
    "cls",
    [cls for cls in _exported_dataclasses() if not inspect.isabstract(cls)],
    ids=lambda cls: cls.__name__,
)
def test_independent_default_objects_compare_with_one_bool(cls):
    first, second = _fresh(cls), _fresh(cls)
    if hasattr(first, "id") and not isinstance(first, SemanticRegion):
        second.id = first.id
    result = first == second
    assert result is True
    assert (first != second) is False
    assert (first == copy.deepcopy(first)) is True


def test_object_ids_still_distinguish_objects():
    assert (Building() == Building()) is False
    assert (Building(id="a") == Building(id="a")) is True


def test_concrete_type_matters():
    assert (Mesh() == PointCloud()) is False
    assert (Grid() == VolumeGrid()) is False
    assert (Mesh() == object()) is False


def test_models_stay_unhashable():
    with pytest.raises(TypeError):
        hash(Mesh())
    with pytest.raises(TypeError):
        hash(SemanticRegion("urn:test:wall"))


@pytest.mark.parametrize(
    "values",
    [np.empty(0), np.array([1.5]), np.arange(12.0).reshape(4, 3)],
    ids=["empty", "singleton", "multi"],
)
def test_equal_and_unequal_arrays_of_any_size(values):
    first = Field(name="t", values=values)
    assert (first == Field(name="t", values=values.copy())) is True
    if values.size:
        changed = values.copy()
        changed.flat[-1] += 1
        assert (first == Field(name="t", values=changed)) is False
    else:
        assert (first == Field(name="t", values=np.array([0.0]))) is False


def test_shape_must_match_exactly():
    # NumPy broadcasting would make these equal; == requires the same shape.
    assert (Field(values=np.array([1.0])) == Field(values=np.array([[1.0]]))) is False
    assert (Mesh(markers=np.array([1.0])) == Mesh(markers=np.array([[1.0]]))) is False
    assert (Raster(data=np.ones((1, 1))) == Raster(data=np.ones((1, 1, 1)))) is False


def test_dtype_must_match():
    ints64 = Mesh(faces=np.array([[0, 1, 2]], dtype=np.int64))
    ints32 = Mesh(faces=np.array([[0, 1, 2]], dtype=np.int32))
    assert (ints64 == ints32) is False
    assert (ints64 == copy.deepcopy(ints64)) is True


def test_values_are_exact_not_tolerant():
    first = PointCloud(points=np.array([[0.1, 0.2, 0.3]]))
    second = PointCloud(points=np.array([[0.1, 0.2, 0.3 + 1e-15]]))
    assert (first == second) is False


def test_nan_equals_nan_in_arrays_and_scalars():
    values = np.array([1.0, np.nan, 3.0])
    assert (Field(values=values) == Field(values=values.copy())) is True
    assert (Field(values=values) == Field(values=np.array([1.0, 2.0, 3.0]))) is False
    assert (Raster() == Raster()) is True
    assert (Raster(nodata=np.nan) == Raster(nodata=-9999.0)) is False


def test_new_subclass_fields_are_compared():
    @dataclasses.dataclass(repr=False)
    class Labelled(Mesh):
        labels: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros(2))

    assert Labelled.__eq__ is model_equal
    assert (Labelled() == Labelled()) is True
    assert (Labelled() == Labelled(labels=np.ones(2))) is False
    assert (Labelled() == Mesh()) is False
    with pytest.raises(TypeError):
        hash(Labelled())


def test_custom_equality_is_kept_and_inherited():
    @dataclasses.dataclass(repr=False)
    class AlwaysEqual(Model):
        payload: np.ndarray = dataclasses.field(default_factory=lambda: np.ones(2))

        def __eq__(self, other):
            return True

    @dataclasses.dataclass(repr=False)
    class Child(AlwaysEqual):
        extra: int = 0

    assert (AlwaysEqual() == AlwaysEqual(payload=np.zeros(2))) is True
    assert Child.__eq__ is AlwaysEqual.__eq__
    assert (Child(extra=1) == Child(extra=2)) is True


def test_shallow_replace_and_deep_copy():
    mesh = Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        faces=np.array([[0, 1, 2]]),
    )
    assert (dataclasses.replace(mesh) == mesh) is True
    deep = copy.deepcopy(mesh)
    assert deep.vertices is not mesh.vertices
    assert (deep == mesh) is True
    deep.vertices[0, 0] = 9.0
    assert (deep == mesh) is False
    assert (dataclasses.replace(mesh, faces=np.array([[2, 1, 0]])) == mesh) is False


def test_coordinate_bounds_cache_is_ignored():
    vertices = np.array([[0.0, 0.0, 0.0], [2.0, 3.0, 0.0]])
    first = PointCloud(points=vertices)
    second = PointCloud(points=vertices.copy())
    first.calculate_bounds()
    assert first.bounds != second._bounds
    assert (first == second) is True

    building = Building(id="b")
    other = Building(id="b")
    building.bounds = Bounds(xmax=5.0, ymax=5.0)
    assert (building == other) is True


def test_grid_domain_is_compared():
    assert (Grid(width=2, height=3) == Grid(width=2, height=3)) is True
    read = Grid(width=2, height=3)
    read.calculate_bounds()
    assert (read == Grid(width=2, height=3)) is True
    moved = Grid(width=2, height=3)
    moved.bounds = Bounds(xmin=10.0, ymin=10.0, xmax=12.0, ymax=13.0)
    assert (moved == Grid(width=2, height=3)) is False

    fresh = VolumeGrid(width=1, height=1, depth=2)
    assert (fresh == VolumeGrid(width=1, height=1, depth=2)) is True
    assert fresh._bounds is None
    deep = VolumeGrid(width=1, height=1, depth=2)
    deep.bounds = Bounds(zmax=5.0)
    assert (deep == fresh) is False


def test_dataset_context_and_other_attributes_are_not_compared():
    first, second = Mesh(), Mesh()
    first.dataset_context = object()
    first.schema_id = "urn:test:schema"
    assert (first == second) is True


def test_named_tuples_and_missing_dates():
    from collections import namedtuple

    Pair = namedtuple("Pair", "x y")
    assert (DatasetValue(Pair(1.0, np.nan)) == (1.0, np.nan)) is True
    assert (DatasetValue([1, 2]) == (1, 2)) is False

    dates = np.array(["NaT", "2020-01-01"], dtype="datetime64[s]")
    assert (DatasetValue(dates) == DatasetValue(dates.copy())) is True
    assert (DatasetValue(dates) == DatasetValue(dates[::-1].copy())) is False


def test_dataset_value_equals_raw_value_when_nested():
    array = np.array([1.0, np.nan])
    wrapped = DatasetValue(array.copy())
    assert ([1, 2] == DatasetValue([1, 2])) is True
    assert (DatasetValue([wrapped, 2]) == [array, 2]) is True
    assert (DatasetValue(DatasetValue(wrapped)) == array) is True
    looped = DatasetValue()
    looped.value = looped
    assert (looped == DatasetValue(looped)) is True
    assert (looped == DatasetValue(1)) is False
    assert (looped == 1) is False
    # Arithmetic through NumPy still treats a DatasetValue as a sequence.
    assert np.sum(DatasetValue([1, 2])) == 3
    assert (DatasetValue({"a": DatasetValue(3)}) == {"a": 3}) is True


def test_record_array_fields_follow_the_same_rules():
    kind = np.dtype([("height", "f8"), ("floors", "i4")])
    records = np.array([(np.nan, 2), (10.0, 3)], dtype=kind)
    assert (Field(values=records) == Field(values=records.copy())) is True
    changed = records.copy()
    changed["floors"][1] = 4
    assert (Field(values=records) == Field(values=changed)) is False

    mixed = np.array([(np.nan, "a")], dtype=[("height", "f8"), ("name", "O")])
    assert (DatasetValue(mixed) == DatasetValue(mixed.copy())) is True


def test_raster_tiles_compare_missing_resolution():
    from datetime import datetime
    from pathlib import Path

    from dtcc_core.model.values.raster_tiles import RasterTile

    def tile(resolution):
        return RasterTile(
            path=Path("tile.tif"),
            id="t1",
            collection="ortho",
            datetime=datetime(2024, 5, 1),
            extent=Bounds(xmax=1.0, ymax=1.0),
            crs="EPSG:3006",
            spektraltyp="RGB",
            resolution=resolution,
            size_bytes=None,
        )

    assert RasterTile.__eq__ is model_equal
    assert (tile(np.nan) == tile(np.nan)) is True
    assert (tile(0.5) == tile(0.25)) is False


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        (complex(np.nan, 1.0), complex(np.nan, 1.0), True),
        (complex(np.nan, 1.0), complex(np.nan, 2.0), False),
        (complex(1.0, np.nan), complex(2.0, np.nan), False),
        (complex(np.nan, 1.0), complex(1.0, np.nan), False),
        (complex(np.nan, np.nan), complex(np.nan, np.nan), True),
        (complex(1.0, 2.0), complex(1.0, 2.0), True),
    ],
)
def test_complex_parts_compare_separately(left, right, expected):
    first = Field(values=np.array([left]))
    second = Field(values=np.array([right]))
    assert (first == second) is expected
    assert (DatasetValue(left) == DatasetValue(right)) is expected


def test_complex_scalar_with_finite_part_differs_from_nan():
    assert (DatasetValue(complex(np.nan, 1.0)) == DatasetValue(float("nan"))) is False
    assert (DatasetValue(complex(np.nan, 0.0)) == DatasetValue(float("nan"))) is True


def test_masked_arrays_compare_masks_then_visible_values():
    def raster(values, mask):
        return Raster(data=np.ma.array([values], mask=[mask]))

    base = raster([1.0, 2.0], [False, True])
    assert (base == raster([1.0, 2.0], [False, True])) is True
    assert (base == raster([1.0, 99.0], [False, True])) is True
    assert (base == raster([1.0, 2.0], [False, False])) is False
    assert (base == raster([3.0, 2.0], [False, True])) is False
    assert (base == Raster(data=np.array([[1.0, 2.0]]))) is False


def test_masked_arrays_keep_nan_rule_and_record_fields():
    nan = Field(values=np.ma.array([np.nan, 1.0], mask=[False, True]))
    assert (nan == Field(values=np.ma.array([np.nan, 2.0], mask=[False, True]))) is True

    dtype = [("a", "f8"), ("b", "i4")]

    def record(mask):
        return Field(values=np.ma.array([(1.0, 2)], dtype=dtype, mask=[mask]))

    assert (record((False, True)) == record((False, True))) is True
    assert (record((False, True)) == record((True, False))) is False


def test_layered_volume_grid_compares_with_unset_and_initialized_domain():
    first = VolumeGrid(width=1, height=1, z_levels=[100, 101, 104])
    second = VolumeGrid(width=1, height=1, z_levels=[100, 101, 104])
    assert (first == second) is True
    assert (second.bounds.zmin, second.bounds.zmax) == (100, 104)
    assert (first == second) is True
    assert first._bounds is None
    np.testing.assert_array_equal(second.coordinates()[:3, 2], [100, 101, 104])


def test_assigned_integer_levels_equal_their_native_round_trip():
    assigned = VolumeGrid(width=1, height=1, depth=2)
    assigned.z_levels = np.array([100, 101, 104], dtype=np.int8)
    restored = exchange.loads(exchange.dumps(copy.deepcopy(assigned)))
    assert (assigned == restored) is True
    assert assigned.z_levels.dtype == np.int8
    assert assigned._bounds is None


@pytest.mark.parametrize(
    "levels",
    [[], [2, 1, 0], [0, np.nan, 1], ["a", "b", "c"], [[0], [1, 2]]],
    ids=["empty", "decreasing", "nan", "text", "ragged"],
)
def test_invalid_assigned_levels_still_compare_to_one_bool(levels):
    first, second = VolumeGrid(depth=2), VolumeGrid(depth=2)
    first.z_levels = levels
    second.z_levels = levels
    assert (first == second) is True
    assert (first == VolumeGrid(depth=2)) is False


def test_levels_that_round_in_float64_stay_distinct():
    first, second = VolumeGrid(depth=1), VolumeGrid(depth=1)
    first.z_levels = np.array([2**53, 2**53 + 2], dtype=np.int64)
    second.z_levels = np.array([2**53 + 1, 2**53 + 2], dtype=np.int64)
    assert (first == second) is False


@pytest.mark.parametrize(
    "make",
    [
        lambda: np.array([1.0, 2.0], dtype=">f8"),
        lambda: np.ma.array([1.0, 2.0]),
    ],
    ids=["big-endian", "masked"],
)
def test_native_round_trip_keeps_values_of_other_representations(make):
    source = Field(name="t", association="sample", values=make())
    restored = exchange.loads(exchange.dumps(source))
    np.testing.assert_array_equal(np.asarray(source.values), restored.values)
    assert type(restored.values) is np.ndarray
    assert (source == restored) is False


@pytest.mark.parametrize(
    "levels",
    [
        np.array(5.0),
        np.ma.masked,
        np.ma.array(5.0, mask=True),
        np.zeros((3, 1)),
        np.array(["2020-01-01"] * 3, dtype="datetime64[D]"),
        np.zeros(3, dtype=[("a", "f8")]),
        np.array([None, 1, 2], dtype=object),
        np.array(["a", "b", "c"]),
        np.array([1j, 2j, 3j]),
    ],
    ids=[
        "0d",
        "masked-constant",
        "0d-masked",
        "2d",
        "datetime",
        "record",
        "object",
        "text",
        "complex",
    ],
)
@pytest.mark.filterwarnings("error:Casting complex values")
def test_invalid_array_levels_still_compare_to_one_bool(levels):
    first, second = VolumeGrid(depth=2), VolumeGrid(depth=2)
    first.z_levels = levels
    second.z_levels = copy.deepcopy(levels)
    assert (first == second) is True
    assert (first == VolumeGrid(depth=2)) is False


def test_masked_levels_compare_their_masks():
    def grid(levels):
        result = VolumeGrid(width=1, height=1, depth=2)
        result.z_levels = levels
        return result

    masked = grid(np.ma.array([0.0, 1.0, 2.0], mask=[False, True, False]))
    assert (masked == grid(np.array([0.0, 1.0, 2.0]))) is False
    other_mask = np.ma.array([0.0, 1.0, 2.0], mask=[True, False, False])
    assert (masked == grid(other_mask)) is False
    hidden = np.ma.array([0.0, 9.0, 2.0], mask=[False, True, False])
    assert (masked == grid(hidden)) is True


def test_mixed_list_levels_that_round_stay_distinct():
    first, second = VolumeGrid(depth=1), VolumeGrid(depth=1)
    first.z_levels = [2**53 + 1, 2.0**53 + 2]
    second.z_levels = [2**53, 2.0**53 + 2]
    assert (first == second) is False


def test_masked_record_with_object_field_compares():
    def value(text):
        array = np.ma.array([(text, 3)], dtype=[("o", "O"), ("k", "i4")])
        return Field(values=array[0])

    assert (value("x") == value("x")) is True
    assert (value("x") == value("y")) is False


@pytest.mark.parametrize("dtype", [[], [("e", []), ("a", "f8")]], ids=["top", "nested"])
@pytest.mark.parametrize("zeros", [np.zeros, np.ma.zeros], ids=["plain", "masked"])
def test_empty_record_dtypes_compare(dtype, zeros):
    assert (
        Field(values=zeros(2, dtype=dtype)) == Field(values=zeros(2, dtype=dtype))
    ) is True


def test_array_like_levels_compare():
    class Levels:
        def __array__(self, dtype=None, copy=None):
            return np.array([0.0, 1.0, 2.0])

    first = VolumeGrid(width=1, height=1, depth=2)
    second = VolumeGrid(width=1, height=1, depth=2)
    first.z_levels, second.z_levels = Levels(), Levels()
    assert (first == second) is True


def test_mixed_list_levels_keep_their_bounds():
    grid = VolumeGrid(width=1, height=1, depth=1)
    grid.z_levels = [2**53 + 1, 2.0**53 + 2]
    assert (grid.bounds.zmin, grid.bounds.zmax) == (2.0**53, 2.0**53 + 2)
    assert len(grid.coordinates()) == grid.num_vertices


@pytest.mark.parametrize("index", [0, 2], ids=["first", "last"])
def test_hidden_end_level_is_ignored(index):
    def grid(hidden):
        values = np.array([0.0, 1.0, 2.0])
        values[index] = hidden
        mask = [i == index for i in range(3)]
        result = VolumeGrid(width=1, height=1, depth=2)
        result.z_levels = np.ma.array(values, mask=mask)
        return result

    assert (grid(-5.0) == grid(9.0)) is True


def test_record_scalar_keeps_nan_rule():
    record = np.array([(np.nan, 1)], dtype=[("a", "f8"), ("b", "i4")])
    assert (Field(values=record[0]) == Field(values=record.copy()[0])) is True
    changed = record.copy()
    changed["b"] = 2
    assert (Field(values=record[0]) == Field(values=changed[0])) is False


@pytest.mark.parametrize("other", [1.0, "x", np.datetime64("NaT")])
def test_record_scalar_differs_from_non_record(other):
    record = np.array([(1.0,)], dtype=[("a", "f8")])[0]
    assert (Field(values=record) == Field(values=other)) is False
    assert (Field(values=other) == Field(values=record)) is False
    assert (DatasetValue(record) == DatasetValue(other)) is False


@pytest.mark.parametrize(
    "nat", [np.datetime64("NaT", "s"), np.timedelta64("NaT", "s")], ids=["M8", "m8"]
)
def test_nat_equals_nat_in_masked_records_and_object_arrays(nat):
    dtype = [("t", nat.dtype), ("x", "f8")]
    field = Field(values=np.ma.array([(nat, 1.0)], dtype=dtype)[0])
    assert (field == copy.deepcopy(field)) is True

    def objects():
        array = np.empty(1, dtype=object)
        array[0] = nat
        return Field(values=array)

    assert (objects() == objects()) is True


def test_datetime_nat_differs_from_timedelta_nat():
    first = DatasetValue(np.datetime64("NaT"))
    second = DatasetValue(np.timedelta64("NaT"))
    assert (first == second) is False
