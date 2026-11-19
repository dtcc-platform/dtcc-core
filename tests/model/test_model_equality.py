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
