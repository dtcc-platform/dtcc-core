"""Validation and lossless state at value serialization boundaries."""

import numpy as np
import pytest
from affine import Affine

from dtcc_core.model import Field, Raster


@pytest.mark.parametrize(
    "values,dim",
    [
        (np.array([1.0, 2.0]), 1),
        (np.array([[1.0], [2.0]]), 1),
        (np.array([[1.0, 2.0, 3.0]]), 3),
        (np.empty((0, 3)), 3),
    ],
)
def test_field_canonical_value_shape(values, dim):
    source = Field(
        name="velocity",
        unit="m/s",
        description="sample",
        values=values,
        dim=dim,
        association="sample",
    )
    restored = Field()
    restored.from_proto(source.to_proto().SerializeToString())
    assert (restored.name, restored.unit, restored.description, restored.dim) == (
        source.name,
        source.unit,
        source.description,
        dim,
    )
    np.testing.assert_array_equal(restored.values, values)


@pytest.mark.parametrize("dim", [0, -1, 1.5, True])
def test_field_rejects_invalid_dimension(dim):
    with pytest.raises(ValueError, match="positive integer"):
        Field(dim=dim).to_proto()


@pytest.mark.parametrize(
    "values,dim",
    [
        (np.ones(6), 3),
        (np.ones((2, 3)), 2),
        (np.ones((1, 2, 3)), 3),
        (np.array(1.0), 1),
        (np.array(["value"]), 1),
        (np.array([1j]), 1),
    ],
)
def test_field_rejects_invalid_value_shape_or_dtype(values, dim):
    with pytest.raises(ValueError, match="Field.values"):
        Field(values=values, dim=dim).to_proto()


def test_empty_raster_has_no_wire_pixels():
    raster = Raster()
    message = raster.to_proto()
    assert message.raster.data.dtype == "<f8"
    assert list(message.raster.data.shape) == []
    restored = Raster(data=np.ones((2, 2)))
    restored.from_proto(message.SerializeToString())
    assert restored.data.shape == ()
    assert np.isnan(restored.data)


def test_raster_bounds_include_rotated_and_sheared_corners():
    raster = Raster(data=np.ones((2, 3)), georef=Affine(1, 2, 10, 3, 4, 20))
    assert raster.bounds.tuple == (10, 20, 17, 37)


@pytest.mark.parametrize(
    "data",
    [
        np.array(42),
        np.ones(2),
        np.ones((1, 2, 3, 4)),
        np.ones((1, 2, 0)),
        np.array([["x"]]),
    ],
)
def test_raster_rejects_invalid_data_shape_or_dtype(data):
    with pytest.raises(ValueError, match="Raster.data"):
        Raster(data=data).to_proto()


@pytest.mark.parametrize("shape", [(0, 3), (0, 0), (2, 3, 1), (2, 3, 2)])
def test_raster_empty_axes_and_canonical_single_channel_shape(shape):
    source = Raster(data=np.zeros(shape, dtype=np.uint8))
    restored = Raster()
    restored.from_proto(source.to_proto().SerializeToString())
    expected = source.data
    np.testing.assert_array_equal(restored.data, expected)
    assert restored.data.dtype == source.data.dtype


def test_field_and_raster_round_trips_compare_equal():
    field = Field(
        name="t", unit="K", values=np.array([1.0, np.nan]), association="sample"
    )
    restored = Field()
    restored.from_proto(field.to_proto().SerializeToString())
    assert (restored == field) is True

    raster = Raster(
        data=np.array([[1.0, np.nan], [3.0, 4.0]]),
        georef=Affine(2.0, 0.0, 10.0, 0.0, -2.0, 20.0),
        nodata=np.nan,
        crs="EPSG:3006",
    )
    restored = Raster()
    restored.from_proto(raster.to_proto().SerializeToString())
    assert (restored == raster) is True
    restored.georef = Affine.identity()
    assert (restored == raster) is False


def test_dataset_value_compares_wrapped_and_raw_values():
    from dtcc_core.model import DatasetValue

    array = np.array([1.0, np.nan])
    wrapped = DatasetValue(array)
    assert (wrapped == DatasetValue(array.copy())) is True
    assert (wrapped == array.copy()) is True
    assert (wrapped == np.array([1.0, 2.0])) is False
    assert (wrapped == np.array([1, 2])) is False
    assert (DatasetValue({"a": [1, 2]}) == {"a": [1, 2]}) is True
    assert (DatasetValue({"a": [1, 2]}) == DatasetValue({"a": [1, 3]})) is False
    assert (DatasetValue(3) == 3) is True
    assert (DatasetValue(None) == DatasetValue(None)) is True
    with pytest.raises(TypeError):
        hash(wrapped)
