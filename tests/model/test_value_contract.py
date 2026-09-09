"""Validation and legacy compatibility at value serialization boundaries."""

import numpy as np
import pytest
from affine import Affine

from dtcc_core.model import Field, Raster, proto


@pytest.mark.parametrize("values,dim", [
    (np.array([1., 2.]), 1),
    (np.array([[1.], [2.]]), 1),
    (np.array([[1., 2., 3.]]), 3),
    (np.empty((0, 3)), 3),
])
def test_field_canonical_value_shape(values, dim):
    source = Field(name="velocity", unit="m/s", description="sample", values=values, dim=dim)
    restored = Field()
    restored.from_proto(source.to_proto().SerializeToString())
    assert (restored.name, restored.unit, restored.description, restored.dim) == (
        source.name, source.unit, source.description, dim
    )
    np.testing.assert_array_equal(restored.values, values.reshape(-1, dim))


@pytest.mark.parametrize("dim", [0, -1, 1.5, True])
def test_field_rejects_invalid_dimension(dim):
    with pytest.raises(ValueError, match="positive integer"):
        Field(dim=dim).to_proto()


@pytest.mark.parametrize("values,dim", [
    (np.ones(6), 3), (np.ones((2, 3)), 2), (np.ones((1, 2, 3)), 3),
    (np.array(1.), 1), (np.array(["value"]), 1), (np.array([1j]), 1),
])
def test_field_rejects_invalid_value_shape_or_dtype(values, dim):
    with pytest.raises(ValueError, match="Field.values"):
        Field(values=values, dim=dim).to_proto()


@pytest.mark.parametrize("message", [proto.Field(dim=0), proto.Field(dim=2, values=[1])])
def test_field_rejects_malformed_wire_values(message):
    with pytest.raises(ValueError, match="Field protobuf"):
        Field().from_proto(message.SerializeToString())


def test_empty_raster_has_no_wire_pixels():
    raster = Raster()
    message = raster.to_proto()
    assert not message.values
    restored = Raster(data=np.ones((2, 2)))
    restored.from_proto(message.SerializeToString())
    assert restored.data.shape == ()
    assert np.isnan(restored.data)


def test_legacy_empty_raster_scalar_is_not_a_pixel():
    # Old Raster() wrote its uninitialized shape-() scalar with zero dimensions.
    message = proto.Raster(values=[12345], dtype="float64", georef=[1, 0, 0, 0, 1, 0])
    restored = Raster(data=np.ones((2, 2)))
    restored.from_proto(message.SerializeToString())
    assert restored.data.shape == ()
    assert np.isnan(restored.data)
    assert not restored.to_proto().values


@pytest.mark.parametrize("dtype", [np.int32, np.uint32, np.int64, np.uint64])
def test_raster_rejects_integer_pixels_that_overflow_their_dtype_on_wire(dtype):
    raster = Raster(data=np.array([[np.iinfo(dtype).max]], dtype=dtype))
    with pytest.raises(ValueError, match="Raster protobuf pixels"):
        raster.to_proto()


def test_raster_rejects_integer_precision_loss_on_wire():
    raster = Raster(data=np.array([[2**24 + 1]], dtype=np.uint32))
    with pytest.raises(ValueError, match="cannot preserve these integer pixels"):
        raster.to_proto()


def test_raster_bounds_include_rotated_and_sheared_corners():
    raster = Raster(data=np.ones((2, 3)), georef=Affine(1, 2, 10, 3, 4, 20))
    assert raster.bounds.tuple == (10, 20, 17, 37)


def test_legacy_raster_grid_and_missing_metadata_reset_receiver():
    message = proto.Raster(values=[1, 2, 3, 4])
    message.grid.height = 2
    message.grid.width = 2
    raster = Raster(georef=Affine.translation(100, 200), crs="EPSG:3006")
    raster.from_proto(message.SerializeToString())
    np.testing.assert_array_equal(raster.data, [[1, 2], [3, 4]])
    assert raster.georef == Affine.identity()
    assert raster.data.dtype == np.float64
    assert raster.crs == ""


@pytest.mark.parametrize("message", [
    proto.Raster(values=[1, 2]),
    proto.Raster(height=-1, width=2, channels=1),
    proto.Raster(height=2, width=2, channels=1, values=[1]),
    proto.Raster(georef=[1, 2]),
    proto.Raster(georef=[1] * 7),
    proto.Raster(georef=[1, 0, float("inf"), 0, 1, 0]),
    proto.Raster(dtype="not-a-dtype"),
    proto.Raster(dtype="object"),
    proto.Raster(height=1, width=1, values=[256], dtype="uint8"),
    proto.Raster(height=1, width=1, values=[1.5], dtype="int16"),
    proto.Raster(height=1, width=1, values=[float("nan")], dtype="int16"),
])
def test_raster_rejects_malformed_wire_data(message):
    with pytest.raises(ValueError, match="[Rr]aster"):
        Raster().from_proto(message.SerializeToString())


@pytest.mark.parametrize("data", [np.array(42), np.ones(2), np.ones((1, 2, 3, 4)), np.ones((1, 2, 0)), np.array([["x"]])])
def test_raster_rejects_invalid_data_shape_or_dtype(data):
    with pytest.raises(ValueError, match="Raster.data"):
        Raster(data=data).to_proto()


@pytest.mark.parametrize("shape", [(0, 3), (0, 0), (2, 3, 1), (2, 3, 2)])
def test_raster_empty_axes_and_canonical_single_channel_shape(shape):
    source = Raster(data=np.zeros(shape, dtype=np.uint8))
    restored = Raster()
    restored.from_proto(source.to_proto().SerializeToString())
    expected = source.data[:, :, 0] if len(shape) == 3 and shape[2] == 1 else source.data
    np.testing.assert_array_equal(restored.data, expected)
    assert restored.data.dtype == source.data.dtype
