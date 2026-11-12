"""Float64 conversion of Field values that refuses to round integers."""

from __future__ import annotations

import numpy as np
import pytest

from dtcc_core.model.values.field import _float64_exact

LIMIT = 2**53


@pytest.mark.parametrize(
    "dtype, value",
    [
        (np.int64, LIMIT + 1),
        (np.int64, -(LIMIT + 1)),
        (np.uint64, 2**64 - 1),
    ],
)
@pytest.mark.parametrize("shape", [(3,), (2, 3)])
def test_an_integer_float64_would_round_raises(dtype, value, shape):
    values = np.zeros(shape, dtype=dtype)
    values[(-1,) * len(shape)] = value

    with pytest.raises(ValueError, match=rf"Field 'speed'.*{value}"):
        _float64_exact(values, "Field 'speed'")


@pytest.mark.parametrize(
    "dtype, value",
    [(np.int64, 2**60), (np.int64, -(2**63)), (np.uint64, 2**63)],
)
@pytest.mark.parametrize("shape", [(3,), (2, 3)])
def test_a_large_integer_float64_holds_exactly_converts(dtype, value, shape):
    values = np.ones(shape, dtype=dtype)
    values[(0,) * len(shape)] = value

    result = _float64_exact(values, "Field 'speed'")

    assert result.dtype == np.float64
    assert result.shape == shape
    assert int(result[(0,) * len(shape)]) == value
    assert result.ravel()[1:].tolist() == [1.0] * (values.size - 1)


@pytest.mark.parametrize("shape", [(4,), (2, 2)])
def test_small_integers_bools_and_floats_pass_through(shape):
    for values in (
        np.arange(4, dtype=np.int64).reshape(shape),
        np.arange(4, dtype=np.uint8).reshape(shape),
        np.array([True, False, True, True]).reshape(shape),
        np.array([0.5, -2.0, np.nan, np.inf]).reshape(shape),
        np.array([2**53, -(2**53), 2**53 - 1, 0], dtype=np.int64).reshape(shape),
    ):
        result = _float64_exact(values, "Field 'x'")

        assert result.dtype == np.float64
        np.testing.assert_array_equal(result, values.astype(float))
