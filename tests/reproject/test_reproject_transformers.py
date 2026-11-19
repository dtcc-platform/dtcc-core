"""reproject_array builds one pyproj transformer per distinct CRS pair."""

from unittest import mock

import numpy as np
import pytest
from pyproj import Transformer

from dtcc_core.reproject.reproject import _transformer, reproject_array

POINTS = np.array([[320000.0, 6399000.0, 5.0], [500000.0, 6500000.0, -2.0]])


@pytest.fixture(autouse=True)
def fresh_cache():
    _transformer.cache_clear()


def test_transformer_built_once_per_crs_pair():
    real = Transformer.from_crs
    with mock.patch.object(Transformer, "from_crs", wraps=real) as from_crs:
        for _ in range(1000):
            reproject_array(POINTS, "EPSG:3006", "EPSG:4326")
        for _ in range(1000):
            reproject_array(POINTS, "EPSG:3006", "EPSG:3857")
    assert from_crs.call_count == 2
    for target in ("EPSG:4326", "EPSG:3857"):
        x, y = real("EPSG:3006", target, always_xy=True).transform(
            POINTS[:, 0], POINTS[:, 1], errcheck=True
        )
        expected = np.column_stack((x, y, POINTS[:, 2]))
        np.testing.assert_array_equal(
            reproject_array(POINTS, "EPSG:3006", target), expected
        )


def test_identity_returns_same_object():
    assert reproject_array(POINTS, "EPSG:3006", "EPSG:3006") is POINTS


def test_validation_precedes_transformer_construction():
    with mock.patch.object(Transformer, "from_crs", side_effect=RuntimeError("built")):
        with pytest.raises(ValueError, match="finite real XYZ"):
            reproject_array(np.zeros((2, 2)), "EPSG:3006", "EPSG:4326")
