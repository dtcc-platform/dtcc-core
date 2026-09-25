"""Opt-in real-size I/O check: DTCC_RUN_LARGE_MODEL_TEST=1 pytest <this file>."""

import os

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.model import Field


@pytest.mark.skipif(os.environ.get('DTCC_RUN_LARGE_MODEL_TEST') != '1',
                    reason='Opt-in: writes >500 MiB and uses several GiB of RAM')
def test_model_over_256_mib_round_trips_through_file_and_package(tmp_path):
    # Incompressible data exercises both the native file and ZIP archive sizes.
    values = np.random.default_rng(0).integers(0, 256, size=260 * 1024 * 1024, dtype=np.uint8)
    field = Field(name='large', association='sample', values=values)
    path = tmp_path / 'large.dtcc'
    io.save_model(field, path)
    assert path.stat().st_size > 256 * 1024 * 1024
    restored = io.load_model(path)
    assert restored.values.dtype == values.dtype
    np.testing.assert_array_equal(restored.values, values)
    del restored

    field.dataset_context = DatasetContext(
        identity={'name': 'large-io-test', 'title': 'Large I/O test'},
        metadata={}, provenance={}, presentation={}, request={'dataset_name': 'large-io-test'},
    )
    package = field.export(tmp_path / 'large.dtccpkg', canonical=True)
    assert package.path.stat().st_size > 256 * 1024 * 1024
    restored = load_model_package(package.path)
    assert restored.dataset_context == field.dataset_context
    assert restored.values.dtype == values.dtype
    np.testing.assert_array_equal(restored.values, values)
    del restored

    # Exercise the archive extraction used by publication without uploading.
    manifest, artifacts = package._extract_archive_package(tmp_path / 'extracted')
    assert manifest.is_file()
    assert artifacts[0].stat().st_size == path.stat().st_size
