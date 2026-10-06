"""Opt-in check of the orthophoto dataset against a configured tile server.

Run only when authorized, from the Chalmers network or VPN::

    DTCC_LIVE_DATASET_TESTS=1 DTCC_ORTHOPHOTO_URL=<service root> \\
        pytest tests/datasets/live/test_orthophoto_live.py --run-live

Each run downloads one original 2025 orthophoto (about 659 MB) into a temporary
cache that is removed afterwards.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.io.data import cache as data_cache
from dtcc_core.io.data import orthophoto as client
from dtcc_core.model import Raster
from dtcc_core.model.values.raster_tiles import RasterTileCollection


pytestmark = pytest.mark.live

# 100 m inside one fully imaged 2025 cell.
BOUNDS = (673_700.0, 6_578_700.0, 673_800.0, 6_578_800.0)
YEAR = 2025
PINNED = [("orto-o2-2025", "o65775_6725_25_mr25")]
TRANSIENT = {"connection", "timeout", "http_5xx"}


def _require_server_url() -> None:
    value = os.environ.get(client.SERVER_URL_ENV)
    if value is None or not value.strip():
        pytest.skip(
            f"Set {client.SERVER_URL_ENV} to the DTCC LM tile server root; it is "
            "reachable from the Chalmers network or VPN only."
        )


def _skip_if_transient(error: DatasetUpstreamError):
    """Skip on failures that say nothing about the code, else re-raise."""
    if error.failure_class in TRANSIENT:
        pytest.skip(
            f"Tile server unreachable ({error.failure_class}); it serves the "
            "Chalmers network or VPN only."
        )
    raise error


def guard_downloads(download_items):
    """Wrap download_items to download only the pinned item.

    The guard checks the items of the call's one manifest before any download,
    so a manifest listing more (or other) items fails without fetching them.
    """

    def guarded(items, **kwargs):
        listed = [(item.collection, item.id) for item in items]
        assert listed == PINNED, f"refusing to download {listed}, expected {PINNED}"
        return download_items(items, **kwargs)

    return guarded


def check_mosaic(raster, bounds, cache_root) -> None:
    """Structural checks of a strict mosaic, independent of the imagery."""
    assert isinstance(raster, Raster)
    assert raster.crs == "EPSG:3006"
    data = raster.data
    assert data.dtype == np.uint8 and data.ndim == 3 and data.shape[2] == 4
    alpha = data[..., 3]
    assert set(np.unique(alpha).tolist()) <= {0, 255}
    assert not data[..., :3][alpha == 0].any()
    context = raster.dataset_context
    health = context.health
    assert health["status"] == "complete" and health["upstream_error_count"] == 0
    assert health["valid_pixels"] > 0
    assert health["valid_pixels"] == int(np.count_nonzero(alpha == 255))
    assert health["total_pixels"] == data.shape[0] * data.shape[1]
    xmin, ymin, xmax, ymax = bounds
    grid = health["bounds"]
    step = health["resolution"]
    assert (grid[0], grid[3]) == (xmin, ymax)
    assert 0 <= grid[2] - xmax < step and 0 <= ymin - grid[1] < step
    sources = [source for source in context.provenance.sources
               if isinstance(source, dict) and source.get("role") == "source_item"]
    known = [source["resolution"] for source in sources
             if source["resolution"] is not None]
    assert sources and (not known or step <= min(known))
    assert context.metadata.collection_period is not None
    assert str(cache_root) not in context.model_dump_json()


def check_tiles(tiles, bounds, cache_root, requests_made) -> None:
    """The tiles call reuses the cache: one manifest request, no file requests."""
    assert isinstance(tiles, RasterTileCollection)
    assert len(tiles) == tiles.dataset_context.health["items_listed"] >= 1
    xmin, ymin, xmax, ymax = bounds
    for tile in tiles:
        path = Path(tile.path)
        assert path.suffix == ".tif" and path.is_file()
        assert Path(cache_root) in path.parents
        if tile.size_bytes is not None:
            assert path.stat().st_size == tile.size_bytes
        assert tile.crs == "EPSG:3006"
        extent = tile.extent
        assert extent.xmin < xmax and xmin < extent.xmax
        assert extent.ymin < ymax and ymin < extent.ymax
    assert [url for url in requests_made if url.endswith("/items")] == [
        requests_made[0]
    ]
    assert [url for url in requests_made if "/files/" in url] == []


def test_orthophoto_live_mosaic_and_tiles(monkeypatch):
    _require_server_url()
    with tempfile.TemporaryDirectory() as directory:
        cache_root = Path(directory)
        monkeypatch.setattr(data_cache, "cache_dir", cache_root)
        monkeypatch.setattr(
            client, "download_items", guard_downloads(client.download_items)
        )
        try:
            raster = datasets.orthophoto(bounds=BOUNDS, year=YEAR, strict_live=True)
        except DatasetUpstreamError as error:
            _skip_if_transient(error)
        check_mosaic(raster, BOUNDS, cache_root)
        requests_made = []
        real_get = client.requests.get

        def recording_get(url, **kwargs):
            requests_made.append(url)
            return real_get(url, **kwargs)

        monkeypatch.setattr(client.requests, "get", recording_get)
        tiles = datasets.orthophoto(
            bounds=BOUNDS, year=YEAR, product="tiles", strict_live=True
        )
        check_tiles(tiles, BOUNDS, cache_root, requests_made)
