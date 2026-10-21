import asyncio
import os
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock

import pytest

from dtcc_core.io.data import geopkg, lidar


@pytest.fixture(params=[(geopkg, "tile.gpkg"), (lidar, "tile.laz")])
def downloader(request, monkeypatch):
    module, filename = request.param
    payload = b"complete tile content"
    response = MagicMock()
    response.status = 200
    response.__aenter__.return_value = response
    response.read = AsyncMock(return_value=payload)

    async def chunks(_size):
        yield payload[:5]
        await asyncio.sleep(0)
        yield payload[5:]

    response.content.iter_chunked.side_effect = chunks
    session = MagicMock()
    session.__aenter__.return_value = session
    session.get.return_value = response
    monkeypatch.setattr(module.aiohttp, "ClientSession", lambda **kwargs: session)
    if module is lidar:
        monkeypatch.setattr(lidar, "_DOWNLOAD_MAX_ATTEMPTS", 1)
    return module, filename, payload, session


def test_concurrent_downloads_publish_complete_tiles(
    downloader, monkeypatch, tmp_path: Path
) -> None:
    module, filename, payload, session = downloader
    replace = os.replace
    barrier = threading.Barrier(2)
    temporary_paths = []

    def replace_together(src, dst):
        assert Path(src).read_bytes() == payload
        temporary_paths.append(Path(src))
        barrier.wait(timeout=5)
        replace(src, dst)

    monkeypatch.setattr(module.os, "replace", replace_together)

    def download():
        module.run_download_files("http://example.test", [filename], tmp_path)

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(download) for _ in range(2)]
        for future in futures:
            future.result(timeout=10)

    assert len(set(temporary_paths)) == 2
    assert all(path.parent == tmp_path for path in temporary_paths)
    assert (tmp_path / filename).read_bytes() == payload
    assert not list(tmp_path.glob("*.part"))
    assert session.get.call_count == 2
    download()
    assert session.get.call_count == 2


def test_download_that_loses_publication_race_keeps_published_tile(
    downloader, monkeypatch, tmp_path: Path
) -> None:
    module, filename, payload, _session = downloader

    # Windows refuses a replace while a concurrent download publishes the tile.
    def lose_race(src, dst):
        Path(dst).write_bytes(payload)
        raise PermissionError("Access is denied")

    monkeypatch.setattr(module.os, "replace", lose_race)
    module.run_download_files("http://example.test", [filename], tmp_path)

    assert (tmp_path / filename).read_bytes() == payload
    assert not list(tmp_path.glob("*.part"))


@pytest.mark.parametrize("error", [OSError, PermissionError])
def test_failed_download_removes_only_its_own_temporary_file(
    downloader, monkeypatch, tmp_path: Path, error
) -> None:
    module, filename, payload, _session = downloader
    other_download = tmp_path / f"{filename}.part"
    other_download.write_bytes(b"another download")

    def fail_replace(src, dst):
        assert Path(src).read_bytes() == payload
        raise error("publication failed")

    monkeypatch.setattr(module.os, "replace", fail_replace)
    expected_error = geopkg.FootprintDownloadError if module is geopkg else OSError
    with pytest.raises(expected_error, match="publication failed"):
        module.run_download_files("http://example.test", [filename], tmp_path)

    assert not (tmp_path / filename).exists()
    assert other_download.read_bytes() == b"another download"
    assert list(tmp_path.glob("*.part")) == [other_download]


@pytest.mark.parametrize("downloader", [(lidar, "tile.laz")], indirect=True)
def test_lidar_retry_cleans_up_partial_download(
    downloader, monkeypatch, tmp_path: Path
) -> None:
    module, filename, payload, session = downloader
    attempts = 0

    async def chunks(_size):
        nonlocal attempts
        attempts += 1
        yield payload[:5]
        if attempts == 1:
            raise module.aiohttp.ClientError("interrupted stream")
        yield payload[5:]

    async def backoff(_delay):
        assert not list(tmp_path.glob("*.part"))
        assert not (tmp_path / filename).exists()

    session.get.return_value.content.iter_chunked.side_effect = chunks
    monkeypatch.setattr(module, "_DOWNLOAD_MAX_ATTEMPTS", 2)
    monkeypatch.setattr(module.asyncio, "sleep", backoff)
    module.run_download_files("http://example.test", [filename], tmp_path)

    assert attempts == 2
    assert (tmp_path / filename).read_bytes() == payload
    assert not list(tmp_path.glob("*.part"))
