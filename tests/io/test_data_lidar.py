from pathlib import Path
from unittest.mock import Mock

import pytest

from dtcc_core.io.data import lidar


def test_run_download_files_reports_all_cached_and_skips_downloader(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    (tmp_path / "tile_ok.laz").write_bytes(b"data")
    info = Mock()
    run = Mock()
    monkeypatch.setattr(lidar, "info", info)
    monkeypatch.setattr(lidar.asyncio, "run", run)

    lidar.run_download_files(
        "http://example.test",
        ["tile_ok.laz"],
        output_dir=str(tmp_path),
    )

    info.assert_called_once_with(
        "Using 1 cached lidar tile file(s) from local cache"
    )
    run.assert_not_called()


def test_run_download_files_reports_partial_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    (tmp_path / "tile_cached.laz").write_bytes(b"data")
    info = Mock()
    captured = {}
    monkeypatch.setattr(lidar, "info", info)

    async def fake_download_all(base_url, filenames, output_dir):
        captured["base_url"] = base_url
        captured["filenames"] = filenames
        captured["output_dir"] = output_dir

    monkeypatch.setattr(lidar, "download_all_lidar_files", fake_download_all)

    lidar.run_download_files(
        "http://example.test",
        ["tile_cached.laz", "tile_missing.laz"],
        output_dir=str(tmp_path),
    )

    info.assert_called_once_with(
        "Using 1 cached lidar tile file(s); downloading 1 missing tile file(s)"
    )
    assert captured == {
        "base_url": "http://example.test",
        "filenames": ["tile_missing.laz"],
        "output_dir": str(tmp_path),
    }
