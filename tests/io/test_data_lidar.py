from pathlib import Path
from unittest.mock import Mock

import pytest

from dtcc_core.io.data import lidar


def test_download_lidar_uses_stable_filename_order(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(lidar, "user_cache_dir", lambda appname: str(tmp_path))
    monkeypatch.setattr(
        lidar,
        "post_lidar_request",
        lambda *args, **kwargs: {
            "tiles": [
                {"filename": "tile_b.laz"},
                {"filename": "tile_a.laz"},
            ]
        },
    )
    plot_bboxes = Mock()
    monkeypatch.setattr(lidar, "plot_bboxes_folium", plot_bboxes)
    captured = {}

    def fake_run(base_url: str, filenames: list[str], output_dir: str) -> None:
        captured["base_url"] = base_url
        captured["filenames"] = filenames
        captured["output_dir"] = output_dir

    monkeypatch.setattr(lidar, "run_download_files", fake_run)

    files = lidar.download_lidar(
        (0, 0, 10, 10), Mock(), base_url="http://example.test"
    )

    assert captured == {
        "base_url": "http://example.test",
        "filenames": ["tile_a.laz", "tile_b.laz"],
        "output_dir": str(tmp_path / "downloaded_laz"),
    }
    assert files == [
        str(tmp_path / "downloaded_laz" / "tile_a.laz"),
        str(tmp_path / "downloaded_laz" / "tile_b.laz"),
    ]
    plot_args, _plot_kwargs = plot_bboxes.call_args
    assert [tile["filename"] for tile in plot_args[1]] == [
        "tile_a.laz",
        "tile_b.laz",
    ]


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
