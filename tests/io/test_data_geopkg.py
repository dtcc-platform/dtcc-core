from pathlib import Path
from unittest.mock import Mock

import pytest

from dtcc_core.io.data import geopkg


def test_download_tiles_raises_when_download_does_not_create_expected_file(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(geopkg, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        geopkg,
        "post_gpkg_request",
        lambda *args, **kwargs: {"tiles": ["tile_missing.gpkg"]},
    )
    monkeypatch.setattr(geopkg, "run_download_files", lambda *args, **kwargs: None)

    with pytest.raises(geopkg.FootprintDownloadError, match="did not produce"):
        geopkg.download_tiles(
            (0, 0, 10, 10), Mock(), server_url="http://example.test"
        )


def test_download_tiles_returns_existing_downloaded_files(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(geopkg, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        geopkg,
        "post_gpkg_request",
        lambda *args, **kwargs: {"tiles": ["tile_ok.gpkg"]},
    )

    def fake_run(_url: str, _filenames: list[str], output_dir: str) -> None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        (Path(output_dir) / "tile_ok.gpkg").write_bytes(b"data")

    monkeypatch.setattr(geopkg, "run_download_files", fake_run)

    files = geopkg.download_tiles(
        (0, 0, 10, 10), Mock(), server_url="http://example.test"
    )

    assert files == [str(tmp_path / "downloaded-gpkg" / "tile_ok.gpkg")]


def test_run_download_files_reports_all_cached_and_skips_downloader(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    (tmp_path / "tile_ok.gpkg").write_bytes(b"data")
    info = Mock()
    run = Mock()
    monkeypatch.setattr(geopkg, "info", info)
    monkeypatch.setattr(geopkg.asyncio, "run", run)

    geopkg.run_download_files(
        "http://example.test",
        ["tile_ok.gpkg"],
        output_dir=str(tmp_path),
    )

    info.assert_called_once_with(
        "Using 1 cached footprint tile file(s) from local cache"
    )
    run.assert_not_called()


def test_run_download_files_reports_partial_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    (tmp_path / "tile_cached.gpkg").write_bytes(b"data")
    info = Mock()
    captured = {}
    monkeypatch.setattr(geopkg, "info", info)

    async def fake_download_all(base_url, filenames, output_dir):
        captured["base_url"] = base_url
        captured["filenames"] = filenames
        captured["output_dir"] = output_dir

    monkeypatch.setattr(geopkg, "download_all_gpkg_files", fake_download_all)

    geopkg.run_download_files(
        "http://example.test",
        ["tile_cached.gpkg", "tile_missing.gpkg"],
        output_dir=str(tmp_path),
    )

    info.assert_called_once_with(
        "Using 1 cached footprint tile file(s); downloading 1 missing tile file(s)"
    )
    assert captured == {
        "base_url": "http://example.test",
        "filenames": ["tile_missing.gpkg"],
        "output_dir": str(tmp_path),
    }
