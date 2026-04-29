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
