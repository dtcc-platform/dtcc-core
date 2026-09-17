"""Cache clearing must never remove files outside its configured root."""

from dtcc_core.io.data import cache, empty_cache


def test_empty_cache_clears_files_and_directories_but_preserves_external_target(tmp_path, monkeypatch):
    root = tmp_path / "cache"
    nested = root / "nested"
    nested.mkdir(parents=True)
    (nested / "tile.bin").write_bytes(b"cached tile")
    (root / "index.json").write_text("{}")
    outside = tmp_path / "keep.txt"
    outside.write_text("user data")
    link = root / "external-link"
    try:
        link.symlink_to(outside)
    except OSError:
        # Symlink creation can require additional privileges on Windows.
        link = None
    monkeypatch.setattr(cache, "cache_dir", root)
    empty_cache()
    empty_cache()  # Repeating the operation is harmless.
    assert root.is_dir()
    assert not nested.exists() and not (root / "index.json").exists()
    assert outside.read_text() == "user data"
    if link is not None:
        assert link.is_symlink()
