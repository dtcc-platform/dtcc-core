"""Source downloads accept live JSON responses without historical byte pins."""

import gzip
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest


@pytest.fixture
def downloader(tmp_path, monkeypatch):
    script = Path(__file__).resolve().parents[1] / 'scripts/download_flagship_sources.py'
    spec = importlib.util.spec_from_file_location('download_flagship_sources', script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, '__file__', str(tmp_path / script.name))
    monkeypatch.setattr(sys, 'argv', [str(script), '--output', str(tmp_path / 'cache')])
    return module


def test_live_sources_download_and_reuse_cache(downloader, tmp_path, monkeypatch):
    water = b'{"type":"FeatureCollection","timeStamp":"2026-09-24T12:00:00Z","features":[]}'
    tile = gzip.compress(b'{"type":"CityJSON"}')
    responses = {'https://example.test/water': water, 'https://example.test/tile': tile}
    (tmp_path / 'flagship-sources.json').write_text(json.dumps({'sources': [
        {'file': 'water.geojson', 'url': 'https://example.test/water'},
        {'file': 'tile.city.json.gz', 'url': 'https://example.test/tile'},
    ]}))
    monkeypatch.setattr(downloader, 'urlopen', lambda url, timeout: io.BytesIO(responses[url]))

    downloader.main()

    assert (tmp_path / 'cache/water.geojson').read_bytes() == water
    assert (tmp_path / 'cache/tile.city.json.gz').read_bytes() == tile

    def unexpected_download(*args, **kwargs):
        pytest.fail('Existing cached sources should not be downloaded again')

    monkeypatch.setattr(downloader, 'urlopen', unexpected_download)
    downloader.main()


def test_malformed_download_is_not_cached(downloader, tmp_path, monkeypatch):
    (tmp_path / 'flagship-sources.json').write_text(json.dumps({'sources': [
        {'file': 'water.geojson', 'url': 'https://example.test/water'},
    ]}))
    monkeypatch.setattr(downloader, 'urlopen', lambda url, timeout: io.BytesIO(b'<html>Error</html>'))

    with pytest.raises(json.JSONDecodeError):
        downloader.main()

    assert list((tmp_path / 'cache').iterdir()) == []
