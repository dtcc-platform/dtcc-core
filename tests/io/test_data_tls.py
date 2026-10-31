import asyncio
import re
import ssl
from pathlib import Path
from unittest.mock import call, patch

import aiohttp
import pytest
import requests

from dtcc_core.io.data import geopkg, lidar
from dtcc_core.io.data.tls import download_ssl_context


def _bundle_certs():
    text = Path(requests.certs.where()).read_text()
    return re.findall(
        r"-----BEGIN CERTIFICATE-----.*?-----END CERTIFICATE-----", text, re.S
    )


@pytest.fixture
def no_ca_settings(monkeypatch):
    monkeypatch.delenv("SSL_CERT_FILE", raising=False)
    monkeypatch.delenv("SSL_CERT_DIR", raising=False)


def test_context_verifies(no_ca_settings):
    context = download_ssl_context()
    assert context.verify_mode == ssl.CERT_REQUIRED
    assert context.check_hostname is True


def test_context_trusts_the_lookup_bundle_without_a_system_store(no_ca_settings):
    with patch.object(ssl.SSLContext, "load_default_certs"):
        context = download_ssl_context()
    assert context.cert_store_stats()["x509_ca"] > 0


def _loaded_bundle(load):
    return call(cafile=requests.certs.where()) in load.call_args_list


@pytest.mark.parametrize("setting", ["SSL_CERT_FILE", "SSL_CERT_DIR"])
def test_ca_settings_keep_the_bundle_out(monkeypatch, no_ca_settings, setting):
    monkeypatch.setenv(setting, "/nonexistent")
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert not _loaded_bundle(load)


def test_empty_ca_setting_keeps_the_bundle(monkeypatch, no_ca_settings):
    monkeypatch.setenv("SSL_CERT_FILE", "")
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _loaded_bundle(load)


def test_ssl_cert_file_trust_matches_the_default_context(
    monkeypatch, tmp_path, no_ca_settings
):
    one_ca = tmp_path / "one.pem"
    one_ca.write_text(_bundle_certs()[0] + "\n")
    monkeypatch.setenv("SSL_CERT_FILE", str(one_ca))
    default = ssl.create_default_context().cert_store_stats()
    assert download_ssl_context().cert_store_stats() == default


@pytest.mark.parametrize(
    ("module", "run"),
    [
        (geopkg, lambda d: geopkg.download_all_gpkg_files("http://x.test", [], d)),
        (lidar, lambda d: lidar.download_all_lidar_files("http://x.test", [], d)),
    ],
)
def test_downloads_use_the_context(module, run, tmp_path):
    context = ssl.create_default_context()
    with (
        patch.object(module, "download_ssl_context", return_value=context),
        patch.object(
            module.aiohttp, "TCPConnector", wraps=aiohttp.TCPConnector
        ) as connector,
    ):
        asyncio.run(run(str(tmp_path)))
    assert connector.call_args.kwargs["ssl"] is context
