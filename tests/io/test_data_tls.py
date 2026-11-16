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

CA_SETTINGS = ["SSL_CERT_FILE", "SSL_CERT_DIR", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"]


def _bundle_certs():
    text = Path(requests.certs.where()).read_text()
    return re.findall(
        r"-----BEGIN CERTIFICATE-----.*?-----END CERTIFICATE-----", text, re.S
    )


@pytest.fixture
def no_ca_settings(monkeypatch):
    for name in CA_SETTINGS:
        monkeypatch.delenv(name, raising=False)


@pytest.fixture
def one_ca(tmp_path):
    path = tmp_path / "one.pem"
    path.write_text(_bundle_certs()[0] + "\n")
    return str(path)


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


def _bundles_loaded(load):
    """The cafile and capath loads; Windows also loads its store as cadata."""
    return [c for c in load.call_args_list if "cadata" not in c.kwargs]


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


@pytest.mark.parametrize("setting", ["REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"])
def test_lookup_bundle_setting_replaces_certifi(
    monkeypatch, no_ca_settings, one_ca, setting
):
    monkeypatch.setenv(setting, one_ca)
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _bundles_loaded(load) == [call(cafile=one_ca)]


def test_lookup_bundle_directory_loads_as_capath(monkeypatch, tmp_path, no_ca_settings):
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", str(tmp_path))
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _bundles_loaded(load) == [call(capath=str(tmp_path))]


def test_lookup_bundle_settings_follow_the_requests_order(
    monkeypatch, tmp_path, no_ca_settings
):
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", str(tmp_path / "a.pem"))
    monkeypatch.setenv("CURL_CA_BUNDLE", str(tmp_path / "b.pem"))
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _bundles_loaded(load) == [call(cafile=str(tmp_path / "a.pem"))]


def test_empty_lookup_bundle_setting_falls_through(
    monkeypatch, tmp_path, no_ca_settings
):
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", "")
    monkeypatch.setenv("CURL_CA_BUNDLE", str(tmp_path / "b.pem"))
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _bundles_loaded(load) == [call(cafile=str(tmp_path / "b.pem"))]


def test_empty_lookup_bundle_settings_keep_certifi(monkeypatch, no_ca_settings):
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", "")
    monkeypatch.setenv("CURL_CA_BUNDLE", "")
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _bundles_loaded(load) == [call(cafile=requests.certs.where())]


@pytest.mark.parametrize("setting", ["SSL_CERT_FILE", "SSL_CERT_DIR"])
def test_lookup_bundle_loads_alongside_ca_settings(
    monkeypatch, no_ca_settings, one_ca, setting
):
    monkeypatch.setenv(setting, "/nonexistent")
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", one_ca)
    with patch.object(ssl.SSLContext, "load_verify_locations") as load:
        download_ssl_context()
    assert _bundles_loaded(load) == [call(cafile=one_ca)]


def test_lookup_bundle_is_trusted(monkeypatch, no_ca_settings, one_ca):
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", one_ca)
    with patch.object(ssl.SSLContext, "load_default_certs"):
        context = download_ssl_context()
    assert context.cert_store_stats()["x509_ca"] == 1
    assert context.verify_mode == ssl.CERT_REQUIRED
    assert context.check_hostname is True


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
