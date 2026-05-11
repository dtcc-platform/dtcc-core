import json

import pytest

from dtcc_core.datasets.publish import (
    DatasetPackageError,
    DatasetPublication,
    DatasetPublishConfigurationError,
    DatasetUploadClient,
    PublishedFile,
    build_publish_idempotency_key,
)


def _write_package(tmp_path):
    file_path = tmp_path / "smoke_slice.geojson"
    manifest_path = tmp_path / "smoke_slice.manifest.json"
    manifest = {
        "name": "smoke_slice",
        "file": "smoke_slice.geojson",
        "format": "geojson",
        "media_type": "application/geo+json",
        "data_kind": "smoke",
    }
    file_path.write_text('{"type":"FeatureCollection","features":[]}', encoding="utf-8")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return manifest_path, (file_path,), manifest


def test_publication_from_response_preserves_known_fields():
    response = {
        "dataset_key": 42,
        "version_id": 123,
        "version_number": "7",
        "status": "published",
        "owner": "dtcc",
        "manifest_sha256": "manifest-sha",
        "file_set_sha256": "files-sha",
        "manifest_url": "https://upload.example/manifests/ver_123",
        "files": [
            {
                "path": "smoke_slice.geojson",
                "original_filename": None,
                "size": "42",
                "sha256": "file-sha",
                "media_type": "application/geo+json",
                "sniffed_media_type": "application/json",
            }
        ],
    }

    publication = DatasetPublication.from_response(
        response, upload_url="https://upload.example/v1/datasets"
    )

    assert publication == DatasetPublication(
        dataset_key="42",
        version_id="123",
        version_number=7,
        status="published",
        owner="dtcc",
        manifest_sha256="manifest-sha",
        file_set_sha256="files-sha",
        files=(
            PublishedFile(
                path="smoke_slice.geojson",
                original_filename=None,
                size=42,
                sha256="file-sha",
                media_type="application/geo+json",
                sniffed_media_type="application/json",
                raw=response["files"][0],
            ),
        ),
        upload_url="https://upload.example/v1/datasets",
        manifest_url="https://upload.example/manifests/ver_123",
        raw=response,
    )


def test_publication_dataclass_optional_defaults_are_independent():
    first = DatasetPublication(
        dataset_key="smoke",
        version_id="ver_123",
        version_number=1,
        status="published",
        owner="dtcc",
        manifest_sha256="manifest-sha",
        file_set_sha256="files-sha",
        files=(),
    )
    second = DatasetPublication(
        dataset_key="smoke",
        version_id="ver_124",
        version_number=2,
        status="published",
        owner="dtcc",
        manifest_sha256="manifest-sha-2",
        file_set_sha256="files-sha-2",
        files=(),
    )

    assert first.upload_url is None
    assert first.manifest_url is None
    assert first.raw == {}
    assert first.raw is not second.raw


@pytest.mark.parametrize(
    ("upload_url", "datasets_url"),
    [
        ("https://upload.example", "https://upload.example/v1/datasets"),
        ("https://upload.example/", "https://upload.example/v1/datasets"),
        ("https://upload.example/v1", "https://upload.example/v1/datasets"),
        ("https://upload.example/v1/", "https://upload.example/v1/datasets"),
    ],
)
def test_client_normalizes_dataset_endpoint_urls(upload_url, datasets_url):
    client = DatasetUploadClient(upload_url, token="secret")

    assert client.datasets_url == datasets_url


def test_client_constructor_preserves_configuration():
    session = object()

    client = DatasetUploadClient(
        "https://upload.example/v1",
        token="secret",
        timeout=12.5,
        session=session,
    )

    assert client.base_url == "https://upload.example/v1"
    assert client.datasets_url == "https://upload.example/v1/datasets"
    assert client.token == "secret"
    assert client.timeout == 12.5
    assert client.session is session


def test_client_constructor_accepts_positional_token():
    client = DatasetUploadClient("https://upload.example/", "secret")

    assert client.base_url == "https://upload.example"
    assert client.datasets_url == "https://upload.example/v1/datasets"
    assert client.token == "secret"


def test_from_config_uses_explicit_url_and_token(monkeypatch):
    monkeypatch.delenv("DTCC_UPLOAD_URL", raising=False)
    monkeypatch.delenv("DTCC_UPLOAD_TOKEN", raising=False)

    client = DatasetUploadClient.from_config(
        upload_url="https://upload.example", token="secret"
    )

    assert client.datasets_url == "https://upload.example/v1/datasets"
    assert client.token == "secret"


def test_from_config_uses_environment_url_and_token(monkeypatch):
    monkeypatch.setenv("DTCC_UPLOAD_URL", "https://upload.example")
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "env-secret")

    client = DatasetUploadClient.from_config()

    assert client.datasets_url == "https://upload.example/v1/datasets"
    assert client.token == "env-secret"


def test_from_config_explicit_values_override_environment(monkeypatch):
    monkeypatch.setenv("DTCC_UPLOAD_URL", "https://env.example")
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "env-secret")

    client = DatasetUploadClient.from_config(
        upload_url="https://explicit.example", token="explicit-secret"
    )

    assert client.datasets_url == "https://explicit.example/v1/datasets"
    assert client.token == "explicit-secret"


def test_from_config_allows_mixed_explicit_and_environment_values(monkeypatch):
    monkeypatch.setenv("DTCC_UPLOAD_URL", "https://env.example")
    monkeypatch.delenv("DTCC_UPLOAD_TOKEN", raising=False)

    client = DatasetUploadClient.from_config(token="explicit-secret")

    assert client.datasets_url == "https://env.example/v1/datasets"
    assert client.token == "explicit-secret"

    monkeypatch.delenv("DTCC_UPLOAD_URL", raising=False)
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "env-secret")

    client = DatasetUploadClient.from_config(upload_url="https://explicit.example")

    assert client.datasets_url == "https://explicit.example/v1/datasets"
    assert client.token == "env-secret"


def test_from_config_missing_url_or_token_raises_configuration_error(monkeypatch):
    monkeypatch.delenv("DTCC_UPLOAD_URL", raising=False)
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "env-secret")

    with pytest.raises(DatasetPublishConfigurationError) as url_error:
        DatasetUploadClient.from_config()

    assert "DTCC_UPLOAD_URL" in str(url_error.value)
    assert url_error.value.failure_class == "configuration"

    monkeypatch.setenv("DTCC_UPLOAD_URL", "https://upload.example")
    monkeypatch.delenv("DTCC_UPLOAD_TOKEN", raising=False)

    with pytest.raises(DatasetPublishConfigurationError) as token_error:
        DatasetUploadClient.from_config()

    assert "DTCC_UPLOAD_TOKEN" in str(token_error.value)
    assert token_error.value.failure_class == "configuration"


def test_from_config_explicit_empty_values_do_not_fall_back_to_environment(monkeypatch):
    monkeypatch.setenv("DTCC_UPLOAD_URL", "https://env.example")
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "env-secret")

    for empty_url in ("", "   "):
        with pytest.raises(DatasetPublishConfigurationError) as url_error:
            DatasetUploadClient.from_config(
                upload_url=empty_url, token="explicit-secret"
            )

        assert "DTCC_UPLOAD_URL" in str(url_error.value)
        assert url_error.value.failure_class == "configuration"

    for empty_token in ("", "   "):
        with pytest.raises(DatasetPublishConfigurationError) as token_error:
            DatasetUploadClient.from_config(
                upload_url="https://explicit.example", token=empty_token
            )

        assert "DTCC_UPLOAD_TOKEN" in str(token_error.value)
        assert token_error.value.failure_class == "configuration"


def test_from_config_whitespace_environment_values_raise_configuration_error(
    monkeypatch,
):
    monkeypatch.setenv("DTCC_UPLOAD_URL", "   ")
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "env-secret")

    with pytest.raises(DatasetPublishConfigurationError) as url_error:
        DatasetUploadClient.from_config()

    assert "DTCC_UPLOAD_URL" in str(url_error.value)
    assert url_error.value.failure_class == "configuration"

    monkeypatch.setenv("DTCC_UPLOAD_URL", "https://upload.example")
    monkeypatch.setenv("DTCC_UPLOAD_TOKEN", "   ")

    with pytest.raises(DatasetPublishConfigurationError) as token_error:
        DatasetUploadClient.from_config()

    assert "DTCC_UPLOAD_TOKEN" in str(token_error.value)
    assert token_error.value.failure_class == "configuration"


def test_build_publish_idempotency_key_is_deterministic(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    manifest_sha256 = "2b74c25b831f55cefe8f159321450bd6a85f940d4cd14dd3b99ae7e35f0a1e70"

    first_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)
    second_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert first_key == second_key
    assert first_key.startswith("dtcc-publish-v1:")
    assert manifest_sha256 not in first_key


def test_build_publish_idempotency_key_changes_when_dataset_key_changes(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)

    first_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)
    second_key = build_publish_idempotency_key("other", manifest_path, files, manifest)

    assert first_key != second_key


def test_build_publish_idempotency_key_changes_when_manifest_content_changes(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)

    first_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)
    manifest["data_kind"] = "air_quality"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    second_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert first_key != second_key


def test_build_publish_idempotency_key_changes_when_file_changes(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)

    first_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)
    files[0].write_text('{"type":"FeatureCollection","features":[{}]}', encoding="utf-8")
    second_key = build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert first_key != second_key


def test_package_validation_rejects_missing_manifest(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    manifest_path.unlink()

    with pytest.raises(DatasetPackageError):
        build_publish_idempotency_key("smoke", manifest_path, files, manifest)


def test_package_validation_rejects_manifest_directory(tmp_path):
    _, files, manifest = _write_package(tmp_path)
    manifest_dir = tmp_path / "manifest_dir"
    manifest_dir.mkdir()

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_dir, files, manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_missing_package_file(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    files[0].unlink()

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_package_file_directory(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    files[0].unlink()
    package_dir = tmp_path / "smoke_slice.geojson"
    package_dir.mkdir(exist_ok=True)

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, (package_dir,), manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_multiple_package_files(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    extra_file = tmp_path / "extra.geojson"
    extra_file.write_text("{}", encoding="utf-8")

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key(
            "smoke", manifest_path, (files[0], extra_file), manifest
        )

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_missing_manifest_file_field(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    manifest.pop("file")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_non_string_manifest_file_field(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    manifest["file"] = ["smoke_slice.geojson"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_non_string_manifest_mapping_file(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    provided_manifest = dict(manifest)
    provided_manifest["file"] = ["smoke_slice.geojson"]

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, files, provided_manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_filename_mismatch(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    mismatched_file = tmp_path / "different.geojson"
    mismatched_file.write_text(files[0].read_text(encoding="utf-8"), encoding="utf-8")

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key(
            "smoke", manifest_path, (mismatched_file,), manifest
        )

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_manifest_mapping_mismatch(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    manifest["file"] = "different.geojson"

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, files, manifest)

    assert error.value.failure_class == "invalid_package"


def test_package_validation_rejects_manifest_metadata_mismatch(tmp_path):
    manifest_path, files, manifest = _write_package(tmp_path)
    provided_manifest = dict(manifest)
    provided_manifest["data_kind"] = "air_quality"

    with pytest.raises(DatasetPackageError) as error:
        build_publish_idempotency_key("smoke", manifest_path, files, provided_manifest)

    assert error.value.failure_class == "invalid_package"
