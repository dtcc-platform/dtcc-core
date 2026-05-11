from __future__ import annotations

from pathlib import Path

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets.publish import DatasetPackageError, DatasetPublication


class RecordingUploader:
    def __init__(self):
        self.calls = []

    def upload_package(
        self,
        dataset_key,
        manifest_path,
        files,
        *,
        manifest=None,
        idempotency_key=None,
    ):
        call = {
            "dataset_key": dataset_key,
            "manifest_path": Path(manifest_path),
            "files": tuple(Path(file) for file in files),
            "manifest": manifest,
            "idempotency_key": idempotency_key,
            "manifest_exists_during_upload": Path(manifest_path).is_file(),
            "files_exist_during_upload": tuple(Path(file).is_file() for file in files),
        }
        self.calls.append(call)
        return DatasetPublication.from_response(
            {
                "dataset_key": dataset_key,
                "version_id": "version-1",
                "version_number": 1,
                "status": "published",
                "owner": "test-owner",
                "manifest_sha256": "manifest-sha",
                "file_set_sha256": "files-sha",
                "files": [
                    {
                        "path": call["manifest"]["file"],
                        "original_filename": call["manifest"]["file"],
                        "size": call["files"][0].stat().st_size,
                        "sha256": "file-sha",
                        "media_type": call["manifest"]["media_type"],
                        "sniffed_media_type": call["manifest"]["media_type"],
                    }
                ],
            }
        )


def test_dataset_descriptor_publish_exports_and_uploads_with_temp_cleanup():
    uploader = RecordingUploader()

    publication = datasets.smoke.publish(
        dataset_key="smoke-slice",
        bounds=[0, 0, 1, 1],
        product="slice",
        resolution=4,
        format="geojson",
        uploader=uploader,
    )

    assert publication.dataset_key == "smoke-slice"
    assert len(uploader.calls) == 1
    call = uploader.calls[0]
    assert call["dataset_key"] == "smoke-slice"
    assert call["manifest_exists_during_upload"] is True
    assert call["files_exist_during_upload"] == (True,)
    assert call["manifest"]["name"] == "smoke"
    assert call["manifest"]["file"] == "smoke_slice.geojson"
    assert call["manifest"]["format"] == "geojson"
    assert call["manifest"]["product"] == "slice"
    assert call["manifest"]["parameters"]["product"] == "slice"
    assert call["manifest"]["parameters"]["resolution"] == 4
    assert call["idempotency_key"] is None
    assert not call["files"][0].exists()
    assert not call["manifest_path"].exists()


def test_dataset_descriptor_publish_keep_export_preserves_files(tmp_path):
    uploader = RecordingUploader()

    datasets.smoke.publish(
        dataset_key="smoke-slice",
        bounds=[0, 0, 1, 1],
        product="slice",
        resolution=4,
        format="geojson",
        output_dir=tmp_path,
        keep_export=True,
        uploader=uploader,
    )

    call = uploader.calls[0]
    assert call["files"][0] == tmp_path / "smoke_slice.geojson"
    assert call["manifest_path"] == tmp_path / "smoke_slice.manifest.json"
    assert call["files"][0].is_file()
    assert call["manifest_path"].is_file()
    assert call["manifest"]["file"] == "smoke_slice.geojson"


def test_dataset_descriptor_publish_uses_explicit_filename(tmp_path):
    uploader = RecordingUploader()

    datasets.smoke.publish(
        dataset_key="smoke-slice",
        bounds=[0, 0, 1, 1],
        product="slice",
        resolution=4,
        filename="custom.geojson",
        output_dir=tmp_path,
        keep_export=True,
        uploader=uploader,
    )

    call = uploader.calls[0]
    assert call["files"][0].name == "custom.geojson"
    assert call["manifest_path"].name == "custom.manifest.json"
    assert call["manifest"]["file"] == "custom.geojson"


@pytest.mark.parametrize(
    "filename",
    [
        "absolute",
        "../escape.geojson",
        "nested/file.geojson",
        "nested\\file.geojson",
        ".hidden.geojson",
    ],
)
def test_dataset_descriptor_publish_rejects_unsafe_explicit_filename(
    tmp_path, filename
):
    bad_filename = tmp_path / "absolute.geojson" if filename == "absolute" else filename
    uploader = RecordingUploader()

    with pytest.raises(ValueError, match="filename|safe|unsafe"):
        datasets.smoke.publish(
            dataset_key="smoke-slice",
            bounds=[0, 0, 1, 1],
            product="slice",
            resolution=4,
            filename=bad_filename,
            format="geojson",
            output_dir=tmp_path / "exports",
            uploader=uploader,
        )

    assert uploader.calls == []
    if isinstance(bad_filename, Path):
        assert not bad_filename.exists()
    if filename == "../escape.geojson":
        assert not (tmp_path / "escape.geojson").exists()


def test_dataset_descriptor_publish_passes_manifest_overrides(tmp_path):
    uploader = RecordingUploader()

    datasets.smoke.publish(
        dataset_key="smoke-slice",
        bounds=[0, 0, 1, 1],
        product="slice",
        resolution=4,
        format="geojson",
        output_dir=tmp_path,
        keep_export=True,
        manifest_id="manifest-1",
        title="Smoke Slice",
        description="A tiny smoke slice.",
        uploader=uploader,
    )

    manifest = uploader.calls[0]["manifest"]
    assert manifest["id"] == "manifest-1"
    assert manifest["title"] == "Smoke Slice"
    assert manifest["description"] == "A tiny smoke slice."


def test_dataset_descriptor_publish_rejects_missing_format_without_filename():
    with pytest.raises(ValueError, match="format"):
        datasets.smoke.publish(
            dataset_key="smoke-slice",
            bounds=[0, 0, 1, 1],
            product="slice",
            resolution=4,
            uploader=RecordingUploader(),
        )


def test_dataset_descriptor_publish_rejects_multifile_format():
    uploader = RecordingUploader()

    with pytest.raises(DatasetPackageError, match="multi-file"):
        datasets.city_volume_mesh.publish(
            dataset_key="volume-mesh",
            bounds=[0, 0, 1, 1],
            format="xdmf",
            uploader=uploader,
        )

    assert uploader.calls == []
