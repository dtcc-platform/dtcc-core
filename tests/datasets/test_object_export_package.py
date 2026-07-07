"""Tests for object-first Dataset Manifest v2 package export."""

from __future__ import annotations

import hashlib
import json
import zipfile

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import DatasetPackage, attach_dataset_context
from dtcc_core.datasets.publish import DatasetPublishConfigurationError
from dtcc_core.model import (
    Building,
    City,
    FootprintCollection,
    GeometryType,
    Mesh,
    Surface,
)


class RecordingUploader:
    def __init__(self):
        self.calls = []

    def upload_package(
        self,
        *,
        dataset_key,
        manifest_path,
        files,
        manifest,
        idempotency_key=None,
    ):
        call = {
            "dataset_key": dataset_key,
            "manifest_path": manifest_path,
            "files": tuple(files),
            "manifest": manifest,
            "idempotency_key": idempotency_key,
            "manifest_exists_during_upload": manifest_path.is_file(),
            "files_exist_during_upload": tuple(path.is_file() for path in files),
        }
        self.calls.append(call)
        return {"published": dataset_key}


def test_city_object_export_writes_manifest_and_artifact(tmp_path):
    city = City()
    attach_dataset_context(city, _context(datasets.city))
    before_export = city.manifest().model_dump(mode="json")

    package = city.export(tmp_path / "city_pkg", format="json")

    manifest_path = tmp_path / "city_pkg" / "manifest.json"
    artifact_path = tmp_path / "city_pkg" / "artifacts" / "city.json"
    manifest = json.loads(manifest_path.read_text())
    artifact = manifest["artifacts"][0]

    assert isinstance(package, DatasetPackage)
    assert package.path == tmp_path / "city_pkg"
    assert package.manifest_path == manifest_path
    assert package.package_format == "directory"
    assert package.artifacts == tuple(package.manifest.artifacts)
    assert package.files == (manifest_path, artifact_path)
    assert artifact_path.exists()

    assert manifest["schema_version"] == "dtcc-dataset-manifest-v2"
    assert manifest["identity"]["name"] == "city"
    assert manifest["request"]["dataset_name"] == "city"
    assert artifact["path"].startswith("artifacts/")
    assert artifact["path"] == "artifacts/city.json"
    assert artifact["role"] == "primary"
    assert artifact["format"] == "json"
    assert artifact["media_type"] == "application/json"
    assert artifact["bounds"] == [0.0, 0.0, 1.0, 1.0]
    assert artifact["size"] == artifact_path.stat().st_size
    assert artifact["sha256"] == _sha256(artifact_path)

    assert city.manifest().model_dump(mode="json") == before_export
    assert city.manifest().artifacts == []


def test_object_export_sanitizes_object_artifact_names(tmp_path):
    city = City()
    city.name = "evil/../City Name"
    attach_dataset_context(city, _context(datasets.city))

    package = city.export(tmp_path / "unsafe_name_pkg", format="json")

    artifact_path = tmp_path / "unsafe_name_pkg" / "artifacts" / "evil_city_name.json"
    artifact = package.manifest.artifacts[0]
    assert artifact_path.exists()
    assert package.manifest.identity.name == "city"
    assert artifact.path == "artifacts/evil_city_name.json"
    assert ".." not in artifact.path
    assert artifact.size == artifact_path.stat().st_size
    assert artifact.sha256 == _sha256(artifact_path)


def test_object_export_without_dataset_context_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="no DatasetContext"):
        City().export(tmp_path / "city_pkg", format="json")


def test_object_publish_exports_manifest_v2_package_and_cleans_temp_files():
    city = City()
    attach_dataset_context(city, _context(datasets.city))
    uploader = RecordingUploader()

    publication = city.publish(
        dataset_key="city-v2",
        format="json",
        uploader=uploader,
        idempotency_key="publish-1",
    )

    assert publication == {"published": "city-v2"}
    assert len(uploader.calls) == 1
    call = uploader.calls[0]
    artifact = call["manifest"]["artifacts"][0]
    assert call["dataset_key"] == "city-v2"
    assert call["idempotency_key"] == "publish-1"
    assert call["manifest"]["schema_version"] == "dtcc-dataset-manifest-v2"
    assert call["manifest"]["identity"]["name"] == "city"
    assert artifact["path"] == "artifacts/city.json"
    assert call["manifest_exists_during_upload"] is True
    assert call["files_exist_during_upload"] == (True,)
    assert call["files"] == (call["manifest_path"].parent / artifact["path"],)
    assert not call["manifest_path"].exists()
    assert not call["files"][0].exists()


def test_object_publish_without_dataset_context_is_rejected():
    with pytest.raises(ValueError, match="no DatasetContext"):
        City().publish(dataset_key="city-v2", uploader=RecordingUploader())


def test_object_publish_requires_upload_configuration_when_no_uploader(monkeypatch):
    city = City()
    attach_dataset_context(city, _context(datasets.city))
    monkeypatch.delenv("DTCC_UPLOAD_URL", raising=False)
    monkeypatch.delenv("DTCC_UPLOAD_TOKEN", raising=False)

    with pytest.raises(DatasetPublishConfigurationError, match="DTCC_UPLOAD_URL"):
        city.publish(dataset_key="city-v2", format="json")


def test_directory_dataset_package_publish_delegates_manifest_v2_files(tmp_path):
    city = City()
    attach_dataset_context(city, _context(datasets.city))
    package = city.export(tmp_path / "city_pkg", format="json")
    uploader = RecordingUploader()

    publication = package.publish(
        dataset_key="city-v2",
        uploader=uploader,
        idempotency_key="publish-2",
    )

    assert publication == {"published": "city-v2"}
    call = uploader.calls[0]
    assert call["manifest_path"] == package.manifest_path
    assert call["files"] == (tmp_path / "city_pkg" / "artifacts" / "city.json",)
    assert call["manifest"] == package.manifest.model_dump(mode="json")
    assert call["idempotency_key"] == "publish-2"


def test_archive_dataset_package_publish_uses_archive_contents(tmp_path):
    city = City()
    attach_dataset_context(city, _context(datasets.city))
    package = city.export(tmp_path / "city.dtccpkg", format="json")
    uploader = RecordingUploader()

    publication = package.publish(dataset_key="city-v2", uploader=uploader)

    assert publication == {"published": "city-v2"}
    call = uploader.calls[0]
    artifact = call["manifest"]["artifacts"][0]
    assert call["manifest"]["schema_version"] == "dtcc-dataset-manifest-v2"
    assert artifact["path"] == "artifacts/city.json"
    assert call["manifest_exists_during_upload"] is True
    assert call["files_exist_during_upload"] == (True,)
    assert call["manifest_path"].name == "manifest.json"
    assert call["files"][0].name == "city.json"
    assert not call["manifest_path"].exists()
    assert not call["files"][0].exists()


def test_mesh_like_object_export_uses_object_serializer(tmp_path):
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            dtype=float,
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
    )
    attach_dataset_context(mesh, _context(datasets.city_surface_mesh))

    package = mesh.export(tmp_path / "mesh_pkg", format="pb")

    artifact_path = tmp_path / "mesh_pkg" / "artifacts" / "city_surface_mesh.pb"
    artifact = package.manifest.artifacts[0]
    assert artifact_path.exists()
    assert artifact.path == "artifacts/city_surface_mesh.pb"
    assert artifact.role == "primary"
    assert artifact.format == "pb"
    assert artifact.size == artifact_path.stat().st_size
    assert artifact.sha256 == _sha256(artifact_path)


def test_footprint_collection_export_defaults_to_geojson(tmp_path):
    building = _building()
    collection = FootprintCollection.from_buildings([building])
    attach_dataset_context(collection, _context(datasets.building_footprints))

    package = collection.export(tmp_path / "footprints_pkg")

    artifact_path = (
        tmp_path / "footprints_pkg" / "artifacts" / "building_footprints.geojson"
    )
    artifact = package.manifest.artifacts[0]
    payload = json.loads(artifact_path.read_text())
    assert payload["type"] == "FeatureCollection"
    assert payload["crs"]["properties"]["name"] == "EPSG:3006"
    assert payload["features"][0]["properties"]["source_id"] == "building-1"
    assert artifact.path == "artifacts/building_footprints.geojson"
    assert artifact.format == "geojson"
    assert artifact.media_type == "application/geo+json"
    assert artifact.crs == "EPSG:3006"
    assert artifact.size == artifact_path.stat().st_size
    assert artifact.sha256 == _sha256(artifact_path)


def test_dtccpkg_export_writes_zip_package(tmp_path):
    city = City()
    attach_dataset_context(city, _context(datasets.city))

    package = city.export(tmp_path / "city.dtccpkg", format="json")

    assert package.package_format == "dtccpkg"
    assert package.path == tmp_path / "city.dtccpkg"
    assert package.files == (tmp_path / "city.dtccpkg",)
    with zipfile.ZipFile(package.path) as archive:
        assert sorted(archive.namelist()) == [
            "artifacts/city.json",
            "manifest.json",
        ]
        manifest = json.loads(archive.read("manifest.json"))
        artifact_bytes = archive.read("artifacts/city.json")

    artifact = manifest["artifacts"][0]
    assert manifest["schema_version"] == "dtcc-dataset-manifest-v2"
    assert artifact["path"] == "artifacts/city.json"
    assert artifact["size"] == len(artifact_bytes)
    assert artifact["sha256"] == hashlib.sha256(artifact_bytes).hexdigest()


def _context(dataset):
    return dataset.create_context(dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)}))


def _building() -> Building:
    building = Building()
    building.id = "building-1"
    building.add_geometry(_surface(), GeometryType.LOD0)
    return building


def _surface() -> Surface:
    surface = Surface()
    surface.vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    return surface


def _sha256(path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
