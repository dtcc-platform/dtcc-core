"""Tests for the synthetic smoke dataset."""

from __future__ import annotations

import json
import importlib
import hashlib
from io import BytesIO

import numpy as np
import pytest
from PIL import Image

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.dataset import DatasetExportResult
from dtcc_core.datasets.smoke import SmokeArgs, SmokeDataset
from dtcc_core.model import FieldSlice, StreamlineCollection, VolumeMesh, proto


def test_smoke_dataset_registered_name():
    assert SmokeDataset().name == "smoke"


def test_smoke_get_dataset_and_module_attribute():
    ds = get_dataset("smoke")
    assert ds is not None
    assert ds.name == "smoke"
    assert hasattr(datasets, "smoke")
    assert callable(datasets.smoke)


def test_smoke_describe_includes_products_and_formats():
    metadata = datasets.smoke.describe()

    assert metadata["name"] == "smoke"
    assert metadata["data_category"] == "simulation"
    assert metadata["result_kind"] == "vector_field"
    assert metadata["python_return_type"] == "dtcc_core.model.VolumeMesh"
    assert set(metadata["supported_formats"]) == {
        "pb",
        "vtu",
        "geojson",
        "png",
        "mp4",
    }

    products = {product["name"]: product for product in metadata["products"]}
    assert products["field"]["formats"] == ["pb", "vtu", "geojson"]
    assert products["field"]["fields"] == ["velocity", "speed", "pressure"]
    assert products["slice"]["python_return_type"] == "dtcc_core.model.FieldSlice"
    assert products["slice"]["formats"] == ["geojson", "png", "mp4"]
    assert products["slice"]["fields"] == ["velocity", "speed", "pressure"]
    assert (
        products["streamlines"]["python_return_type"]
        == "dtcc_core.model.StreamlineCollection"
    )
    assert products["streamlines"]["formats"] == ["geojson", "png", "mp4"]
    assert products["streamlines"]["fields"] == ["velocity", "speed", "pressure"]
    assert metadata["field_names"] == ["velocity", "speed", "pressure"]


def test_smoke_default_returns_volume_mesh_with_fields():
    mesh = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        zmax=40.0,
        resolution=5,
    )

    assert isinstance(mesh, VolumeMesh)
    assert len(mesh.vertices) == 5**3
    assert len(mesh.cells) == 5 * 4**3
    assert np.allclose(mesh.vertices[0], [0.0, 0.0, 0.0])
    assert np.allclose(mesh.vertices[-1], [10.0, 20.0, 40.0])

    fields = {field.name: field for field in mesh.fields}
    assert fields["velocity"].dim == 3
    assert fields["velocity"].values.shape == (5**3, 3)
    assert fields["speed"].dim == 1
    assert fields["speed"].values.shape == (5**3, 1)
    assert fields["pressure"].dim == 1
    assert fields["pressure"].values.shape == (5**3, 1)

    center = np.where(np.all(np.isclose(mesh.vertices, [5.0, 10.0, 20.0]), axis=1))[0]
    assert len(center) == 1
    assert np.allclose(fields["velocity"].values[center[0]], [0.0, 0.0, 0.0])
    assert np.allclose(fields["speed"].values[center[0]], [0.0])


def test_smoke_protobuf_format_returns_volume_mesh_bytes():
    payload = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        format="pb",
    )

    assert isinstance(payload, bytes)
    pb = proto.Geometry.FromString(payload)
    assert len(pb.volume_mesh.vertices) == 4**3 * 3
    assert len(pb.volume_mesh.cells) == 5 * 3**3 * 4
    assert [field.name for field in pb.fields] == ["velocity", "speed", "pressure"]


def test_smoke_field_geojson_format_returns_points():
    payload = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        format="geojson",
    )

    assert isinstance(payload, bytes)
    data = json.loads(payload.decode("utf-8"))
    assert data["type"] == "FeatureCollection"
    assert data["metadata"]["product"] == "field"
    assert data["metadata"]["crs"] == "EPSG:3006"
    assert data["crs"]["properties"]["name"] == "EPSG:3006"
    assert len(data["features"]) == 4**3
    assert data["features"][0]["geometry"]["type"] == "Point"
    assert len(data["features"][0]["geometry"]["coordinates"]) == 3
    assert {"u", "v", "w", "speed", "pressure"} <= set(
        data["features"][0]["properties"]
    )


def test_smoke_slice_geojson_format_returns_plane_points():
    payload = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=6,
        product="slice",
        slice_axis="y",
        slice_position=0.25,
        format="geojson",
    )

    data = json.loads(payload.decode("utf-8"))
    assert data["metadata"]["product"] == "slice"
    assert data["metadata"]["slice_axis"] == "y"
    assert data["metadata"]["slice_position"] == 0.25
    assert len(data["features"]) == 6**2
    assert {feature["geometry"]["type"] for feature in data["features"]} == {"Point"}
    assert {"velocity", "speed", "pressure"} <= set(data["metadata"]["fields"])
    assert "pressure" in data["features"][0]["properties"]


def test_smoke_slice_without_format_returns_field_slice():
    result = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        product="slice",
    )

    assert isinstance(result, FieldSlice)
    assert result.dataset_context is not None
    assert result.metadata is not None
    assert result.provenance is not None
    assert result.presentation is not None
    assert len(result.points) == 4**2
    fields = {field.name: field for field in result.fields}
    assert set(fields) == {"velocity", "speed", "pressure"}
    assert fields["velocity"].dim == 3
    assert fields["speed"].dim == 1
    assert fields["pressure"].dim == 1
    assert fields["velocity"].values.shape == (len(result.points), 3)
    assert fields["speed"].values.shape == (len(result.points), 1)
    assert fields["pressure"].values.shape == (len(result.points), 1)
    assert result.to_python()["type"] == "FeatureCollection"
    assert result.to_python()["metadata"]["product"] == "slice"
    json.dumps(result.manifest().model_dump(mode="json"))


def test_smoke_streamlines_without_format_returns_streamline_collection():
    result = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="streamlines",
        streamline_count=9,
        streamline_steps=40,
    )

    assert isinstance(result, StreamlineCollection)
    assert result.dataset_context is not None
    assert len(result) > 0
    assert len(result) <= 9
    for line in result:
        assert len(line.vertices) >= 2
        fields = {field.name: field for field in line.fields}
        assert set(fields) == {"velocity", "speed", "pressure"}
        assert fields["velocity"].dim == 3
        assert fields["speed"].dim == 1
        assert fields["pressure"].dim == 1
        assert fields["velocity"].values.shape == (len(line.vertices), 3)
        assert fields["speed"].values.shape == (len(line.vertices), 1)
        assert fields["pressure"].values.shape == (len(line.vertices), 1)


def test_smoke_slice_snapshots_are_time_dependent():
    first = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=6,
        product="slice",
        time=0.0,
    )
    second = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=6,
        product="slice",
        time=2.0,
        period=8.0,
    )

    first_speed = first.field("speed").values[0, 0]
    second_speed = second.field("speed").values[0, 0]
    assert first.time == 0.0
    assert second.time == 2.0
    assert first.period == 8.0
    assert not np.isclose(first_speed, second_speed)


def test_smoke_geojson_can_omit_crs_and_z_coordinates():
    payload = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        product="slice",
        crs=None,
        include_z=False,
        format="geojson",
    )

    data = json.loads(payload.decode("utf-8"))
    assert "crs" not in data
    assert "crs" not in data["metadata"]
    assert len(data["features"][0]["geometry"]["coordinates"]) == 2


def test_smoke_streamlines_geojson_format_returns_lines():
    payload = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="streamlines",
        streamline_count=9,
        streamline_steps=40,
        format="geojson",
    )

    data = json.loads(payload.decode("utf-8"))
    assert data["metadata"]["product"] == "streamlines"
    assert len(data["features"]) > 0
    assert len(data["features"]) <= 9
    assert {feature["geometry"]["type"] for feature in data["features"]} == {
        "LineString"
    }
    assert all(len(feature["geometry"]["coordinates"]) >= 2 for feature in data["features"])
    assert {"velocity", "speed", "pressure"} <= set(data["metadata"]["fields"])
    assert {"velocity", "speed", "pressure"} <= set(
        data["features"][0]["properties"]
    )


def test_smoke_slice_png_format_returns_exact_size_image():
    payload = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=128,
        product="slice",
        format="png",
        width=320,
        height=240,
    )

    assert isinstance(payload, bytes)
    assert payload.startswith(b"\x89PNG")
    with Image.open(BytesIO(payload)) as image:
        assert image.size == (320, 240)
        assert image.format == "PNG"
        assert np.asarray(image.convert("RGB")).std() > 0


def test_smoke_field_slice_object_export_defaults_to_png_package(tmp_path):
    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=16,
        product="slice",
        width=160,
        height=90,
    )

    package = field_slice.export(tmp_path / "smoke_slice_pkg")

    artifact_path = _assert_smoke_visual_package(
        package,
        package_dir=tmp_path / "smoke_slice_pkg",
        product="slice",
        expected_artifact_path="artifacts/smoke_slice.png",
    )
    with Image.open(artifact_path) as image:
        assert image.size == (160, 90)
        assert image.format == "PNG"
        assert np.asarray(image.convert("RGB")).std() > 0


def test_smoke_field_slice_plot_shows_presentation_by_default():
    import matplotlib

    matplotlib.use("Agg", force=True)

    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=8,
        product="slice",
        width=160,
        height=90,
    )

    ax = field_slice.plot(show=False)

    try:
        assert len(ax.figure.axes) >= 2
        panel_ax = ax.figure.axes[1]
        assert panel_ax.patch.get_visible()
        assert panel_ax.get_facecolor()[:3] == pytest.approx(
            (247 / 255, 247 / 255, 247 / 255)
        )
        assert any(
            "Synthetic smoke field" in text.get_text()
            for figure_ax in ax.figure.axes
            for text in figure_ax.texts
        )
    finally:
        import matplotlib.pyplot as plt

        plt.close(ax.figure)


def test_smoke_field_slice_plot_can_hide_presentation():
    import matplotlib

    matplotlib.use("Agg", force=True)

    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=8,
        product="slice",
        width=160,
        height=90,
    )

    ax = field_slice.plot(show=False, presentation=False)

    try:
        assert len(ax.figure.axes) == 1
        assert not any(
            "Synthetic smoke field" in text.get_text()
            for figure_ax in ax.figure.axes
            for text in figure_ax.texts
        )
    finally:
        import matplotlib.pyplot as plt

        plt.close(ax.figure)


def test_smoke_descriptor_plot_shows_presentation_by_default():
    import matplotlib

    matplotlib.use("Agg", force=True)

    ax = datasets.smoke.plot(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=8,
        product="slice",
        width=160,
        height=90,
        show=False,
    )

    try:
        assert len(ax.figure.axes) >= 2
        panel_ax = ax.figure.axes[1]
        assert panel_ax.patch.get_visible()
        assert panel_ax.get_facecolor()[:3] == pytest.approx(
            (247 / 255, 247 / 255, 247 / 255)
        )
        assert any(
            "Synthetic smoke field" in text.get_text()
            for figure_ax in ax.figure.axes
            for text in figure_ax.texts
        )
    finally:
        import matplotlib.pyplot as plt

        plt.close(ax.figure)


def test_smoke_streamline_collection_object_export_defaults_to_png_package(tmp_path):
    streamlines = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="streamlines",
        streamline_count=9,
        streamline_steps=40,
        width=160,
        height=90,
    )

    package = streamlines.export(tmp_path / "smoke_streamlines_pkg")

    artifact_path = _assert_smoke_visual_package(
        package,
        package_dir=tmp_path / "smoke_streamlines_pkg",
        product="streamlines",
        expected_artifact_path="artifacts/smoke_streamlines.png",
    )
    with Image.open(artifact_path) as image:
        assert image.size == (160, 90)
        assert image.format == "PNG"
        assert np.asarray(image.convert("RGB")).std() > 0


def test_smoke_object_export_supports_explicit_geojson(tmp_path):
    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        product="slice",
    )

    package = field_slice.export(tmp_path / "smoke_slice_geojson", format="geojson")

    artifact = package.manifest.artifacts[0]
    artifact_path = tmp_path / "smoke_slice_geojson" / artifact.path
    payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert artifact.format == "geojson"
    assert artifact.media_type == "application/geo+json"
    assert payload["metadata"]["product"] == "slice"


def test_smoke_object_export_rejects_mp4_with_dataset_level_hint(tmp_path):
    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        product="slice",
    )

    with pytest.raises(ValueError, match="dataset-level MP4 path"):
        field_slice.export(tmp_path / "smoke_slice_mp4", format="mp4")


def test_smoke_field_product_keeps_resolution_guardrail():
    with pytest.raises(ValueError, match="resolution <= 64"):
        datasets.smoke(
            bounds=(0.0, 0.0, 10.0, 20.0),
            product="field",
            resolution=128,
        )


def test_smoke_field_product_rejects_video_format():
    with pytest.raises(ValueError, match="product='field' does not support"):
        datasets.smoke(
            bounds=(0.0, 0.0, 10.0, 20.0),
            product="field",
            format="mp4",
        )


def test_smoke_mp4_rejects_transparent_background():
    with pytest.raises(ValueError, match="transparent"):
        datasets.smoke(
            bounds=(0.0, 0.0, 10.0, 20.0),
            product="slice",
            format="mp4",
            transparent=True,
        )


def test_smoke_streamlines_png_export_writes_visualization_manifest(tmp_path):
    path = tmp_path / "smoke_streamlines.png"

    result = datasets.smoke.export(
        path,
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="streamlines",
        streamline_count=9,
        streamline_steps=40,
        width=320,
        height=240,
        cmap="inferno",
    )

    assert result.path == path
    assert result.manifest_path == tmp_path / "smoke_streamlines.manifest.json"
    with Image.open(path) as image:
        assert image.size == (320, 240)
        assert image.format == "PNG"
        assert np.asarray(image.convert("RGB")).std() > 0

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["format"] == "png"
    assert manifest["media_type"] == "image/png"
    assert manifest["data_kind"] == "raster"
    assert manifest["product"] == "streamlines"
    assert manifest["fields"] == ["velocity", "speed", "pressure"]

    visualization = manifest["visualization"]
    assert visualization["profile"] == "table"
    assert visualization["width"] == 320
    assert visualization["height"] == 240
    assert visualization["colormap"] == "inferno"
    assert visualization["product_kind"] == "streamlines"
    assert visualization["visual_axes"] == ["x", "y"]
    assert visualization["extent"] == [0.0, 0.0, 10.0, 20.0]
    assert visualization["origin"] == "lower"
    assert visualization["requested_line_count"] == 9
    assert visualization["fields"] == ["velocity", "speed", "pressure"]


def test_smoke_streamlines_mp4_export_writes_video_manifest(tmp_path, monkeypatch):
    path = tmp_path / "smoke_streamlines.mp4"
    smoke_module = importlib.import_module("dtcc_core.datasets.smoke")

    monkeypatch.setattr(
        smoke_module,
        "render_product_mp4",
        lambda product_factory, options: b"fake mp4",
    )

    result = datasets.smoke.export(
        path,
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="streamlines",
        streamline_count=9,
        streamline_steps=40,
        width=320,
        height=240,
        fps=12,
        duration=2.0,
        period=4.0,
        time=1.0,
    )

    assert result.path == path
    assert path.read_bytes() == b"fake mp4"

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["format"] == "mp4"
    assert manifest["media_type"] == "video/mp4"
    assert manifest["data_kind"] == "video"
    assert manifest["product"] == "streamlines"

    visualization = manifest["visualization"]
    assert visualization["profile"] == "table"
    assert visualization["width"] == 320
    assert visualization["height"] == 240
    assert visualization["fps"] == 12
    assert visualization["duration"] == 2.0
    assert visualization["frame_count"] == 24
    assert visualization["start_time"] == 1.0
    assert visualization["end_time"] == 3.0
    assert visualization["time_period"] == 4.0
    assert visualization["codec"] == "h264"
    assert visualization["container"] == "mp4"
    assert visualization["loop"] is True
    assert visualization["product_kind"] == "streamlines"


def test_smoke_export_writes_payload_and_manifest(tmp_path):
    path = tmp_path / "smoke_slice.geojson"

    result = datasets.smoke.export(
        path,
        bounds=(319720, 6397660, 320220, 6398160),
        resolution=4,
        product="slice",
    )

    assert result.path == path
    assert result.manifest_path == tmp_path / "smoke_slice.manifest.json"
    assert path.exists()
    assert result.manifest_path.exists()
    assert result.files == (path,)
    assert result.format == "geojson"

    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["metadata"]["product"] == "slice"

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest == result.manifest
    assert manifest["manifest_schema_version"] == "dtcc-dataset-manifest-v1"
    assert manifest["created_by"]["package"] == "dtcc-core"
    assert isinstance(manifest["created_by"]["version"], str)
    assert manifest["created_by"]["version"]
    descriptor = datasets.smoke.describe()
    for key, value in descriptor.items():
        assert manifest[key] == value
    assert "id" not in manifest
    assert "projectionBbox" not in manifest
    assert manifest["file"] == "smoke_slice.geojson"
    assert manifest["format"] == "geojson"
    assert manifest["product"] == "slice"
    assert manifest["bounds"] == [319720, 6397660, 320220, 6398160]
    assert manifest["parameters"]["bounds"] == [319720, 6397660, 320220, 6398160]
    assert manifest["parameters"]["product"] == "slice"
    assert manifest["parameters"]["format"] == "geojson"
    assert manifest["parameters"]["resolution"] == 4


def test_smoke_export_manifest_can_be_disabled_or_overridden(tmp_path):
    path = tmp_path / "slice.geojson"
    manifest_path = tmp_path / "atlas.json"

    result = datasets.smoke.export(
        path,
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
        manifest=False,
    )

    assert path.exists()
    assert result.manifest is None
    assert result.manifest_path is None
    assert result.files == (path,)
    assert result.format == "geojson"
    assert not (tmp_path / "slice.manifest.json").exists()

    result = datasets.smoke.export(
        path,
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
        manifest_path=manifest_path,
        manifest_id="atlas-smoke",
        title="Atlas Smoke",
        description="Smoke data prepared for Atlas++.",
    )

    assert result.manifest_path == manifest_path
    assert result.manifest["id"] == "atlas-smoke"
    assert result.manifest["name"] == "smoke"
    assert result.manifest["title"] == "Atlas Smoke"
    assert result.manifest["description"] == "Smoke data prepared for Atlas++."
    assert result.manifest["bounds"] == [0.0, 0.0, 10.0, 20.0]
    assert "projectionBbox" not in result.manifest
    assert json.loads(manifest_path.read_text(encoding="utf-8")) == result.manifest


def test_dataset_export_result_supports_legacy_construction(tmp_path):
    path = tmp_path / "slice.geojson"
    manifest_path = tmp_path / "slice.manifest.json"
    manifest = {"name": "smoke"}

    result = DatasetExportResult(
        path=path,
        manifest_path=manifest_path,
        manifest=manifest,
    )

    assert result.path == path
    assert result.manifest_path == manifest_path
    assert result.manifest == manifest
    assert result.files == (path,)
    assert result.format == "geojson"


def test_smoke_export_result_publish_delegates_to_uploader(tmp_path):
    path = tmp_path / "slice.geojson"
    publication = object()
    calls = []

    class FakeUploader:
        def upload_package(
            self,
            *,
            dataset_key,
            manifest_path,
            files,
            manifest,
            idempotency_key,
        ):
            calls.append(
                {
                    "dataset_key": dataset_key,
                    "manifest_path": manifest_path,
                    "files": files,
                    "manifest": manifest,
                    "idempotency_key": idempotency_key,
                }
            )
            return publication

    result = datasets.smoke.export(
        path,
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
    )

    assert (
        result.publish(
            dataset_key="atlas/smoke",
            uploader=FakeUploader(),
            idempotency_key="upload-1",
        )
        is publication
    )
    assert calls == [
        {
            "dataset_key": "atlas/smoke",
            "manifest_path": result.manifest_path,
            "files": result.files,
            "manifest": result.manifest,
            "idempotency_key": "upload-1",
        }
    ]


def test_smoke_publish_methods_are_keyword_only(tmp_path):
    result = datasets.smoke.export(
        tmp_path / "slice.geojson",
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
    )

    with pytest.raises(TypeError):
        result.publish("atlas/smoke", uploader=object())
    with pytest.raises(TypeError):
        datasets.smoke.publish("atlas/smoke")


def test_smoke_export_result_publish_requires_manifest(tmp_path):
    from dtcc_core.datasets.publish import DatasetPackageError

    result = datasets.smoke.export(
        tmp_path / "slice.geojson",
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
        manifest=False,
    )

    with pytest.raises(DatasetPackageError, match="without a manifest") as excinfo:
        result.publish(dataset_key="atlas/smoke", uploader=object())

    assert excinfo.value.failure_class == "invalid_package"


def test_smoke_rejects_non_geojson_visualization_formats():
    with pytest.raises(ValueError, match="product='slice' only supports"):
        SmokeDataset().build(
            SmokeArgs(
                bounds=(0.0, 0.0, 10.0, 20.0),
                product="slice",
                format="pb",
            )
        )


def _assert_smoke_visual_package(
    package,
    *,
    package_dir,
    product: str,
    expected_artifact_path: str,
):
    manifest_path = package_dir / "manifest.json"
    artifact = package.manifest.artifacts[0]
    artifact_path = package_dir / artifact.path
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert package.path == package_dir
    assert package.manifest_path == manifest_path
    assert manifest["schema_version"] == "dtcc-dataset-manifest-v2"
    assert manifest["identity"]["name"] == "smoke"
    assert manifest["request"]["parameters"]["product"] == product
    assert artifact.path.startswith("artifacts/")
    assert artifact.path == expected_artifact_path
    assert artifact.format == "png"
    assert artifact.media_type == "image/png"
    assert artifact.data_kind == "raster"
    assert artifact.role == "primary"
    assert artifact.size == artifact_path.stat().st_size
    assert artifact.sha256 == _sha256(artifact_path)
    return artifact_path


def _sha256(path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
