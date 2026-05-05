"""Tests for the synthetic smoke dataset."""

from __future__ import annotations

import json

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.smoke import SmokeArgs, SmokeDataset
from dtcc_core.model import VolumeMesh, proto


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
    assert set(metadata["supported_formats"]) == {"pb", "vtu", "geojson"}

    products = {product["name"]: product for product in metadata["products"]}
    assert products["field"]["formats"] == ["pb", "vtu", "geojson"]
    assert products["slice"]["formats"] == ["geojson"]
    assert products["streamlines"]["formats"] == ["geojson"]


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
    assert [field.name for field in pb.fields] == ["velocity", "speed"]


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
    assert {"u", "v", "w", "speed"} <= set(data["features"][0]["properties"])


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


def test_smoke_slice_without_format_returns_geojson_dict():
    result = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
        product="slice",
    )

    assert isinstance(result, dict)
    assert result["type"] == "FeatureCollection"
    assert result["metadata"]["product"] == "slice"


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

    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["metadata"]["product"] == "slice"

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest == result.manifest
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


def test_smoke_rejects_non_geojson_visualization_formats():
    with pytest.raises(ValueError, match="product='slice' only supports"):
        SmokeDataset().build(
            SmokeArgs(
                bounds=(0.0, 0.0, 10.0, 20.0),
                product="slice",
                format="pb",
            )
        )
