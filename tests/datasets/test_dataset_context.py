"""Tests for initial Dataset v2 context attachment."""

from __future__ import annotations

import json

import pytest

import dtcc_core.datasets as datasets
import dtcc_core.model as model
from dtcc_core.datasets import (
    Dataset,
    DatasetContext,
    DatasetManifest,
    attach_dataset_context,
)
from dtcc_core.model import DatasetCollection, DatasetValue, VolumeMesh


def test_dataset_alias_keeps_descriptor_compatibility():
    assert isinstance(datasets.smoke, Dataset)


def test_smoke_returns_native_model_with_dataset_context():
    mesh = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        zmax=40.0,
        resolution=4,
    )

    assert isinstance(mesh, VolumeMesh)
    assert isinstance(mesh.dataset_context, DatasetContext)
    assert mesh.metadata is mesh.dataset_context.metadata
    assert mesh.provenance is mesh.dataset_context.provenance
    assert mesh.presentation is mesh.dataset_context.presentation

    manifest = mesh.manifest()
    assert isinstance(manifest, DatasetManifest)
    assert manifest.schema_version == "dtcc-dataset-manifest-v2"
    assert manifest.identity.name == "smoke"
    assert manifest.identity.title == "Smoke"
    assert manifest.metadata.description == datasets.smoke.description
    assert manifest.metadata.data_category == "simulation"
    assert manifest.metadata.result_kind == "vector_field"
    assert manifest.metadata.formats == ["pb", "vtu", "geojson", "png", "mp4"]
    assert manifest.provenance.generated_by["package"] == "dtcc-core"
    assert manifest.presentation.headline == "Smoke"
    assert manifest.request.dataset_name == "smoke"
    assert manifest.request.bounds == [0.0, 0.0, 10.0, 20.0]
    assert manifest.request.parameters["resolution"] == 4
    assert manifest.artifacts == []


def test_dataset_object_info_includes_presentation_tables():
    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
        resolution=4,
        width=80,
        height=64,
    )

    text = field_slice.info(print=False)

    assert "Tangible Table Metadata" in text
    assert "C-P1 Provider / Source" in text
    assert "DTCC Platform (generator)" in text
    assert "C-S7 Processing / Methodology" in text
    assert "Presentation" in text
    assert "No live data dependency" in text
    assert "velocity" in text
    assert "speed" in text
    assert "pressure" in text


def test_dataset_object_info_can_hide_presentation_tables():
    field_slice = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        product="slice",
        resolution=4,
        width=80,
        height=64,
    )

    text = field_slice.info(print=False, presentation=False)

    assert "Tangible Table Metadata" not in text
    assert "Presentation" not in text


def test_descriptor_metadata_flows_to_context_manifest():
    context = datasets.smoke.create_context(
        datasets.smoke.validate({"bounds": (0.0, 0.0, 10.0, 20.0)})
    )
    manifest = context.manifest()

    assert manifest.metadata.provider == [
        {"name": "DTCC Platform", "role": "generator"}
    ]
    assert manifest.metadata.source == [
        "Synthetic analytical velocity, speed, and pressure fields"
    ]
    assert manifest.metadata.license == "MIT"
    assert manifest.metadata.update_frequency == "generated on demand"
    assert manifest.provenance.processing_steps[:3] == [
        "Map requested bounds to the normalized smoke domain",
        "Evaluate deterministic analytical smoke fields",
        "Generate requested field, slice, or streamline product",
    ]
    assert manifest.provenance.processing_steps[-1] == "Build dataset 'smoke'"
    assert manifest.presentation.summary == datasets.smoke.presentation_summary
    assert manifest.presentation.view_hints["preferred_media_types"] == [
        "image/png",
        "video/mp4",
        "application/geo+json",
    ]


def test_public_dataset_context_metadata_audit():
    missing = {}
    for name, dataset in datasets.list().items():
        if not dataset.__class__.__module__.startswith("dtcc_core.datasets."):
            continue
        context = dataset.create_context(
            dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)})
        )
        failures = []
        if not context.metadata.description:
            failures.append("description")
        if not context.metadata.provider:
            failures.append("provider")
        if not context.metadata.source:
            failures.append("source")
        if not context.metadata.license:
            failures.append("license")
        if not context.metadata.update_frequency:
            failures.append("update_frequency")
        if not context.provenance.processing_steps:
            failures.append("processing_steps")
        if not context.provenance.generated_by:
            failures.append("generated_by")
        if not context.presentation.headline:
            failures.append("presentation.headline")
        if not context.presentation.summary:
            failures.append("presentation.summary")
        if failures:
            missing[name] = failures

    assert missing == {}


def test_known_epsg3006_datasets_have_context_crs_without_request_crs():
    dataset_names = [
        "point_cloud",
        "building_footprints",
        "buildings",
        "city",
        "terrain_surface_mesh",
        "city_flat_mesh",
        "city_surface_mesh",
        "city_volume_mesh",
        "trees",
        "roads",
        "space_syntax",
    ]

    for name in dataset_names:
        dataset = datasets.get_dataset(name)
        context = dataset.create_context(
            dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)})
        )
        assert "EPSG:3006" in context.metadata.crs


def test_dataset_context_serializes_to_json_safe_data():
    mesh = datasets.smoke(
        bounds=(0.0, 0.0, 10.0, 20.0),
        resolution=4,
    )

    context_data = mesh.dataset_context.model_dump(mode="json")
    manifest_data = mesh.manifest().model_dump(mode="json")

    json.dumps(context_data)
    json.dumps(manifest_data)
    assert json.loads(mesh.manifest().model_dump_json()) == manifest_data


def test_model_without_dataset_context_has_no_contextual_metadata():
    mesh = VolumeMesh()

    assert mesh.dataset_context is None
    assert mesh.metadata is None
    assert mesh.provenance is None
    assert mesh.presentation is None
    with pytest.raises(ValueError, match="no DatasetContext"):
        mesh.manifest()


def test_attach_dataset_context_leaves_bare_containers_unchanged():
    context = datasets.smoke.create_context(
        datasets.smoke.validate({"bounds": (0.0, 0.0, 10.0, 20.0)})
    )
    values = [{"value": 1}]
    mapping = {"value": 1}

    assert attach_dataset_context(values, context) is values
    assert attach_dataset_context(mapping, context) is mapping
    assert not hasattr(values, "dataset_context")
    assert not hasattr(mapping, "dataset_context")


def test_dataset_collection_can_carry_dataset_context():
    context = datasets.smoke.create_context(
        datasets.smoke.validate({"bounds": (0.0, 0.0, 10.0, 20.0)})
    )
    collection = DatasetCollection(items=["a", "b"])

    result = attach_dataset_context(collection, context)

    assert result is collection
    assert len(collection) == 2
    assert list(collection) == ["a", "b"]
    assert collection[0] == "a"
    assert collection.to_list() == ["a", "b"]
    assert collection.dataset_context is context
    assert collection.metadata is context.metadata
    assert isinstance(collection.manifest(), DatasetManifest)


def test_dataset_value_can_carry_dataset_context():
    context = datasets.smoke.create_context(
        datasets.smoke.validate({"bounds": (0.0, 0.0, 10.0, 20.0)})
    )
    value = DatasetValue({"type": "FeatureCollection", "features": []})

    result = attach_dataset_context(value, context)

    assert result is value
    assert value["type"] == "FeatureCollection"
    assert value.get("missing", "fallback") == "fallback"
    assert list(value.keys()) == ["type", "features"]
    assert list(value.items())[0] == ("type", "FeatureCollection")
    assert value.to_python() == {"type": "FeatureCollection", "features": []}
    assert value.dataset_context is context
    assert value.provenance is context.provenance
    assert isinstance(value.manifest(), DatasetManifest)


def test_no_dataset_result_or_run_api_is_introduced():
    assert not hasattr(datasets, "DatasetResult")
    assert not hasattr(datasets, "VectorLayer")
    assert not hasattr(model, "VectorLayer")
    for dataset in (
        datasets.smoke,
        datasets.buildings,
        datasets.building_footprints,
        datasets.trees,
        datasets.calibration_grid,
    ):
        assert not hasattr(dataset, "run")
