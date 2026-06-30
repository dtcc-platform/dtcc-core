"""Tests for initial Dataset v2 context attachment."""

from __future__ import annotations

import json

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import (
    Dataset,
    DatasetContext,
    DatasetManifest,
    attach_dataset_context,
)
from dtcc_core.model import VolumeMesh


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
