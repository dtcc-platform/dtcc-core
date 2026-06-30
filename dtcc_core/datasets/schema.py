"""Dataset v2 schema models.

These models are the initial dtcc-core representation of the Dataset v2
concepts described in ``docs/design/datasets-v2.md``.
"""

from __future__ import annotations

from typing import Any, ClassVar

from pydantic import BaseModel, ConfigDict, Field


JsonObject = dict[str, Any]
JsonListItem = str | JsonObject


class DatasetSchemaModel(BaseModel):
    """Base model for Dataset v2 schema objects."""

    model_config = ConfigDict(extra="forbid")


class DatasetIdentity(DatasetSchemaModel):
    """Stable dataset naming and addressing."""

    name: str
    title: str
    version: str | None = None


class DatasetMetadata(DatasetSchemaModel):
    """Concise factual discovery and catalog metadata."""

    description: str = ""
    provider: list[JsonListItem] = Field(default_factory=list)
    source: list[JsonListItem] = Field(default_factory=list)
    license: str | JsonObject | None = None
    collection_period: str | JsonObject | None = None
    crs: list[str] = Field(default_factory=list)
    lod: str | None = None
    data_types: list[str] = Field(default_factory=list)
    formats: list[str] = Field(default_factory=list)
    geographic_coverage: str | None = None
    update_frequency: str | None = None
    data_category: str | None = None
    result_kind: str | None = None
    python_return_type: str | None = None


class DatasetProvenance(DatasetSchemaModel):
    """Lineage and reproducibility information."""

    sources: list[JsonListItem] = Field(default_factory=list)
    processing_steps: list[JsonListItem] = Field(default_factory=list)
    generated_by: str | JsonObject | None = None
    generated_at: str | None = None
    derived_from: list[JsonListItem] = Field(default_factory=list)


class DatasetPresentation(DatasetSchemaModel):
    """Human explanation and display guidance."""

    headline: str | None = None
    summary: str | None = None
    narrative: list[JsonListItem] = Field(default_factory=list)
    key_points: list[str] = Field(default_factory=list)
    legend: str | JsonObject | None = None
    annotations: list[JsonListItem] = Field(default_factory=list)
    view_hints: JsonObject | None = None
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class DatasetRequest(DatasetSchemaModel):
    """Concrete parameter set used to produce a dataset object."""

    dataset_name: str
    parameters: JsonObject = Field(default_factory=dict)
    bounds: list[float] | None = None


class DatasetArtifact(DatasetSchemaModel):
    """Concrete file inside an exported or published dataset package."""

    path: str
    role: str
    format: str
    media_type: str
    data_kind: str
    crs: str | None = None
    bounds: list[float] | None = None
    size: int | None = None
    sha256: str | None = None


class DatasetManifest(DatasetSchemaModel):
    """Machine-readable Dataset Manifest v2 package contract."""

    schema_version: str = "dtcc-dataset-manifest-v2"
    identity: DatasetIdentity
    metadata: DatasetMetadata
    provenance: DatasetProvenance
    presentation: DatasetPresentation
    request: DatasetRequest
    artifacts: list[DatasetArtifact] = Field(default_factory=list)


class DatasetContext(DatasetSchemaModel):
    """Dataset context attached to a native DTCC model object."""

    manifest_schema_version: ClassVar[str] = "dtcc-dataset-manifest-v2"

    identity: DatasetIdentity
    metadata: DatasetMetadata
    provenance: DatasetProvenance
    presentation: DatasetPresentation
    request: DatasetRequest
    health: JsonObject | None = None
    warnings: list[str] = Field(default_factory=list)

    def manifest(
        self,
        artifacts: list[DatasetArtifact] | None = None,
    ) -> DatasetManifest:
        """Return a Dataset Manifest v2 snapshot for this context."""
        return DatasetManifest(
            schema_version=self.manifest_schema_version,
            identity=self.identity,
            metadata=self.metadata,
            provenance=self.provenance,
            presentation=self.presentation,
            request=self.request,
            artifacts=list(artifacts or ()),
        )
