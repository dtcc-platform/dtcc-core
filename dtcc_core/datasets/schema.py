"""Dataset v2 schema models.

These models implement the dataset concepts described in ``DESIGN.md``.
"""

from __future__ import annotations

from typing import Annotated, Any, ClassVar, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    FiniteFloat,
    NonNegativeInt,
    PositiveInt,
    field_validator,
)

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


_Box = Annotated[list[FiniteFloat], Field(min_length=4, max_length=4)]


class RasterPlacement(DatasetSchemaModel):
    """Grid of a map-placed raster artifact.

    ``transform`` holds the affine coefficients ``a, b, c, d, e, f`` mapping
    pixel column/row to ``crs`` coordinates, ``bounds`` is the grid envelope in
    ``crs`` and ``bounds_wgs84`` its longitude/latitude envelope, both
    ``[min x, min y, max x, max y]``. ``nodata`` is the source raster's nodata
    as declared, written as the string ``"NaN"`` when it is not a number; a NaN
    declaration is recorded as ``"NaN"`` even when a GeoTIFF of an integer
    raster carries no nodata tag.
    """

    crs: str = Field(min_length=1)
    transform: Annotated[list[FiniteFloat], Field(min_length=6, max_length=6)]
    width: PositiveInt
    height: PositiveInt
    bounds: _Box
    bounds_wgs84: _Box
    nodata: FiniteFloat | Literal["NaN"]

    @field_validator("bounds", "bounds_wgs84")
    @classmethod
    def _ordered(cls, value: list[float]) -> list[float]:
        if value[0] > value[2] or value[1] > value[3]:
            raise ValueError("bounds must be ordered [min x, min y, max x, max y]")
        return value


class RasterRendering(DatasetSchemaModel):
    """How a raster artifact's pixels were produced from the model.

    ``bands`` lists the source bands written as they are; ``field`` names the
    field a ``to_raster`` hook rendered. ``band``, ``value_range`` and
    ``colormap`` describe a colour-mapped single band. ``transparent_nodata``
    is true when invalid pixels are transparent in the image.
    """

    bands: list[NonNegativeInt] | None = None
    field: str | None = None
    band: NonNegativeInt | None = None
    unit: str | None = None
    value_range: (
        Annotated[list[FiniteFloat], Field(min_length=2, max_length=2)] | None
    ) = None
    colormap: str | None = None
    transparent_nodata: bool
    components: list[str] | None = Field(
        default=None, exclude_if=lambda value: value is None
    )

    @field_validator("value_range")
    @classmethod
    def _ordered(cls, value: list[float] | None) -> list[float] | None:
        if value is not None and value[0] > value[1]:
            raise ValueError("value_range must be ordered [low, high]")
        return value


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
    model_type: str | None = Field(default=None, exclude_if=lambda value: value is None)
    model_schema_version: int | None = Field(
        default=None, exclude_if=lambda value: value is None
    )
    derived_from: str | None = Field(
        default=None, exclude_if=lambda value: value is None
    )
    placement: RasterPlacement | None = Field(
        default=None, exclude_if=lambda value: value is None
    )
    rendering: RasterRendering | None = Field(
        default=None, exclude_if=lambda value: value is None
    )


class DatasetManifest(DatasetSchemaModel):
    """Machine-readable legacy v2 and canonical v3 package manifest."""

    schema_version: str = "dtcc-dataset-manifest-v2"
    identity: DatasetIdentity
    metadata: DatasetMetadata
    provenance: DatasetProvenance
    presentation: DatasetPresentation
    request: DatasetRequest
    artifacts: list[DatasetArtifact] = Field(default_factory=list)
    # Canonical package v3 carries Context health; legacy v2 emission is unchanged.
    health: JsonObject | None = Field(
        default=None, exclude_if=lambda value: value is None
    )
    warnings: list[str] | None = Field(
        default=None, exclude_if=lambda value: value is None
    )


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
