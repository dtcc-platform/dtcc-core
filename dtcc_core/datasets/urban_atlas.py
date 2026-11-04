import tempfile
from pathlib import Path
from typing import Literal

from pydantic import Field

import dtcc_core
from dtcc_core.model import Landuse

from ..common import warning
from .dataset import DatasetBaseArgs, DatasetDescriptor, DatasetUpstreamError
from .providers import provider_entry


class UrbanAtlasArgs(DatasetBaseArgs):
    """Arguments for the Copernicus Urban Atlas land use dataset."""

    year: Literal[2018] = Field(
        2018,
        description=(
            "Urban Atlas reference year. Only 2018 is available without a "
            "Copernicus account."
        ),
    )
    format: Literal["pb", "geojson", "gpkg"] | None = Field(
        None,
        description="Output format (pb for protobuf bytes, geojson/gpkg for vector bytes)",
    )


class UrbanAtlasDataset(DatasetDescriptor):
    """Urban Atlas dataset, registered as ``urban_atlas``.

    Copernicus Urban Atlas land use/land cover polygons for the requested
    bounds, returned as a Landuse object in EPSG:3006.

    Call it with keyword arguments defined by ``UrbanAtlasArgs``.
    """

    name = "urban_atlas"
    title = "Urban Atlas Land Use"
    description = (
        "Copernicus Urban Atlas land use/land cover polygons for the requested "
        "bounds, returned as a Landuse object in EPSG:3006."
    )
    ArgsModel = UrbanAtlasArgs
    data_category = "raw"
    result_kind = "landuse"
    python_return_type = "dtcc_core.model.Landuse"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("copernicus-land", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "EEA ArcGIS REST service UA_UrbanAtlas_2018",
            "role": "source_provider",
            "service": "image.discomap.eea.europa.eu",
            "supported_years": [2018],
            "doi": "10.2909/fb4dffa1-6ceb-4cc0-8372-1ed354c285e6",
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "Copernicus full, free and open data policy (Regulation (EU) 2021/696). "
        "Attribution: \"Generated using European Union's Copernicus Land "
        "Monitoring Service information; https://doi.org/10.2909/"
        'fb4dffa1-6ceb-4cc0-8372-1ed354c285e6". State any modification and do '
        "not imply EU endorsement. Requires review before redistribution."
    )
    collection_period = "Urban Atlas 2018 reference year (imagery 2017-2019)"
    data_types = ["vector", "polygons", "landuse"]
    geographic_coverage = "European Functional Urban Areas (12 in Sweden), constrained by requested bounds"
    update_frequency = "Urban Atlas releases every three to six years"
    processing_steps = [
        "Query the anonymous EEA ArcGIS REST service by bounding box in EPSG:3006",
        (
            "Fetch the Functional-Urban-Area-wide 'Other roads' polygon separately "
            "with 0.5 m server-side simplification"
        ),
        "Repair and clip polygons to the requested bounds",
        "Map the 27 Urban Atlas classes onto dtcc-core LanduseClasses",
        "Keep the original Urban Atlas code and class as surface-aligned attributes",
        (
            "On service failure, raise a typed upstream error with strict_live, "
            "otherwise warn and return an empty, uncached result marked partial"
        ),
        "Return a Landuse object or serialize the selected vector/protobuf format",
    ]
    presentation_headline = "Urban Atlas Land Use"
    presentation_summary = (
        "Block-scale land use polygons from the EU Copernicus Urban Atlas, "
        "harmonised across European cities."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each surface is an Urban Atlas polygon clipped to the requested "
                "bounds. Housing is graded by how much of the ground is sealed."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Polygons are city blocks, not individual buildings or parcels. "
                "The minimum mapping unit is 0.25 ha in urban areas, with a target "
                "accuracy of at least 85%."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Only Functional Urban Areas are covered, the 2018 data may be out "
                "of date, and the coarse LanduseClasses mapping loses detail; use "
                "the ua_code attribute for the original 27 classes."
            ),
        },
    ]
    key_points = [
        "Covers 100% of the requested area inside a Functional Urban Area, including roads",
        "Housing classes are graded by sealing degree (from <10% to >80%)",
        "ua_code and ua_class keep the original Urban Atlas nomenclature",
        "Only 2018 is available without a Copernicus account",
    ]
    presentation_legend = {
        "title": "Urban Atlas attributes",
        "entries": [
            {
                "label": "ua_code",
                "meaning": "Urban Atlas class code, such as 11210 or 14100",
            },
            {"label": "ua_class", "meaning": "Urban Atlas class name"},
            {"label": "landuses", "meaning": "Mapped dtcc-core LanduseClasses value"},
        ],
    }
    view_hints = {
        "preferred_geometry": "polygons",
        "default_crs": "EPSG:3006",
        "default_color_attribute": "ua_code",
    }
    presentation_warnings = [
        "Copernicus attribution and DOI must accompany redistributed data.",
        "The EEA REST service is not documented for programmatic use and has no SLA.",
    ]
    presentation_limitations = [
        "Areas outside Functional Urban Areas return no polygons.",
        "Blocks do not resolve individual buildings, streets inside blocks, or small green areas.",
        "Urban Atlas 2021 needs a Copernicus account and is not supported yet.",
    ]

    def build(self, args: UrbanAtlasArgs):
        bounds = self.parse_bounds(args.bounds)
        urban_atlas = dtcc_core.io.data.urban_atlas
        errors: list[DatasetUpstreamError] = []
        if args.format in ("geojson", "gpkg"):
            try:
                gdf = urban_atlas.download_urban_atlas_geodataframe(
                    bounds=bounds, year=args.year
                )
            except DatasetUpstreamError as exc:
                self._degrade(exc, args)
                gdf = urban_atlas.empty_urban_atlas_geodataframe()
            return self._export_gdf_to_bytes(gdf, args.format)

        try:
            landuse: Landuse = dtcc_core.io.data.download_urban_atlas(
                bounds=bounds, year=args.year
            )
        except DatasetUpstreamError as exc:
            self._degrade(exc, args)
            errors.append(exc)
            landuse = urban_atlas.landuse_from_urban_atlas(
                urban_atlas.empty_urban_atlas_geodataframe(), year=args.year
            )
            landuse.attributes["requested_bounds"] = list(bounds.tuple)
        self.apply_result_health_metadata(landuse.attributes, upstream_errors=errors)
        if args.format == "pb":
            return landuse.to_proto().SerializeToString()
        return landuse

    @staticmethod
    def _degrade(exc: DatasetUpstreamError, args: UrbanAtlasArgs) -> None:
        """Re-raise ``exc`` with ``strict_live``; otherwise warn about the empty result."""
        if args.strict_live:
            raise exc
        warning(
            f"{exc} Returning an empty, uncached Urban Atlas result; "
            "pass strict_live=True to raise instead."
        )

    @staticmethod
    def _export_gdf_to_bytes(gdf, format: str) -> bytes:
        with tempfile.TemporaryDirectory() as tmpdir:
            if format == "geojson":
                path = Path(tmpdir) / "urban_atlas.geojson"
                gdf.to_file(path, driver="GeoJSON")
            else:
                path = Path(tmpdir) / "urban_atlas.gpkg"
                gdf.to_file(path, layer="urban_atlas", driver="GPKG")
            return path.read_bytes()


__all__ = ["UrbanAtlasArgs", "UrbanAtlasDataset"]
