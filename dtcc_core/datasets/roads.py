import dtcc_core
from dtcc_core.model import RoadNetwork
from pydantic import Field
from typing import Literal, Optional

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_entry


class RoadsArgs(DatasetBaseArgs):
    source: Literal["OSM"] = Field(
        "OSM",
        description=(
            "Road-network source. OSM uses OpenStreetMap highway ways through "
            "the Overpass/cache download path."
        ),
    )
    format: Optional[Literal["pb"]] = Field(
        None, description="Output format (pb for protobuf bytes)"
    )


class RoadsDataset(DatasetDescriptor):
    """OpenStreetMap road-network dataset."""

    name = "roads"
    title = "OpenStreetMap Roads"
    description = (
        "Road network data from OpenStreetMap via the dtcc-core Overpass/cache "
        "download path, returned as a RoadNetwork in EPSG:3006."
    )
    ArgsModel = RoadsArgs
    data_category = "raw"
    result_kind = "road_network"
    python_return_type = "dtcc_core.model.RoadNetwork"
    provider = [provider_entry("openstreetmap", role="source_provider")]
    source = [
        {
            "name": "OpenStreetMap road network data via Overpass",
            "role": "source_provider",
            "service": "overpass",
            "selected_when": 'source="OSM"',
            "source_terms_status": "requires_review",
            "license": "ODbL",
        }
    ]
    license = (
        "Requires review: OpenStreetMap data is licensed under ODbL; verify "
        "attribution, share-alike, and redistribution requirements before "
        "publishing derived packages."
    )
    collection_period = (
        "Current OpenStreetMap/Overpass response or local cache hit at request "
        "time; no fixed historical snapshot is recorded by this dataset."
    )
    default_crs = "EPSG:3006"
    data_types = ["road_network", "lines", "graph", "openstreetmap"]
    geographic_coverage = (
        "Global OpenStreetMap coverage, constrained by requested bounds"
    )
    update_frequency = "depends on OpenStreetMap edits and Overpass availability"
    processing_steps = [
        "Validate requested bounds and road source",
        "Query the dtcc-core road download wrapper for source='OSM'",
        "Use cached Overpass data when a suitable cached bbox is available",
        (
            "Otherwise query Overpass for ways tagged with highway inside the "
            "requested bbox"
        ),
        "Convert provider geometries to a DTCC RoadNetwork in EPSG:3006",
        (
            "Preserve edge-aligned OpenStreetMap tags such as highway and "
            "oneway when available"
        ),
        "Return a RoadNetwork or serialize protobuf bytes when format='pb'",
    ]
    presentation_headline = "OpenStreetMap Road Network"
    presentation_summary = (
        "Line-segment road graph from OpenStreetMap highway ways for the "
        "requested bounds."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "The dataset returns road centerline segments as a RoadNetwork. "
                "Vertices and edges describe graph topology; edge attributes "
                "carry provider tags when the download path supplies them."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "The highway tag is an OpenStreetMap road-class label, not a "
                "traffic capacity model. Segment lengths are geometric lengths "
                "in the output CRS."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Coverage, classification, turn restrictions, private access, "
                "and completeness depend on current OpenStreetMap contributors "
                "and Overpass/cache availability."
            ),
        },
    ]
    key_points = [
        "Currently source='OSM' is the only public roads source",
        "Road classes come from OpenStreetMap highway tags when present",
        "Coordinates are requested and returned through the EPSG:3006 path",
        "The dataset is a raw network source for derived analyses such as space_syntax",
    ]
    presentation_legend = {
        "title": "Road network attributes",
        "entries": [
            {
                "label": "highway",
                "meaning": (
                    "OpenStreetMap road-class tag, such as residential or "
                    "primary"
                ),
            },
            {
                "label": "oneway",
                "meaning": "OpenStreetMap one-way tag when present in provider data",
            },
            {
                "label": "length",
                "meaning": "Road segment length in the output CRS",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "lines",
        "default_crs": "EPSG:3006",
        "table_role": "network_context",
        "default_color_attribute": "highway",
    }
    presentation_warnings = [
        (
            "OpenStreetMap/ODbL source terms require review before redistributing "
            "generated packages."
        ),
        (
            "Overpass availability and rate limits can affect live downloads; "
            "cached data may not represent the latest OpenStreetMap edits."
        ),
        (
            "Road classes and access semantics are provider tags and should not "
            "be read as validated engineering classifications."
        ),
    ]
    presentation_limitations = [
        "Only the OpenStreetMap/Overpass path is exposed by the public dataset.",
        (
            "The dataset does not validate turn restrictions, lane counts, "
            "speeds, or traffic capacity."
        ),
        "Completeness and tagging quality vary by area and contributor activity.",
    ]

    def build(self, args: RoadsArgs):
        bounds = self.parse_bounds(args.bounds)
        roads: RoadNetwork = dtcc_core.io.data.download_roadnetwork(
            bounds=bounds,
            provider=args.source,
            epsg="3006",
        )
        if args.format == "pb":
            return roads.to_proto().SerializeToString()
        return roads


__all__ = ["RoadsArgs", "RoadsDataset"]
