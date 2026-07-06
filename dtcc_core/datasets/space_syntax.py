"""Road-network segment space syntax dataset."""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import Field, model_validator
from dtcc_core.builder.roadnetwork.space_syntax import (
    DEFAULT_SPACE_SYNTAX_MEASURES,
    SpaceSyntaxCost,
    SpaceSyntaxMeasure,
    analyze_space_syntax,
)

from .dataset import DatasetBaseArgs, DatasetDescriptor
from .providers import provider_entry


_RADIUS_UNITS = {
    "topological": "steps",
    "metric": "meters",
    "angular": "degrees",
}


class SpaceSyntaxArgs(DatasetBaseArgs):
    """Arguments for the road-network space syntax dataset."""

    source: Literal["OSM"] = Field("OSM", description="Road data source")
    cost: SpaceSyntaxCost = Field(
        "topological",
        description=(
            "Segment transition cost: topological steps, metric center-to-center "
            "distance, or angular deflection in degrees."
        ),
    )
    radius: Optional[float] = Field(
        None,
        ge=0.0,
        description=(
            "Optional local search radius. Units follow radius_unit, or the "
            "selected cost model if radius_unit is omitted."
        ),
    )
    radius_unit: Optional[Literal["steps", "meters", "degrees"]] = Field(
        None,
        description=(
            "Optional radius unit. Must match cost: steps for topological, "
            "meters for metric, degrees for angular."
        ),
    )
    measures: tuple[SpaceSyntaxMeasure, ...] = Field(
        DEFAULT_SPACE_SYNTAX_MEASURES,
        description=(
            "Segment measures to compute: connectivity, reach, mean_depth, "
            "integration, and/or choice."
        ),
    )
    include_disconnected: bool = Field(
        True,
        description=(
            "Keep all disconnected road components in the analysis. If false, "
            "only the largest connected segment component is analysed."
        ),
    )
    normalize: bool = Field(
        True,
        description=(
            "Normalize integration and choice to approximately unitless values."
        ),
    )
    format: Optional[Literal["pb"]] = Field(
        None,
        description="Output format (pb for protobuf bytes)",
    )

    @model_validator(mode="after")
    def validate_space_syntax_options(self):
        if len(self.measures) == 0:
            raise ValueError("At least one space syntax measure must be requested.")
        expected_unit = _RADIUS_UNITS[self.cost]
        if self.radius_unit is not None and self.radius_unit != expected_unit:
            raise ValueError(
                f"radius_unit must be {expected_unit!r} when cost={self.cost!r}"
            )
        return self


class SpaceSyntaxDataset(DatasetDescriptor):
    """Segment-based road-network space syntax measures."""

    name = "space_syntax"
    title = "Road Space Syntax"
    description = (
        "Derived segment-based space syntax measures on a DTCC RoadNetwork. "
        "Downloads roads for the requested bounds, builds the dual segment graph, "
        "and attaches connectivity, reach, mean depth, integration, choice, and "
        "component attributes to the returned RoadNetwork."
    )
    ArgsModel = SpaceSyntaxArgs
    data_category = "derived"
    result_kind = "road_network"
    python_return_type = "dtcc_core.model.RoadNetwork"
    timeout_hint = 300
    provider = [provider_entry("dtcc-platform", role="processor")]
    source = [
        {
            "name": "roads dataset",
            "role": "upstream_dataset",
            "dataset": "roads",
            "source_terms_status": "inherits_requires_review",
        },
        {
            "name": "dtcc-core segment space syntax analyzer",
            "role": "derived_processor",
            "module": "dtcc_core.builder.roadnetwork.space_syntax",
        },
    ]
    derived_from = [
        {
            "name": "roads",
            "role": "upstream_dataset",
            "source": "OpenStreetMap/Overpass when source='OSM'",
        }
    ]
    license = (
        "Derived from OpenStreetMap data; review ODbL attribution, share-alike, "
        "and redistribution requirements before publishing derived packages."
    )
    collection_period = (
        "Derived on demand from the roads dataset response for the requested "
        "bounds; inherits the upstream OpenStreetMap/Overpass or cache timing."
    )
    default_crs = "EPSG:3006"
    data_types = [
        "road_network",
        "lines",
        "graph",
        "space_syntax_measures",
        "derived_analysis",
    ]
    geographic_coverage = (
        "Global OpenStreetMap coverage, constrained by requested bounds"
    )
    update_frequency = "derived on demand from OpenStreetMap source data"
    processing_steps = [
        "Fetch the upstream roads dataset for the requested bounds and source",
        "Treat each road segment as a node in a dual segment graph",
        "Connect segment-nodes when their road segments share an endpoint",
        (
            "Assign transition costs from the selected model: topological "
            "steps, metric center-to-center distance, or angular deflection "
            "in degrees"
        ),
        (
            "Optionally restrict shortest-path searches to the requested "
            "radius in the selected cost units"
        ),
        "Optionally drop disconnected components except the largest segment component",
        (
            "Compute requested connectivity, reach, mean_depth, integration, "
            "and/or choice measures"
        ),
        "Normalize integration and choice when normalize=True",
        "Attach edge-aligned space_syntax_* attributes to the returned RoadNetwork",
        "Return a RoadNetwork or serialize protobuf bytes when format='pb'",
    ]
    presentation_headline = "Road-Network Space Syntax"
    presentation_summary = (
        "Road network enriched with segment-based connectivity, reach, depth, "
        "integration, and choice measures."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each road segment receives graph measures computed on the dual "
                "segment graph. Adjacent segments in the source network become "
                "neighboring nodes in the analysis graph."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Connectivity counts immediate segment neighbors. Reach counts "
                "segments reachable within the optional radius. Mean depth is "
                "average shortest-path cost. Integration highlights segments "
                "that are shallow to many others, while choice highlights "
                "segments that lie on many shortest paths."
            ),
        },
        {
            "heading": "Cost models",
            "body": (
                "Topological cost uses segment-to-segment steps, metric cost "
                "uses center-to-center distance in meters, and angular cost "
                "uses deflection in degrees."
            ),
        },
    ]
    key_points = [
        "Measures are edge-aligned RoadNetwork attributes prefixed with space_syntax_",
        "radius units are steps, meters, or degrees depending on cost",
        (
            "include_disconnected=False assigns zero-valued measures outside "
            "the largest component"
        ),
        "normalize=True scales integration and choice for easier relative comparison",
        (
            "The analysis is descriptive graph analysis, not a traffic "
            "assignment or routing model"
        ),
    ]
    presentation_legend = {
        "title": "Space syntax measures",
        "entries": [
            {
                "label": "space_syntax_connectivity",
                "meaning": "Number of adjacent road segments in the dual graph",
            },
            {
                "label": "space_syntax_reach",
                "meaning": "Number of reachable segments within the optional radius",
            },
            {
                "label": "space_syntax_mean_depth",
                "meaning": "Average shortest-path cost to reachable segments",
            },
            {
                "label": "space_syntax_integration",
                "meaning": "Relative closeness-like accessibility measure",
            },
            {
                "label": "space_syntax_choice",
                "meaning": "Betweenness-like shortest-path choice measure",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "lines",
        "default_crs": "EPSG:3006",
        "table_role": "network_analysis",
        "default_color_attribute": "space_syntax_integration",
    }
    presentation_warnings = [
        (
            "The output inherits OpenStreetMap/ODbL source-term review "
            "requirements from the upstream roads dataset."
        ),
        (
            "Space syntax measures depend on graph construction choices and "
            "should not be interpreted as validated traffic, safety, or demand "
            "models."
        ),
        (
            "Small changes in OSM topology, disconnected components, radius, or "
            "cost model can materially change the resulting measures."
        ),
    ]
    presentation_limitations = [
        "Turn restrictions, lane counts, speeds, and capacities are not modeled.",
        (
            "Metric and angular costs are computed from segment geometry, not "
            "observed travel behavior."
        ),
        (
            "Disconnected components require careful interpretation, especially "
            "when include_disconnected=True."
        ),
        (
            "Domain validation for local planning use remains required before "
            "table-ready publication."
        ),
    ]

    def build(self, args: SpaceSyntaxArgs):
        import dtcc_core.datasets as datasets

        bounds = self.parse_bounds(args.bounds)
        roads = datasets.roads(bounds=bounds, source=args.source)
        result = analyze_space_syntax(
            roads,
            cost=args.cost,
            radius=args.radius,
            measures=args.measures,
            include_disconnected=args.include_disconnected,
            normalize=args.normalize,
        )
        if args.format == "pb":
            return result.to_proto().SerializeToString()
        return result


__all__ = ["SpaceSyntaxArgs", "SpaceSyntaxDataset"]
