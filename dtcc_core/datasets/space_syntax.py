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
        description="Normalize integration and choice to approximately unitless values.",
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
