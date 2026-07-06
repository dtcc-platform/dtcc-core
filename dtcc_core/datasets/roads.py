import dtcc_core
from dtcc_core.model import RoadNetwork
from pydantic import Field
from typing import Literal, Optional

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_entry


class RoadsArgs(DatasetBaseArgs):
    source: Literal["OSM"] = Field("OSM", description="Data source")
    format: Optional[Literal["pb"]] = Field(
        None, description="Output format (pb for protobuf bytes)"
    )


class RoadsDataset(DatasetDescriptor):
    name = "roads"
    description = "Road network data from OSM/Overpass."
    ArgsModel = RoadsArgs
    data_category = "raw"
    result_kind = "road_network"
    python_return_type = "dtcc_core.model.RoadNetwork"
    provider = [provider_entry("openstreetmap")]
    source = ["OpenStreetMap road network data via Overpass"]
    license = "Review OpenStreetMap/ODbL terms before redistribution."
    default_crs = "EPSG:3006"
    geographic_coverage = "Global OpenStreetMap coverage, constrained by requested bounds"
    update_frequency = "depends on OpenStreetMap edits and Overpass availability"
    processing_steps = ["Download road network data for requested bounds"]
    presentation_summary = "Road network for the requested bounds from OpenStreetMap."

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
