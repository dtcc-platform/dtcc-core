from __future__ import annotations

from typing import Literal, Optional

from pydantic import Field

from dtcc_core.common.progress import ProgressTracker
from dtcc_core.model import City

from ._city_mesh_common import (
    condition_city_meshing_footprints,
    prepare_footprint_city_from_bounds,
)
from .dataset import DatasetBaseArgs, DatasetDescriptor


class CityFootprintsArgs(DatasetBaseArgs):
    max_mesh_size: Optional[float] = Field(
        10.0,
        description=(
            "Maximum target edge length in meters used by the meshing-footprint "
            "stage for size-dependent metadata and downstream contracts"
        ),
    )
    min_building_detail: float = Field(
        0.5, description="Minimum building feature size to resolve in meters"
    )
    min_building_area: float = Field(
        15.0, description="Smallest building footprint area to include (m²)"
    )
    merge_buildings: bool = Field(
        True, description="Whether to merge adjacent building footprints"
    )
    merge_tolerance: float = Field(
        0.5, description="Distance tolerance for merging building footprints"
    )
    cleaning_diagnostics: bool = Field(
        True,
        description="Whether to emit detailed footprint-conditioning diagnostics",
    )
    show_footprints: bool = Field(
        False,
        description="Whether to show a live matplotlib view of raw vs conditioned footprints",
    )
    footprint_cleaning_plot_block: bool = Field(
        True,
        description="Whether the optional footprint cleaning plot should block until the window is closed",
    )
    pipeline_mode: Literal["strict"] = Field(
        "strict",
        description="Meshing pipeline mode",
    )


class CityFootprintsDataset(DatasetDescriptor):
    name = "city_footprints"
    description = (
        "Meshing-ready conditioned building footprints prepared from a city tile."
    )
    ArgsModel = CityFootprintsArgs
    data_category = "derived"
    result_kind = "city_meshing_footprints"
    python_return_type = "dtcc_core.datasets._city_mesh_common.CityMeshingFootprints"

    def _build_footprints_from_city(self, city: City, args: CityFootprintsArgs):
        return condition_city_meshing_footprints(
            city,
            min_building_detail=args.min_building_detail,
            min_building_area=args.min_building_area,
            merge_tolerance=args.merge_tolerance,
            merge_buildings=args.merge_buildings,
            max_mesh_size=args.max_mesh_size,
            cleaning_diagnostics=args.cleaning_diagnostics,
            show_footprints=args.show_footprints,
            footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
            pipeline_mode=args.pipeline_mode,
        )

    def build_from_city(self, city: City, **kwargs):
        args = self.validate(kwargs)
        return self._build_footprints_from_city(city, args)

    def build(self, args: CityFootprintsArgs):
        progress_phases = {
            "download_footprints": 0.20,
            "condition_footprints": 0.80,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            city = prepare_footprint_city_from_bounds(
                bounds,
                progress=progress,
            )
            with progress.phase(
                "condition_footprints",
                "Conditioning city footprints...",
            ):
                return self._build_footprints_from_city(city, args)
