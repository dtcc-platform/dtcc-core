from __future__ import annotations

from pydantic import Field

from dtcc_core.common.progress import ProgressTracker
from dtcc_core.model import City

from ._city_mesh_common import prepare_footprint_city_from_bounds
from .dataset import DatasetBaseArgs, DatasetDescriptor


class CityFootprintsArgs(DatasetBaseArgs):
    """Arguments for the cleaned city footprints dataset."""

    min_building_detail: float = Field(
        0.5,
        description=(
            "Scale of footprint cleaning in meters: minimum feature size, "
            "fidelity scale and merge distance"
        ),
    )
    min_building_area: float = Field(
        15.0, description="Smallest cleaned building footprint area to include (m²)"
    )
    show_footprints: bool = Field(
        False,
        description="Whether to show a live matplotlib view of raw vs cleaned footprints",
    )


class CityFootprintsDataset(DatasetDescriptor, register=False):
    """Internal city footprints dataset; not registered publicly.

    A city tile whose building footprints are cleaned in unison with
    ``clean_footprints``, ready for meshing.

    Call it with keyword arguments defined by ``CityFootprintsArgs``.
    """

    name = "city_footprints"
    description = "City tile with cleaned building footprints, ready for meshing."
    ArgsModel = CityFootprintsArgs
    data_category = "derived"
    result_kind = "city_model"
    python_return_type = "dtcc_core.model.City"

    def _build_footprints_from_city(self, city: City, args: CityFootprintsArgs):
        return city.clean_footprints(
            min_building_detail=args.min_building_detail,
            min_area=args.min_building_area,
            show=args.show_footprints,
        )

    def build_from_city(self, city: City, **kwargs):
        args = self.validate(kwargs)
        return self._build_footprints_from_city(city, args)

    def build(self, args: CityFootprintsArgs):
        progress_phases = {
            "download_footprints": 0.20,
            "clean_footprints": 0.80,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            city = prepare_footprint_city_from_bounds(bounds, progress=progress)
            with progress.phase("clean_footprints", "Cleaning city footprints..."):
                return self._build_footprints_from_city(city, args)
