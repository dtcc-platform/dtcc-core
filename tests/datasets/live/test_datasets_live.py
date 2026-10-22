"""Live checks that each downloading dataset still returns data upstream.

Each case builds one dataset over a small Gothenburg bbox, so a failure points
at the dataset and the upstream service behind it. SMHI, roads and transit
datasets have their own live test files.

Run with ``DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.model import (
    BuildingCollection,
    City,
    DeSO,
    FootprintCollection,
    Mesh,
    PointCloud,
    RoadNetwork,
    TreeCollection,
)


pytestmark = pytest.mark.live

# 200 m x 200 m in central Gothenburg (EPSG:3006).
SMALL_BBOX = (319_900.0, 6_397_800.0, 320_100.0, 6_398_000.0)


@dataclass(frozen=True)
class LiveDatasetCase:
    """One dataset call and the upstream service it exercises."""

    dataset: str
    upstream: str
    expected_type: type
    is_populated: Callable[[Any], bool]
    kwargs: dict[str, Any] = field(default_factory=dict)

    @property
    def id(self) -> str:
        source = self.kwargs.get("source")
        return f"{self.dataset}-{source}" if source else self.dataset


CASES = [
    LiveDatasetCase(
        dataset="point_cloud",
        upstream="DTCC LiDAR server (Lantmäteriet)",
        expected_type=PointCloud,
        is_populated=lambda pc: len(pc.points) > 0,
    ),
    LiveDatasetCase(
        dataset="building_footprints",
        upstream="DTCC GeoPackage server (Lantmäteriet)",
        expected_type=FootprintCollection,
        is_populated=lambda fc: len(fc) > 0,
        kwargs={"source": "LM"},
    ),
    LiveDatasetCase(
        dataset="building_footprints",
        upstream="OpenStreetMap Overpass",
        expected_type=FootprintCollection,
        is_populated=lambda fc: len(fc) > 0,
        kwargs={"source": "OSM"},
    ),
    LiveDatasetCase(
        dataset="buildings",
        upstream="DTCC LiDAR + GeoPackage servers (Lantmäteriet)",
        expected_type=BuildingCollection,
        is_populated=lambda bc: len(bc) > 0,
    ),
    LiveDatasetCase(
        dataset="city",
        upstream="DTCC LiDAR + GeoPackage servers (Lantmäteriet)",
        expected_type=City,
        is_populated=lambda city: len(city.buildings) > 0,
    ),
    LiveDatasetCase(
        dataset="terrain_surface_mesh",
        upstream="DTCC LiDAR server (Lantmäteriet)",
        expected_type=Mesh,
        is_populated=lambda mesh: len(mesh.faces) > 0,
    ),
    LiveDatasetCase(
        dataset="trees",
        upstream="DTCC LiDAR server (Lantmäteriet)",
        expected_type=TreeCollection,
        is_populated=lambda trees: trees is not None,
    ),
    LiveDatasetCase(
        dataset="space_syntax",
        upstream="OpenStreetMap Overpass",
        expected_type=RoadNetwork,
        is_populated=lambda roads: len(roads.edges) > 0,
    ),
    LiveDatasetCase(
        dataset="deso",
        upstream="SCB DeSO",
        expected_type=DeSO,
        is_populated=lambda deso: len(deso) > 0,
    ),
]


@pytest.mark.parametrize("case", CASES, ids=[case.id for case in CASES])
def test_dataset_live(case: LiveDatasetCase) -> None:
    """Build the dataset from live upstream data and check it is not empty."""
    build = getattr(datasets, case.dataset)
    result = build(bounds=SMALL_BBOX, strict_live=True, **case.kwargs)

    assert isinstance(result, case.expected_type), (
        f"{case.id} ({case.upstream}) returned {type(result).__name__}, "
        f"expected {case.expected_type.__name__}"
    )
    assert case.is_populated(result), (
        f"{case.id} ({case.upstream}) returned an empty result"
    )
