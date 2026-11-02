"""Live checks that each downloading dataset still returns data upstream.

Each case builds one dataset over a small Gothenburg bbox, so a failure points
at the dataset and the upstream service behind it. SMHI, roads, orthophoto and
transit datasets have their own live test files. Transient upstream failures
skip, as described in ``live_policy``.

Downloads are cached, so a warm local cache can answer a case without
contacting its upstream. CI starts with an empty cache.

Run with ``DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

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
from tests.datasets.live.live_policy import skip_if_transient

pytestmark = pytest.mark.live

# 200 m x 200 m in central Gothenburg (EPSG:3006), inside one LiDAR tile.
SMALL_BBOX = (319_700.0, 6_397_800.0, 319_900.0, 6_398_000.0)


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
        is_populated=lambda trees: len(trees) > 0,
    ),
    LiveDatasetCase(
        dataset="space_syntax",
        upstream="OpenStreetMap Overpass",
        expected_type=RoadNetwork,
        is_populated=lambda roads: len(roads.edges) > 0,
    ),
    LiveDatasetCase(
        dataset="deso",
        upstream="SCB DeSO geometry + statistics API",
        expected_type=DeSO,
        # SCB suppresses some area values (NaN, stored as None), so one
        # published value is enough to show the statistics API answered.
        is_populated=lambda deso: any(
            area.attributes.get("population_total") is not None for area in deso.areas
        ),
        kwargs={"statistics": ["population"]},
    ),
]


@pytest.mark.parametrize("case", CASES, ids=[case.id for case in CASES])
def test_dataset_live(case: LiveDatasetCase) -> None:
    """Build the dataset from live upstream data and check it is not empty."""
    build = getattr(datasets, case.dataset)
    # strict_live makes datasets that support it raise instead of returning a
    # partial result. The others ignore it and raise on upstream failures anyway.
    try:
        result = build(bounds=SMALL_BBOX, strict_live=True, **case.kwargs)
    except Exception as exc:
        skip_if_transient(exc, f"{case.id} ({case.upstream})")
        raise

    assert isinstance(result, case.expected_type), (
        f"{case.id} ({case.upstream}) returned {type(result).__name__}, "
        f"expected {case.expected_type.__name__}"
    )
    assert case.is_populated(result), (
        f"{case.id} ({case.upstream}) returned an empty result"
    )
