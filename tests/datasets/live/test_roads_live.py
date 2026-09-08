"""Opt-in live checks for the OpenStreetMap roads dataset.

Run with ``DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k roads``.
"""

from __future__ import annotations

import pytest

import dtcc_core.datasets as datasets
from dtcc_core.model import RoadNetwork


GOTHENBURG_CORE_BBOX_EPSG3006 = (
    318_000.0,
    6_397_500.0,
    320_000.0,
    6_399_500.0,
)


@pytest.mark.live
def test_roads_live_overpass_structural_smoke() -> None:
    """Validate a small live Overpass road-network response when available."""
    try:
        roads = datasets.roads(bounds=GOTHENBURG_CORE_BBOX_EPSG3006)
    except RuntimeError as exc:
        if "Overpass" in str(exc) or "endpoint" in str(exc).lower():
            pytest.skip(f"Overpass roads live case unavailable: {exc}")
        raise

    assert isinstance(roads, RoadNetwork)
    assert len(roads.vertices) > 0
    assert len(roads.edges) > 0
    assert len(roads.length) == len(roads.edges)
    if "highway" in roads.attributes:
        assert len(roads.attributes["highway"]) == len(roads.edges)
