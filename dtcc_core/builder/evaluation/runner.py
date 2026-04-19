from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, List, Optional

from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings
from dtcc_core.builder.geometry_builders.surface import extrude_surface
from dtcc_core.builder.evaluation.dataset import EvalDataset, EvalSample
from dtcc_core.builder.evaluation.metrics import (
    stage_outcome,
    roof_type_prediction,
    roof_type_correct,
    detected_plane_count,
    plane_count_delta,
    point_coverage_ratio,
    top2_area_share,
    slope_symmetry_error,
    is_watertight,
    ridge_height_error,
    eave_height_error,
    semantic_overlap,
)


@dataclass
class BuildingResult:
    building_id: str
    stage_outcome: str
    roof_type_predicted: Optional[str]
    roof_type_truth: str
    roof_type_correct: bool
    confidence: Optional[float]
    classifier_source: Optional[str]
    fallback_reason_attr: Optional[str]
    plane_count: int
    plane_count_delta: Optional[int]
    coverage_ratio: Optional[float]
    top2_area_share: Optional[float]
    symmetry_error_deg: Optional[float]
    watertight: bool
    ridge_height_error: Optional[float]
    eave_height_error: Optional[float]
    semantic_iou: Dict[str, float]
    timings: Dict[str, float]
    total_ms: float


class Runner:
    """Runs the LoD2 pipeline across an EvalDataset and computes per-building metrics."""

    def __init__(self, build_lod1_if_missing: bool = True, config=None):
        self._build_lod1_if_missing = build_lod1_if_missing
        self._config = config

    def run(self, dataset: EvalDataset) -> List[BuildingResult]:
        results: List[BuildingResult] = []
        for sample in dataset:
            results.append(self._run_one(sample))
        return results

    def _run_one(self, sample: EvalSample) -> BuildingResult:
        building = self._sample_to_building(sample)
        diagnostics: Dict[str, Dict] = {}
        total_start = time.perf_counter()
        build_lod2_buildings([building], diagnostics_out=diagnostics, config=self._config)
        total_ms = (time.perf_counter() - total_start) * 1000.0

        rec = diagnostics.get(building.id, {})
        timings = rec.get("timings", {})
        filtered_n = rec.get("filtered_point_count", 0)

        iou: Dict[str, float] = {}
        if sample.ground_truth.geometry is not None and building.lod2 is not None:
            raw_iou = semantic_overlap(building.lod2, sample.ground_truth.geometry)
            iou = {k.name: v for k, v in raw_iou.items()}

        return BuildingResult(
            building_id=sample.building_id,
            stage_outcome=stage_outcome(building),
            roof_type_predicted=roof_type_prediction(building),
            roof_type_truth=sample.ground_truth.roof_type,
            roof_type_correct=roof_type_correct(building, sample.ground_truth),
            confidence=building.attributes.get("roof_confidence"),
            classifier_source=building.attributes.get("classifier_source"),
            fallback_reason_attr=building.attributes.get("fallback_reason"),
            plane_count=detected_plane_count(rec),
            plane_count_delta=plane_count_delta(rec, sample.ground_truth),
            coverage_ratio=point_coverage_ratio(rec, filtered_n) if filtered_n else None,
            top2_area_share=top2_area_share(rec),
            symmetry_error_deg=slope_symmetry_error(rec),
            watertight=is_watertight(building),
            ridge_height_error=ridge_height_error(rec, sample.ground_truth),
            eave_height_error=eave_height_error(rec, sample.ground_truth),
            semantic_iou=iou,
            timings=dict(timings),
            total_ms=total_ms,
        )

    def _sample_to_building(self, sample: EvalSample) -> Building:
        b = Building(id=sample.building_id)
        b.add_geometry(sample.footprint, GeometryType.LOD0)
        wall_height = sample.ground_truth.eave_height or 5.0
        if self._build_lod1_if_missing:
            lod1 = extrude_surface(sample.footprint, wall_height)
            b.add_geometry(lod1, GeometryType.LOD1)
        b.add_geometry(sample.point_cloud, GeometryType.POINT_CLOUD)
        b.attributes["ground_height"] = sample.ground_truth.ground_height
        b.attributes["height"] = wall_height
        return b
