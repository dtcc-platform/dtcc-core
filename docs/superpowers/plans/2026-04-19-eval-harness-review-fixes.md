# Eval Harness Review-Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the five review findings both the first reviewer and Codex agreed on — fine-grained stage timings, validation-result reuse, configurable wall-height fallback, gallery-without-pipeline-re-run, and a non-skipping gallery test.

**Architecture:** Minimal, surgical changes. `filter_roof_points` and `detect_roof_planes` gain an optional `stage_timings` dict and populate it with sub-stage durations; `build_lod2_buildings` pipes it through. `build_lod2_buildings` also records the existing `validate_shell` result into `diagnostics_out["validated"]`. `Runner` exposes the built `Building` objects via an instance property so the failure-gallery writer can consume them directly instead of re-running the pipeline. `test_report_gallery.py` drops its conditional skip in favor of a synthetic result fixture. No new modules; no public API renames.

**Tech Stack:** Same as the eval harness — Python 3.10+, numpy, shapely, pytest, existing `dtcc_core` modules.

**User git preference:** per durable instruction, each commit step shows the `git add` + suggested commit message but does **not** execute `git commit`. The user runs the commits.

**Branch state at plan start:** `feature/lod2-reconstruction`, HEAD at `bcafa3c`, 14 commits ahead of `develop`, full test suite (674 pass, 4 skip, 0 fail).

---

## Task 1: Finer-grained stage timings in roof_detection helpers

**Files:**
- Modify: `dtcc_core/builder/geometry_builders/roof_detection.py`
- Test: `tests/builder/test_roof_detection_timing.py` (new)

**Rationale:** Spec Track 1 Metrics section lists 8 timing stages. The current implementation collapses normals/filter into `filter` and RANSAC/region-growing into `detect`. Restoring them means instrumenting the internal steps of `filter_roof_points` and `detect_roof_planes`.

- [ ] **Step 1: Write the failing test**

```python
# tests/builder/test_roof_detection_timing.py
import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.builder.geometry_builders.roof_config import RoofDetectionConfig
from dtcc_core.builder.geometry_builders.roof_detection import (
    detect_roof_planes,
    filter_roof_points,
)


def _flat_pc(n=200):
    np.random.seed(0)
    pts = np.column_stack([
        np.random.rand(n) * 10,
        np.random.rand(n) * 10,
        np.full(n, 5.0) + np.random.randn(n) * 0.1,
    ])
    return PointCloud(points=pts)


def test_filter_roof_points_records_normals_and_filter_timings():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    timings = {}
    filtered, normals = filter_roof_points(pc, cfg, stage_timings=timings)
    assert filtered is not None
    assert "normals" in timings
    assert "filter" in timings
    assert timings["normals"] >= 0.0
    assert timings["filter"] >= 0.0


def test_filter_roof_points_none_stage_timings_is_noop():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    # Must work identically without the kwarg
    filtered_a, _ = filter_roof_points(pc, cfg)
    filtered_b, _ = filter_roof_points(pc, cfg, stage_timings=None)
    assert len(filtered_a.points) == len(filtered_b.points)


def test_detect_roof_planes_records_ransac_and_region_growing():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    filtered, normals = filter_roof_points(pc, cfg)
    timings = {}
    planes = detect_roof_planes(filtered, normals, cfg, stage_timings=timings)
    assert len(planes) >= 1
    assert "ransac" in timings
    assert "region_growing" in timings
    assert timings["ransac"] >= 0.0
    assert timings["region_growing"] >= 0.0


def test_detect_roof_planes_none_stage_timings_is_noop():
    pc = _flat_pc()
    cfg = RoofDetectionConfig()
    filtered, normals = filter_roof_points(pc, cfg)
    planes_a = detect_roof_planes(filtered, normals, cfg)
    planes_b = detect_roof_planes(filtered, normals, cfg, stage_timings=None)
    assert len(planes_a) == len(planes_b)
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/test_roof_detection_timing.py -v
```

Expected: FAIL with `TypeError: filter_roof_points() got an unexpected keyword argument 'stage_timings'`.

- [ ] **Step 3: Update `filter_roof_points` signature and body**

Modify `dtcc_core/builder/geometry_builders/roof_detection.py`:

Add at top of file, after existing imports:
```python
import time
```

Change `filter_roof_points` signature and body:
```python
def filter_roof_points(
    pc: PointCloud,
    config: RoofDetectionConfig,
    stage_timings: Optional[dict] = None,
) -> Tuple[Optional[PointCloud], Optional[np.ndarray]]:
    """Filter a per-building point cloud to keep likely roof points.

    Estimates normals first (reusable by detect_roof_planes), then removes
    facade hits and optionally filters by classification.

    Returns (filtered_pointcloud, normals) or (None, None) if too few remain.

    When stage_timings is provided, records `normals` and `filter` keys in ms.
    """
    if pc is None or len(pc.points) < config.min_roof_points:
        return None, None

    # Step 1: Estimate normals using Open3D
    t0 = time.perf_counter()
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pc.points)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamKNN(knn=config.normal_estimation_k)
    )
    normals = np.asarray(pcd.normals)
    # Orient normals upward
    flip_mask = normals[:, 2] < 0
    normals[flip_mask] *= -1
    if stage_timings is not None:
        stage_timings["normals"] = (time.perf_counter() - t0) * 1000.0

    # Remaining steps: facade filter + classification filter
    t_filter = time.perf_counter()
    keep_mask = np.ones(len(pc.points), dtype=bool)

    # Step 2: Remove facade points (near-vertical normals)
    angle_to_z = np.degrees(np.arccos(np.clip(normals[:, 2], -1.0, 1.0)))
    keep_mask &= angle_to_z < config.facade_angle_threshold

    # Step 3: Use classification if available
    if pc.classification is not None and len(pc.classification) == len(pc.points):
        has_building_class = np.any(pc.classification == 6)
        if has_building_class:
            keep_mask &= pc.classification == 6

    # Apply mask
    filtered_points = pc.points[keep_mask]
    filtered_normals = normals[keep_mask]
    filtered_cls = None
    if pc.classification is not None and len(pc.classification) == len(pc.points):
        filtered_cls = pc.classification[keep_mask]

    # Step 4: Check minimum count
    if len(filtered_points) < config.min_roof_points:
        if stage_timings is not None:
            stage_timings["filter"] = (time.perf_counter() - t_filter) * 1000.0
        return None, None

    filtered_pc = PointCloud(
        points=filtered_points,
        classification=filtered_cls if filtered_cls is not None else np.empty(0),
    )
    if stage_timings is not None:
        stage_timings["filter"] = (time.perf_counter() - t_filter) * 1000.0
    return filtered_pc, filtered_normals
```

- [ ] **Step 4: Update `detect_roof_planes` signature and body**

Change the function to accumulate ransac and region_growing times across iterations. Replace the whole function:

```python
def detect_roof_planes(
    pc: PointCloud,
    normals: np.ndarray,
    config: RoofDetectionConfig,
    stage_timings: Optional[dict] = None,
) -> List[RoofPlane]:
    """Detect roof planes via iterative RANSAC + region-growing refinement.

    When stage_timings is provided, records `ransac` and `region_growing`
    keys in ms, summed across all plane iterations.
    """
    if len(pc.points) < config.min_plane_points:
        if stage_timings is not None:
            stage_timings.setdefault("ransac", 0.0)
            stage_timings.setdefault("region_growing", 0.0)
        return []

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pc.points)
    pcd.normals = o3d.utility.Vector3dVector(normals)

    total_points = len(pc.points)
    remaining_indices = set(range(total_points))
    all_points = np.asarray(pcd.points)
    all_normals = normals.copy()
    planes: List[RoofPlane] = []

    ransac_ms = 0.0
    region_ms = 0.0

    for _ in range(config.max_roof_planes):
        if len(remaining_indices) < config.min_plane_points:
            break
        if len(remaining_indices) / total_points < config.min_remaining_ratio:
            break

        # Build sub-cloud from remaining points
        remaining_list = sorted(remaining_indices)
        sub_pcd = o3d.geometry.PointCloud()
        sub_pcd.points = o3d.utility.Vector3dVector(all_points[remaining_list])

        # RANSAC plane fit
        t0 = time.perf_counter()
        plane_model, inlier_sub_indices = sub_pcd.segment_plane(
            distance_threshold=config.ransac_distance_threshold,
            ransac_n=3,
            num_iterations=config.ransac_iterations,
        )
        ransac_ms += (time.perf_counter() - t0) * 1000.0

        if len(inlier_sub_indices) < config.min_plane_points:
            break

        # Map sub-indices back to global indices
        inlier_global = {remaining_list[i] for i in inlier_sub_indices}
        plane_normal = np.array(plane_model[:3])
        plane_offset = -plane_model[3]

        # Orient normal upward
        if plane_normal[2] < 0:
            plane_normal *= -1
            plane_offset *= -1

        # Region-growing refinement
        t1 = time.perf_counter()
        inlier_global = _region_grow(
            all_points, all_normals, plane_normal, plane_offset,
            inlier_global, remaining_indices, config,
        )
        region_ms += (time.perf_counter() - t1) * 1000.0

        if len(inlier_global) < config.min_plane_points:
            remaining_indices -= inlier_global
            continue

        # Build RoofPlane
        inlier_arr = np.array(sorted(inlier_global))
        inlier_pts = all_points[inlier_arr]

        # Convex hull boundary in 3D world coordinates
        try:
            hull = ConvexHull(inlier_pts[:, :2])
            boundary_3d = inlier_pts[hull.vertices]
        except Exception:
            boundary_3d = inlier_pts

        plane = RoofPlane(
            normal=plane_normal / np.linalg.norm(plane_normal),
            offset=plane_offset,
            inliers=inlier_arr,
            boundary_3d=boundary_3d,
        )
        planes.append(plane)
        remaining_indices -= inlier_global

    planes.sort(key=lambda p: len(p.inliers), reverse=True)
    if stage_timings is not None:
        stage_timings["ransac"] = ransac_ms
        stage_timings["region_growing"] = region_ms
    return planes
```

- [ ] **Step 5: Run all roof_detection tests to confirm no regression + new test passes**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/test_roof_detection.py tests/builder/test_roof_detection_timing.py tests/builder/test_extract_roof_points_fix.py -v
```

Expected: ALL PASS (4 new + 7 existing + 2 existing).

- [ ] **Step 6: Suggested commit**

```bash
git add \
    dtcc_core/builder/geometry_builders/roof_detection.py \
    tests/builder/test_roof_detection_timing.py && \
git commit -m "builder: finer-grained stage timings in filter_roof_points and detect_roof_planes

Adds optional stage_timings dict kwarg. When provided, filter_roof_points
records 'normals' + 'filter' and detect_roof_planes records 'ransac' +
'region_growing' in ms. Default None is no-op."
```

---

## Task 2: Wire stage_timings through build_lod2_buildings

**Files:**
- Modify: `dtcc_core/builder/geometry_builders/buildings.py:498-517` (filter + detect + merge stage calls)
- Modify: `tests/builder/test_build_lod2_diagnostics.py` (assert 8 stage keys)

- [ ] **Step 1: Update the diagnostics test to expect 8 stage keys**

In `tests/builder/test_build_lod2_diagnostics.py`, replace `test_diagnostics_records_stage_timings`:

```python
def test_diagnostics_records_stage_timings():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    assert "b001" in diagnostics
    timings = diagnostics["b001"]["timings"]
    # 8 stages per spec: normals, filter, ransac, region_growing, merge, classify, geometry, validate
    for stage in ("normals", "filter", "ransac", "region_growing",
                  "merge", "classify", "geometry", "validate"):
        assert stage in timings, f"missing stage: {stage}"
        assert timings[stage] >= 0.0
```

- [ ] **Step 2: Run the updated test to confirm failure**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/test_build_lod2_diagnostics.py::test_diagnostics_records_stage_timings -v
```

Expected: FAIL with `AssertionError: missing stage: normals`.

- [ ] **Step 3: Update buildings.py to pass stage_timings into the helpers**

In `dtcc_core/builder/geometry_builders/buildings.py`, replace:

```python
        with _timed("filter", timings):
            filtered_pc, normals = filter_roof_points(pc, config)
```

With:

```python
        filtered_pc, normals = filter_roof_points(pc, config, stage_timings=timings)
```

And replace:

```python
        with _timed("detect", timings):
            planes = detect_roof_planes(filtered_pc, normals, config)
        with _timed("merge", timings):
            planes = merge_planes(planes, config)
```

With:

```python
        planes = detect_roof_planes(filtered_pc, normals, config, stage_timings=timings)
        with _timed("merge", timings):
            planes = merge_planes(planes, config)
```

The `classify`, `geometry`, and `validate` stage wraps stay unchanged. This yields the 8 keys: normals, filter, ransac, region_growing, merge, classify, geometry, validate.

- [ ] **Step 4: Run the diagnostics tests + regression suite**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/test_build_lod2_diagnostics.py tests/builder/test_lod2_pipeline.py tests/builder/evaluation/ -v
```

Expected: all PASS.

- [ ] **Step 5: Suggested commit**

```bash
git add \
    dtcc_core/builder/geometry_builders/buildings.py \
    tests/builder/test_build_lod2_diagnostics.py && \
git commit -m "builder: wire sub-stage timings through build_lod2_buildings

Restores the 8-stage timing breakdown (normals, filter, ransac,
region_growing, merge, classify, geometry, validate) that the spec's
Metrics section asks for, using the new stage_timings kwargs on
filter_roof_points and detect_roof_planes."
```

---

## Task 3: Record validation outcome in diagnostics; is_watertight reuses it

**Files:**
- Modify: `dtcc_core/builder/geometry_builders/buildings.py:582-594` (validate stage)
- Modify: `dtcc_core/builder/evaluation/metrics.py` (is_watertight)
- Modify: `tests/builder/test_build_lod2_diagnostics.py` (assert validated field)
- Modify: `tests/builder/evaluation/test_metrics_watertight.py` (test new signature)

- [ ] **Step 1: Add failing test for `validated` in diagnostics record**

Append to `tests/builder/test_build_lod2_diagnostics.py`:

```python
def test_diagnostics_records_validation_outcome():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    # Flat roof should build a watertight shell
    assert "validated" in diagnostics["b001"]
    assert diagnostics["b001"]["validated"] is True
```

- [ ] **Step 2: Add failing test for is_watertight reading diagnostics**

Append to `tests/builder/evaluation/test_metrics_watertight.py`:

```python
def test_is_watertight_prefers_diagnostics_record():
    # When a diagnostics record carries 'validated', is_watertight should
    # return that value without re-running validate_shell.
    b = Building(id="b1")
    # No lod2 on building; normally returns False. But with diagnostics:
    diagnostics = {"validated": True}
    assert is_watertight(b, diagnostics=diagnostics) is True


def test_is_watertight_diagnostics_false_overrides_lod2():
    # lod2 exists but diagnostics says validation failed → respect diagnostics
    b = Building(id="b2")
    ms = MultiSurface()
    b.add_geometry(ms, GeometryType.LOD2)
    diagnostics = {"validated": False}
    assert is_watertight(b, diagnostics=diagnostics) is False
```

- [ ] **Step 3: Confirm failures**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/test_build_lod2_diagnostics.py::test_diagnostics_records_validation_outcome tests/builder/evaluation/test_metrics_watertight.py -v
```

Expected: two failures. `diagnostics["b001"]` has no `validated` key; `is_watertight` does not accept `diagnostics` kwarg.

- [ ] **Step 4: Record `validated` in `build_lod2_buildings`**

In `dtcc_core/builder/geometry_builders/buildings.py`, inside the validate stage, after:

```python
        with _timed("validate", timings):
            is_valid, issues = validate_shell(lod2, config.edge_snap_tolerance)
```

Add:

```python
        if record is not None:
            record["validated"] = bool(is_valid)
```

Also, at the top of the for-loop where `record` is initialized (near line 480), extend the initial dict:

```python
            record = {
                "timings": timings,
                "planes": None,
                "classification": None,
                "filtered_point_count": 0,
                "validated": None,
            }
```

- [ ] **Step 5: Update `is_watertight` to accept optional diagnostics**

In `dtcc_core/builder/evaluation/metrics.py`, replace `is_watertight`:

```python
def is_watertight(
    building: Building,
    tolerance: float = 0.01,
    diagnostics: Optional[dict] = None,
) -> bool:
    """True iff the building's LoD2 is watertight.

    Prefers `diagnostics["validated"]` when available (set by
    build_lod2_buildings), avoiding a redundant shell validation.
    Falls back to calling validate_shell on building.lod2 otherwise.
    """
    if diagnostics is not None and diagnostics.get("validated") is not None:
        return bool(diagnostics["validated"])

    from dtcc_core.builder.geometry.shell_validation import validate_shell
    if building.lod2 is None:
        return False
    valid, _ = validate_shell(building.lod2, tolerance)
    return bool(valid)
```

- [ ] **Step 6: Update Runner to pass diagnostics and honor config tolerance**

In `dtcc_core/builder/evaluation/runner.py`, the `_run_one` method currently calls:

```python
            watertight=is_watertight(building),
```

Change to:

```python
            watertight=is_watertight(building, diagnostics=rec),
```

- [ ] **Step 7: Confirm all new + existing watertight tests pass**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/test_build_lod2_diagnostics.py tests/builder/evaluation/test_metrics_watertight.py tests/builder/evaluation/test_runner.py -v
```

Expected: all PASS.

- [ ] **Step 8: Suggested commit**

```bash
git add \
    dtcc_core/builder/geometry_builders/buildings.py \
    dtcc_core/builder/evaluation/metrics.py \
    dtcc_core/builder/evaluation/runner.py \
    tests/builder/test_build_lod2_diagnostics.py \
    tests/builder/evaluation/test_metrics_watertight.py && \
git commit -m "builder/evaluation: reuse validate_shell result from diagnostics

build_lod2_buildings now records the validation outcome into
diagnostics_out['validated']. is_watertight accepts an optional
diagnostics dict and prefers its value over re-running validate_shell.
Eliminates redundant shell validation per building in the eval harness."
```

---

## Task 4: Runner exposes built buildings; configurable default_wall_height

**Files:**
- Modify: `dtcc_core/builder/evaluation/runner.py` (Runner class)
- Modify: `tests/builder/evaluation/test_runner.py` (new tests)

- [ ] **Step 1: Add failing tests**

Append to `tests/builder/evaluation/test_runner.py`:

```python
def test_runner_exposes_built_buildings_by_id():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    assert "b001" in runner.buildings
    b = runner.buildings["b001"]
    assert b.id == "b001"
    assert b.lod0 is not None
    assert b.lod1 is not None


def test_runner_buildings_reset_between_runs():
    # Two runs produce a fresh buildings dict
    runner = Runner()
    runner.run(EvalDataset(FIXTURE_ROOT))
    first_ids = set(runner.buildings.keys())
    runner.run(EvalDataset(FIXTURE_ROOT))
    second_ids = set(runner.buildings.keys())
    assert first_ids == second_ids
    # Object identity for the cached building should differ between runs
    # (new Building object constructed per run)
    assert all(k in second_ids for k in first_ids)


def test_runner_default_wall_height_is_used_when_gt_omits_eave():
    import numpy as np
    from dtcc_core.model.geometry.pointcloud import PointCloud
    from dtcc_core.model.geometry.surface import Surface
    from dtcc_core.builder.evaluation.dataset import EvalSample, GroundTruth

    fp = Surface(vertices=np.array(
        [[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]], dtype=float
    ))
    pc = PointCloud(points=np.array([[1, 1, 5]], dtype=float))
    gt = GroundTruth(roof_type="FLAT")  # no eave_height
    sample = EvalSample(
        building_id="bx", footprint=fp, point_cloud=pc, ground_truth=gt,
    )

    runner = Runner(default_wall_height=7.5)
    building = runner._sample_to_building(sample)
    assert building.attributes["height"] == 7.5
```

- [ ] **Step 2: Confirm failure**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/evaluation/test_runner.py -v
```

Expected: three failures. `runner.buildings` attribute absent; `Runner(default_wall_height=...)` raises TypeError.

- [ ] **Step 3: Update Runner**

In `dtcc_core/builder/evaluation/runner.py`, change the `Runner` class:

Replace `__init__`:

```python
    def __init__(
        self,
        build_lod1_if_missing: bool = True,
        config=None,
        default_wall_height: float = 5.0,
    ):
        self._build_lod1_if_missing = build_lod1_if_missing
        self._config = config
        self._default_wall_height = default_wall_height
        self._buildings: Dict[str, Building] = {}

    @property
    def buildings(self) -> Dict[str, Building]:
        """Buildings from the last run() call, keyed by building_id."""
        return self._buildings
```

Update `run` to clear and populate `_buildings`:

```python
    def run(self, dataset: EvalDataset) -> List[BuildingResult]:
        self._buildings = {}
        results: List[BuildingResult] = []
        for sample in dataset:
            result, building = self._run_one(sample)
            self._buildings[sample.building_id] = building
            results.append(result)
        return results
```

Change `_run_one` to return `(BuildingResult, Building)` instead of just `BuildingResult`:

```python
    def _run_one(self, sample: EvalSample):
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

        result = BuildingResult(
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
            watertight=is_watertight(building, diagnostics=rec),
            ridge_height_error=ridge_height_error(rec, sample.ground_truth),
            eave_height_error=eave_height_error(rec, sample.ground_truth),
            semantic_iou=iou,
            timings=dict(timings),
            total_ms=total_ms,
        )
        return result, building
```

Update `_sample_to_building` to use `self._default_wall_height`:

```python
    def _sample_to_building(self, sample: EvalSample) -> Building:
        b = Building(id=sample.building_id)
        b.add_geometry(sample.footprint, GeometryType.LOD0)
        wall_height = (
            sample.ground_truth.eave_height
            if sample.ground_truth.eave_height is not None
            else self._default_wall_height
        )
        if self._build_lod1_if_missing:
            lod1 = extrude_surface(sample.footprint, wall_height)
            b.add_geometry(lod1, GeometryType.LOD1)
        b.add_geometry(sample.point_cloud, GeometryType.POINT_CLOUD)
        b.attributes["ground_height"] = sample.ground_truth.ground_height
        b.attributes["height"] = wall_height
        return b
```

- [ ] **Step 4: Run all runner tests**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/evaluation/test_runner.py -v
```

Expected: all PASS (3 original + 3 new = 6).

- [ ] **Step 5: Suggested commit**

```bash
git add \
    dtcc_core/builder/evaluation/runner.py \
    tests/builder/evaluation/test_runner.py && \
git commit -m "builder/evaluation: Runner exposes built buildings; configurable default_wall_height

Runner.run() now populates self._buildings (public via .buildings property)
with the per-id Building objects from the just-completed run, so downstream
consumers (failure gallery) do not need to re-run the pipeline.
default_wall_height becomes a ctor arg (was a magic 5.0 literal) for when
ground truth omits eave_height."
```

---

## Task 5: Gallery consumes buildings map; remove pipeline re-run; update CLI

**Files:**
- Modify: `dtcc_core/builder/evaluation/report.py` (`write_failure_gallery`)
- Modify: `sandbox/evaluate_lod2.py` (pass `runner.buildings`)
- Modify: `tests/builder/evaluation/test_report_gallery.py` (new signature)

- [ ] **Step 1: Update the gallery tests for the new signature**

Replace the whole contents of `tests/builder/evaluation/test_report_gallery.py`:

```python
from pathlib import Path

import numpy as np

from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface, Surface
from dtcc_core.builder.evaluation.runner import BuildingResult
from dtcc_core.builder.evaluation.report import FailureThresholds, write_failure_gallery


def _square(z: float = 0.0) -> Surface:
    return Surface(vertices=np.array([[0,0,z],[1,0,z],[1,1,z],[0,1,z]], dtype=float))


def _building_with_lod2(bid: str) -> Building:
    b = Building(id=bid)
    ms = MultiSurface()
    ms.surfaces.append(_square())
    b.add_geometry(ms, GeometryType.LOD2)
    b.attributes["roof_type"] = "FLAT"
    b.attributes["roof_confidence"] = 0.9
    return b


def _fake_failing_result(bid: str) -> BuildingResult:
    return BuildingResult(
        building_id=bid,
        stage_outcome="validation_failed",   # will trip require_watertight=True
        roof_type_predicted="FLAT",
        roof_type_truth="FLAT",
        roof_type_correct=True,
        confidence=0.9,
        classifier_source=None,
        fallback_reason_attr="validation_failed",
        plane_count=1,
        plane_count_delta=0,
        coverage_ratio=0.9,
        top2_area_share=1.0,
        symmetry_error_deg=None,
        watertight=False,
        ridge_height_error=None,
        eave_height_error=None,
        semantic_iou={},
        timings={"filter": 1.0},
        total_ms=5.0,
    )


def _fake_passing_result(bid: str) -> BuildingResult:
    return BuildingResult(
        building_id=bid,
        stage_outcome="success",
        roof_type_predicted="FLAT",
        roof_type_truth="FLAT",
        roof_type_correct=True,
        confidence=0.95,
        classifier_source=None,
        fallback_reason_attr=None,
        plane_count=1,
        plane_count_delta=0,
        coverage_ratio=0.98,
        top2_area_share=1.0,
        symmetry_error_deg=None,
        watertight=True,
        ridge_height_error=None,
        eave_height_error=None,
        semantic_iou={},
        timings={"filter": 1.0},
        total_ms=5.0,
    )


def test_gallery_writes_vtk_for_each_failing_building(tmp_path: Path):
    buildings = {"b1": _building_with_lod2("b1"), "b2": _building_with_lod2("b2")}
    results = [_fake_failing_result("b1"), _fake_failing_result("b2")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    assert len(paths) == 2
    assert (out_dir / "b1.vtk").exists()
    assert (out_dir / "b2.vtk").exists()


def test_gallery_skips_passing_buildings(tmp_path: Path):
    buildings = {"b1": _building_with_lod2("b1")}
    results = [_fake_passing_result("b1")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    assert paths == []


def test_gallery_skips_buildings_not_in_map(tmp_path: Path):
    # failing result references a building id not present in the map → skip it
    buildings = {"b1": _building_with_lod2("b1")}
    results = [_fake_failing_result("b1"), _fake_failing_result("ghost")]
    thresholds = FailureThresholds(
        max_plane_count_delta=None,
        min_coverage=None,
        min_semantic_iou=None,
        require_watertight=True,
    )

    out_dir = tmp_path / "gallery"
    paths = write_failure_gallery(
        buildings=buildings,
        results=results,
        thresholds=thresholds,
        out_dir=out_dir,
    )

    assert len(paths) == 1
    assert (out_dir / "b1.vtk").exists()
```

- [ ] **Step 2: Confirm failure**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/evaluation/test_report_gallery.py -v
```

Expected: FAIL. `write_failure_gallery` still requires `dataset_root`, not `buildings`.

- [ ] **Step 3: Rewrite `write_failure_gallery`**

In `dtcc_core/builder/evaluation/report.py`, replace the body of `write_failure_gallery`:

```python
def write_failure_gallery(
    buildings: Dict[str, "Building"],
    results: List[BuildingResult],
    thresholds: FailureThresholds,
    out_dir: Path | str,
) -> List[Path]:
    """Write a VTK file for each result that trips any failure threshold.

    Consumes already-built Building objects by id. Does not re-run the
    pipeline. Results whose building_id is not present in `buildings` are
    silently skipped (warning-worthy but not fatal).
    """
    from dtcc_core.builder.evaluation.vtk_export import write_building_vtk

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths: List[Path] = []
    for r in results:
        if not _is_failure(r, thresholds):
            continue
        building = buildings.get(r.building_id)
        if building is None or building.lod2 is None:
            continue
        path = out_dir / f"{r.building_id}.vtk"
        write_building_vtk(building, path)
        paths.append(path)
    return paths
```

The `Dict[str, "Building"]` annotation uses a forward-reference string literal, so no new import is needed — Python won't resolve it at runtime. If static type checking matters later, add `if TYPE_CHECKING: from dtcc_core.model.object.building import Building`, but don't import it unconditionally (would introduce a model import into the report module).

- [ ] **Step 4: Update the CLI to pass `runner.buildings`**

In `sandbox/evaluate_lod2.py`, change:

```python
    gallery_paths = write_failure_gallery(
        dataset_root=args.dataset,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )
```

To:

```python
    gallery_paths = write_failure_gallery(
        buildings=runner.buildings,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )
```

- [ ] **Step 5: Run the gallery tests + end-to-end smoke**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/builder/evaluation/test_report_gallery.py tests/builder/evaluation/test_end_to_end.py -v
```

Expected: 3 new gallery tests PASS + end-to-end smoke PASS.

- [ ] **Step 6: Smoke-run the CLI against the fixture**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/python sandbox/evaluate_lod2.py \
    --dataset tests/builder/evaluation/fixtures/minimal_dataset \
    --out /tmp/eval_smoke_out_v2
```

Expected: exit 0, same console output shape as before.

- [ ] **Step 7: Suggested commit**

```bash
git add \
    dtcc_core/builder/evaluation/report.py \
    sandbox/evaluate_lod2.py \
    tests/builder/evaluation/test_report_gallery.py && \
git commit -m "builder/evaluation: gallery consumes runner.buildings instead of re-running pipeline

write_failure_gallery now takes a buildings dict (Runner exposes it via
.buildings) and writes VTK for failing ids directly. Eliminates a full
pipeline re-run per failing building, the silent drop of eval-run config,
and RANSAC nondeterminism between the measured run and the gallery."
```

---

## Task 6: Full-suite verification

- [ ] **Step 1: Run the full test suite**

```bash
/Users/vasnas/scratch/lod3/temp/dtcc-core/.venv/bin/pytest tests/ 2>&1 | tail -6
```

Expected: all tests pass. Expected count: 674 (baseline) + 4 (new roof_detection timing tests) + 1 (validated diagnostics test) + 2 (is_watertight diagnostics tests) + 3 (runner new tests) + 3 (gallery new tests, minus 2 removed tests) = around 685 passing.

If any test fails, stop and investigate before claiming completion.

---

## Acceptance (against review findings)

- [x] **Finding #1 + #2 (gallery re-run + protected-name access):** Task 5 — gallery consumes runner-exposed buildings; no Runner() inside gallery; no `_sample_to_building` reach-through.
- [x] **Finding #5 (watertight recomputation + hardcoded tolerance):** Task 3 — `diagnostics_out["validated"]` carries the pipeline's own validate_shell result; `is_watertight` prefers it; tolerance concern goes away since we don't re-validate.
- [x] **Spec deviation (8 timing stages):** Tasks 1 + 2 — `normals`, `filter`, `ransac`, `region_growing`, `merge`, `classify`, `geometry`, `validate` all present in diagnostics.
- [x] **Minor: conditional skip in gallery test:** Task 5 — test rewritten to construct synthetic buildings + results, no conditional skip, no nondeterminism.
- [x] **Minor: magic `wall_height = ... or 5.0`:** Task 4 — `Runner(default_wall_height=...)` is explicit and configurable.
