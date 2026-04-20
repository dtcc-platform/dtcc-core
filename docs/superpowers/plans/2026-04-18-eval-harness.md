# Eval Harness + Baseline Commit — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Commit the currently-uncommitted rule-based LoD2 baseline as a coherent git history, then build an evaluation harness that runs the baseline across a labeled dataset and produces per-building metrics (stage outcomes, detection diagnostics, roof-type accuracy, watertightness, timings) plus a failure gallery.

**Architecture:** Two components land in sequence. First, the existing working-tree changes (roof detection/classification/geometry/validation modules + model changes + tests + sandbox scripts) are committed as a frozen baseline. Second, a new `dtcc_core/builder/evaluation/` package provides a pure-Python pipeline runner driven by a small dataset convention, with metrics as isolated pure functions and a reporter that writes CSV + summary JSON + VTK failure exports. One minor additive extension to `build_lod2_buildings()` — an opt-in `diagnostics_out` kwarg — makes per-stage timings and intermediate artifacts (`RoofPlane` list, `ClassificationResult`) available to the runner without duplicating pipeline orchestration.

**Tech Stack:** Python 3.10+, numpy, shapely, pytest, existing `dtcc_core` modules. No new runtime dependencies (reuses `meshio` and existing VTK export in `sandbox/export_vtk.py`).

**Scope note — spec → plans mapping:** The approved design doc (`docs/superpowers/specs/2026-04-18-ml-for-lod2-plus-design.md`) covers three independently-shippable pieces. This plan covers pieces 1 + 2 from the "Rollout sequencing" section (baseline commit + Track 1 eval harness). Track 2 Phase 1 (classifier parity) and Track 2 Phase 2 (Tier 2 expansion) will be written as separate plans once this lands and we have baseline metrics in hand.

**User git preference:** per user's durable instruction ("leave all git commits to me"), each commit step shows the exact `git add` and a suggested `git commit -m "..."` message but does **not** execute the commit. The user runs the commits at their discretion.

---

## Phase 0 — Commit the rule-based baseline

The uncommitted working tree has ~25 files changed or added for the rule-based LoD2 pipeline. This phase carves that into a coherent commit sequence so subsequent work has a frozen, versioned baseline to measure against.

### Task 0.1: Verify the uncommitted baseline is green

**Files:**
- No file edits; diagnostics only.

- [ ] **Step 1: Run the LoD2-related test suite against the working tree**

```bash
cd /Users/vasnas/scratch/lod3/temp/dtcc-core
pytest tests/builder/test_lod2_pipeline.py tests/builder/test_lod1_regression.py tests/builder/test_roof_config.py tests/builder/test_roof_detection.py tests/builder/test_roof_classification.py tests/builder/test_roof_geometry.py tests/builder/test_shell_validation.py tests/builder/test_extract_roof_points_fix.py tests/model/test_multisurface_semantics.py tests/model/test_pointcloud_normals.py -v
```

Expected: all tests PASS. If any fail, stop and fix before proceeding with commits.

- [ ] **Step 2: Sanity-run one validation script to confirm pipeline works end-to-end**

```bash
cd /Users/vasnas/scratch/lod3/temp/dtcc-core
python sandbox/validate_lod2.py
```

Expected: script runs to completion, no tracebacks. Output artifacts appear under `sandbox/output/`. Look at the console summary to confirm at least one building reached `lod2` successfully.

If either step fails, do not proceed. Report findings for repair.

### Task 0.2: Stage and propose commit breakdown

**Files:**
- No new files; suggested staging only. The USER runs each commit.

**Proposed commit breakdown (five commits):**

1. **Model changes** (enums, surface semantics, pointcloud normals, proto/generated)
2. **Builder — shell validation helper**
3. **Builder — roof pipeline (config, detection, classification, geometry, buildings integration)**
4. **IO — CityJSON per-surface semantics**
5. **Tests + sandbox scripts**

- [ ] **Step 1: Suggested commit 1 — model changes**

```bash
git add dtcc_core/model/enums.py \
        dtcc_core/model/geometry/surface.py \
        dtcc_core/model/geometry/pointcloud.py \
        dtcc_core/model/object/object.py \
        dtcc_core/model/mixins/city/builder_mixin.py \
        dtcc_core/model/dtcc_pb2.py \
        dtcc_core/proto/dtcc.proto
```

Suggested commit message:
```
model: add RoofType + SurfaceSemantic enums, wire semantic labels through model layer

- new RoofType, SurfaceSemantic enums in dtcc_core/model/enums.py
- Surface/MultiSurface carry per-surface semantic labels
- PointCloud gains normals field
- proto schema + generated code updated for new fields
- city builder mixin gains LoD2 build helper
```

User: `git commit -m "..."` when ready.

- [ ] **Step 2: Suggested commit 2 — shell validation helper**

```bash
git add dtcc_core/builder/geometry/shell_validation.py
```

Suggested commit message:
```
builder/geometry: add shell_validation for watertightness checks
```

- [ ] **Step 3: Suggested commit 3 — rule-based LoD2 roof pipeline**

```bash
git add dtcc_core/builder/geometry_builders/roof_config.py \
        dtcc_core/builder/geometry_builders/roof_detection.py \
        dtcc_core/builder/geometry_builders/roof_classification.py \
        dtcc_core/builder/geometry_builders/roof_geometry.py \
        dtcc_core/builder/geometry_builders/buildings.py
```

Suggested commit message:
```
builder: add rule-based LoD2 roof reconstruction (Tier 1)

- roof_config: RoofDetectionConfig + RoofPlane dataclass
- roof_detection: facade filtering, normal estimation (Open3D),
  RANSAC plane detection, region-growing refinement
- roof_classification: plane merging + flat/gabled/hipped classifier
- roof_geometry: watertight MultiSurface construction per roof type
- buildings: build_lod2_buildings() wires the full pipeline with
  fallback chain (insufficient_points, classification_failed,
  complex_footprint, geometry_construction_failed, validation_failed)
```

- [ ] **Step 4: Suggested commit 4 — CityJSON semantic mapping**

```bash
git add dtcc_core/io/cityjson/converters.py
```

Suggested commit message:
```
io/cityjson: emit per-surface RoofSurface/WallSurface/GroundSurface semantics
```

- [ ] **Step 5: Suggested commit 5 — tests and sandbox**

```bash
git add tests/builder/test_lod1_regression.py \
        tests/builder/test_lod2_pipeline.py \
        tests/builder/test_roof_classification.py \
        tests/builder/test_roof_config.py \
        tests/builder/test_roof_detection.py \
        tests/builder/test_roof_geometry.py \
        tests/builder/test_shell_validation.py \
        tests/builder/test_extract_roof_points_fix.py \
        tests/model/test_multisurface_semantics.py \
        tests/model/test_pointcloud_normals.py \
        sandbox/validate_lod2.py \
        sandbox/validate_lod2_real.py \
        sandbox/export_vtk.py
```

Suggested commit message:
```
tests + sandbox: coverage and validators for rule-based LoD2 pipeline
```

- [ ] **Step 6: Confirm a clean tree (apart from the brainstorm-dump .txt and sandbox/output/)**

```bash
git status
```

Expected: the `2026-04-13-...txt` transcript dump and `sandbox/output/` may remain untracked — these are ephemeral and should not be committed. All source/test/sandbox additions from the baseline are either committed or explicitly unstaged.

At this point, the rule-based baseline is a committed, reproducible history. Proceed to Phase 1.

### Task 0.3: Update pyproject.toml dependencies unrelated to Phase 0 (if any emerged)

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Check pyproject.toml diff**

```bash
git diff HEAD pyproject.toml
```

If the diff is empty, skip this task.
If the diff touches only a dependency needed by the baseline (e.g. `open3d` pinning), stage and commit it with the baseline:

```bash
git add pyproject.toml
```

Suggested commit message:
```
build: pin dependency versions required by rule-based LoD2 pipeline
```

---

## Phase 1 — Dataset convention + loader

### Task 1.1: Define the dataset data model

**Files:**
- Create: `dtcc_core/builder/evaluation/__init__.py`
- Create: `dtcc_core/builder/evaluation/dataset.py`
- Test: `tests/builder/evaluation/test_dataset_model.py` (new test directory)
- Create: `tests/builder/evaluation/__init__.py` (empty)

- [ ] **Step 1: Write the failing test**

```python
# tests/builder/evaluation/test_dataset_model.py
import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface
from dtcc_core.builder.evaluation.dataset import EvalSample, GroundTruth


def test_ground_truth_supports_minimum_fields():
    gt = GroundTruth(roof_type="FLAT")
    assert gt.roof_type == "FLAT"
    assert gt.ridge_height is None
    assert gt.eave_height is None
    assert gt.geometry is None


def test_ground_truth_supports_optional_heights():
    gt = GroundTruth(roof_type="GABLED", ridge_height=8.0, eave_height=5.0)
    assert gt.ridge_height == 8.0
    assert gt.eave_height == 5.0


def test_eval_sample_bundles_inputs_and_truth():
    fp = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float))
    pc = PointCloud(points=np.zeros((10, 3)))
    gt = GroundTruth(roof_type="FLAT")

    sample = EvalSample(
        building_id="b001",
        footprint=fp,
        point_cloud=pc,
        ground_truth=gt,
    )

    assert sample.building_id == "b001"
    assert sample.footprint is fp
    assert sample.point_cloud is pc
    assert sample.ground_truth.roof_type == "FLAT"
```

- [ ] **Step 2: Run the test and confirm it fails**

```bash
pytest tests/builder/evaluation/test_dataset_model.py -v
```

Expected: `ModuleNotFoundError: dtcc_core.builder.evaluation.dataset` (module does not exist yet).

- [ ] **Step 3: Implement the minimal module**

```python
# dtcc_core/builder/evaluation/__init__.py
"""Evaluation harness for the LoD2 reconstruction pipeline."""
```

```python
# dtcc_core/builder/evaluation/dataset.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import MultiSurface, Surface


@dataclass
class GroundTruth:
    """Ground-truth attributes for a single building in the eval dataset."""

    roof_type: str
    ground_height: float = 0.0           # Z of the footprint in world coordinates
    ridge_height: Optional[float] = None
    eave_height: Optional[float] = None
    geometry: Optional[MultiSurface] = None
    expected_plane_count: Optional[int] = None
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class EvalSample:
    """A single building (input + ground truth) for the eval harness to process."""

    building_id: str
    footprint: Surface
    point_cloud: PointCloud
    ground_truth: GroundTruth
```

Create empty `tests/builder/evaluation/__init__.py`.

- [ ] **Step 4: Run the test and confirm it passes**

```bash
pytest tests/builder/evaluation/test_dataset_model.py -v
```

Expected: all three tests PASS.

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/__init__.py \
        dtcc_core/builder/evaluation/dataset.py \
        tests/builder/evaluation/__init__.py \
        tests/builder/evaluation/test_dataset_model.py
```

Suggested commit message:
```
builder/evaluation: add dataset data model (EvalSample, GroundTruth)
```

### Task 1.2: Implement EvalDataset loader with a fixture directory

**Files:**
- Modify: `dtcc_core/builder/evaluation/dataset.py`
- Test: `tests/builder/evaluation/test_dataset_loader.py`
- Create fixture: `tests/builder/evaluation/fixtures/minimal_dataset/b001/footprint.geojson`
- Create fixture: `tests/builder/evaluation/fixtures/minimal_dataset/b001/points.npy` (numpy-saved for simplicity in tests; the real loader must also handle LAS)
- Create fixture: `tests/builder/evaluation/fixtures/minimal_dataset/b001/ground_truth.json`

- [ ] **Step 1: Create the fixture dataset**

```bash
mkdir -p tests/builder/evaluation/fixtures/minimal_dataset/b001
```

Write `tests/builder/evaluation/fixtures/minimal_dataset/b001/footprint.geojson`:

```json
{
  "type": "Feature",
  "geometry": {
    "type": "Polygon",
    "coordinates": [[[0,0],[10,0],[10,10],[0,10],[0,0]]]
  },
  "properties": {}
}
```

Write `tests/builder/evaluation/fixtures/minimal_dataset/b001/ground_truth.json`:

```json
{
  "roof_type": "FLAT",
  "ridge_height": null,
  "eave_height": 5.0,
  "expected_plane_count": 1
}
```

Generate the point cloud fixture by running (one-shot script, output committed as a binary fixture):

```python
# helper: generate tests/builder/evaluation/fixtures/minimal_dataset/b001/points.npy
import numpy as np
np.random.seed(0)
pts = np.column_stack([
    np.random.rand(200) * 10.0,
    np.random.rand(200) * 10.0,
    np.full(200, 5.0) + np.random.randn(200) * 0.05,
])
np.save("tests/builder/evaluation/fixtures/minimal_dataset/b001/points.npy", pts)
```

- [ ] **Step 2: Write the failing loader test**

```python
# tests/builder/evaluation/test_dataset_loader.py
from pathlib import Path

from dtcc_core.builder.evaluation.dataset import EvalDataset, EvalSample


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_eval_dataset_discovers_buildings():
    ds = EvalDataset(FIXTURE_ROOT)
    ids = [sample.building_id for sample in ds]
    assert ids == ["b001"]


def test_eval_dataset_len_matches_iter():
    ds = EvalDataset(FIXTURE_ROOT)
    assert len(ds) == 1


def test_eval_dataset_loads_sample_fields():
    ds = EvalDataset(FIXTURE_ROOT)
    sample = next(iter(ds))
    assert isinstance(sample, EvalSample)
    assert sample.building_id == "b001"
    assert sample.footprint.vertices.shape[1] == 3   # (N, 3), z from ground_truth.ground_height
    assert sample.point_cloud.points.shape == (200, 3)
    assert sample.ground_truth.roof_type == "FLAT"
    assert sample.ground_truth.eave_height == 5.0
    assert sample.ground_truth.expected_plane_count == 1
```

- [ ] **Step 3: Run tests and confirm failure**

```bash
pytest tests/builder/evaluation/test_dataset_loader.py -v
```

Expected: `ImportError` / `AttributeError`: `EvalDataset` not defined.

- [ ] **Step 4: Implement the loader**

Append to `dtcc_core/builder/evaluation/dataset.py`:

```python
import json
from pathlib import Path
from typing import Iterator

import numpy as np
import fiona
from shapely.geometry import shape as shapely_shape


class EvalDataset:
    """Iterates over buildings on disk using a simple per-building folder convention.

    Directory layout:
        <root>/
          <building_id>/
            footprint.geojson   OR   footprint.shp   OR   footprint.gpkg
            points.npy          OR   points.las      OR   points.laz
            ground_truth.json
    """

    def __init__(self, root: Path | str):
        self._root = Path(root)
        if not self._root.is_dir():
            raise FileNotFoundError(f"Dataset root not found: {self._root}")
        self._building_dirs = sorted(
            d for d in self._root.iterdir() if d.is_dir()
        )

    def __len__(self) -> int:
        return len(self._building_dirs)

    def __iter__(self) -> Iterator[EvalSample]:
        for bdir in self._building_dirs:
            yield self._load(bdir)

    def _load(self, bdir: Path) -> EvalSample:
        gt = self._load_ground_truth(bdir / "ground_truth.json")
        footprint = self._load_footprint(bdir, gt)
        points = self._load_points(bdir)
        pc = PointCloud(points=points)
        return EvalSample(
            building_id=bdir.name,
            footprint=footprint,
            point_cloud=pc,
            ground_truth=gt,
        )

    @staticmethod
    def _load_ground_truth(path: Path) -> GroundTruth:
        with open(path) as f:
            raw = json.load(f)
        known = {
            "roof_type", "ground_height", "ridge_height",
            "eave_height", "expected_plane_count",
        }
        return GroundTruth(
            roof_type=raw["roof_type"],
            ground_height=float(raw.get("ground_height", 0.0)),
            ridge_height=raw.get("ridge_height"),
            eave_height=raw.get("eave_height"),
            expected_plane_count=raw.get("expected_plane_count"),
            extra={k: v for k, v in raw.items() if k not in known},
        )

    @staticmethod
    def _load_footprint(bdir: Path, gt: GroundTruth) -> Surface:
        for candidate in ("footprint.geojson", "footprint.shp", "footprint.gpkg"):
            fp_path = bdir / candidate
            if fp_path.exists():
                break
        else:
            raise FileNotFoundError(f"No footprint found in {bdir}")
        z = gt.ground_height
        with fiona.open(fp_path) as src:
            feat = next(iter(src))
            geom = shapely_shape(feat["geometry"])
            coords = list(geom.exterior.coords)
            if coords[0] == coords[-1]:
                coords = coords[:-1]
            verts = np.array([[x, y, z] for x, y in coords], dtype=float)
        return Surface(vertices=verts)

    @staticmethod
    def _load_points(bdir: Path) -> np.ndarray:
        npy = bdir / "points.npy"
        if npy.exists():
            return np.load(npy)
        for ext in ("points.las", "points.laz"):
            las_path = bdir / ext
            if las_path.exists():
                import laspy
                with laspy.open(las_path) as src:
                    las = src.read()
                    return np.column_stack([las.x, las.y, las.z])
        raise FileNotFoundError(f"No point cloud found in {bdir}")
```

- [ ] **Step 5: Run tests and confirm pass**

```bash
pytest tests/builder/evaluation/test_dataset_loader.py -v
```

Expected: all three tests PASS.

- [ ] **Step 6: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/dataset.py \
        tests/builder/evaluation/test_dataset_loader.py \
        tests/builder/evaluation/fixtures/minimal_dataset/b001/footprint.geojson \
        tests/builder/evaluation/fixtures/minimal_dataset/b001/ground_truth.json \
        tests/builder/evaluation/fixtures/minimal_dataset/b001/points.npy
```

Suggested commit message:
```
builder/evaluation: add EvalDataset loader with folder-per-building convention

Supports footprint as geojson/shp/gpkg and point cloud as npy/las/laz.
Fixture dataset added under tests/ with one flat-roof sample.
```

---

## Phase 2 — Instrumentation hook on build_lod2_buildings

The existing pipeline already records stage outcomes into `building.attributes["fallback_reason"]` (values: `insufficient_points`, `classification_failed`, `complex_footprint`, `geometry_construction_failed`, `validation_failed`). The eval harness reads those directly.

It also needs:
1. Per-stage timing (not currently recorded).
2. The `RoofPlane` list from detection (used by detection-diagnostic metrics; not preserved after the pipeline returns).
3. The full `ClassificationResult` (confidence, ridge_line — useful for metrics and diagnostic output).

This phase adds one opt-in kwarg `diagnostics_out` to `build_lod2_buildings()` that populates a user-supplied dict when provided. Default `None` → zero overhead, no behavior change.

### Task 2.1: Add `diagnostics_out` kwarg to build_lod2_buildings

**Files:**
- Modify: `dtcc_core/builder/geometry_builders/buildings.py` (around `build_lod2_buildings`)
- Test: `tests/builder/test_build_lod2_diagnostics.py` (new)

- [ ] **Step 1: Write the failing test**

```python
# tests/builder/test_build_lod2_diagnostics.py
import numpy as np
import pytest

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface
from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings
from dtcc_core.builder.geometry_builders.surface import extrude_surface


def _flat_roof_building(bid="b001"):
    b = Building(id=bid)
    fp = Surface(vertices=np.array([[0,0,0],[10,0,0],[10,10,0],[0,10,0]], dtype=float))
    b.add_geometry(fp, GeometryType.LOD0)
    b.attributes["ground_height"] = 0.0
    b.attributes["height"] = 5.0
    b.add_geometry(extrude_surface(fp, 5.0), GeometryType.LOD1)
    np.random.seed(0)
    pts = np.column_stack([
        np.random.rand(200) * 10,
        np.random.rand(200) * 10,
        np.full(200, 5.0) + np.random.randn(200) * 0.1,
    ])
    b.add_geometry(PointCloud(points=pts), GeometryType.POINT_CLOUD)
    return b


def test_diagnostics_records_stage_timings():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    assert "b001" in diagnostics
    timings = diagnostics["b001"]["timings"]
    for stage in ("filter", "detect", "merge", "classify", "geometry", "validate"):
        assert stage in timings
        assert timings[stage] >= 0.0


def test_diagnostics_records_planes_for_successful_build():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    planes = diagnostics["b001"]["planes"]
    assert planes is not None
    assert len(planes) >= 1


def test_diagnostics_records_classification_result():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    clsr = diagnostics["b001"]["classification"]
    assert clsr is not None
    assert clsr.roof_type.name == "FLAT"


def test_diagnostics_none_is_noop():
    # baseline behavior unchanged when diagnostics_out=None
    b = _flat_roof_building()
    build_lod2_buildings([b])  # no diagnostics kwarg
    assert b.lod2 is not None
```

- [ ] **Step 2: Run test and confirm failure**

```bash
pytest tests/builder/test_build_lod2_diagnostics.py -v
```

Expected: `TypeError: build_lod2_buildings() got an unexpected keyword argument 'diagnostics_out'`.

- [ ] **Step 3: Implement the instrumentation hook**

In `dtcc_core/builder/geometry_builders/buildings.py`, add at module scope (above `build_lod2_buildings`):

```python
import time
from contextlib import contextmanager


@contextmanager
def _timed(stage: str, timings: dict | None):
    """Records elapsed ms into timings[stage] if dict provided; no-op otherwise."""
    if timings is None:
        yield
        return
    start = time.perf_counter()
    try:
        yield
    finally:
        timings[stage] = (time.perf_counter() - start) * 1000.0
```

Modify the signature of `build_lod2_buildings` to accept `diagnostics_out: Optional[Dict[str, Dict[str, Any]]] = None`.

Wrap each pipeline stage call with the timing context:

```python
# Inside build_lod2_buildings:
for building in buildings:
    timings: dict | None = None
    record: dict | None = None
    if diagnostics_out is not None:
        timings = {}
        record = {"timings": timings, "planes": None, "classification": None}
        diagnostics_out[building.id] = record

    # existing rebuild / lod1-null guards unchanged

    pc = building.point_cloud
    if pc is None or len(pc.points) < config.min_roof_points:
        _set_lod2_fallback(building, "insufficient_points", 0.0, fallback_to_flat, build_fallback_geometry)
        continue

    with _timed("filter", timings):
        filtered_pc, normals = filter_roof_points(pc, config)
    if filtered_pc is None:
        _set_lod2_fallback(building, "insufficient_points", 0.0, fallback_to_flat, build_fallback_geometry)
        continue

    with _timed("detect", timings):
        planes = detect_roof_planes(filtered_pc, normals, config)
    with _timed("merge", timings):
        planes = merge_planes(planes, config)
    if record is not None:
        record["planes"] = planes

    # footprint/area unchanged

    with _timed("classify", timings):
        result = classify_roof(
            planes, footprint_area, len(filtered_pc.points), config,
            footprint_vertices=fp_verts,
            all_points=filtered_pc.points,
        )
    if record is not None:
        record["classification"] = result

    # existing attribute writes + classification_failed / complex_footprint gates unchanged

    with _timed("geometry", timings):
        # existing geometry construction unchanged
        ...

    with _timed("validate", timings):
        is_valid, issues = validate_shell(lod2, config.edge_snap_tolerance)
    # existing validation handling unchanged
```

**Important:** every existing `continue` in the fallback branches must still execute. Make sure each fallback path populates `record["timings"]` with at least the stages that did run (those before the fallback), leaving later stages absent rather than zero — this is a signal the stage never ran.

Add `typing` imports at the top if not already present: `from typing import Any, Dict, Optional`.

- [ ] **Step 4: Run the new test and confirm pass**

```bash
pytest tests/builder/test_build_lod2_diagnostics.py -v
```

Expected: all four tests PASS.

- [ ] **Step 5: Run the entire pre-existing LoD2 test suite to confirm no regression**

```bash
pytest tests/builder/test_lod2_pipeline.py tests/builder/test_lod1_regression.py tests/builder/test_extract_roof_points_fix.py -v
```

Expected: all pre-existing tests still PASS unchanged.

- [ ] **Step 6: Suggested commit**

```bash
git add dtcc_core/builder/geometry_builders/buildings.py \
        tests/builder/test_build_lod2_diagnostics.py
```

Suggested commit message:
```
builder: add opt-in diagnostics_out kwarg to build_lod2_buildings

Populates a per-building dict with per-stage timings, detected planes,
and ClassificationResult when supplied. Default None is no-op; no
behavior change for existing callers.
```

---

## Phase 3 — Metric functions

Each metric is a pure function in `dtcc_core/builder/evaluation/metrics.py`. They consume either the processed `Building`, the `diagnostics_out[bid]` record, or a `GroundTruth`. Each gets its own test.

### Task 3.1: Stage outcome extraction

**Files:**
- Create: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_stage_outcome.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/builder/evaluation/test_metrics_stage_outcome.py
from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface
from dtcc_core.builder.evaluation.metrics import stage_outcome


def test_stage_outcome_success_when_lod2_present_and_no_fallback():
    b = Building(id="b1")
    b.add_geometry(MultiSurface(), GeometryType.LOD2)
    assert stage_outcome(b) == "success"


def test_stage_outcome_reads_fallback_reason():
    b = Building(id="b2")
    b.attributes["fallback_reason"] = "insufficient_points"
    assert stage_outcome(b) == "insufficient_points"


def test_stage_outcome_unknown_when_no_lod2_and_no_reason():
    b = Building(id="b3")
    assert stage_outcome(b) == "unknown"
```

- [ ] **Step 2: Confirm test failure**

```bash
pytest tests/builder/evaluation/test_metrics_stage_outcome.py -v
```

Expected: `ImportError: cannot import name 'stage_outcome'`.

- [ ] **Step 3: Implement the function**

```python
# dtcc_core/builder/evaluation/metrics.py
"""Per-building evaluation metrics. Each function is a pure, standalone measurement."""
from __future__ import annotations

from typing import Optional

from dtcc_core.model.object.building import Building


STAGE_OUTCOMES = (
    "success",
    "insufficient_points",
    "classification_failed",
    "complex_footprint",
    "geometry_construction_failed",
    "validation_failed",
    "unknown",
)


def stage_outcome(building: Building) -> str:
    """Return the terminal pipeline outcome for one building.

    Reads building.attributes['fallback_reason'] (written by build_lod2_buildings).
    Returns 'success' if an LoD2 geometry is present and no fallback_reason set,
    'unknown' if neither is present.
    """
    reason = building.attributes.get("fallback_reason")
    if reason:
        return str(reason)
    if building.lod2 is not None:
        return "success"
    return "unknown"
```

- [ ] **Step 4: Run and confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_stage_outcome.py -v
```

Expected: all three PASS.

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_stage_outcome.py
```

Suggested commit message:
```
builder/evaluation: add stage_outcome metric
```

### Task 3.2: Roof type accuracy metric

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py` (append)
- Test: `tests/builder/evaluation/test_metrics_roof_type.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/builder/evaluation/test_metrics_roof_type.py
from dtcc_core.model.object.building import Building
from dtcc_core.builder.evaluation.dataset import GroundTruth
from dtcc_core.builder.evaluation.metrics import roof_type_prediction, roof_type_correct


def test_roof_type_prediction_reads_from_attributes():
    b = Building(id="b1")
    b.attributes["roof_type"] = "GABLED"
    assert roof_type_prediction(b) == "GABLED"


def test_roof_type_prediction_missing_returns_none():
    b = Building(id="b2")
    assert roof_type_prediction(b) is None


def test_roof_type_correct_matches_ground_truth():
    b = Building(id="b3")
    b.attributes["roof_type"] = "FLAT"
    assert roof_type_correct(b, GroundTruth(roof_type="FLAT")) is True
    assert roof_type_correct(b, GroundTruth(roof_type="GABLED")) is False


def test_roof_type_correct_none_prediction_is_false():
    b = Building(id="b4")
    assert roof_type_correct(b, GroundTruth(roof_type="FLAT")) is False
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_metrics_roof_type.py -v
```

- [ ] **Step 3: Implement**

Append to `dtcc_core/builder/evaluation/metrics.py`:

```python
from dtcc_core.builder.evaluation.dataset import GroundTruth


def roof_type_prediction(building: Building) -> Optional[str]:
    """Return the predicted roof type name written to building.attributes, or None."""
    pred = building.attributes.get("roof_type")
    return str(pred) if pred is not None else None


def roof_type_correct(building: Building, gt: GroundTruth) -> bool:
    """True if predicted roof type matches ground-truth label (case-sensitive name match)."""
    pred = roof_type_prediction(building)
    if pred is None:
        return False
    return pred == gt.roof_type
```

- [ ] **Step 4: Run and confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_roof_type.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_roof_type.py
```

Suggested commit message:
```
builder/evaluation: add roof_type_prediction and roof_type_correct metrics
```

### Task 3.3: Plane count delta metric

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_plane_count.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_plane_count.py
from dtcc_core.builder.evaluation.dataset import GroundTruth
from dtcc_core.builder.evaluation.metrics import (
    detected_plane_count,
    expected_plane_count,
    plane_count_delta,
)


def test_detected_plane_count_from_diagnostics():
    rec = {"planes": [object(), object(), object()]}
    assert detected_plane_count(rec) == 3


def test_detected_plane_count_none_planes():
    assert detected_plane_count({"planes": None}) == 0
    assert detected_plane_count({}) == 0


def test_expected_plane_count_explicit_overrides_type_default():
    gt = GroundTruth(roof_type="GABLED", expected_plane_count=5)
    assert expected_plane_count(gt) == 5


def test_expected_plane_count_defaults_per_type():
    assert expected_plane_count(GroundTruth(roof_type="FLAT")) == 1
    assert expected_plane_count(GroundTruth(roof_type="GABLED")) == 2
    assert expected_plane_count(GroundTruth(roof_type="HIPPED")) == 4


def test_expected_plane_count_unknown_type_returns_none():
    assert expected_plane_count(GroundTruth(roof_type="MANSARD")) is None


def test_plane_count_delta_absolute():
    rec = {"planes": [1, 2, 3, 4]}
    assert plane_count_delta(rec, GroundTruth(roof_type="GABLED")) == 2


def test_plane_count_delta_none_when_expected_unknown():
    rec = {"planes": [1, 2]}
    assert plane_count_delta(rec, GroundTruth(roof_type="MANSARD")) is None
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_metrics_plane_count.py -v
```

- [ ] **Step 3: Implement**

Append to `dtcc_core/builder/evaluation/metrics.py`:

```python
_DEFAULT_PLANE_COUNT = {
    "FLAT": 1,
    "GABLED": 2,
    "HIPPED": 4,
}


def detected_plane_count(diagnostics: dict) -> int:
    planes = diagnostics.get("planes") if diagnostics else None
    return 0 if planes is None else len(planes)


def expected_plane_count(gt: GroundTruth) -> Optional[int]:
    if gt.expected_plane_count is not None:
        return int(gt.expected_plane_count)
    return _DEFAULT_PLANE_COUNT.get(gt.roof_type)


def plane_count_delta(diagnostics: dict, gt: GroundTruth) -> Optional[int]:
    expected = expected_plane_count(gt)
    if expected is None:
        return None
    return abs(detected_plane_count(diagnostics) - expected)
```

- [ ] **Step 4: Run and pass**

```bash
pytest tests/builder/evaluation/test_metrics_plane_count.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_plane_count.py
```

Suggested commit message:
```
builder/evaluation: add plane count metrics (detected, expected, delta)
```

### Task 3.4: Point coverage ratio

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_coverage.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_coverage.py
import numpy as np

from dtcc_core.builder.geometry_builders.roof_config import RoofPlane
from dtcc_core.builder.evaluation.metrics import point_coverage_ratio


def _plane_with_inliers(n):
    # RoofPlane requires boundary_3d; area property is computed from it, but
    # we do not use area here, so any valid triangle works.
    boundary = np.array([[0.0, 0.0, 5.0], [1.0, 0.0, 5.0], [0.0, 1.0, 5.0]])
    return RoofPlane(
        normal=np.array([0.0, 0.0, 1.0]),
        offset=0.0,
        inliers=np.arange(n),
        boundary_3d=boundary,
    )


def test_coverage_ratio_is_sum_inliers_over_filtered():
    planes = [_plane_with_inliers(40), _plane_with_inliers(30)]
    rec = {"planes": planes, "classification": None}
    filtered_n = 100
    assert point_coverage_ratio(rec, filtered_n) == 0.70


def test_coverage_ratio_empty_planes():
    rec = {"planes": [], "classification": None}
    assert point_coverage_ratio(rec, 100) == 0.0


def test_coverage_ratio_zero_filtered_returns_none():
    rec = {"planes": [_plane_with_inliers(10)], "classification": None}
    assert point_coverage_ratio(rec, 0) is None


def test_coverage_ratio_clamps_at_one():
    # inlier sum can exceed filtered if planes overlap — clamp
    planes = [_plane_with_inliers(80), _plane_with_inliers(80)]
    rec = {"planes": planes, "classification": None}
    assert point_coverage_ratio(rec, 100) == 1.0
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_metrics_coverage.py -v
```

- [ ] **Step 3: Implement**

Append:

```python
def point_coverage_ratio(diagnostics: dict, filtered_point_count: int) -> Optional[float]:
    """Fraction of filtered roof points assigned to a detected plane.

    Mirrors the `coverage_factor` term used in the rule-based classifier's confidences.
    Returns None if filtered_point_count is zero.
    """
    if filtered_point_count <= 0:
        return None
    planes = diagnostics.get("planes") or []
    inlier_total = sum(len(p.inliers) for p in planes)
    return min(1.0, inlier_total / filtered_point_count)
```

Note: the runner is responsible for passing `filtered_point_count` separately, because the filtered PointCloud is not preserved in `diagnostics` by default. The runner captures it via `len(filtered_pc.points)` at the time of the run.

To support this, extend Task 2.1's instrumentation to also record `filtered_point_count` in the diagnostics record:

Modify `dtcc_core/builder/geometry_builders/buildings.py` inside `build_lod2_buildings` after the filter stage (the existing successful branch):

```python
if record is not None:
    record["filtered_point_count"] = len(filtered_pc.points)
```

Then update Task 2.1's test file `test_build_lod2_diagnostics.py` to add:

```python
def test_diagnostics_records_filtered_point_count():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    assert diagnostics["b001"]["filtered_point_count"] > 0
```

Run both test files after the modification:

```bash
pytest tests/builder/test_build_lod2_diagnostics.py tests/builder/evaluation/test_metrics_coverage.py -v
```

- [ ] **Step 4: Confirm pass**

Expected: all pass.

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_coverage.py \
        dtcc_core/builder/geometry_builders/buildings.py \
        tests/builder/test_build_lod2_diagnostics.py
```

Suggested commit message:
```
builder/evaluation: add point_coverage_ratio metric; expose filtered_point_count in diagnostics
```

### Task 3.5: Top-2 area share

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_area_share.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_area_share.py
import numpy as np

from dtcc_core.builder.geometry_builders.roof_config import RoofPlane
from dtcc_core.builder.evaluation.metrics import top2_area_share, plane_areas


def _square_plane(side: float) -> RoofPlane:
    """A horizontal plane whose boundary_3d is an axis-aligned square of
    side length `side`, so the computed `area` property equals side**2."""
    boundary = np.array([
        [0.0, 0.0, 5.0],
        [side, 0.0, 5.0],
        [side, side, 5.0],
        [0.0, side, 5.0],
    ])
    return RoofPlane(
        normal=np.array([0.0, 0.0, 1.0]),
        offset=0.0,
        inliers=np.arange(10),
        boundary_3d=boundary,
    )


def _plane(area: float) -> RoofPlane:
    import math
    return _square_plane(math.sqrt(area) if area > 0 else 0.0)


def test_plane_areas_extracts_areas():
    planes = [_plane(3.0), _plane(5.0), _plane(2.0)]
    areas = plane_areas(planes)
    assert len(areas) == 3
    assert areas[0] == pytest.approx(3.0, abs=1e-9)
    assert areas[1] == pytest.approx(5.0, abs=1e-9)
    assert areas[2] == pytest.approx(2.0, abs=1e-9)


def test_top2_area_share_with_three_planes():
    planes = [_plane(3.0), _plane(5.0), _plane(2.0)]  # total=10, top2=8
    rec = {"planes": planes}
    assert top2_area_share(rec) == pytest.approx(0.8, abs=1e-9)


def test_top2_area_share_single_plane_is_one():
    planes = [_plane(4.0)]
    rec = {"planes": planes}
    assert top2_area_share(rec) == pytest.approx(1.0, abs=1e-9)


def test_top2_area_share_no_planes_returns_none():
    rec = {"planes": []}
    assert top2_area_share(rec) is None


def test_top2_area_share_zero_total_area_returns_none():
    # Degenerate planes with zero boundary area
    planes = [_plane(0.0), _plane(0.0)]
    rec = {"planes": planes}
    assert top2_area_share(rec) is None
```

(include `import pytest` at the top).

**Design note:** `RoofPlane.area` is a computed `@property` (it runs a 2D `ConvexHull` over `boundary_3d[:, :2]`). Test helpers must construct a valid `boundary_3d` that yields the target area; they cannot assign to `.area` directly.

- [ ] **Step 2: Confirm test failure**

```bash
pytest tests/builder/evaluation/test_metrics_area_share.py -v
```

- [ ] **Step 3: Implement**

Append:

```python
from typing import List


def plane_areas(planes) -> List[float]:
    """Return the list of plane areas in the order given."""
    return [float(getattr(p, "area", 0.0)) for p in planes]


def top2_area_share(diagnostics: dict) -> Optional[float]:
    """Area of the two largest detected planes / total detected plane area.

    Parallels the hipped-classification area-ratio rule. Returns None when
    no planes or when all areas are zero.
    """
    planes = diagnostics.get("planes") or []
    areas = plane_areas(planes)
    total = sum(areas)
    if not areas or total <= 0.0:
        return None
    top2 = sum(sorted(areas, reverse=True)[:2])
    return top2 / total
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_area_share.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_area_share.py
```

Suggested commit message:
```
builder/evaluation: add top2_area_share metric
```

### Task 3.6: Slope symmetry error

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_symmetry.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_symmetry.py
import math

import numpy as np
import pytest

from dtcc_core.builder.geometry_builders.roof_config import RoofPlane
from dtcc_core.builder.evaluation.metrics import slope_symmetry_error


def _plane_with_tilt(tilt_deg: float, side: float) -> RoofPlane:
    """A square plane tilted by tilt_deg about the Y axis, with an approximate
    boundary_3d that gives an area of side**2 (good enough for area ordering)."""
    rad = math.radians(tilt_deg)
    n = np.array([math.sin(rad), 0.0, math.cos(rad)])
    boundary = np.array([
        [0.0, 0.0, 5.0],
        [side, 0.0, 5.0],
        [side, side, 5.0],
        [0.0, side, 5.0],
    ])
    return RoofPlane(
        normal=n,
        offset=0.0,
        inliers=np.arange(10),
        boundary_3d=boundary,
    )


def test_perfectly_symmetric_returns_zero():
    planes = [_plane_with_tilt(30.0, 3.0), _plane_with_tilt(30.0, 3.0)]
    assert slope_symmetry_error({"planes": planes}) == pytest.approx(0.0, abs=1e-9)


def test_asymmetric_returns_abs_difference_in_degrees():
    planes = [_plane_with_tilt(30.0, 3.0), _plane_with_tilt(20.0, 2.5)]
    err = slope_symmetry_error({"planes": planes})
    assert err == pytest.approx(10.0, abs=0.01)


def test_single_plane_returns_none():
    planes = [_plane_with_tilt(30.0, 3.0)]
    assert slope_symmetry_error({"planes": planes}) is None


def test_no_planes_returns_none():
    assert slope_symmetry_error({"planes": []}) is None
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_metrics_symmetry.py -v
```

- [ ] **Step 3: Implement**

Append:

```python
def slope_symmetry_error(diagnostics: dict) -> Optional[float]:
    """Absolute tilt difference (degrees) between the two largest detected planes.

    Parallels the symmetric-pair test in the rule-based gabled classifier.
    Returns None if fewer than two planes. Uses RoofPlane.slope_deg (angle between
    plane normal and vertical).
    """
    planes = diagnostics.get("planes") or []
    if len(planes) < 2:
        return None
    sorted_by_area = sorted(planes, key=lambda p: p.area, reverse=True)
    top2 = sorted_by_area[:2]
    return abs(top2[0].slope_deg - top2[1].slope_deg)
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_symmetry.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_symmetry.py
```

Suggested commit message:
```
builder/evaluation: add slope_symmetry_error metric
```

### Task 3.7: Watertightness pass-through

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_watertight.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_watertight.py
import numpy as np

from dtcc_core.model.object.building import Building
from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.object.object import GeometryType
from dtcc_core.builder.evaluation.metrics import is_watertight


def test_is_watertight_none_lod2_returns_false():
    b = Building(id="b1")
    assert is_watertight(b) is False


def test_is_watertight_simple_box_returns_true():
    # construct a cube MultiSurface and verify shell_validation accepts it
    ms = MultiSurface()
    verts = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float)
    faces = [
        [0, 1, 2, 3],  # bottom
        [4, 5, 6, 7],  # top
        [0, 1, 5, 4],  # -Y
        [1, 2, 6, 5],  # +X
        [2, 3, 7, 6],  # +Y
        [3, 0, 4, 7],  # -X
    ]
    for face in faces:
        ms.surfaces.append(Surface(vertices=verts[face]))
    b = Building(id="b2")
    b.add_geometry(ms, GeometryType.LOD2)
    assert is_watertight(b) is True
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_metrics_watertight.py -v
```

- [ ] **Step 3: Implement**

Append:

```python
from dtcc_core.builder.geometry.shell_validation import validate_shell


def is_watertight(building: Building, tolerance: float = 0.01) -> bool:
    """True iff building.lod2 exists and passes validate_shell at tolerance."""
    if building.lod2 is None:
        return False
    valid, _ = validate_shell(building.lod2, tolerance)
    return bool(valid)
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_watertight.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_watertight.py
```

Suggested commit message:
```
builder/evaluation: add is_watertight metric
```

### Task 3.8: Ridge/eave height error

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_heights.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_heights.py
from dtcc_core.builder.evaluation.dataset import GroundTruth
from dtcc_core.builder.evaluation.metrics import ridge_height_error, eave_height_error


class _StubResult:
    def __init__(self, ridge_height=None, eave_height=None):
        self.ridge_height = ridge_height
        self.eave_height = eave_height


def test_ridge_height_error_abs_difference():
    rec = {"classification": _StubResult(ridge_height=8.0)}
    gt = GroundTruth(roof_type="GABLED", ridge_height=7.5)
    assert ridge_height_error(rec, gt) == 0.5


def test_ridge_height_error_no_gt_returns_none():
    rec = {"classification": _StubResult(ridge_height=8.0)}
    gt = GroundTruth(roof_type="FLAT")
    assert ridge_height_error(rec, gt) is None


def test_ridge_height_error_no_classification_returns_none():
    gt = GroundTruth(roof_type="GABLED", ridge_height=7.5)
    assert ridge_height_error({}, gt) is None


def test_eave_height_error_behaves_symmetrically():
    rec = {"classification": _StubResult(eave_height=5.1)}
    gt = GroundTruth(roof_type="GABLED", eave_height=5.0)
    assert eave_height_error(rec, gt) == pytest.approx(0.1, abs=1e-9)
```

(include `import pytest` at the top).

- [ ] **Step 2: Confirm failure**

- [ ] **Step 3: Implement**

Append:

```python
def _predicted_height(diagnostics: dict, attr: str) -> Optional[float]:
    clsr = diagnostics.get("classification") if diagnostics else None
    if clsr is None:
        return None
    v = getattr(clsr, attr, None)
    return None if v is None else float(v)


def ridge_height_error(diagnostics: dict, gt: GroundTruth) -> Optional[float]:
    if gt.ridge_height is None:
        return None
    pred = _predicted_height(diagnostics, "ridge_height")
    if pred is None:
        return None
    return abs(pred - gt.ridge_height)


def eave_height_error(diagnostics: dict, gt: GroundTruth) -> Optional[float]:
    if gt.eave_height is None:
        return None
    pred = _predicted_height(diagnostics, "eave_height")
    if pred is None:
        return None
    return abs(pred - gt.eave_height)
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_heights.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_heights.py
```

Suggested commit message:
```
builder/evaluation: add ridge/eave height error metrics
```

### Task 3.9: Per-surface semantic overlap

**Files:**
- Modify: `dtcc_core/builder/evaluation/metrics.py`
- Test: `tests/builder/evaluation/test_metrics_semantic_overlap.py`

Implementation approach: for each semantic class in the predicted and ground-truth MultiSurfaces, project the member surfaces to the XY plane, union them into a shapely `MultiPolygon`, then compute `|intersection| / |union|`. This gives a first-pass footprint-level overlap; works well for ROOF and GROUND (horizontal or near-horizontal) and provides a usable WALL-coverage proxy for Tier 1 walls.

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_metrics_semantic_overlap.py
import numpy as np

from dtcc_core.model.geometry.surface import MultiSurface, Surface
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.evaluation.metrics import semantic_overlap


def _unit_square_at_z(z: float) -> Surface:
    return Surface(vertices=np.array([[0,0,z],[1,0,z],[1,1,z],[0,1,z]], dtype=float))


def test_identical_shapes_overlap_is_one():
    pred = MultiSurface()
    pred.surfaces.append(_unit_square_at_z(0))
    pred.semantics = [SurfaceSemantic.ROOF]

    truth = MultiSurface()
    truth.surfaces.append(_unit_square_at_z(0))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == pytest.approx(1.0, abs=1e-9)


def test_no_overlap_is_zero():
    pred = MultiSurface()
    pred.surfaces.append(Surface(vertices=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]], dtype=float)))
    pred.semantics = [SurfaceSemantic.ROOF]

    truth = MultiSurface()
    truth.surfaces.append(Surface(vertices=np.array([[10,10,0],[11,10,0],[11,11,0],[10,11,0]], dtype=float)))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == 0.0


def test_half_overlap():
    # pred: [0,0]-[2,1]; truth: [1,0]-[3,1]; IoU = 1 / 3
    pred = MultiSurface()
    pred.surfaces.append(Surface(vertices=np.array([[0,0,0],[2,0,0],[2,1,0],[0,1,0]], dtype=float)))
    pred.semantics = [SurfaceSemantic.ROOF]

    truth = MultiSurface()
    truth.surfaces.append(Surface(vertices=np.array([[1,0,0],[3,0,0],[3,1,0],[1,1,0]], dtype=float)))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == pytest.approx(1.0/3.0, abs=1e-9)


def test_missing_semantic_in_pred_is_zero():
    pred = MultiSurface()
    pred.surfaces.append(_unit_square_at_z(0))
    pred.semantics = [SurfaceSemantic.WALL]

    truth = MultiSurface()
    truth.surfaces.append(_unit_square_at_z(0))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == 0.0
```

(include `import pytest`).

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_metrics_semantic_overlap.py -v
```

- [ ] **Step 3: Implement**

Append:

```python
from collections import defaultdict
from typing import Dict

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union

from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.model.geometry.surface import MultiSurface


def _project_xy(surface) -> Optional[ShapelyPolygon]:
    verts = surface.vertices
    if verts is None or len(verts) < 3:
        return None
    ring = [(float(x), float(y)) for x, y, *_ in verts]
    poly = ShapelyPolygon(ring)
    if not poly.is_valid or poly.area <= 0.0:
        poly = poly.buffer(0)
    return poly if poly.is_valid and poly.area > 0.0 else None


def _union_by_semantic(ms: MultiSurface) -> Dict[SurfaceSemantic, "BaseGeometry"]:
    if ms is None or not ms.surfaces:
        return {}
    semantics = ms.semantics or [None] * len(ms.surfaces)
    bag = defaultdict(list)
    for surf, sem in zip(ms.surfaces, semantics):
        poly = _project_xy(surf)
        if poly is None:
            continue
        bag[sem].append(poly)
    return {k: unary_union(v) for k, v in bag.items() if v}


def semantic_overlap(
    predicted: MultiSurface,
    truth: MultiSurface,
) -> Dict[SurfaceSemantic, float]:
    """Per-semantic intersection-over-union of XY-projected surfaces.

    Returns a dict keyed by SurfaceSemantic (drawn from predicted∪truth) with
    values in [0, 1]. Missing semantic in either side → 0.0.
    """
    pred_unions = _union_by_semantic(predicted)
    truth_unions = _union_by_semantic(truth)

    keys = set(pred_unions) | set(truth_unions)
    out: Dict[SurfaceSemantic, float] = {}
    for k in keys:
        if k is None:
            continue
        p = pred_unions.get(k)
        t = truth_unions.get(k)
        if p is None or t is None:
            out[k] = 0.0
            continue
        inter = p.intersection(t).area
        union = p.union(t).area
        out[k] = float(inter / union) if union > 0 else 0.0
    return out
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_metrics_semantic_overlap.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/metrics.py \
        tests/builder/evaluation/test_metrics_semantic_overlap.py
```

Suggested commit message:
```
builder/evaluation: add semantic_overlap (XY-projected IoU per semantic class)
```

---

## Phase 4 — Runner

The runner iterates a dataset, invokes the pipeline with diagnostics enabled, evaluates each metric, and returns a structured results list.

### Task 4.1: Results data model + runner skeleton

**Files:**
- Create: `dtcc_core/builder/evaluation/runner.py`
- Test: `tests/builder/evaluation/test_runner.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_runner.py
from pathlib import Path

from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.runner import Runner, BuildingResult


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_runner_produces_one_result_per_sample():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    assert len(results) == 1
    assert isinstance(results[0], BuildingResult)


def test_runner_result_has_expected_fields_for_flat_roof_fixture():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    r = results[0]
    assert r.building_id == "b001"
    assert r.stage_outcome in ("success", "validation_failed")  # flat fixture may validate differently
    assert r.roof_type_predicted is not None
    assert r.roof_type_truth == "FLAT"
    assert r.timings is not None
    assert "filter" in r.timings
    # detection-stage diagnostics may be present (or None if fallback before detect)
    assert r.plane_count >= 0


def test_runner_records_total_wall_clock():
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))
    r = results[0]
    assert r.total_ms >= 0.0
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_runner.py -v
```

- [ ] **Step 3: Implement**

```python
# dtcc_core/builder/evaluation/runner.py
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings
from dtcc_core.builder.geometry_builders.surface import extrude_surface
from dtcc_core.builder.evaluation.dataset import EvalDataset, EvalSample
from dtcc_core.builder.evaluation.metrics import (
    stage_outcome,
    roof_type_prediction,
    roof_type_correct,
    detected_plane_count,
    expected_plane_count,
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
    classifier_source: Optional[str]        # populated when Track 2 lands; None for now
    fallback_reason_attr: Optional[str]
    plane_count: int
    plane_count_delta: Optional[int]
    coverage_ratio: Optional[float]
    top2_area_share: Optional[float]
    symmetry_error_deg: Optional[float]
    watertight: bool
    ridge_height_error: Optional[float]
    eave_height_error: Optional[float]
    semantic_iou: Dict[str, float]          # keyed by SurfaceSemantic.name for serializability
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
            classifier_source=building.attributes.get("classifier_source"),  # None until Track 2
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
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_runner.py -v
```

Expected: all three tests PASS.

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/runner.py \
        tests/builder/evaluation/test_runner.py
```

Suggested commit message:
```
builder/evaluation: add Runner + BuildingResult; drives pipeline and collects metrics
```

---

## Phase 5 — Reporter

### Task 5.1: CSV writer

**Files:**
- Create: `dtcc_core/builder/evaluation/report.py`
- Test: `tests/builder/evaluation/test_report_csv.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_report_csv.py
import csv
from pathlib import Path

from dtcc_core.builder.evaluation.runner import BuildingResult
from dtcc_core.builder.evaluation.report import write_per_building_csv


def _fake_result(bid="b001", outcome="success"):
    return BuildingResult(
        building_id=bid,
        stage_outcome=outcome,
        roof_type_predicted="FLAT",
        roof_type_truth="FLAT",
        roof_type_correct=True,
        confidence=0.9,
        classifier_source=None,
        fallback_reason_attr=None,
        plane_count=1,
        plane_count_delta=0,
        coverage_ratio=0.95,
        top2_area_share=1.0,
        symmetry_error_deg=None,
        watertight=True,
        ridge_height_error=None,
        eave_height_error=None,
        semantic_iou={"ROOF": 0.98, "GROUND": 0.99},
        timings={"filter": 3.0, "detect": 12.0},
        total_ms=20.5,
    )


def test_csv_one_row_per_result(tmp_path: Path):
    out = tmp_path / "per_building.csv"
    write_per_building_csv([_fake_result(), _fake_result(bid="b002", outcome="validation_failed")], out)
    with out.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 2
    assert rows[0]["building_id"] == "b001"
    assert rows[0]["stage_outcome"] == "success"
    assert rows[1]["stage_outcome"] == "validation_failed"


def test_csv_flattens_semantic_iou_columns(tmp_path: Path):
    out = tmp_path / "per_building.csv"
    write_per_building_csv([_fake_result()], out)
    with out.open() as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
    assert "iou_ROOF" in headers
    assert "iou_GROUND" in headers


def test_csv_flattens_stage_timings(tmp_path: Path):
    out = tmp_path / "per_building.csv"
    write_per_building_csv([_fake_result()], out)
    with out.open() as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
    assert "timing_filter_ms" in headers
    assert "timing_detect_ms" in headers
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_report_csv.py -v
```

- [ ] **Step 3: Implement**

```python
# dtcc_core/builder/evaluation/report.py
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable, List

from dtcc_core.builder.evaluation.runner import BuildingResult


def write_per_building_csv(results: Iterable[BuildingResult], path: Path | str) -> None:
    """Write one CSV row per result with flattened timing + semantic-IoU columns."""
    results = list(results)
    if not results:
        path_obj = Path(path)
        path_obj.parent.mkdir(parents=True, exist_ok=True)
        path_obj.write_text("")
        return

    iou_keys = sorted({k for r in results for k in r.semantic_iou})
    timing_keys = sorted({k for r in results for k in r.timings})

    base_headers = [
        "building_id", "stage_outcome", "roof_type_predicted", "roof_type_truth",
        "roof_type_correct", "confidence", "classifier_source", "fallback_reason_attr",
        "plane_count", "plane_count_delta", "coverage_ratio", "top2_area_share",
        "symmetry_error_deg", "watertight", "ridge_height_error", "eave_height_error",
        "total_ms",
    ]
    headers = base_headers + [f"iou_{k}" for k in iou_keys] + [f"timing_{k}_ms" for k in timing_keys]

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        for r in results:
            row = {h: getattr(r, h, None) for h in base_headers}
            for k in iou_keys:
                row[f"iou_{k}"] = r.semantic_iou.get(k)
            for k in timing_keys:
                row[f"timing_{k}_ms"] = r.timings.get(k)
            w.writerow(row)
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_report_csv.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/report.py \
        tests/builder/evaluation/test_report_csv.py
```

Suggested commit message:
```
builder/evaluation: add per-building CSV writer
```

### Task 5.2: Aggregate stats + summary JSON

**Files:**
- Modify: `dtcc_core/builder/evaluation/report.py`
- Test: `tests/builder/evaluation/test_report_summary.py`

- [ ] **Step 1: Write failing test**

```python
# tests/builder/evaluation/test_report_summary.py
import json
from pathlib import Path

from dtcc_core.builder.evaluation.report import summarize, write_summary_json
from tests.builder.evaluation.test_report_csv import _fake_result


def test_summarize_counts_stage_outcomes():
    results = [
        _fake_result(bid="a", outcome="success"),
        _fake_result(bid="b", outcome="success"),
        _fake_result(bid="c", outcome="insufficient_points"),
    ]
    summary = summarize(results)
    assert summary["count"] == 3
    assert summary["stage_outcome_counts"]["success"] == 2
    assert summary["stage_outcome_counts"]["insufficient_points"] == 1


def test_summarize_reports_roof_type_accuracy_on_success_subset():
    a = _fake_result(bid="a", outcome="success")
    a.roof_type_correct = True
    b = _fake_result(bid="b", outcome="success")
    b.roof_type_correct = False
    c = _fake_result(bid="c", outcome="insufficient_points")
    c.roof_type_correct = False   # should not count in accuracy
    summary = summarize([a, b, c])
    assert summary["roof_type_accuracy_on_success"] == 0.5


def test_summarize_reports_watertight_rate():
    a = _fake_result(bid="a", outcome="success")
    a.watertight = True
    b = _fake_result(bid="b", outcome="success")
    b.watertight = False
    summary = summarize([a, b])
    assert summary["watertight_rate"] == 0.5


def test_summarize_timing_percentiles():
    results = []
    for i, ms in enumerate([1.0, 2.0, 3.0, 4.0, 100.0]):
        r = _fake_result(bid=f"b{i}")
        r.total_ms = ms
        results.append(r)
    summary = summarize(results)
    assert summary["total_ms"]["median"] == 3.0
    assert summary["total_ms"]["max"] == 100.0


def test_summary_json_is_written(tmp_path: Path):
    out = tmp_path / "summary.json"
    write_summary_json([_fake_result()], out)
    with out.open() as f:
        data = json.load(f)
    assert data["count"] == 1
```

- [ ] **Step 2: Confirm failure**

- [ ] **Step 3: Implement**

Append to `dtcc_core/builder/evaluation/report.py`:

```python
from collections import Counter
from statistics import median
from typing import Dict, Any


def summarize(results: List[BuildingResult]) -> Dict[str, Any]:
    n = len(results)
    if n == 0:
        return {"count": 0}

    stage_counts = Counter(r.stage_outcome for r in results)
    success = [r for r in results if r.stage_outcome == "success"]
    rt_acc = (
        sum(r.roof_type_correct for r in success) / len(success)
        if success else None
    )
    wt_rate = sum(1 for r in results if r.watertight) / n

    def _pct(vals: List[float], q: float) -> float:
        s = sorted(vals)
        k = max(0, min(len(s) - 1, int(round(q * (len(s) - 1)))))
        return s[k]

    total_times = [r.total_ms for r in results]
    timing = {
        "min": min(total_times),
        "median": median(total_times),
        "p95": _pct(total_times, 0.95),
        "max": max(total_times),
    }

    pc_deltas = [r.plane_count_delta for r in results if r.plane_count_delta is not None]
    coverage = [r.coverage_ratio for r in results if r.coverage_ratio is not None]
    sym = [r.symmetry_error_deg for r in results if r.symmetry_error_deg is not None]

    return {
        "count": n,
        "stage_outcome_counts": dict(stage_counts),
        "roof_type_accuracy_on_success": rt_acc,
        "watertight_rate": wt_rate,
        "total_ms": timing,
        "plane_count_delta": {
            "median": median(pc_deltas) if pc_deltas else None,
            "max": max(pc_deltas) if pc_deltas else None,
        },
        "coverage_ratio": {
            "median": median(coverage) if coverage else None,
        },
        "slope_symmetry_deg": {
            "median": median(sym) if sym else None,
        },
    }


def write_summary_json(results: List[BuildingResult], path: Path | str) -> None:
    summary = summarize(results)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w") as f:
        json.dump(summary, f, indent=2, default=str)
```

- [ ] **Step 4: Confirm pass**

```bash
pytest tests/builder/evaluation/test_report_summary.py -v
```

- [ ] **Step 5: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/report.py \
        tests/builder/evaluation/test_report_summary.py
```

Suggested commit message:
```
builder/evaluation: add summarize() + write_summary_json()
```

### Task 5.3: Failure gallery (VTK export for threshold-crossing buildings)

**Files:**
- Modify: `dtcc_core/builder/evaluation/report.py`
- Test: `tests/builder/evaluation/test_report_gallery.py`

**Design:** the gallery inspects each result against configurable thresholds; for every building that fails any of them, writes `{building_id}.vtk` under a specified directory, re-invoking the pipeline for that sample so we have a live Building to hand `sandbox/export_vtk.py`'s writer. Since `export_vtk.py` currently lives in `sandbox/`, this task also promotes its writer function into `dtcc_core/builder/evaluation/` as a reusable helper.

- [ ] **Step 1: Write the failing vtk_export test**

```python
# tests/builder/evaluation/test_vtk_export.py
from pathlib import Path

import numpy as np

from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface, Surface
from dtcc_core.builder.evaluation.vtk_export import write_building_vtk


def _square(z: float = 0.0) -> Surface:
    return Surface(vertices=np.array([[0,0,z],[1,0,z],[1,1,z],[0,1,z]], dtype=float))


def test_write_building_vtk_creates_file(tmp_path: Path):
    b = Building(id="b1")
    ms = MultiSurface()
    ms.surfaces.append(_square())
    b.add_geometry(ms, GeometryType.LOD2)
    out = tmp_path / "b1.vtk"
    write_building_vtk(b, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_write_building_vtk_no_lod2_raises(tmp_path: Path):
    import pytest
    b = Building(id="b2")
    out = tmp_path / "b2.vtk"
    with pytest.raises(ValueError):
        write_building_vtk(b, out)
```

- [ ] **Step 2: Confirm failure**

```bash
pytest tests/builder/evaluation/test_vtk_export.py -v
```

Expected: `ModuleNotFoundError: dtcc_core.builder.evaluation.vtk_export`.

- [ ] **Step 3: Create `dtcc_core/builder/evaluation/vtk_export.py` by lifting the existing fan-triangulation logic**

The existing `sandbox/export_vtk.py` has a `multisurface_to_meshio(ms, building_id)` helper that fan-triangulates a MultiSurface and returns `(vertices, triangles, semantics, building_ids)`. Lift that helper into the package and add a single-building wrapper.

Create `dtcc_core/builder/evaluation/vtk_export.py`:

```python
"""Reusable VTK writer for Building/MultiSurface objects. Fan-triangulates
each surface and writes a meshio Mesh with per-triangle cell data."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import meshio

from dtcc_core.model.geometry.surface import MultiSurface
from dtcc_core.model.object.building import Building


def multisurface_to_meshio(ms: MultiSurface, building_id: int = 0):
    """Convert a MultiSurface to meshio-compatible vertices, triangles, and cell data.

    Lifted from sandbox/export_vtk.py so it can be shared with the eval harness.
    """
    vertices: list = []
    triangles: list = []
    semantic_values: list = []
    building_ids: list = []
    v_offset = 0

    for i, surface in enumerate(ms.surfaces):
        verts = surface.vertices
        n = len(verts)
        for v in verts:
            vertices.append(v)

        sem_val = 1  # default WALL
        if ms.semantics is not None and i < len(ms.semantics):
            sem_val = ms.semantics[i].value

        for j in range(1, n - 1):
            triangles.append([v_offset, v_offset + j, v_offset + j + 1])
            semantic_values.append(sem_val)
            building_ids.append(building_id)

        v_offset += n

    return (
        np.array(vertices),
        np.array(triangles),
        np.array(semantic_values),
        np.array(building_ids),
    )


def write_building_vtk(
    building: Building,
    path: Path | str,
    lod: str = "lod2",
) -> Path:
    """Write a single Building's LoD1 or LoD2 MultiSurface to a VTK file.

    Raises ValueError if the requested LoD is not present.
    """
    ms = building.lod2 if lod == "lod2" else building.lod1
    if ms is None:
        raise ValueError(f"Building {building.id} has no {lod} geometry")

    verts, tris, sems, bids = multisurface_to_meshio(ms, building_id=0)
    if len(tris) == 0:
        raise ValueError(f"Building {building.id}: no triangles to export")

    roof_type = building.attributes.get("roof_type", "NONE")
    rt_code = {"FLAT": 0, "GABLED": 1, "HIPPED": 2, "UNKNOWN": 3}.get(roof_type, 4)
    confidence = float(building.attributes.get("roof_confidence", 0.0))

    cell_data = {
        "building_id": [bids],
        "roof_type": [np.full(len(tris), rt_code)],
        "confidence": [np.full(len(tris), confidence)],
    }
    if lod == "lod2":
        cell_data["semantic"] = [sems]

    mesh = meshio.Mesh(
        points=verts,
        cells=[("triangle", tris)],
        cell_data=cell_data,
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(path, mesh)
    return path
```

- [ ] **Step 4: Rewire `sandbox/export_vtk.py` to import the lifted helper**

Edit `sandbox/export_vtk.py` to replace the local `multisurface_to_meshio` definition with an import:

```python
# sandbox/export_vtk.py — change:
from dtcc_core.builder.evaluation.vtk_export import multisurface_to_meshio
```

Delete the local `def multisurface_to_meshio(...)` function body from the file. Leave `export_buildings_vtk(...)` and `main()` intact — they now call the imported helper.

- [ ] **Step 5: Confirm tests pass and sandbox script still imports cleanly**

```bash
pytest tests/builder/evaluation/test_vtk_export.py -v
python -c "import sandbox.export_vtk"    # smoke: the import should not error
```

Expected: both succeed.

- [ ] **Step 6: Write the failing gallery test**

```python
# tests/builder/evaluation/test_report_gallery.py
from pathlib import Path

from dtcc_core.builder.evaluation.runner import Runner
from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.report import FailureThresholds, write_failure_gallery


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_gallery_exports_vtk_for_threshold_crossing(tmp_path: Path):
    # very aggressive thresholds so the fixture will always "fail" and be exported
    thresholds = FailureThresholds(max_plane_count_delta=0, min_coverage=0.999, require_watertight=True)
    runner = Runner()
    results = runner.run(EvalDataset(FIXTURE_ROOT))

    gallery_dir = tmp_path / "gallery"
    write_failure_gallery(
        dataset_root=FIXTURE_ROOT,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )
    assert any(gallery_dir.iterdir())
```

- [ ] **Step 7: Implement**

The gallery needs live `Building` objects. Two options:
(a) Modify `Runner.run` to return `(results, buildings)` tuple.
(b) Re-run the dataset inside `write_failure_gallery` to reconstruct the buildings for failing ids.

Option (b) keeps `Runner.run`'s return shape simple. Use it.

Append to `dtcc_core/builder/evaluation/report.py`:

```python
from dataclasses import dataclass
from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.vtk_export import write_building_vtk


@dataclass
class FailureThresholds:
    max_plane_count_delta: Optional[int] = 2
    min_coverage: Optional[float] = 0.7
    min_semantic_iou: Optional[float] = 0.7
    require_watertight: bool = True


def _is_failure(r: BuildingResult, t: FailureThresholds) -> bool:
    if t.max_plane_count_delta is not None and r.plane_count_delta is not None:
        if r.plane_count_delta > t.max_plane_count_delta:
            return True
    if t.min_coverage is not None and r.coverage_ratio is not None:
        if r.coverage_ratio < t.min_coverage:
            return True
    if t.min_semantic_iou is not None and r.semantic_iou:
        if any(v < t.min_semantic_iou for v in r.semantic_iou.values()):
            return True
    if t.require_watertight and not r.watertight:
        return True
    if r.stage_outcome != "success":
        return True
    return False


def write_failure_gallery(
    dataset_root,
    results: List[BuildingResult],
    thresholds: FailureThresholds,
    out_dir: Path | str,
) -> List[Path]:
    """For each failing building, re-run the pipeline and write a VTK file."""
    from dtcc_core.builder.evaluation.runner import Runner
    from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings

    failing_ids = {r.building_id for r in results if _is_failure(r, thresholds)}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths: List[Path] = []
    ds = EvalDataset(dataset_root)
    runner = Runner()   # reuse sample→building conversion
    for sample in ds:
        if sample.building_id not in failing_ids:
            continue
        building = runner._sample_to_building(sample)  # noqa: protected, intentional reuse
        build_lod2_buildings([building])
        path = out_dir / f"{sample.building_id}.vtk"
        write_building_vtk(building, path)
        paths.append(path)
    return paths
```

- [ ] **Step 8: Confirm pass**

```bash
pytest tests/builder/evaluation/test_report_gallery.py -v
```

- [ ] **Step 9: Suggested commit**

```bash
git add dtcc_core/builder/evaluation/report.py \
        dtcc_core/builder/evaluation/vtk_export.py \
        sandbox/export_vtk.py \
        tests/builder/evaluation/test_report_gallery.py \
        tests/builder/evaluation/test_vtk_export.py
```

Suggested commit message:
```
builder/evaluation: add failure gallery + promote VTK writer to package
```

---

## Phase 6 — CLI driver and end-to-end smoke test

### Task 6.1: sandbox/evaluate_lod2.py

**Files:**
- Create: `sandbox/evaluate_lod2.py`

- [ ] **Step 1: Implement the CLI driver**

```python
#!/usr/bin/env python3
"""
Evaluate the rule-based LoD2 pipeline against a labeled dataset.

Usage:
    python sandbox/evaluate_lod2.py --dataset path/to/dataset --out path/to/output
"""
import argparse
import sys
from pathlib import Path

from dtcc_core.builder.evaluation.dataset import EvalDataset
from dtcc_core.builder.evaluation.runner import Runner
from dtcc_core.builder.evaluation.report import (
    FailureThresholds,
    summarize,
    write_failure_gallery,
    write_per_building_csv,
    write_summary_json,
)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--max-plane-count-delta", type=int, default=2)
    p.add_argument("--min-coverage", type=float, default=0.7)
    p.add_argument("--min-semantic-iou", type=float, default=0.7)
    p.add_argument("--no-watertight-required", action="store_true")
    args = p.parse_args(argv)

    ds = EvalDataset(args.dataset)
    runner = Runner()
    results = runner.run(ds)

    args.out.mkdir(parents=True, exist_ok=True)
    csv_path = args.out / "per_building.csv"
    summary_path = args.out / "summary.json"
    gallery_dir = args.out / "failures"

    write_per_building_csv(results, csv_path)
    write_summary_json(results, summary_path)
    thresholds = FailureThresholds(
        max_plane_count_delta=args.max_plane_count_delta,
        min_coverage=args.min_coverage,
        min_semantic_iou=args.min_semantic_iou,
        require_watertight=not args.no_watertight_required,
    )
    gallery_paths = write_failure_gallery(
        dataset_root=args.dataset,
        results=results,
        thresholds=thresholds,
        out_dir=gallery_dir,
    )

    summary = summarize(results)
    print(f"Evaluated {summary['count']} buildings")
    for reason, count in sorted(summary["stage_outcome_counts"].items()):
        print(f"  {reason}: {count}")
    if summary.get("roof_type_accuracy_on_success") is not None:
        print(f"Roof-type accuracy on success subset: "
              f"{summary['roof_type_accuracy_on_success']:.3f}")
    print(f"Watertight rate: {summary['watertight_rate']:.3f}")
    print(f"Wrote: {csv_path}, {summary_path}")
    print(f"Failure gallery: {len(gallery_paths)} VTK files under {gallery_dir}")


if __name__ == "__main__":
    main(sys.argv[1:])
```

- [ ] **Step 2: Run the CLI against the fixture dataset as a smoke test**

```bash
cd /Users/vasnas/scratch/lod3/temp/dtcc-core
python sandbox/evaluate_lod2.py \
    --dataset tests/builder/evaluation/fixtures/minimal_dataset \
    --out /tmp/eval_smoke_out
```

Expected: console output reports 1 building, CSV and summary.json files written, gallery directory created. Exit code 0.

- [ ] **Step 3: Suggested commit**

```bash
git add sandbox/evaluate_lod2.py
```

Suggested commit message:
```
sandbox: add evaluate_lod2.py CLI driver for the eval harness
```

### Task 6.2: End-to-end smoke test in pytest

**Files:**
- Test: `tests/builder/evaluation/test_end_to_end.py`

- [ ] **Step 1: Write the e2e test**

```python
# tests/builder/evaluation/test_end_to_end.py
from pathlib import Path
from sandbox import evaluate_lod2


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"


def test_cli_runs_end_to_end(tmp_path: Path):
    # Calling main() directly — avoid shelling out for test speed + cross-platform.
    evaluate_lod2.main([
        "--dataset", str(FIXTURE_ROOT),
        "--out", str(tmp_path),
    ])
    assert (tmp_path / "per_building.csv").exists()
    assert (tmp_path / "summary.json").exists()
```

**Note:** importing `sandbox.evaluate_lod2` requires `sandbox/__init__.py`. If one does not already exist in the repo, add one (empty file) in a separate commit alongside this test.

- [ ] **Step 2: Confirm pass**

```bash
pytest tests/builder/evaluation/test_end_to_end.py -v
```

- [ ] **Step 3: Suggested commit**

```bash
git add tests/builder/evaluation/test_end_to_end.py \
        sandbox/__init__.py    # only if added
```

Suggested commit message:
```
tests: add end-to-end smoke test for the eval harness CLI
```

### Task 6.3: Final full-suite verification

**Files:** none.

- [ ] **Step 1: Run the full test suite**

```bash
cd /Users/vasnas/scratch/lod3/temp/dtcc-core
pytest tests/ -v
```

Expected: all tests PASS. If any fail, stop and fix before declaring the plan complete.

- [ ] **Step 2: Run the eval harness against one real dataset (user provides)**

If the user has a Swedish / Lantmäteriet sample dataset or another labeled set, run:

```bash
python sandbox/evaluate_lod2.py --dataset <path-to-real-dataset> --out sandbox/output/eval_baseline
```

Expected: script completes; `summary.json` contains non-trivial counts; `per_building.csv` has meaningful rows; at least some failures in the gallery (if any buildings cross thresholds).

This produces the baseline metrics for the rule-based LoD2 pipeline. These numbers become the benchmark for Track 2.

---

## Acceptance (against spec)

Check the spec's "Acceptance for Track 1" section:

- [x] Runs on at least one real dataset end-to-end — Task 6.3 step 2
- [x] Produces a CSV with all metrics for every building — Tasks 5.1, 6.1, 6.3
- [x] Produces aggregate stats written to console and a summary JSON — Task 5.2, 6.1
- [x] Failure gallery exported for threshold-crossing buildings — Task 5.3, 6.1
- [x] Runs in under 10 minutes for a dataset of 1,000 buildings — verify empirically at Task 6.3; if violated, profile with the captured per-stage timings before declaring failure

Against the spec's Metrics section:

- [x] Stage outcome (Task 3.1)
- [x] Per-stage timing (Task 2.1)
- [x] Roof type accuracy (Task 3.2)
- [x] Detection diagnostics: plane count delta (3.3), coverage ratio (3.4), top-2 area share (3.5), slope symmetry (3.6)
- [x] Per-surface semantic overlap (Task 3.9)
- [x] Watertightness (Task 3.7)
- [x] Ridge/eave height error (Task 3.8)
- [x] Classifier provenance — captured as `classifier_source` and `fallback_reason_attr` columns in the result; these will be `None` / rule-fallback until Track 2 lands. Harness is ready to display them from day one.

Phase 0 (baseline commit) is the first rollout-sequencing step per the spec; Track 1 is the second. Both are covered.
