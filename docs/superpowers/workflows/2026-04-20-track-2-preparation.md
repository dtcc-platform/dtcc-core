# Track 2 Preparation Workflow

**Date:** 2026-04-20
**Purpose:** What you need to do between "Track 1 eval harness is committed" and "write the Track 2 classifier-parity plan."
**Reference spec:** `docs/superpowers/specs/2026-04-18-ml-for-lod2-plus-design.md`

---

## What Track 2 needs from Track 1

Track 2 Phase 1 is a learned roof-type classifier that must **match rule-based accuracy on Tier 1** on a held-out subset of your real data. To write and validate that plan, we need, in order:

1. **A labeled dataset in the harness format** — a directory of buildings with footprint + point cloud + `ground_truth.json`.
2. **Baseline metrics from the rule-based pipeline on that dataset** — `summary.json` + `per_building.csv` produced by `sandbox/evaluate_lod2.py`.
3. **Three decisions** (see "Decision points" below) that shape the Track 2 plan.

Without (1) and (2), Track 2 Phase 1 acceptance criteria are abstract. Without (3), the plan will have placeholders. All three are your side of the handoff.

---

## Prerequisites checklist

Before starting the workflow steps, confirm you can tick these:

- [ ] A python environment where `sandbox/evaluate_lod2.py` runs against the fixture dataset (`tests/builder/evaluation/fixtures/minimal_dataset`) without error.
- [ ] Access to a real dataset candidate — at minimum: building footprints as shapefile/GeoPackage/GeoJSON, and a LiDAR point cloud (LAS/LAZ) covering the same area. Lantmäteriet data is the obvious first choice for Sweden.
- [ ] A list of building IDs you care about, or a bounding box, or a municipality — whatever scope feels right for "enough to be representative, small enough to iterate on." Target: **100–2,000 buildings** for the first pass. Fewer than 100 and statistics are noisy; more than 2,000 and iteration is slow.
- [ ] Either (a) existing roof-type labels for those buildings, (b) a plan to hand-label a sample, or (c) acceptance that the first run will measure geometry quality and stage outcomes without roof-type accuracy.

---

## Step-by-step workflow

### Step 1 — Prepare the dataset directory

The harness loader expects this layout under `<dataset-root>/`:

```
<dataset-root>/
  <building_id>/
    footprint.shp       (or .gpkg / .geojson)
    points.las          (or .laz / .npy)
    ground_truth.json
```

Each `ground_truth.json` needs at minimum:

```json
{
  "roof_type": "FLAT",
  "ground_height": 0.0
}
```

Allowed `roof_type` values at Track 1 scope: `FLAT`, `GABLED`, `HIPPED`, `UNKNOWN`.

Optional fields that unlock more metrics when present:

```json
{
  "roof_type": "GABLED",
  "ground_height": 42.1,
  "eave_height": 47.1,
  "ridge_height": 49.8,
  "expected_plane_count": 2
}
```

### Step 2 — One-shot slicing script

If your data is "one shapefile + one point-cloud tile," write a throwaway script that produces the directory layout above. Suggested outline (adapt to your CRS, attribute names, and IO conventions):

```python
# scripts/prepare_eval_dataset.py  (one-shot, not committed)
from pathlib import Path
import json
import geopandas as gpd
import laspy
import numpy as np
from shapely.strtree import STRtree
from shapely.geometry import Point

ROOT = Path("eval_datasets/sthlm_01")
ROOT.mkdir(parents=True, exist_ok=True)

fps = gpd.read_file("data/footprints.shp").to_crs("EPSG:3006")
las = laspy.read("data/area.laz")
pts_xyz = np.column_stack([las.x, las.y, las.z])
pts_cls = np.array(las.classification, dtype=int)

# Build spatial index over point-cloud XY
tree = STRtree([Point(x, y) for x, y, _ in pts_xyz])

for i, row in fps.iterrows():
    bid = f"b{i:06d}"
    bdir = ROOT / bid
    bdir.mkdir(exist_ok=True)

    # footprint: single feature in its own file
    gpd.GeoDataFrame([row], crs=fps.crs).to_file(bdir / "footprint.shp")

    # clip point cloud to this footprint
    poly = row.geometry
    buf = poly.buffer(0.5)           # small buffer catches roof overhangs
    idx = [i for i, pt in enumerate(tree.geometries)
           if buf.contains_properly(pt)]
    pts_here = pts_xyz[idx]
    cls_here = pts_cls[idx]
    if len(pts_here) < 5:
        # skip buildings with no LiDAR coverage
        continue

    np.save(bdir / "points.npy", pts_here)   # npy is fast + tiny; or write a .laz

    # ground truth
    gt = {
        "roof_type": str(row.get("roof_type_label", "UNKNOWN")).upper(),
        "ground_height": float(row.get("ground_height", 0.0)),
    }
    if row.get("eave_height") is not None:
        gt["eave_height"] = float(row["eave_height"])
    (bdir / "ground_truth.json").write_text(json.dumps(gt, indent=2))
```

Run it once. Verify the output with a quick sanity check:

```bash
ls eval_datasets/sthlm_01 | wc -l           # expect 100–2000
ls eval_datasets/sthlm_01/b000000            # expect 3 files
cat eval_datasets/sthlm_01/b000000/ground_truth.json
```

### Step 3 — Run the eval harness on a small subset first

Before running on the full dataset, smoke-test with a handful of buildings. Copy the first 5–10 into a sub-folder and run against that:

```bash
mkdir -p eval_datasets/sthlm_01_smoke
cp -r eval_datasets/sthlm_01/b00000{0,1,2,3,4} eval_datasets/sthlm_01_smoke/

.venv/bin/python sandbox/evaluate_lod2.py \
    --dataset eval_datasets/sthlm_01_smoke \
    --out sandbox/output/eval_smoke_run
```

Expected: exit 0, `per_building.csv` has 5 rows, `summary.json` exists, `failures/` may or may not have VTKs. If anything blows up, fix the dataset or the slicing script before running the full set.

### Step 4 — Run against the full dataset

```bash
.venv/bin/python sandbox/evaluate_lod2.py \
    --dataset eval_datasets/sthlm_01 \
    --out sandbox/output/eval_baseline_2026_04_20
```

Expect this to take under 10 minutes for ~1000 buildings (spec acceptance criterion). If it's slower, per-stage timings in the CSV will tell you which stage to blame.

### Step 5 — Inspect results

Three questions to answer from the outputs. Note the answers somewhere (the Track 2 plan will reference them):

**Q1 — Where does the pipeline bleed?** Open `summary.json`, find `stage_outcome_counts`. Example:

```json
"stage_outcome_counts": {
  "success": 812,
  "insufficient_points": 54,
  "classification_failed": 18,
  "complex_footprint": 91,
  "geometry_construction_failed": 7,
  "validation_failed": 22
}
```

- `success` dominant, `classification_failed` sizable → Track 2 classifier has room to help.
- `insufficient_points` or `complex_footprint` dominant → upstream data issue; the ML classifier won't help until you improve the input side.
- `geometry_construction_failed` dominant → geometry synthesis bug, not classifier scope.

**Q2 — How accurate is the rule-based classifier on clean cases?** Same file, `roof_type_accuracy_on_success`. Example: `0.847`. That becomes Phase 1's target for the learned classifier to **match or exceed**. If this number is `null`, you didn't supply roof-type labels — see "If you don't have labels" below.

**Q3 — How long does it take?** `total_ms.p95`. Phase 2 benchmarking will care.

### Step 6 — Open failures in ParaView

```bash
# any VTK in sandbox/output/eval_baseline_2026_04_20/failures/
open sandbox/output/eval_baseline_2026_04_20/failures/b000123.vtk
```

Scan 10–20 of them. Patterns you'll spot:

- Misclassified gabled (reported as hipped or vice-versa) — learnable.
- Roof with a dormer/chimney confusing RANSAC — Tier 3 problem, out of scope for Phase 1.
- Complex footprint (L-shape) dropping to fallback — expected, not a failure mode to fix here.
- Too-sparse LiDAR (missing roof entirely) — upstream data problem.

Note the rough fraction of each category — that informs whether Phase 2 should push hard on RoofN3D training or whether something else is the bigger lever.

---

## Decision points before the Track 2 plan can be written

Answer these three before I write the plan. Defaults in **bold**:

### Decision 1 — Labels strategy

How are roof-type labels obtained for the training set in Phase 1?

- **(a) Pseudo-labels from the rule-based classifier** (spec default) — cheap, easy; caps ML at rule-based performance. Good for proving plumbing.
- (b) Hand-labeled subset — a few hundred buildings manually labeled. Buys the option of ML actually exceeding rules on Tier 1, but only useful with (a) as a fallback for the rest.
- (c) A public dataset like RoofN3D pre-trained, then fine-tuned on your data — skips straight to Phase 2 effectively.

### Decision 2 — Training/test split

What data goes into training vs. held-out evaluation?

- **(a) Split the single dataset 80/20 randomly** — simple, fine for a parity demo.
- (b) Train on the synthetic/public set, evaluate only on your real set — cleaner separation, more work.
- (c) Two separate real datasets (one for training, one for eval) — ideal but usually impractical.

### Decision 3 — Where does the trained model file live?

- **(a) Not checked into git; users train it with `sandbox/train_classifier.py`** (spec default) — keeps repo lean.
- (b) Checked in under `dtcc_core/builder/ml/models/` as part of the package — only feasible if the model stays small (<5 MB) and infrequently retrained.

---

## If you don't have roof-type labels at all

Two usable paths:

**Path A — Geometry-only baseline.** Run the harness with `"roof_type": "UNKNOWN"` in every `ground_truth.json`. You get: stage outcomes, plane-count stats, coverage ratios, top-2 area share, watertight rate, per-stage timings. You don't get: roof-type accuracy, confusion matrix. This still tells us where the pipeline bleeds. Good enough to start writing Track 2.

**Path B — Generate pseudo-labels.** Run the current rule-based pipeline once, accept its outputs as ground truth, and use those for Track 2 Phase 1 parity measurement. Phase 1 becomes "the learned classifier agrees with rules 95%+ of the time" — honest but self-referential. This is what the spec describes as Phase 1 training; you can also use it as the eval ground truth, understanding the accuracy number is meaningless in absolute terms (always 1.0 at parity) and only useful as a plumbing check.

My recommendation: **Path A first.** It gives you a real "where does the pipeline bleed" signal immediately. Follow up with hand-labeling a 100-building sample when you have an hour.

---

## What to hand to me for the Track 2 plan

When you're ready for me to write the plan, paste these into our chat:

1. `summary.json` from the full baseline run.
2. Dataset size: number of buildings, number per roof type (if labeled).
3. Answers to Decisions 1, 2, 3.
4. Any patterns you noticed in the failure gallery (optional but helps).

With those in hand, I'll write `docs/superpowers/plans/<date>-learned-classifier-parity.md` scoped to your actual numbers.

---

## Parallel housekeeping worth doing

Not gated on Track 2; can happen any time:

- **Resolve the protobuf version pin.** Current `.venv` has protobuf 6.33.6; `pyproject.toml` still pins `< 6.0.0`. Either bump the pin to `< 7.0.0` or regenerate `dtcc_core/model/dtcc_pb2.py` with a protoc 5.x compiler.
- **Decide GroundTruth.geometry serialization format.** If you do have per-surface LoD2 ground truth anywhere, tell me the format and I'll wire it into the loader so `semantic_iou` columns start populating.
- **Push the branch, open a PR, or merge** when Track 1 is at a point you're happy to publish.
