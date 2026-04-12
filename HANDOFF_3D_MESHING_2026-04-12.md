## Current state

Branch: `dev/cleaning`

This commit captures the upstream footprint-conditioning work done to improve the 3D meshing pipeline without adding downstream PLC/TetGen hacks.

The main source changes are in:

- `dtcc_core/builder/cleaning/footprints.py`
- `dtcc_core/builder/geometry_builders/meshes.py`
- `sandbox/compare_stockholm_flat_mesh.py`
- `sandbox/mesh_quality_survey_3d.py`
- `tests/builder/test_cleaning_footprints.py`
- `tests/builder/test_cleaning_meshing_integration.py`

## What changed

### Stage contracts and diagnostics

- Added stage-audit support to the 3D builder path and to `sandbox/mesh_quality_survey_3d.py`.
- Added explicit stage contracts for:
  - conditioned footprints
  - 2D triangle meshes
  - PLC / TetGen precheck
- Added a mesher-ready coverage revalidation step in `meshes.py` so the stage-1 output is rechecked before stage 2 if polygon normalization reintroduces defects.
- Added a guard that rejects revalidation candidates with insufficient clearance gain.

### Sandbox harness

- Fixed `sandbox/compare_stockholm_flat_mesh.py` so it can inspect stage 2 meshing from already-conditioned footprints without accidentally sending those footprints back through stage 1 again.

### Cleaner architecture changes

- Unified point-touch and close-pair contact repair around one connector-family operator.
- Moved contact merging upstream into normal conditioning rather than relying on late `coverage_meshing_regularization`.
- Fixed the source-coordinate-recovery contract so it does not degrade exact point touches into tiny close gaps.
- Added cluster-aware handling for point-touch clusters.
- Added upstream ring-contact repair and short-edge angle-opening operators.
- Added post-contact local repair acceptance logic so sub-grid but real clearance improvements are kept when they are not worse on the coverage contract.

## Key findings

### Case 62

`case 62` is no longer the main blocker.

Current stage-1 output:

- conditioned minimum clearance: `0.36253060910424234 m`
- no pair issues
- no ring contacts
- no short edges

Current 3D output is materially improved and usable:

- `num_cells = 294445`
- `ARmax = 292.5824`
- `edge p01 = 2.3689`
- `vol p01 = 2.4876`
- `Q<0.05 = 114`

The remaining residual in `case 62` is a single-polygon `exterior <-> hole` clearance issue, but it is not currently the highest-priority problem.

### Case 23

`case 23` remains the real stage-1 blocker.

Current stage-1 output:

- conditioned minimum clearance improved from `0.014600733992761013 m` to `0.03125 m`
- no pair issues
- no ring contacts
- no short edges

However the remaining worst polygon is now:

- conditioned polygon id: `29`
- source ids:
  `132,133,134,135,137,138,139,172,173,176,177,178,180,181,182,183,185,187,188,189,190,192,193`
- clearance class: `hole_3 <-> hole_2`
- minimum clearance: `0.03125 m`

This is the important conclusion: the remaining `case 23` problem is no longer a pair-contact or ring-contact problem. It is a pure single-polygon `hole-hole` self-clearance defect inside one merged footprint.

Current 3D result for `case 23` is still bad:

- `num_cells = 287897`
- `ARmax = 1004.3823`
- `edge p01 = 0.1635`
- `vol p01 = 0.000633`
- `Q<0.05 = 467`

So the next principled fix should be in stage 1:

- detect subscale `hole-hole` passages / walls inside a single conditioned polygon
- remove or merge them upstream
- do not try to solve this in PLC preparation or TetGen

## Recommended next step

Focus on the remaining `case 23` polygon-29 self-clearance defect.

The likely correct policy is:

- if two interior holes are separated only by a subscale wall, merge the holes / remove the wall
- keep the resulting polygon contract-monotone:
  - no pair issues
  - no ring contacts
  - no short edges
  - improved minimum clearance
  - bounded fidelity drift

In other words, treat this the same way the work treated point-touch and courtyard-contact classes: remove subscale topology rather than carrying it downstream.

## Useful commands

Run the flat comparison harness on the representative cases:

```bash
PYTHONPATH=/Users/logg/scratch/dtcc/dtcc-core \
/Users/logg/scratch/dtcc/venv/bin/python \
/Users/logg/scratch/dtcc/dtcc-core/sandbox/compare_stockholm_flat_mesh.py \
  --cases 23 55 62 83 \
  --mode auto \
  --label laptop_resume
```

Run a single 3D survey case:

```bash
PYTHONPATH=/Users/logg/scratch/dtcc/dtcc-core \
/Users/logg/scratch/dtcc/venv/bin/python \
/Users/logg/scratch/dtcc/dtcc-core/sandbox/mesh_quality_survey_3d.py \
  23 \
  --no-plots \
  --output-dir /Users/logg/scratch/dtcc/dtcc-core/sandbox/output_3d_case23_resume
```

Run the hotspot audit bucket:

```bash
PYTHONPATH=/Users/logg/scratch/dtcc/dtcc-core \
/Users/logg/scratch/dtcc/venv/bin/python \
/Users/logg/scratch/dtcc/dtcc-core/sandbox/mesh_quality_survey_3d.py \
  16 23 55 62 83 84 \
  --stage-audit \
  --no-plots
```

## Verification slice used before handoff

```bash
/Users/logg/scratch/dtcc/venv/bin/python -m py_compile \
  /Users/logg/scratch/dtcc/dtcc-core/dtcc_core/builder/cleaning/footprints.py \
  /Users/logg/scratch/dtcc/dtcc-core/dtcc_core/builder/geometry_builders/meshes.py \
  /Users/logg/scratch/dtcc/dtcc-core/tests/builder/test_cleaning_footprints.py \
  /Users/logg/scratch/dtcc/dtcc-core/tests/builder/test_cleaning_meshing_integration.py

/Users/logg/scratch/dtcc/venv/bin/python -m pytest -q \
  /Users/logg/scratch/dtcc/dtcc-core/tests/builder/test_cleaning_footprints.py \
  -k "post_contact_local_repair_accepts_subgrid_clearance_gain or regularize_coverage_contacts_records_post_contact_local_repair_metrics or regularize_coverage_contacts_accepts_local_point_bridge_on_case23_pair or regularize_coverage_contacts_uses_point_cluster_bridge_for_case62_pair or self_clearance_connector_resolves_case23_style_vertex_edge_junction or apply_local_polygon_repairs_resolves_case62_style_self_clearance_junction"

/Users/logg/scratch/dtcc/venv/bin/python -m pytest -q \
  /Users/logg/scratch/dtcc/dtcc-core/tests/builder/test_cleaning_meshing_integration.py \
  -k "normalize_mesher_ready_coverage or conditioned_footprint_contract_audit or triangle_mesh_contract_from_audit or tetgen_plc_contract_from_audit or stage_audit_records_stage_contracts or uses_shared_surface_pipeline"
```
