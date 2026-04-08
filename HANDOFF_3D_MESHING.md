## 3D Meshing Hand-Off

### Repos

- `dtcc-core`
- `dtcc-tetgen-wrapper`

### Environment

Use the same editable install pattern in `~/scratch/dtcc/venv`:

```sh
/Users/logg/scratch/dtcc/venv/bin/pip install -e /Users/logg/scratch/dtcc/dtcc-tetgen-wrapper
CMAKE_ARGS='-DDTCC_USE_TRIANGLE=ON -DDTCC_TRIANGLE_DIR=/Users/logg/scratch/dtcc/dtcc-core/dtcc_core/cpp/external/triangle' \
  /Users/logg/scratch/dtcc/venv/bin/pip install -e /Users/logg/scratch/dtcc/dtcc-core
```

### What Changed

- Flat, surface, and volume meshing now share the same Python-side cleaning / coverage preparation path.
- The 3D survey is split into:
  - `sandbox/mesh_quality_survey_2d.py`
  - `sandbox/mesh_quality_survey_3d.py`
  - `sandbox/mesh_quality_survey.py` as a 2D compatibility wrapper
- TetGen PLC prechecks were added in `dtcc_core/builder/meshing/tetgen_utils.py`.
- The box top for 3D closure is now meshed independently instead of copying the full ground triangulation.
- Cleaning now treats hole/exterior ring contacts as a real defect and enforces `min_feature_size` globally at the end of conditioning.
- `compare_stockholm_flat_mesh.py` has a `--show` flag for live Matplotlib inspection.
- `dtcc-tetgen-wrapper` was cleaned up so the API/docs agree, effective switches are reported correctly, and failure dumps include a native `.poly` repro.

### Current 10x10 3D Survey Status

Settings:

- `mesher=dtcc_mesher`
- `max_mesh_size=10`
- `quality=q1.6/25`
- `domain_height=100`
- `max_added_points=50000`
- merged buildings

Reconciled result:

- `74/100` success
- `18/100` failed
- `8/100` timeout

Failure buckets:

- `12` `dtcc_mesher` coverage region seed conflicts
- `4` TetGen self-intersection failures
- `2` TetGen code-2 failures
- `8` timeouts

Worst successful quality outliers:

- case `54`
- case `46`
- case `28`
- case `20`

### Most Important Next Step

Fix the `dtcc_mesher` coverage seed-conflict bucket first. That is now the largest deterministic blocker and it happens before TetGen. After that:

1. work through the timeout / TetGen-failure cases
2. improve the successful-but-bad quality tail, starting with case `54`

### Useful Commands

Single 3D case:

```sh
cd ~/scratch/dtcc/dtcc-core
~/scratch/dtcc/venv/bin/python sandbox/mesh_quality_survey_3d.py 55 --mesher dtcc_mesher
```

Interactive footprint comparison:

```sh
cd ~/scratch/dtcc/dtcc-core
MPLBACKEND=MacOSX ~/scratch/dtcc/venv/bin/python sandbox/compare_stockholm_flat_mesh.py --cases 34 --mode new --label live_case34 --delay 0 --show
```

Core regression slices:

```sh
cd ~/scratch/dtcc/dtcc-core
~/scratch/dtcc/venv/bin/pytest tests/builder/test_meshing.py -q
~/scratch/dtcc/venv/bin/pytest tests/builder/test_cleaning_meshing_integration.py -q
~/scratch/dtcc/venv/bin/pytest tests/builder/test_cleaning_footprints.py -k "regularizes_hole_touching_exterior or residual_scale_polygon_simplifier or should_attempt_meshing_local_candidate or removes_meshing_hostile_short_edges" -q
```

Wrapper regression:

```sh
cd ~/scratch/dtcc/dtcc-tetgen-wrapper
~/scratch/dtcc/venv/bin/pytest tests -q
```
