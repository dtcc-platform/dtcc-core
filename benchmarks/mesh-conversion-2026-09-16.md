# Mesh conversion checks and timing — 16 September 2026

## Scope

Fix `Mesh.to_cpp()` and `VolumeMesh.to_cpp()` to pass the model object to their
converters. Validate numerical inputs at the native constructors before casting
or indexing: real finite coordinates, triangular/tetrahedral integer connectivity
within the vertex range, and optional integer markers aligned with faces/cells
and representable by native `int`. Nonempty float indices/markers now fail instead
of silently truncating. Empty models remain supported.

Volume export uses owning NumPy arrays, following the existing Mesh export.
There is no borrowed storage or change to mesh algorithms. Transform, CRS,
normals and other metadata preservation are outside this change; the round-trip
checks cover vertices, connectivity and markers.

## Benchmark

Run `uv run python benchmarks/benchmark_mesh_conversion.py`. It uses a fixed RNG
seed, one warm-up and seven timed repetitions per case. Times are medians in
seconds in the JSON output. Model creation, meshing and result destruction are
outside the timed region; the converters' own model construction is included.

Measured on macOS ARM64 / Python 3.11.14, in an isolated uv environment. Baseline
native extension and Python conversion code were from `3da764e`. The same script
was run after a fresh native rebuild. Each large case has 100,000 vertices and
100,000 faces or tetrahedra, with float64 coordinates, int64 connectivity and
int32 markers. The script also measures 1,000-element cases.

| Model | Direction | Before (ms) | After (ms) |
| --- | --- | ---: | ---: |
| Mesh | `to_cpp` | 1.203 | 1.160 |
| Mesh | `from_cpp` | 1.584 | 0.179 |
| VolumeMesh | `to_cpp` | 1.261 | 1.016 |
| VolumeMesh | `from_cpp` | 171.286 | 0.175 |

The volume export improvement removes per-element Python object construction.
Both exports fill their new arrays directly, removing per-value bounds checks
from Mesh export. Input conversion reserves native vector capacity and checks
values while copying. These are synthetic conversion results, not end-to-end
meshing speedups; small timing differences can vary between runs.
The benchmark is a manual diagnostic, not a timing threshold in CI.

## Verification

- Fresh isolated native build passed.
- 29 focused conversion tests passed: numerical roundtrips, owned/copy lifetime,
  empty meshes, strided/non-native-endian and unaligned arrays, invalid shape,
  nonfinite coordinates, invalid index types/ranges, and marker alignment/range.
- Full local suite against the final native build: **1,979 passed, 6 skipped,
  66 opt-in live tests deselected**. No live dataset or cross-platform checks ran.
- Independent code review found no blocking issue; whitespace checks passed.

## Follow-up: coordinate frames and face normals

Native Mesh conversion now copies optional face normals in both directions.
Direct Mesh/VolumeMesh conversion rejects nonidentity transforms, SRS, fields,
semantic regions, Dataset Context and schema declarations that native numerical
objects cannot carry. Python merge/snap operations retain a copy of the frame;
merging requires identical affine transforms and SRS. Distances remain in local
units. Merge preserves supplied face normals and computes missing ones when
needed; snapping recomputes supplied normals and rejects degenerate faces.

The same 100,000-element benchmark after this follow-up measured Mesh input/output
at 0.786/0.138 ms and VolumeMesh input/output at 0.806/0.158 ms. These cases have
no normals; they check that metadata admission does not materially regress the
existing numerical path. Normal transfer adds a linear copy when normals exist.

The rebuilt extension passed 46 focused conversion/metadata tests and the full
local suite: **1,996 passed, 6 skipped, 66 opt-in live tests deselected**.
