# DTCC Model performance on the complete 3DBAG tile

Status: completed, 10 September 2026
Authority: DESIGN.md, docs/design/geometry-representations.md and the user's
request to continue with performance profiling and optimization.

## Acceptance boundary

Measure the existing public real-data workflow and optimize its largest measured
avoidable cost. Keep the single native model, v3 wire bytes, v1/v2 readers,
independent mutable arrays and all boundary validation. No production dependencies,
new wire version, cache registry or polygon-storage redesign in this slice.
Preserve all previous uncommitted milestones and unrelated work.

## Checkpoints

- [x] Record fresh-process CPU and memory measurements per workflow stage, plus
      call profiles and numerical-storage accounting on the unchanged 3DBAG tile.
- [x] Implement the smallest measured optimization and check byte-exact exchange,
      independence of mutable data and meaningful malformed-input failures.
- [x] Repeat the measurements and complete public file/package/CityJSON workflow;
      run affected checks and document evidence and remaining costs.

## Verification

Use the existing checksum-pinned source and v3 artifacts under /private/tmp.
Keep profiled timings separate from ordinary timings. Record process baseline,
preparation and peak RSS so interpreter/native libraries and retained allocations
are not confused with numerical payload size. Use focused canonical/representation
checks first, then relevant model/I/O/dataset checks. Do not remove validation to
improve timing. Run the existing full real-data acceptance example after changes.

## Findings and completion evidence

Fresh-process, three-run stage medians are recorded in
`docs/design/model-performance.md`; raw samples are in
`/private/tmp/dtcc-model-perf/{before,after}/`. The numerical buffers total
19,127,848 bytes and are unchanged. Process imports use roughly 200 MiB.

Call profiles identified 1,201,024 scalar transform checks per validation and
147,906 Transform constructions during decode, for only 75,064 actual transforms.
The optimized codec relies on validated float dtype for exact wire-double
representation, retains scalar checks for integers/wider floats, initializes
native metadata directly and fills Protobuf object/geometry destinations in place.
No constructor bypass, cached validity, schema change or shared mutable storage.

Encode: 1.908 -> 1.195 s; decode: 3.045 -> 2.200 s. Encode peak RSS:
590.2 -> 517.6 MiB. Import/package/CityJSON timings also improved; native storage
and decode memory are essentially unchanged. Instrumented validation cost fell
from 3.01 to 1.36 seconds and decode constructs exactly 75,064 Transforms.

The complete pre-optimization canonical tile re-encodes byte-for-byte unchanged.
42 focused canonical/representation checks passed, including two new tests for
nonfinite float64 transforms and independent mutable decoded polygon data.
The full public representation example passed: exact file/package round trips,
Context preservation and all CityJSON semantic/topological/coordinate checks.
Its CityJSON output is byte-identical to the prior official-schema-validated
artifact, so that external schema check was not needlessly repeated. The unchanged
whole acceptance command peaked at 902.0 MiB versus the prior 1,192.8 MiB; this is
single-run workflow evidence, not the fresh-process three-run stage measurement.

Final affected regression run: **1,301 passed, 1 skipped, 53 deselected** (35.64 s):

```bash
../venv/bin/python -m pytest tests/model tests/io tests/reproject tests/datasets tests/builder/test_meshing.py tests/builder/test_semantic_meshing.py -q --maxfail=3
```

`git diff --check` and compilation of the codec/profiler passed. Self-review found
no remaining blocking issue or needless production abstraction. All earlier
uncommitted work is preserved.
No production dependency, network request, commit or publication.

## Independent-agent handoff

Implement `.agent/plans/2026-09-10-model-performance.md` through its checkpoints.
Keep the plan updated as material decisions or status change, preserve unrelated
work and earlier uncommitted milestones, and run the specified verification.
Follow DESIGN.md and docs/design/geometry-representations.md; use measured costs
and preserve wire bytes, model facts, mutable-array independence and validation.
