# DTCC Model live-memory investigation

Status: completed, 10 September 2026
Authority: DESIGN.md and the user-approved next step following the completed
model-performance milestone.

## Acceptance boundary

Measure live native-model allocations after loading the complete 3DBAG tile,
separating data buffers, object/array overhead, caches and retained loader memory.
Check scaling with a controlled larger workload. Optimize only an avoidable cost
that can be removed without changing Python access, independent mutable arrays,
wire bytes or validation. A storage redesign requires a separate design decision;
do not introduce one to force a memory improvement. Preserve all prior work.

## Checkpoints

- [x] Profile live memory, ownership and post-release retention for one and three
      independent copies of the recorded tile.
- [x] Remove a measured avoidable allocation if it fits the existing contract;
      otherwise document the boundary and recommendation with concrete evidence.
- [x] Verify file/package/CityJSON behavior and relevant regressions, document
      measured results and keep reproduction commands in the existing sandbox.

## Verification

Separate instrumented live-allocation accounting from uninstrumented RSS/time.
Use the existing recorded canonical tile and public load entry point. Keep source
and canonical bytes unchanged. Exercise mutations and meaningful failure paths
if code changes affect those boundaries. The larger workload is explicitly a
synthetic collection of independent copies, not another real city.

## Findings and completion evidence

The existing native layout owns 228,864 arrays (18.24 MiB of buffers), 68,399
Surfaces, 75,064 Transforms and 210,574 empty sequences per tile. The decoded
model has no bounds caches to remove. Its 101,045 string objects have only 3,707
distinct values. A bounded decode-local string table removes this duplication;
no global cache, mutable sharing, class-layout or wire change was introduced.

Live traced allocations: 84.88 -> 79.22 MiB for one tile; 254.45 -> 237.46 MiB for
three independently decoded copies (6.7% less). After release, traced live memory
returns to about 0.085 MiB for both workloads. Uninstrumented RSS and instrumented
live memory are recorded separately; the diagnostic inventory is excluded from
uninstrumented measurements. See `docs/design/model-live-memory.md` for full
methods/limits and `/private/tmp/dtcc-model-memory/{before,after}/` for raw data.

Three-run warmed alternating-order codec comparison: 2.056 -> 2.060 s decode,
with no material latency change observed. The complete canonical tile re-encodes
byte-for-byte unchanged. The full public file/package/CityJSON workflow passed, preserving all tile
representations, shells, regions, source defects and extent evidence. Its CityJSON
output is byte-identical to the previous official-schema-validated artifact;
that schema check was not repeated on identical bytes. Whole-workflow outputs
are in `/private/tmp/dtcc-model-memory/workflow-after/`.

Final affected regression run: **1,302 passed, 1 skipped, 53 deselected** (40.82 s):

```bash
../venv/bin/python -m pytest tests/model tests/io tests/reproject tests/datasets tests/builder/test_meshing.py tests/builder/test_semantic_meshing.py -q --maxfail=3
```

`git diff --check` and compilation of the reader/profiler passed. Self-review found
no blocking issue; string reuse is bounded and private, and existing mutable
containers/arrays, validation and wire contracts are preserved.

Remaining substantial costs are structural. A lazy-state or polygon-storage
redesign is explicitly deferred; next domain work can return to openings and
boundary relationships. Previous uncommitted work is preserved.

## Independent-agent handoff

Implement `.agent/plans/2026-09-10-model-live-memory.md` through its checkpoints.
Keep the plan updated as material findings/status change, preserve unrelated work
and earlier uncommitted milestones, and run the specified verification. Follow
DESIGN.md and the existing representation contract. Do not introduce shared
mutable defaults, a second geometry store or a new public API to force a saving.
