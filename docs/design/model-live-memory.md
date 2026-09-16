# DTCC Model live memory

Status: bounded memory pass, 10 September 2026. This follows the
[codec performance work](model-performance.md) and preserves the
[representation contract](geometry-representations.md) under [Core Design](../../DESIGN.md).

## Result

Reusing equal immutable strings within one canonical decode reduces live traced
allocations by **6.7%** on the recorded 3DBAG tile. No native class, numerical
buffer, public accessor, canonical wire byte or validation rule changed.

| Controlled workload | Before | After | Saving |
| --- | ---: | ---: | ---: |
| One decoded tile | 84.88 MiB | 79.22 MiB | 5.66 MiB |
| Three independently decoded copies | 254.45 MiB | 237.46 MiB | 16.99 MiB |

These are live allocations tracked by `tracemalloc` after loading and garbage
collection. They include tracked Python/NumPy allocations, not every native
allocation or the interpreter's pre-existing memory. They are not process RSS.
Three copies form a synthetic collection of independent models; this is not a
larger real geographic dataset. It contains 3,330 Buildings, 3,333 Parts and
205,197 polygon Surfaces. Each model retains its own ID/reference scope.

The separate uninstrumented runs observed loaded-process RSS of 382.5 → 375.1 MiB
for one tile and 652.9 → 606.8 MiB for three. These single-run RSS measurements
include roughly 185 MiB of imports and retained loader/allocator memory, and are
not a precise measure of model size or the saving attributable to string reuse.
The live-allocation table is the evidence for the memory reduction.

A separate warmed, alternating-order comparison of the old/new codec measured
three-run median decode times of **2.056 and 2.060 seconds**. This shows no material
latency change in that comparison. Instrumented timings were not used to assess
speed. Both revisions consumed the same 24,504,507-byte canonical artifact.

## Where the memory goes

The one-tile inventory contains:

| Stored state | Count or size |
| --- | ---: |
| Numerical array buffers | 19,127,848 bytes (18.24 MiB) |
| NumPy arrays, including empty arrays | 228,864 |
| Empty NumPy arrays | 70,581 |
| Polygon Surface objects | 68,399 |
| Transform objects | 75,064 |
| Empty lists/tuples | 210,574 |
| Empty dictionaries | 11,112 |
| Distinct string contents | 3,707 |

The two largest initial live-allocation sites were decoded arrays and affine
array copies: about 25.4 and 16.0 MiB respectively, including their tracked storage
and overhead. Bounds caches are absent immediately after canonical loading;
removing such caches cannot reduce this measured load footprint. Instance
`__dict__` mappings were not materialized by the inventory, because inspection
itself can otherwise change allocation behavior.

Before the change, 101,045 string objects occupied 5.90 MiB although their distinct
contents needed only 0.23 MiB. Repeated attribute names, semantic URIs, CRS strings,
LoD labels and representation IDs explain this avoidable duplication. Afterward,
this tile has 3,707 stored string objects. Strings remain case-sensitive and are
never normalized, interpreted or fetched from a schema by the reader.

After deleting the loaded models and collecting garbage, tracked live allocations
returned to about 0.085 MiB above the tracing start for both one and three copies.
The check therefore found no retained Python/NumPy model allocation leak. RSS did
not return to the import baseline; tracked liveness alone cannot attribute all
remaining native/allocator memory. Instrumentation and the inventory themselves
also consume memory, so their RSS peaks are not application peaks.

## Implementation boundary

Only the canonical reader in [exchange.py](../../dtcc_core/model/exchange.py)
changes. A decode-local dictionary reuses equal strings and is discarded after
the operation. Its bookkeeping is capped at 4,096 entries; once full, unseen
strings are decoded normally and already-known strings can still be reused.
Data with many unique strings remains valid and lossless. This is an internal
allocation choice, not a profile vocabulary, public registry or process-wide cache.

Mutable dictionaries, lists, transforms and numerical arrays remain independent.
Editing one object's metadata cannot edit another object's metadata. All existing
boundary validation and file/package integrity checks remain in place. The writer
and the wire schema are unchanged; v1/v2/v3 readers retain their previous coverage.

The remaining overhead belongs largely to the current per-polygon native object
layout. Meaningfully larger reductions would need a separate evaluation of lazy
optional state or packed polygon storage, with explicit treatment of mutation,
copying and Python access. This pass introduces neither. The next domain milestone
can return to building openings and their relationships to boundary surfaces.

## Reproduce

Use the canonical output of the existing full-tile acceptance example. The
recorded artifact SHA-256 for this comparison is
`e8d9b7a88754bab5b1c87ba1f43bcb211498e539e4ff4d2dc871294af9d3e909`.
Its source tile and source checksum are recorded in the
[representation report](geometry-representations.md).

Run each command in a fresh process using the existing DTCC environment:

```bash
../venv/bin/python sandbox/model_profiles/profile_model_memory.py /path/to/model.dtcc /tmp/memory-1-rss.json
../venv/bin/python sandbox/model_profiles/profile_model_memory.py /path/to/model.dtcc /tmp/memory-1-trace.json --trace
../venv/bin/python sandbox/model_profiles/profile_model_memory.py /path/to/model.dtcc /tmp/memory-3-rss.json --copies 3
../venv/bin/python sandbox/model_profiles/profile_model_memory.py /path/to/model.dtcc /tmp/memory-3-trace.json --copies 3 --trace
```

The instrumented mode records live allocation sites and a native-state inventory;
the uninstrumented mode omits that inventory so its post-release RSS is not
confounded by a diagnostic graph walk. Traced post-release RSS includes tracing
and inventory effects; only its tracked-live value is used for the release check.
Reports include artifact hash, platform, baseline, loaded and released snapshots.
Raw evidence from this run is in `/private/tmp/dtcc-model-memory/{before,after}/`.
The timing comparison is in `/private/tmp/dtcc-model-memory/decode-timing.json`.
Python 3.12.12 and macOS 26.6.2 arm64 were used; no production dependency or network
operation was added. Profiling remains development tooling.

## Verification

The full recorded canonical tile re-encodes byte-for-byte identically. Focused
checks cover independently editable repeated/nested metadata, Unicode/empty
strings, 6,000 unique string values, numerical-array independence, malformed data
and the previous version fixtures. The full file/package/CityJSON workflow passed;
its CityJSON output is byte-identical to the prior official-schema-validated
artifact, so the schema check was not repeated on identical bytes. The affected
suite passed **1,302 tests**, with 1 skipped and 53 deselected. Diff whitespace
checks and Python compilation also passed. Commands and completion evidence
are recorded in the [completed plan](../../.agent/plans/2026-09-10-model-live-memory.md).
