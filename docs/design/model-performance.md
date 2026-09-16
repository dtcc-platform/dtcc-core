# DTCC Model performance: complete 3DBAG tile

Status: first measured codec optimization, 10 September 2026. This continues the
[representation milestone](geometry-representations.md) under [Core Design](../../DESIGN.md).
The native model, semantic schema and canonical wire version/bytes are unchanged.

## Workload and method

The same unchanged [3DBAG tile](real-buildings.md) contains 1,110 Buildings,
1,111 Parts, 4,443 representations and 68,399 polygon surfaces. Source SHA-256:
`2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba`.
Its canonical v3 artifact is 24,504,507 bytes.

Measurements use Python 3.12.12 on macOS 26.6.2 arm64 and the existing DTCC virtual
environment. Each stage runs three times in fresh processes, sequentially, with
no simultaneous benchmark/test work. The table reports medians. OS file caches
are not reset. Stage timings exclude Python startup, imports and input preparation;
import/export/package timings include their public I/O operation. Encode/decode
use the ordinary codec entry points. Source import explicitly recomputes extents
and retains their discrepancies in Dataset Context, as in the prior acceptance run.

Peak RSS includes process imports, preparation and the result while it is still
alive. Each sample also records RSS and peak RSS after imports and preparation.
Peaks are **not** incremental allocations or native model size. The environment
loads roughly 200 MiB before these public workflows run. Call profiling is a
separate, instrumented measurement; its timings are not mixed into the table.

## Measured change

| Operation | Before | After | Less elapsed time | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| Strict CityJSON import | 3.156 s | 2.693 s | 14.7% | 380.2 → 380.8 MiB |
| Canonical encode | 1.908 s | 1.195 s | 37.4% | 590.2 → 517.6 MiB |
| Canonical decode | 3.045 s | 2.200 s | 27.7% | 397.8 → 397.7 MiB |
| Canonical package write | 2.082 s | 1.446 s | 30.6% | 591.5 → 518.5 MiB |
| Canonical package read | 3.095 s | 2.267 s | 26.8% | 423.0 → 424.4 MiB |
| Strict CityJSON export | 2.919 s | 2.330 s | 20.2% | 381.6 → 382.9 MiB |

Small RSS differences outside encode/package-write are not treated as improvements.
The source/native geometry buffers are unchanged: 228,864 arrays containing
19,127,848 bytes (18.24 MiB), including transforms, shell indices and regions.
That buffer count excludes Python objects, dictionaries, lists, bounds caches,
array headers, native libraries and temporary conversion structures. A model has
2,222 Objects and 72,842 Geometry instances including nested polygon Surfaces;
these per-polygon structures explain why raw coordinate size alone is not a useful
estimate of process memory. The earlier 1,192.8 MiB peak covered a whole repeated
verification workflow, not one model or one codec call.

## Why these changes

The initial call profile attributed about 3.01 of 4.11 instrumented encode seconds
to native validation, with 2.31 seconds in transform validation. Each validation
pass invoked scalar precision checks 1,201,024 times. Decoding constructed 147,906
Transforms even though the model only needs 75,064: most geometries first created
a default Transform and immediately replaced it with the decoded Transform.

Afterward, the separate call profile records 75,064 Transform constructions,
with no scalar transform precision calls for this float64-only tile. Instrumented
encode validation fell from 3.01 to 1.36 seconds. These call counts substantiate
the eliminated work; the uninstrumented table above is the performance result.

The implementation changes only [the canonical codec](../../dtcc_core/model/exchange.py):

- After validating shape, real numeric dtype, finite coefficients and the final
  affine row, float16/32/64 coefficients need no scalar exactness loop: every such
  value is exactly representable as a wire double. Integer and wider-float
  coefficients retain the original exactness checks. Mutable data is validated
  anew at every exchange boundary; there is no cached validity flag.
- Decode supplies transforms, fields, regions and the absent bounds cache to each
  native constructor directly. It no longer allocates default values that it will
  immediately replace. Arrays remain independent and writable.
- Encode fills nested Protobuf destinations directly. It no longer constructs
  complete intermediate geometry/object trees and copies those subtrees into their
  parents. Numerical arrays, attribute mappings and the wire schema are unchanged.

Unknown-field/version checks, containment/reference validation, field association,
array dtype/shape/connectivity, Solid shell partitions, semantic regions, size
limits and atomic persistence remain in place. No production dependency, alternate
model representation, global cache, shared mutable default or new wire version was
introduced.

## Reproduce

First produce the existing complete-tile acceptance artifacts:

```bash
../venv/bin/python sandbox/model_profiles/representations_example.py /path/to/9-284-556.city.json /tmp/dtcc-representations
```

Then run the stage measurements. `psutil` is used only by this development script
and was already installed in the measured environment.

```bash
../venv/bin/python sandbox/model_profiles/profile_model.py /path/to/9-284-556.city.json /tmp/dtcc-representations /tmp/dtcc-performance
```

The output directory contains all per-run JSON measurements, operation logs and
one summary `report.json`. Use a different output directory for each revision.
For a separate call profile, add `--profile --repeats 1`; this writes `.prof` files
readable with Python's `pstats`. Instrumented times are not comparable to the
ordinary timing table. The script performs no network operations.

Before/after measurements from this run are retained locally in
`/private/tmp/dtcc-model-perf/before/` and `/private/tmp/dtcc-model-perf/after/`.
The measured source file was reused unchanged for both runs.

## Verification and remaining work

The complete pre-optimization 3DBAG canonical artifact decodes and re-encodes
byte-for-byte identically. Focused checks retain malformed-input and fail-before-
replacement behavior, frozen v1/v2 compatibility, and add explicit checks for
nonfinite float64 transforms and independent mutable decoded polygon data.

The complete public file/package/CityJSON acceptance command also passed after
the optimization, including every representation, shell, region, source defect
and recorded extent discrepancy. Its CityJSON export is byte-for-byte identical
to the previously official-schema-validated export (24,116,694 bytes); schema
validation was not rerun on those identical bytes.

The unchanged complete verification command peaked at **902.0 MiB**, compared
with the prior **1,192.8 MiB** (24.4% lower). These are single whole-workflow runs,
separate from the three-sample stage medians. The repeated codec medians inside
that run were 1.163 s encode and 2.238 s decode; construction, disk and LinkML
remain outside those two timings.

The final affected regression suite passed: **1,301 passed, 1 skipped,
53 deselected** across model, I/O, reprojection, datasets and affected meshing
tests. `git diff --check` and Python compilation also passed. Commands and
completion evidence are in the completed plan:
[model performance](../../.agent/plans/2026-09-10-model-performance.md).

This is one measured optimization, not a city-scale performance guarantee. Decode
memory and native storage are essentially unchanged. Future storage changes should
first distinguish live Python objects, NumPy array headers, Protobuf allocations
and allocator retention using representative larger inputs. Redesigning polygon
storage or changing wire defaults would need a separate contract decision; the
current measurements do not justify introducing a second geometry authority.


The following [live-memory pass](model-live-memory.md) separates native allocations
from process RSS and reduces repeated immutable metadata while preserving this
codec's bytes, validation and numerical arrays.
