# 3DBAG semantic mapping, then terrain and transportation crosswalk

Status: complete, 12 September 2026. Authority: user-approved slice 1–4 followed directly by crosswalk;
DESIGN.md and docs/design/model-contract.md. No repository plan template exists.

## Acceptance boundary

One explicit loader for the audited 3DBAG b3_h_dak attribute convention. Preserve
supplier attributes, owners and geometry; map elevations with an explicit vertical
reference and physical roof codes. No inferred building height, acquisition date,
release identification or automatic code translation. Extend the versioned schema
using existing nested records, without native classes, wire changes or dependencies.
Then document the property-level terrain/transport comparison and a bounded next
implementation slice; do not implement the entire CityGML hierarchy.

## Checkpoints

- [x] Audit primary 3DBAG definitions and cached tile; document mapping decisions.
- [x] Implement schema 0.7.0 and explicit public loader, preserving source evidence.
- [x] Exercise public load/access/native/package/CityJSON, meaningful failure,
      representative performance and affected tests; installed-wheel smoke.
- [x] Cross-check terrain and transportation against CityGML 3.0, CityJSON 2.0.2
      and actual native/strict adapter code. Document gaps and next acceptance slice.

## Verification

Focused synthetic tests cover missing values, negative elevations, collisions,
wrong CRS, malformed source numbers, qualified record validation and persistence.
The pinned 1,110-building real tile proves unchanged source attributes/geometries,
5,550 elevation records, native/package preservation and simple Python access.
Measure schema evaluation before/after enrichment and default native I/O, separating
warm timings from initialization. Run model/I/O checks and wheel smoke. Crosswalk
claims distinguish native storage from implemented standards conversion.

## Completion and limitations

Implemented `io.load_3dbag`, standard 0.7.0 ElevationMeasurement and the public
real-data example. Existing base classes, wire v6, height convenience properties
and generic record validator are unchanged. Explicit source mapping precedes the
single semantic evaluation. Source CRS, collisions and malformed values fail;
unknown physical classification meanings stay supplier metadata.

Focused mapping/record checks: 7 passed in 4.08 s; final focused check after review:
7 passed in 3.87 s. Model/I/O regression: 686 passed, 1 skipped in 46.63 s.
Real tile: all 1,110 buildings, 5,550 elevations and 1,110 qualified roof codes;
removing enrichment recovers the complete original native payload. Native and
canonical package round trips are exact, including package Dataset Context.
Strict CityJSON retains attributes, IDs, representations, shells and regions;
maximum coordinate error 0.000375 m at 0.001 m quantization. The first example
attempt incorrectly demanded byte equality across CityJSON; diagnosis isolated
coordinate quantization, and the example now reuses the existing representation
comparison with its 0.00050001 m tolerance. No production geometry change was made.
Removing the required vertical reference rejects a save without replacing the file.

Warm medians: unchanged/enriched schema evaluation 0.2894/0.6926 s (five alternating
runs); enriched native encode/decode 1.9317/2.9157 s (three runs). Enriched native
payload 26,854,163 bytes, previously 24,504,572. Loader with initialized schema
3.4185 s (one run). No extra optimization or production evidence machinery added.
Artifacts and raw samples: `/private/tmp/dtcc-3dbag-mapped/report.json`.

Built and installed the wheel outside the checkout. Its loader/schema match final
source exactly; only standard 0.7.0 is bundled. The archived 0.6.0 matches the
previous milestone wheel byte-for-byte. Installed public loader on the real tile,
native/package equivalence, version, missing-reference path, atomic rejection and
explicit bypass all passed. Existing dependencies reused without resolution;
no new platform/browser consumer certification. Relative document links and
`git diff --check` pass. Unrelated dirty-tree work is preserved; no publication.

Crosswalk completed in `docs/design/terrain-transportation-crosswalk.md`, with
Relief/Transportation property tables, conceptual-vs-XSD cardinality differences,
CityJSON encoding differences and actual native/strict/permissive adapter gaps.
Primary XML snapshots in `/private/tmp`:
Relief SHA256 `e07afafa959e7e4c5e02c63bf96b21683639972af0731676413e307c4c0a8c64`;
Transportation SHA256 `ce6621639f1263c8596d73553da87bce7f6d48e0dd3ea0b45b68c8d762961a7d`.
The recommended next implementation is strict triangular TINRelief exchange with
existing Terrain/Mesh, retaining all features/representations; it is not implemented
by this crosswalk. DEM construction and physical-road/routing semantics remain
explicit future workflows, with concrete risks recorded rather than silently fixed.

Cached source release is unspecified. The mapping supports its named attribute
convention, not every 3DBAG release, and retains source owners without assigning
roof-statistic membership to Parts. No inferred height, capture status/date, datum
conversion, code-list service or new raster unit semantics. Source validity codes
remain evidence, not certification. Dataset Context requires canonical packaging;
standalone native files carry the qualified record source and original attributes.

Logs: `/private/tmp/dtcc-3dbag-{focused,final-focused,regression,example,build,install,installed}.log`.
Wheel: `/private/tmp/dtcc-3dbag-dist`.
Installed smoke: `/private/tmp/dtcc-3dbag-installed-smoke.py`.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-3dbag-terrain-transport.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Finish the 3DBAG slice and then the
terrain/transport crosswalk. Keep one schema validator and existing generic data
carriers; do not infer height differences or expand to full standards conformance.
