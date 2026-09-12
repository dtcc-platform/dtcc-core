# Geometry representations and Solid exchange

Status: completed, 10 September 2026
Authority: DESIGN.md and the user-approved contract in
`docs/design/geometry-representations.md`. Earlier milestones are completed.

## Acceptance boundary

One authoritative `Object.geometry` mapping of local IDs to representation
records; exact LoD/role selection, explicit ambiguity and native convenience
properties. Solid preserves polygon surfaces, shell membership and regions.
Complete unchanged 3DBAG import -> native file/package -> CityJSON with explicit
source precision/validity handling. Preserve prior uncommitted work.

## Checkpoints

- [x] Implement records/selectors and migrate in-repository geometry consumers.
- [x] Add Solid and ordered canonical v3 exchange, retaining v1/v2 readers and
      explicit fail-on-loss legacy exports.
- [x] Extend strict CityJSON to exact LoDs, repeated representations and Solids.
- [x] Exercise ambiguity, interior shells, compatibility and real tile workflow;
      measure codec bytes/time/peak memory and run affected regressions.
- [x] Update contract, examples and completion evidence; no replacement stores,
      production dependencies or publication.

## Verification

Run narrow representation/exchange tests first, then model/I/O/dataset/reprojection
and affected builder tests. Migrate storage assertions to the approved record
contract, retaining their behavioral invariants. Exercise public file and package
APIs on the local checksum-pinned 3DBAG tile. Check source extents against its
coordinate precision without inventing missing semantic facts. Keep malformed
geometry distinct from unsupported topology and avoid source repair.

## Material decisions

- V3 uses repeated representation entries; old map field 5 remains v1/v2-only.
  Old readers migrate sorted slots. Legacy export rejects loss of new descriptors,
  Solid or insertion order. Direct enum-key dictionary consumers migrate to the
  single record store; downstream repositories still need a coordinated update.
- Strict source extent validation remains the default. Explicit recomputation
  retains all discrepant original summaries in root Dataset Context; packages
  preserve Context. Half the source coordinate grid is the tolerance per axis.
- Preserve pre-existing collapsed rings without repair. Reject new collapse below
  three distinct vertices caused by export quantization. Meshing still rejects
  unusable polygons. Solid meshing/reprojection and semantic MultiSurface
  reprojection remain unsupported, with explicit errors at affected entry points.
- Source bounds exposed a mutable union alias: parent aggregation now copies its
  first bounds value, preventing corruption of child/surface bounds.

## Completion evidence

Final affected regression run: **1,299 passed, 1 skipped, 53 deselected** (31.54 s):

```bash
../venv/bin/python -m pytest tests/model tests/io tests/reproject tests/datasets tests/builder/test_meshing.py tests/builder/test_semantic_meshing.py -q --maxfail=3
```

Earlier focused checks also passed: 90 object/representation checks and 65
reprojection checks. Meaningful failure coverage includes ambiguity, malformed
shells, extent discrepancies, unsupported reprojection, legacy loss and export
quantization failure without overwriting an existing destination.

The checksum-pinned real workflow in `sandbox/model_profiles/representations_example.py`
passed for the complete unchanged 3DBAG tile: 1,110 Buildings, 1,111 Parts, 4,443
representations, 3,333 Solids and 68,399 polygon surfaces. Canonical file/package
round trips are exact, including package Context. CityJSON round trip checks all
feature IDs, attributes, attachments, shells and regions; coordinate error is at
most 0.000375 m. Four existing collapsed rings survive and 1,111 discrepant source
extent summaries are retained. The exported tile passed the official CityJSON
2.0.2 JSON Schema (101.5 s). Local outputs and full measurement report are in
`/private/tmp/dtcc-representations-evidence/`; durable measurements and reproduction
instructions are in `docs/design/geometry-representations.md`.

Three-run medians: 24,504,507 bytes, encode 1.77 s, decode 3.01 s. Import 3.05 s;
whole verification process peak 1,192.8 MiB, not isolated codec allocation. This
is measured evidence for further profiling, not an efficiency completion claim.

Both public profile examples passed with v3 file/package artifacts. The isolated
LinkML evaluators passed all 7 buildings-profile checks and all 18 original
experiment outcomes. No LinkML production dependency or schema vocabulary change.
`git diff --check` and Python compilation passed. No publication or commit.

## Independent-agent handoff

Implement `.agent/plans/2026-09-10-representations-solid.md` through its checkpoints.
Keep it updated as material decisions/status change, preserve unrelated work and
the previous uncommitted milestones, and run the specified verification. Follow
`docs/design/geometry-representations.md` and DESIGN.md. Do not add production
dependencies, publish data or introduce a second geometry store.
