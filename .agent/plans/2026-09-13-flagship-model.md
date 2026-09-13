# Flagship native model fixture

Status: complete locally, 13 September 2026. Authority: the user requested a rich development `.dtcc`
test case combining detailed CityJSON with real or synthetic enrichment.
No applicable PLAN_TEMPLATE.md exists in this workspace.

## Acceptance boundary

Create one reproducible, moderately sized neighbourhood fixture from the existing
checksum-pinned 3DBAG CityJSON tile. Preserve selected source buildings and all
their LoDs and parts. Add an explicitly synthetic demonstration area and numerical
data using the existing model and bundled schema. Include provenance/licensing
inside the standalone file, a canonical package, a readable inventory and preview,
and examples of ordinary Python access. Do not extend the native format/schema,
add dependencies, invent real observations, or publish/deploy the new artifacts.

## Checkpoints

- [x] Select and attribute real source data; author deterministic enrichment.
- [x] Generate native file, canonical package, source subset and visual inventory.
- [x] Verify default validated save/load, exact native/package preservation,
      deterministic native output, and a meaningful rejected edit.
- [x] Document access/reproduction and review the generated preview and code.

## Verification

Run the generator through public I/O with the pinned real tile, then reload its
outputs. One focused offline regression protects synthetic fixture invariants and
field/reference failure behavior. Verify unchanged source building geometry and
metadata, exact native re-encoding, package context, and representative field and
opening access. Run the generator twice and compare native bytes. Keep large
generated binaries in ignored `data/flagship`, with source/generator/docs in Git.

## Completion evidence

Generated `data/flagship/flagship.dtcc`: 1,725,527 bytes; 163 objects, 309
representations, 19 fields, all eleven native geometry classes plus Raster and all
six field associations. Real data: 64 buildings and 64 parts, selected from the
existing checksum-pinned 3DBAG sample; original source facts remain exact after
removing only the added boundary-mesh representation. Synthetic objects/fields
are explicitly labelled, with attribution and changes embedded in the native file.

The companion package preserves native bytes and Dataset Context exactly. The
compact source CityJSON, preview, inventory and README accompany it. The source
tile is retained in ignored `data/flagship-source/3dbag.city.json`.
Two complete final generations produced identical native SHA-256
`02abb1a57c552816f55c35acaabfdfaf199a57719690265cfe244f4fa603eadc`.
The second took 11.115 seconds in the current environment; this is not a benchmark
guarantee. Generated ZIP timestamps and inventory timing are not deterministic.

29 focused tests passed (flagship regression and existing canonical exchange).
The copy-and-run Python documentation passed against the real artifact, including
49-vertex footprint access, 315 velocity vectors, and a 42×55 DEM. Invalid storey
counts and dangling references fail without replacing the valid file. The preview
was visually inspected; sorting all faces in one collection fixes terrain
occlusion, and the pavilion detail triangulates the opening holes correctly.

The generator now runs without arguments from either the checkout root or scripts/,
using checkout-relative source/output defaults. That real invocation regenerated
identical native bytes; a missing source produces an actionable CLI error.

No runtime/schema/dependency changes or remote publication for this fixture. The new
generator, guide and focused test are repository files; binaries are ignored local
development artifacts. This does not certify source geometric validity, full
CityJSON export of enriched state, or solver/observation accuracy.

## Independent-agent handoff

Implement `/Users/logg/scratch/dtcc/dtcc-core/.agent/plans/2026-09-13-flagship-model.md`
through its checkpoints, keep the plan updated as material decisions or status
change, preserve unrelated work, and run the specified verification. Keep existing
model/schema authorities and label all synthetic additions explicitly.
