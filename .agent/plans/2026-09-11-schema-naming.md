# Consistent DTCC property names

Status: complete, 11 September 2026.
Authority: the user's decision for PascalCase classes and snake_case attributes
and relationships throughout DTCC; DESIGN.md and the model contract.

## Acceptance boundary

Publish standard schema 0.3.0 with snake_case properties, update native projection,
active example profiles and callers, and translate recognized CityJSON names only
at import/export. Preserve arbitrary metadata spelling and reject reserved-name
collisions instead of overwriting data. No change to height meaning, units, schema
coverage, native class hierarchy, dependencies or the Protobuf definition.

Earlier standard schema definitions are preserved byte-for-byte in Git history,
outside the checkout and runtime bundle. No compatibility aliases or
automatic old-attribute migration. Earlier optional profiles remain historical;
current examples use new snake_case versions. Unsupported standard versions fail
default validation; the existing explicit semantic bypass remains available.

## Checkpoints

- [x] Rename the six standard domain properties and native projection; version
  current optional profiles and update executable callers.
- [x] Add explicit, feature-specific CityJSON mappings at both strict and ordinary
  boundaries; reject ambiguous/reserved destination names without touching unknown
  metadata keys or nested values.
- [x] Verify CityJSON → native → .dtcc → CityJSON, default semantic failure,
  intrinsic Tree rules, package version selection and browser example consistency.
- [x] Update current documentation and record observed results.

## Verification

Run focused model profile and affected I/O tests, then the model/I/O suite if changes
cross several adapters. Add a representative naming round trip and collision tests;
preserve existing semantic and failure-atomicity checks. Run the real browser
example with its renamed map keys. Check packaged schema resources because the
runtime version changes. Preserve unrelated edits; do not add broad naming lint
that rejects arbitrary user metadata or external schemas.

## Completion and evidence

Standard schema 0.3.0 renames measured_height, roof_type, crown_radius,
crown_diameter, trunk_diameter and native_landuses; current optional profiles are
city 0.2.0 and buildings 0.3.0. Native projections, tests and executable examples
use the same spellings. Updated current contract/inventory/I/O/browser documentation
and added docs/design/model-naming.md. Historical crosswalk evidence retains its
original names with a prominent current-naming note.

The CityJSON adapter has one explicit, feature-specific mapping in both directions.
It preserves unrelated keys and nested values and rejects destination-name collisions,
including destination-only keys whose spelling could not round-trip. The permissive
loader propagates this specific naming error instead of silently skipping a building.
No global metadata renaming, compatibility alias, height change or new dependency.

Observed verification:

- Initial affected profile/I/O checks: 58 passed; focused new naming checks: 4 passed.
- Full model and I/O suite: **663 passed, 1 skipped** in 32.37 s.
- Real Chromium 152 browser: untouched/edited Python–browser–Python round trips
  passed with schema 0.3.0, exact numerical values and snake_case attribute keys;
  negative measurement, malformed array and unknown wire field were rejected with
  receiver preservation. Current returns were 4,956 and 4,957 bytes.
- Built wheel contains only model/0.3.0 as its standard schema. Archived 0.1.0 and
  0.2.0 files are byte-identical to those in the previous wheel. Installed outside
  the checkout and verified import location, schema resources, default numerical
  measurement validation, failed-write atomicity and CityJSON name translation.
  The smoke reused existing dependency directories; it is not clean cross-platform
  dependency-resolution certification.
- Building, mixed-city, openings and native schema-extension examples passed, along
  with their updated portable evaluators. No stale current-profile spelling remains.
- JavaScript syntax, Python compile checks, current schema naming inspection, local
  documentation links and git diff --check passed.

Evidence: /private/tmp/dtcc-schema-naming-tests.log,
/private/tmp/dtcc-schema-naming-build.log, /private/tmp/dtcc-schema-naming-dist,
/private/tmp/dtcc-browser-naming/report.json and /private/tmp/dtcc-naming-examples.
Earlier standard versions require explicit semantic bypass; there is no migration.
Remaining height/units/expanded-schema decisions are unchanged by this milestone.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-schema-naming.md` through its checkpoints, keep
the plan updated as material decisions or status change, preserve unrelated work,
and run the specified verification. Keep external spellings in explicit adapters;
do not introduce old-name aliases, case conversion of arbitrary metadata or new
height semantics.
