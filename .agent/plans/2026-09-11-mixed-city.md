# Mixed-city semantic coverage

Status: completed, 11 September 2026
Authority: DESIGN.md, model-contract.md and the user's approved mixed-city step.

## Acceptance boundary

Map DTCC's existing semantic/geometry types against relevant CityGML 3.0 and
CityJSON 2.0.2 concepts. Implement one useful city slice: building, park, tree and
bench, with intuitive existing Python access, optional profile validation and
faithful native file/package exchange. Extend strict CityJSON where the mapping
preserves meaning. No new production dependency, generated class hierarchy,
terrain/road redesign or full CityGML conformance claim.

## Decisions

- Preserve native Tree and Landuse concrete types and their extra fields using
  canonical v5 typed payloads. Existing semantic-only Object still carries bench
  classifications without a new Python class. Retain readers v1–v4 and freeze v4.
- Tree.position/height/crown_radius remain their existing native facts. Landuse
  codes remain their explicit native enum list. Do not mirror these into mutable
  attributes, infer missing measurements, equate land-use codes with source codes,
  or use a float32 legacy serializer for canonical data.
- A new self-contained city profile 0.1.0 incorporates the buildings vocabulary
  and adds LandUse, SolitaryVegetationObject, native Tree and CityFurniture.
  Retain the buildings profile versions. No live schema fetching or imports.
- CityJSON LandUse maps to Landuse with an empty native code list; original
  classifications stay in attributes. SolitaryVegetationObject maps to Object,
  because the source type includes bushes. CityFurniture maps to Object. A
  singleton MultiPoint maps to Point; larger point collections remain unsupported.
- Strict export rejects native Tree's extra state and nonempty Landuse code lists
  until their external mapping is explicitly defined. The native typed example
  and the standards-based mixed source example each round-trip through canonical
  files/packages; the latter also round-trips through strict CityJSON. Do not
  silently downcast Tree, guess vegetation subtypes, or duplicate location data.

## Checkpoints

- [x] Record standards/native coverage and concrete gaps in one coverage document.
- [x] Add canonical typed-state support, strict selected-feature mappings and
      the portable city schema without changing existing class behavior.
- [x] Run the mixed native/source examples, schema-only subtype experiment,
      version coexistence and focused malformed-state/unsupported-export checks.
- [x] Verify affected regressions, official CityJSON schema and unchanged-building
      compatibility; update this plan and docs with evidence and remaining gaps.

## Verification

Use the existing DTCC Python environment for public workflows and tests. Use the
isolated LinkML environment and recorded official CityJSON schema for development
checks. Start with new mixed-city tests, then model/IO/dataset/reproject and affected
meshing regressions. Exercise old wire fixtures and current building examples.
No production tooling or dependency is added just to produce evidence.

## Completion evidence

`docs/design/model-semantic-coverage.md` records the standards/native coverage map,
exact mappings, Python access and remaining gaps. City profile 0.1.0 is a single
self-contained file. The new CityFeature base unifies building, vegetation,
land-use and furniture schema objects; BuildingPart retains its own containment
rules. Existing buildings profiles were not edited. Native Tree/Landuse behavior
is unchanged; canonical v5 now preserves all their existing typed fields.

Both mixed source/native examples passed public file/package workflows. The source
example passed strict CityJSON round trips, preserving each feature's canonical
facts after reimport. Footprint area is 96 m², park area is 360 m², and tree/bench
attributes and point access work. The explicitly authored native Tree retains
float64 position and double measurements; native enum list identity/order remain.
The source plant stays generic and absent measurements remain absent.

The input/exported mixed fixtures passed the recorded official CityJSON 2.0.2 JSON
Schema. The portable city evaluator passed source/native records, optional and
invalid measurements, wrong containment, version coexistence, and a schema-only
ChargingBench with a required attribute. Existing buildings (7 checks) and
openings evaluators and public examples also passed.

Affected regression command:

```bash
../venv/bin/python -m pytest tests/model tests/io tests/reproject tests/datasets tests/builder/test_meshing.py tests/builder/test_semantic_meshing.py -q --maxfail=3
```

Result: **1,332 passed, 1 skipped, 53 deselected** (41.39 s). A subsequent focused
run passed **65 checks** (2.25 s), covering the final malformed-type guard and
format-mapping cleanup. Compilation
of changed Python modules and `git diff --check` passed. Frozen v1–v4 fixtures
remain readable; the v4 fixture rewrites identically after restoring its header.
The complete 24,504,507-byte 3DBAG v3 tile also rewrites identically after restoring
only the old envelope version. No new latency benchmark is claimed for this slice.

Evidence artifacts are in `/private/tmp/dtcc-mixed-city-evidence/`, with prior-profile
reruns under `/private/tmp/dtcc-mixed-city-buildings.*` and
`/private/tmp/dtcc-mixed-city-openings/`. Self-review found no blocking issue or
second native data authority. Native Tree extra state and nonempty Landuse codes
still require explicit external mappings; strict export fails rather than losing
them. General terrain/road/point-collection coverage, Tree position-based bounds,
production profile integration and downstream consumers remain separate work.
All previous uncommitted milestones and unrelated files are preserved.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-mixed-city.md` through its checkpoints. Keep the
plan updated as material decisions or status change, preserve all previous
uncommitted milestones and unrelated files, and run the specified verification.
Follow DESIGN.md; keep native field and geometry authorities singular and avoid
inferring domain meaning from incomplete external classifications.
