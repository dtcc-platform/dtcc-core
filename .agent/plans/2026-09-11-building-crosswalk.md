# Building semantic crosswalk

Status: complete, 11 September 2026.
Authority: DESIGN.md, docs/design/model-contract.md and the user's request to
continue the CityGML/CityJSON crosswalk after the browser interchange milestone.

## Acceptance boundary

Document the building/part and boundary/opening slice against CityGML 3.0's
conceptual model, CityJSON 2.0.2 and the actual native/schema/adapter behavior.
For relevant properties and relationships, identify source type/cardinality/units,
current DTCC representation and validation, gaps and the recommended next decision.
Include inherited construction facts and distinguish generic metadata preservation
from declared semantic validation. Give deliberate omissions, not a claim of full
CityGML conformance. Keep native access simple and future semantic edits in YAML
where existing native data suffices.

This milestone is an evidence-backed design crosswalk, not an automatic expansion
of the runtime or schema. Do not mutate released semantic version files, invent
permanent vocabulary URIs, add dependencies or implement all building interiors.

## Checkpoints

- [x] Inspect authoritative standards and current schema, model and strict adapter.
- [x] Write a source-linked property/relationship crosswalk and prioritized decisions.
- [x] Exercise representative public entry points to verify material gaps; update
  inventory and plan with observed evidence and a bounded next implementation slice.

## Verification

Use the actual local schema and public save/load APIs to check any behavioral
findings, including height access versus stored measuredHeight and default semantic
validation versus strict CityJSON admission. Preserve unrelated work. Documentation
changes need link/whitespace checks, not a new production test matrix.

## Completion

Added `docs/design/building-semantic-crosswalk.md` and updated the native inventory.
The crosswalk covers identity, all non-ADE AbstractBuilding attributes, inherited
construction properties, relevant lifespan/space facts, building ownership and
surface/opening relationships, plus geometry/access restrictions. It separates
actual validation from generic preservation and lists deferred structured/interior
work. This completes the agreed building slice, not all thematic modules.

Inspected the relevant sections of the full downloaded OGC 20-010 HTML, the
Building/Construction/Core 3.0.0 XSDs, GML basic types and CityJSON 2.0.2. Recorded
the discrepancy between conceptual catalogue association lower bounds and published
XSD occurrence constraints rather than silently selecting one as universally true.
The source files remain in `/private/tmp/dtcc-citygml-*`; conceptual HTML SHA-256:
`ed55c7146e9dbb25f0bb6c89ab23b7368a763776e236b68b93ee972bb28890f3`.

Observed public-API probe: measuredHeight 12.5 survives save/load but height access
returns geometry span 3.0 or separate height attribute 8.0; undeclared negative
storeys/invalid date/structured height survive as metadata; declared scalar name
and measuredHeight reject structured values; strict CityJSON accepts/exports a
negative measuredHeight that default canonical save rejects; strict CityJSON rejects
unmapped address; ClosureSurface receives generic admission. Confirmed the five
induced Building slots and the backend restriction against inlined class slots.
Evidence: `/private/tmp/dtcc-building-crosswalk-probe.py` and
`/private/tmp/dtcc-building-crosswalk-evidence/report.json`.
The completed probe exited successfully. Local documentation links, cited OGC
section anchors, new-file whitespace checks and git diff --check passed.

The next bounded implementation slice is height authority and ordinary Python
access, including existing builder producers/consumers. Basic building attributes
and remaining surface declarations follow in a new immutable semantic version;
default strict-adapter admission and generic inlined records are subsequent steps.
No production code, schema, dependency, wire format or sibling repository changed.
No tests were added to freeze the current gaps as desired behavior.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-building-crosswalk.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Separate verified standards facts,
observed DTCC behavior and recommendations. Do not claim the full CityGML ontology
is covered or change runtime/schema behavior during this crosswalk milestone.
