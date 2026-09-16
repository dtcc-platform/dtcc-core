# Basic building attributes and boundary vocabulary

Status: complete, 11 September 2026.
Authority: user-approved next step; DESIGN.md, model contract and the reviewed
building semantic crosswalk. Prior height milestone is complete.

## Acceptance boundary

Publish immutable standard schema 0.5.0: optional description/class strings,
function/usage string lists and nonnegative integer storeys_above_ground and
storeys_below_ground. No code coercion, invented enums, default counts or inferred
purposes. Use the existing LinkML/JSON Schema number semantics; booleans and
fractional counts must fail. Preserve open undeclared metadata and native values.

Declare ClosureSurface, OuterCeilingSurface, OuterFloorSurface,
InteriorWallSurface, CeilingSurface and FloorSurface as semantic regions. Express
Building and BuildingPart as sibling schema classes under abstract AbstractBuilding;
retain current object ownership and opening-host restrictions. No new Python
classes, validator, dependency, geometry certification, qualified records or wire
layout. Add only explicit CityJSON storey-name mappings. Default semantic admission
in external adapters is the next milestone, not part of this one.

## Checkpoints

- [x] Add schema 0.5.0 and archive 0.4.0 unchanged; align current version selection.
- [x] Map CityJSON storey spellings without changing values or generic metadata.
- [x] Exercise public CityJSON → .dtcc → CityJSON with Building/Part attributes
  and all six surface declarations; verify new type, count and owner failures.
- [x] Update current contract/inventory and document scope, sources and version
  consequences. Record focused tests and installed-wheel resource/workflow proof.

## Verification

Run focused I/O/schema tests first, including boolean versus integer counts, scalar
versus list codes, missing values, failed-save preservation, optional bypass and
known-surface owner admission. Then run model/I/O tests for shared inheritance
and persistence effects. Build and inspect the wheel and smoke default I/O outside
the checkout. No browser rerun for an unchanged wire layout; check its fixture
version and syntax. Preserve unrelated work; do not run live datasets.

## Completion and evidence

Implemented standard schema 0.5.0 with the six optional attributes and six
surface declarations. Building and BuildingPart now share AbstractBuilding as
sibling schema classes; existing CityFeature name/parent rules are inherited
instead of duplicated. Part containment and surface ownership still accept both
building kinds. No native class, backend validator, dependency or wire change.
The existing CityJSON name adapter maps the two storey spellings bidirectionally
and rejects ambiguous keys; generic/nested metadata remains intact.

Added docs/design/building-attributes.md with exact attribute rules, code-space
limits, JSON Schema integrality, surface meaning, version consequences and links
to the reviewed primary sources. Updated current contract, inventory, naming,
height and schema-I/O documentation. The historical crosswalk retains original
observations with an updated implementation note. The browser fixture selects
0.5.0 and its earlier observed browser evidence remains labelled 0.3.0.

Observed verification:

- Initial focused run: 21 passed and one fixture failure because semantic regions
  require MultiSurface/Mesh/Solid rather than Surface. Corrected that fixture to
  MultiSurface; the targeted owner and generic-metadata checks then both passed.
- Full model/I/O suite: **672 passed, 1 skipped in 35.22 s**. This includes the
  public CityJSON → .dtcc → CityJSON workflow with Building/Part metadata and all
  six added region classes; numeric/code/list rejection, boolean counts, explicit
  semantic bypass and failed-write preservation; standalone geometry owners and
  abstract-type rejection; existing opening-host and height rules.
- Wheel build passed. It bundles only standard schema 0.5.0. Archived 0.4.0 is
  byte-identical to the previous height-milestone wheel.
- Installed outside the checkout: confirmed package import/resource location and
  default schema selection, repeated the Building/Part attribute and six-surface
  external/native round trip, and rejected a boolean storey count while preserving
  the existing file. This reused existing dependencies; it is not a fresh platform
  dependency-resolution test.
- JavaScript syntax, local documentation links, schema diff self-review and
  git diff --check passed. No browser rerun, live dataset run or new dependency.

Evidence: /private/tmp/dtcc-building-attributes-focused.log,
/private/tmp/dtcc-building-attributes-tests.log,
/private/tmp/dtcc-building-attributes-build.log,
/private/tmp/dtcc-building-attributes-dist,
/private/tmp/dtcc-building-attributes-installed-smoke.py and
/private/tmp/dtcc-building-attributes-installed.log.

Declared integral values retain JSON Schema semantics (2.0 accepted, True and 2.5
rejected) without native conversion. The scope does not certify orientation,
watertightness, rooms, code-space records or CityGML conformance. Strict CityJSON
still admits native structure rather than automatically evaluating the selected
standard schema; applying that semantic boundary, including explicit bypass, is
next. Existing representation-ID/role export restrictions remain unchanged.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-building-attributes.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Keep semantic rules in YAML and external
spelling conversion in the existing adapter; do not broaden to qualified records,
new Python classes, CityGML conformance or external default schema admission.
