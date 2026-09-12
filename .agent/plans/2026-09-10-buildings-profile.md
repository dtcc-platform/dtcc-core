# Buildings profile and faithful CityJSON surfaces

Status: completed, 10 September 2026
Authority: DESIGN.md, docs/design/model-contract.md, and the user's requirements
for conceptual simplicity, numerical efficiency, direct Python access and portable
versioned schemas. Earlier canonical/profile plans are completed milestones.

## Acceptance boundary

One complete local workflow: CityJSON 2.0 Building/BuildingPart with MultiSurface
geometry at distinct integer LoDs, polygon holes and roof/wall/ground assignments
-> native DTCC -> canonical file/package -> native DTCC -> CityJSON. Preserve IDs,
attributes, containment, CRS, coordinates and semantic grouping. Unsupported
required input must fail explicitly in the strict workflow. Keep current legacy
routes explicit during migration. No uploads or new production dependencies.

## Decisions

- Record CityGML conceptual mappings and CityJSON differences before declaring
  conformance. The first profile is experimental, not full CityGML certification.
- A self-contained, versioned LinkML YAML is the authority for semantic rules.
  Keep schema evaluation explicit and outside the numerical codec; reuse the
  installed isolated LinkML environment. No global mutable profile or generated
  DTCC class hierarchy. Prove an edited profile accepts a new semantic type with
  no native class or Protobuf change.
- Add one generic SemanticRegion value type: semantic URI, optional local ID,
  typed attributes and a NumPy array of element indices. Geometry.regions stores
  these regions; MultiSurface elements are surfaces, Mesh elements are faces.
  No per-triangle Python semantic objects or copied coordinate arrays.
- Add Surface/MultiSurface and regions to versioned canonical exchange, retaining
  the v1 reader. Preserve polygon rings/holes and nested geometry metadata.
- Keep ordinary access through building.id, building.attributes, building.lod0,
  building.footprint().to_polygon(), geometry.regions and region indices. Do not
  infer physical surface roles from normals or mesh markers.
- Implement strict CityJSON admission as an explicit option on the existing
  loader/exporter, sharing their geometry conversion. Preserve legacy defaults
  until broader type coverage is established; do not hide partial imports behind
  the strict mode.

## Checkpoints

- [x] Standards mapping and precise scope/limitations.
- [x] Self-contained buildings schema and generic native region representation.
- [x] Surface/region canonical exchange with explicit reader compatibility.
- [x] Strict CityJSON workflow and direct Python access example.
- [x] Schema evolution, focused failures, affected regressions and measured codec
      throughput on a repeatable local model; record evidence here.

## Verification

Exercise the public file and package APIs, exact canonical equality, grouping of
multiple polygons under one semantic entity, holes, multiple LoDs, part identity,
missing assignments and invalid indices/references. Check that schema changes do
not enter codec timing. Use an explicit local synthetic benchmark, reporting data
size, run count and hardware/runtime context without claiming city-scale results.
Run affected model/I/O/dataset tests, the existing 18 profile experiment cases and
the new building profile workflow. No broad test matrix or remote data calls.

## Completion evidence

- Mapping, decisions, ordinary Python example and limitations are recorded in
  `docs/design/buildings-profile.md`. The portable profile is
  `schemas/profiles/buildings/0.1.0/schema.yaml`.
- `buildings_example.py` passed public strict CityJSON import/export, native file
  and canonical package round trips: footprint 96 m² with a courtyard hole, part
  identity, two roof polygons sharing one region, and unclassified surfaces.
- The isolated LinkML evaluator passed seven profile checks: portable single-file
  schema, optional/invalid height, a schema-only SolarRoofSurface subtype with a
  required efficiency property, and coexistence of two profile versions. The
  original 18-case experiment also passed all expected outcomes.
- Affected suite: `../venv/bin/python -m pytest tests/model tests/io tests/datasets
  tests/reproject tests/builder/test_meshing.py -q` passed **1,278 tests**, with
  1 skipped and 53 deselected. Includes v1 reader/package compatibility, invalid
  indices/references, strict unsupported-content failures, region reindexing,
  mutable bounds, and rejection at the unsupported C++ conversion boundary.
- Fixture and exported example passed the official bundled CityJSON 2.0.2 JSON
  Schema, retrieved from
  `https://3d.bk.tudelft.nl/schemas/cityjson/2.0.2/cityjson.min.schema.json`
  (SHA-256 `74e128f4505429775fb97da4132a06e2fb3b88116332d7f6893529841fb3505e`).
  The same schema rejected empty MultiSurface boundaries; strict DTCC import and
  export now reject that case too. Empty native MultiSurfaces remain canonical.
  This check is structural schema validation, not full geometric certification.
- Generated protobuf binding matched a fresh generation byte for byte;
  `git diff --check` passed. Existing legacy protobuf schema/binding are unchanged.
- Synthetic benchmark: 500 buildings / 4,000 polygons, Python 3.12.12 on macOS
  26.6.2 arm64, three warmed runs. Median with regions: 1,776,445 bytes,
  154.8 ms encode / 234.3 ms decode; same geometry without regions: 1,616,945 bytes,
  138.8 ms / 210.6 ms. No disk I/O or LinkML in codec timing. This is a baseline,
  not a claim about real-city throughput or peak memory.
- Remaining scope is explicit: solid/shell topology, repeated/fractional LoDs,
  opening relations, richer source metadata, and region-aware triangulation.
  LinkML is still an isolated development dependency; no production profile loader
  or public namespace release was introduced. No data was uploaded or published.

## Independent-agent handoff

Implement `.agent/plans/2026-09-10-buildings-profile.md` through its checkpoints.
Keep the plan updated when material decisions/status change; preserve unrelated
work and the previous uncommitted milestones. Follow DESIGN.md and the user's
simplicity, performance, Python usability and schema portability constraints.
Run the specified verification. Do not add production dependencies, publish data,
or claim unsupported CityGML/CityJSON coverage.
