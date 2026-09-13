# Complete the first exterior-city profile

Status: complete (12 September 2026). Authority: user requested a finite boundary, a plan first, subagents,
and implementation through completion in this task. DESIGN.md and the model contract
control scope. No applicable AGENTS.md or PLAN_TEMPLATE.md was found in this repo.

## Completion boundary

Deliver an explicitly bounded exterior-city profile on the existing generic native
model, with one new standard schema version and one current coverage document.
Retain implemented buildings/parts, land use, furniture, individual vegetation and
triangular TINs. Add physical transportation features (Road, Railway, Waterway,
TransportSquare), WaterBody and PlantCover, with meaningful surface classifications
and common qualified class/function/usage metadata. Support the line and polygon
representations needed by those features through strict CityJSON and canonical
native/package exchange, preserving source geometry distinctions explicitly.

Complete a native terrain-raster interpretation workflow with explicit elevation
unit, vertical reference and sampling meaning, and address concrete ground-filter,
nodata and CRS risks in its builder. A raster stays a raster in native exchange;
no automatic CityJSON conversion is claimed. Geometry numbers remain NumPy arrays.

Consolidate the schema namespace/identity, examples, coverage and current contracts.
Exercise a mixed exterior city through Python access, strict CityJSON, native files,
canonical packages and an installed wheel. Use focused failure checks and a
representative performance check. This milestone does not publish a release.

## Explicit exclusions

No full CityGML conformance or XML adapter; interiors, installations, bridges,
tunnels, OtherConstruction and CityObjectGroup; appearance/textures/instancing;
MultiSolid/CompositeSolid; advanced routing/lane topology, occupancy or temporal
models; arbitrary supplier code conversion, full address model or remote code-list
resolution. No generated Python feature hierarchy, new production dependency,
service, wire replacement, sibling-repository edits or publication. Unrelated dirty
work and historical schema snapshots must remain intact.

## Work and ownership

1. Coordinator: one standard schema revision and final identifiers; integrate schema,
   authoritative coverage/usage documentation, end-to-end evidence and final review.
2. Adapter subagent: bounded transportation/water/cover strict CityJSON feature,
   geometry and semantic-region mappings; focused adapter tests. Coordinate names
   and contracts before changing shared files.
3. Terrain subagent: native DEM interpretation and a truthful terrain-raster builder
   workflow, focused tests, source/schema requirements reported to coordinator.
4. Independent assessment subagent: verify primary standard definitions and audit
   scope/acceptance risks; then review integrated implementation. No competing plan.

## Checkpoints

- [x] Confirm exact feature/geometry/metadata and raster contracts from code and
      primary specifications; fold material decisions into this plan.
- [x] Implement schema, adapter and terrain workflows within the boundary above.
- [x] Demonstrate a mixed exterior city and native raster through public access and
      all claimed persistence boundaries; unsupported data fails without loss.
- [x] Run focused then affected regression checks, representative timing and
      installed-wheel smoke; resolve independent review findings.
- [x] Finalize one current coverage matrix with implemented/restricted/deferred
      status, reproducible examples, limitations and completion evidence.

## Verification and evidence

Keep tests proportional: one representative mixed workflow, relevant shape/semantic
failures and a raster ground-filter/nodata/CRS regression. Reuse existing wire,
schema and package tests. Native arrays and metadata must round-trip exactly;
CityJSON coordinates use its explicit quantization tolerance. Schema bypass cannot
bypass independent format or numerical integrity boundaries. Measure a representative
mixed model without concurrent heavy jobs. Installed smoke reuses dependencies;
no unrun browser/other-language/platform certification is claimed.

## Decisions and completion

Implemented decisions:

- Standard schema 0.9.0, logical namespace `https://github.com/dtcc-platform/dtcc-core/schemas/model#`; previous 0.8.0 archived unchanged. No runtime URI aliases or legacy migration.
- Six new generic features; direct top-level containment only. Transportation admits lines and polygon aggregates; WaterBody also Solid, PlantCover polygon aggregates/Solid. Waterway semantic regions conservatively excluded following CityJSON 2.0.2 section 3.3.
- CompositeSurface is preserved by a representation role on MultiSurface; TIN retains its triangular Mesh mapping. No connectedness/watertightness certification.
- Feature classifications reuse qualified scalar class and list function/usage. Region codes are scalar, water_level categorical; PlantCover average_height uses metres. Existing LandUse/Furniture scalar metadata remains generic.
- Same-theme source region grouping is preserved and schema-checked without claiming geometric hosting.
- Terrain elevation_rasters references a local Raster by geometry_id. A narrow schema annotation uses existing native bindings and semantic projection to validate this reference; no second metadata validator or numerical copying.
- DEM construction records explicit units/reference and cell-center IDW, preserves valid zero/negative samples and CRS, requires actual class 2 for ground-only, rejects ambiguous local transforms, and makes gap filling explicit. Qualified wrapper defaults to no window/gap filling; lower-level defaults retained.

The milestone is complete when the named exterior themes
and the native DEM workflow work as documented, not when every CityGML concept has
been added. Further themes require a separate user-visible scope decision.

### Completion evidence

- Implementation is in this repository, with three collaborating subagents and independent reviews of adapter/schema and terrain paths. Review found and resolved missing semantic host rules and unsupported annotation owners; the latter now fail schema construction rather than leaving a rule unevaluated.
- Current authority: `docs/design/exterior-city-profile.md`; schema `dtcc_core/schemas/dtcc.yaml`; DEM contract `docs/design/terrain-dem.md`. Wire remains v6; no new production dependency or native feature hierarchy.
- Exterior example: all 12 external types, 13 objects, 17 representations, 18,693 native bytes; native/package/context exact, strict CityJSON geometry/roles/regions/attributes verified, default-invalid and bypass-independent failures preserve prior files. Source/export passed the cached official CityJSON 2.0.2 JSON schema.
- Native DEM example preserves zero/negative samples, 14 missing cells, CRS and explicit interpretation across native/package boundaries; invalid metadata is rejected without overwriting. Owner-local representation tests cover deleted/replaced targets, bypass, schema-only annotation configuration and source region-host preservation.
- Final affected regression: `python -m pytest tests/model tests/io tests/builder/test_terrain_dem.py tests/builder/test_flatten_terrain_raster.py tests/builder/test_terrain_meshing.py tests/builder/test_semantic_meshing.py tests/reproject/test_reproject_object.py -q`: **717 passed, 1 skipped**, 33.46s. After final annotation-owner tightening, the 3 focused metadata tests passed (5.32s). Focused schema-only editing example also passed. No broad test failures remain.
- Real pinned 3DBAG tile: 1,110 buildings, 5,550 elevation records, exact native/package, strict coordinate error <=0.000375m, missing-reference save rejection. Run `python -m sandbox.model_profiles.three_d_bag_example SOURCE OUTPUT`; current artifacts/report `/private/tmp/dtcc-exterior-3dbag/`.
- Sequential warm timings: synthetic 1,300-feature/1,700-representation model validation0.063s, validated encode0.166s/decode0.231s; bypass encode0.102s/decode0.170s. Real3DBAG enriched semantics0.687s, encode1.977s/decode3.064s. DEM250,000points/251,001cells rawbackend0.0054s vs hardened0.0151s including support pass. Timings exclude disk/cold schema and are machine-specific. Scripts/reports `/private/tmp/dtcc-exterior-performance.py`, `/private/tmp/dtcc-exterior-performance.json`, `/private/tmp/dtcc-dem-performance.py`, `/private/tmp/dtcc-exterior-dem-performance.json`.
- Final wheel built successfully to `/private/tmp/dtcc-exterior-dist/`, installed without dependency changes in `/private/tmp/dtcc-standard-schema-env/`. Both exterior and DEM examples passed from `/private/tmp/dtcc-exterior-wheel-smoke/` outside the source checkout; artifacts `/private/tmp/dtcc-exterior-wheel-final/`. Installed final backend/projection/schema bytes were verified equal to source. No other-platform/language or release certification is claimed.
- Archived schema0.8 is byte-exact with the preceding wheel (SHA256 `e5a4ac09195d962bdad9338fbfb7e8a23f9291fbc7cf7e4f1594bd6c9b4e3854`). No old URI migration, publication, commit or unrelated cleanup performed.

The finite milestone is closed. Wider CityGML concepts and producer-specific
adapters remain explicit exclusions rather than implied unfinished work.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-exterior-city-profile.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Complete the bounded exterior-city and
native DEM workflows using generic carriers and one schema authority. Coordinate
shared files with the other agents and do not expand into the explicit exclusions.
