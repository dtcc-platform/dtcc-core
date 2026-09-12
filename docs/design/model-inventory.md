# Native model, wire and semantic inventory

For the current strict CityJSON feature and geometry coverage, see the
[exterior-city profile](exterior-city-profile.md). This inventory separately
tracks native carriers, which support more numerical state than CityJSON exchange.

Status: 12 September 2026, after the exterior-city milestone; wire **6**,
standard schema **0.9.0**.
This is the current inventory; the earlier phase-1 serialization tables are
historical. It covers all 36 exported Model subclasses: 27 supported concrete
roots, eight explicitly unsupported result/container types and abstract Geometry.
Model itself and module-visible Polygon are noted separately below.

The wire authority is [dtcc.proto](../../dtcc_core/proto/dtcc.proto), generated as
[dtcc_pb2.py](../../dtcc_core/model/dtcc_pb2.py). The default semantic contract is
[model/0.9.0/schema.yaml](../../dtcc_core/schemas/model/0.9.0/schema.yaml).
Earlier standard schemas are retained unchanged under `schemas/archive/model` as
historical definitions, outside the runtime bundle. Default evaluation requires an
available schema version; it does not migrate old property names. Every admitted
concrete native type has an explicit 0.9.0 schema binding.
A binding to a generic class is not validation of every possible domain meaning.
The [terrain/transportation crosswalk](terrain-transportation-crosswalk.md) audits
their properties and distinguishes historical gaps from the subsequent implementation.
The [exterior-city profile](exterior-city-profile.md) now admits 12 strict CityJSON
feature types, including triangular TINs, physical transportation, water and vegetation
cover. The [native DEM workflow](terrain-dem.md) adds explicit raster interpretation
without claiming a CityJSON raster conversion. The current namespace is
`https://github.com/dtcc-platform/dtcc-core/schemas/model#`.

## Objects

All supported Objects preserve concrete type, identity, attributes, children,
representations, transforms, semantic/profile labels and named references.
Derived bounds caches are reconstructed, not persisted; root schema selection is
retained. Public mutation still requires explicit refresh for an up-to-date cached
view. The [spatial contract](model-spatial-contract.md) describes local-coordinate
bounds and the remaining global transform-composition limits.

| Native class | Additional authoritative state in wire | Schema class / semantic coverage | CityGML / CityJSON assessment |
| --- | --- | --- | --- |
| Object | None beyond shared state | Object; declared URI can select a known semantic class, unfamiliar URIs stay generic | Carrier for the strict individual-vegetation/furniture and six new exterior feature classes; no automatic equivalence for arbitrary URIs |
| City | Shared state | City / generic containment | Aggregate; CityJSON places the aggregate at document level |
| CityObject | Shared state | CityObject / generic Object rules | No automatic equivalence to every CityGML CityObject restriction |
| Building | Shared state | Building under AbstractBuilding; [optional attributes](building-attributes.md) and containment | Bounded mapping implemented; [property crosswalk](building-semantic-crosswalk.md) identifies the remaining rules; [height authority](building-height.md) is implemented |
| BuildingPart | Shared state | BuildingPart; optional owner must be Building/BuildingPart | Bounded part mapping implemented; standalone native roots allowed |
| Terrain | Shared state, including Mesh/Raster representations | Terrain with optional elevation_rasters interpretation records; explicit TINRelief specialization | Strict triangular TINRelief ↔ Mesh and native DEM workflow implemented; other relief conversions remain unmapped |
| Tree | Position array, height, crown_radius | Tree / individual-vegetation rules | Native Tree state is narrower than SolitaryVegetationObject and has no strict CityJSON mapping; no automatic conversion from a generic plant |
| Landuse | Ordered native LanduseClasses names | LandUse; native codes distinct from source class/function/usage | Bounded LandUse mapping implemented |
| RoadNetwork | Vertices, edges and length arrays | RoadNetwork / generic Object rules | Routing graph remains separate from physical Road/Railway/Waterway/TransportSquare, now mapped through generic Object; no inferred graph conversion |
| SensorCollection | Shared state; child stations retain concrete type | SensorCollection / generic Object rules | Observation metadata not covered by the current city vocabulary |
| VehicleCollection | Shared state; child vehicles retain concrete type | VehicleCollection / generic Object rules | Live snapshot; temporal/moving-object semantics remain open |
| DeSO | Shared state; area attributes and explicitly associated fields | DeSO / generic Object rules | Statistical area collection, not automatically LandUse |
| BuildingCollection | **Unsupported**: buildings list | No binding | Decide whether this is a result wrapper or durable model root |
| FootprintCollection | **Unsupported**: surfaces, source IDs/indices | No binding | Preserve source identity if admitting it; footprints are representations, not a new city class |
| TreeCollection | **Unsupported**: trees list | No binding | Result-wrapper decision remains |
| CalibrationGrid | **Unsupported**: bounds/divisions/CRS/features/metadata | No binding | Computational result; not a thematic CityGML feature by default |

## Geometry, values and support

Concrete Geometry roots preserve transforms, attached fields and supported semantic
regions. Arrays preserve dtype, shape and values with little-endian row-major bytes.
Raster/Bounds can also be representation values. Grid bounds define the numerical
domain, including explicitly zero-area domains. The
[spatial contract](model-spatial-contract.md) records refresh behavior and raw-envelope
aggregation limits. Other geometry bounds are
derived caches.

| Native class | Authoritative state / wire coverage | Semantic coverage and standards assessment |
| --- | --- | --- |
| Geometry | Abstract base | Abstract Geometry vocabulary; no standalone wire value |
| Point | XYZ doubles | Native Point; bounded singleton CityJSON MultiPoint conversion |
| LineString | Vertex array, 2D or 3D | Generic Geometry; strict 3D line member within supported MultiLineString; a standalone LineString representation has no automatic external mapping |
| MultiLineString | Ordered LineString children with their metadata | Strict CityJSON transportation/WaterBody line mapping implemented; requires nonempty 3D lines, rejects fields and other unmapped metadata |
| Surface | Exterior vertices, holes, normal | Native Surface; polygon mapping in implemented city slice |
| MultiSurface | Ordered Surface children | Bounded strict polygon mapping; role `cityjson:CompositeSurface` preserves a source CompositeSurface declaration without connectivity certification |
| Solid | Surface children and indexed shells | Native Solid; strict building/part, furniture, individual vegetation, WaterBody and PlantCover mapping; `.mesh()` triangulates all boundary surfaces with regions, retaining shell partition on the source; no watertightness claim |
| Mesh | Vertices, faces, markers, normals | Native Mesh; numerical connectivity checks; generic mesh is not automatically TINRelief |
| VolumeMesh | Vertices, tetrahedral cells, markers | Generic Geometry; simulation mesh is not a CityGML Solid shell |
| PointCloud | Points, classification, intensity, return_number, num_returns | Generic Geometry; classification numbers do not establish city-feature semantics |
| Grid | Width, height, domain bounds | Generic Geometry; numerical discretization; explicit domain retained independently of cell resolution |
| VolumeGrid | Width, height, depth, domain bounds | Generic Geometry; numerical discretization; explicit domain retained independently of cell resolution |
| Field | Values, name, unit, description, dim, association | Field metadata projected to schema; arrays and counts checked natively; general time/scenario axes missing |
| Raster | Data array, six affine coefficients, nodata, CRS | Generic Entity; Terrain can associate schema-declared elevation interpretation by representation ID. Native DEM support is not a CityGML RasterRelief or CityJSON conversion |
| Bounds | Six finite ordered bounds | Generic Entity; standalone bounds or explicit representation, not cached Object extent |
| Transform | SRS and affine 4x4 matrix | Generic Entity; local/global transform composition remains a separate model concern |
| FieldSlice | **Unsupported**: PointCloud state plus slice/time/domain/metadata | No binding; must preserve all slice context before admission |
| StreamlineCollection | **Unsupported**: lines, seeds/integration parameters/time/domain/metadata | No binding; do not serialize only its lines |
| DatasetCollection | **Unsupported**: arbitrary items | No binding; define durable collection meaning first |
| DatasetValue | **Unsupported**: arbitrary Python value | No binding; do not invent a universal Python-object serializer |

The base Model provides the shared serialization methods but is not a supported
root itself. Polygon is module-visible, not exported, and remains abstract because
it has no calculate_bounds implementation. GeometryRepresentation and SemanticRegion
are contained dataclasses, not standalone Model roots. Representation ID/LoD/role
are wire metadata. SemanticRegion has its own schema binding, known building,
transportation and water boundary rules, membership indices and optional local host index. Strict external owner
vocabularies and scalar semantic-attribute restrictions are documented in the
exterior profile. GeometryType, RoadType and LanduseClasses are enums, not
standalone models.

Field associations are explicit: vertex, edge, face, cell, sample or geometry.
The last means one value for the entire geometry, used for a DeSO area statistic.
PointCloud sample fields have one row per point; MultiLineString sample fields have
one row per line. Other associations use the corresponding supported owner element
counts. A standalone Field has no owner count check. None remains valid while
editing but cannot be serialized; no association is inferred from array length.

## What the standards comparison establishes

There has **not** been a complete property-by-property CityGML conformance audit.
The [building semantic crosswalk](building-semantic-crosswalk.md) now reviews
building/part attributes, relevant inherited construction facts and boundary/opening
relationships against the conceptual model, published XSDs and actual DTCC behavior.
It records source multiplicity discrepancies, units, deliberate restrictions and
recommendations initially assessed against schema 0.2.0. Naming, qualified codes,
height/elevation records and the explicit 3DBAG mapping were subsequently implemented.
The current standard is 0.9.0; historical recommendations are not all outstanding work.
Neither the crosswalk nor later milestones claim full CityGML conformance.
Earlier work inspected the building/part, boundary/opening, LandUse, individual
vegetation and furniture slice, implemented selected rules and tested a strict
CityJSON subset. The present inventory completes the native class classification,
not a complete standards conformance audit. Generic schema bindings introduced for numerical
and collection types make their admission explicit without claiming richer semantics.

CityJSON explicitly implements a subset of the CityGML 3.0 conceptual model. Its
feature/geometry restrictions are useful adapter requirements; its flexible
attributes do not establish that DTCC implements every CityGML property.
[CityJSON 2.0.2, sections 2, 3 and 10](https://www.cityjson.org/specs/2.0.2/).
CityGML separates conceptual modeling from encodings, so DTCC can use the concepts
with Protobuf while retaining its native computational types.
[OGC CityGML overview](https://www.ogc.org/standards/citygml/).

| Standards area | Current assessment | Next decision |
| --- | --- | --- |
| Building/part/boundary/opening | Implemented bounded slice; property/relationship crosswalk documented | Height authority, basic exterior attributes, eleven surface types, [strict schema admission](cityjson-schema-io.md) and [qualified height/code records](qualified-values.md) are implemented. Other qualified records, rooms/installations and addresses remain deferred. |
| Land use/individual vegetation/furniture | Implemented bounded slice | Extend vocabulary and mappings only for concrete DTCC workflows |
| Relief, transportation, water, vegetation cover | The exterior-city profile implements physical transportation, WaterBody, PlantCover and triangular TIN exchange; native DEM interpretation is explicit | Keep documented exclusions: Waterway regions, transport Solid, additional relief conversions and richer CityGML relationships require separate scope |
| Bridges, tunnels, other constructions, city-object groups | Explicitly outside the first exterior-city profile | Revisit only for an agreed workflow and separately bounded mapping |
| Appearance, templates/instancing | No canonical representation | Assess actual visualization/storage requirements |
| Time, observations, scenarios and provenance | Snapshot fields/context exist; full relationships and axes missing | Define computational contract first; assess relevant standards modules separately |

## Acceptance evidence and separately scoped follow-ups

The [issue-closeout assessment](model-issue-85-closeout.md) is the authority for
current #85 acceptance and its exact verification. This inventory defines the
supported native boundary; its follow-ups do not imply that every wrapper,
CityGML theme or external adapter must be implemented in this milestone.

1. Python → browser JavaScript → Python has a
   [reproducible browser example](../../sandbox/protobuf_interop/README.md) using
   the public dtcc.proto. Its versioned evidence covers IDs, footprint access,
   exact uint64 values, float64 coordinates, fields, regions, unknown-type policy
   and browser edits, with invalid semantics/arrays and preserved unknown wire
   fields rejected by Python. The schema 0.9.0 / wire 6 run passed in Chromium 152;
   the closeout assessment records its exact evidence. Production TypeScript/downstream integration and C++ evidence
   at actual consumer boundaries remain separate; this is not a general browser validator.
2. The [exterior-city profile](exterior-city-profile.md) establishes a finite
   standards boundary, using the current repository-owned semantic namespace.
   Additional themes require an explicit workflow and acceptance scope. New
   semantic subtypes can use YAML where the generic native carriers suffice;
   no CityGML-shaped Python hierarchy is required.
3. The eight unsupported wrappers/results fail explicitly. A future durable
   workflow must decide whether to unwrap an existing model, retain a convenience
   result, or admit additional authoritative state. Downcasting merely to claim
   serialization coverage is not acceptable.
4. Canonical v3 packages work through `export(canonical=True)`; default artifact-only
   v2 packages, downstream consumers and other external I/O adapters require
   separately bounded adoption work. Native and strict CityJSON default validation
   does not imply that every external adapter now evaluates the standard schema.
   External service `format="pb"` selectors emit ModelFile bytes; no legacy native
   model codec remains.
5. [Local coordinates and explicit bounds refresh](model-spatial-contract.md) are
   the current spatial contract. Global transform composition, general temporal/
   scenario axes and per-entity provenance associations need their own concrete
   workflow before extending the contract. Further performance work must follow
   measured consumer needs rather than weakening admission for presumed speed.

Executable native coverage is in
[test_serialization_matrix.py](../../tests/model/test_serialization_matrix.py).
The unified-format workflow and boundary checks are in
[test_unified_protobuf.py](../../tests/io/test_unified_protobuf.py).
The [unified-format plan](../../.agent/plans/2026-09-11-unified-protobuf.md) records
that historical milestone. Current exterior evidence is in the
[exterior-city plan](../../.agent/plans/2026-09-12-exterior-city-profile.md); current
issue acceptance and follow-ups are recorded in the closeout assessment linked above.
