# Terrain and transportation: CityGML, CityJSON and DTCC

Current implementation: the [exterior-city profile](exterior-city-profile.md)
and [native DEM workflow](terrain-dem.md), using schema 0.9.0, complete the
bounded transportation, water, cover, TIN and native raster work selected from
this comparison. The tables below retain the 0.7.0 audit and are historical gap
analysis, not the current supported-format matrix.

Reviewed 12 September 2026 against standard DTCC schema 0.7.0 and Protobuf v6.
This extends the [building crosswalk](building-semantic-crosswalk.md). It is a
design comparison with concrete implementation gaps, not a claim of full standards
conformance or an implemented terrain/transport adapter.

The important distinction is between a computational carrier and a domain claim.
DTCC already preserves Terrain, Mesh, Raster, PointCloud, MultiLineString and
RoadNetwork natively. Their numerical admission does not establish relief meaning,
road classification, traffic direction, or a faithful CityJSON conversion.

## Authorities inspected

- [CityGML 3.0 conceptual model, OGC 20-010](https://docs.ogc.org/is/20-010/20-010.html):
  Relief §8.10, tables 365–383; Transportation §8.11, tables 390–430 and enums.
- Published [Relief 3.0 XSD](https://schemas.opengis.net/citygml/relief/3.0/relief.xsd)
  and [Transportation 3.0 XSD](https://schemas.opengis.net/citygml/transportation/3.0/transportation.xsd):
  encoded property cardinalities and enum values, inspected directly as XML.
- [CityJSON 2.0.2](https://www.cityjson.org/specs/2.0.2/): §§2.11, 2.12, 3.3,
  plus the linked machine schemas. Its payload version is still `2.0`.
- DTCC [native inventory](model-inventory.md), [standard schema](https://github.com/dtcc-platform/dtcc-core/blob/349cdfa8ae3eabceb65a7131aa6a1e8c847c0e46/schemas/archive/model/0.7.0/schema.yaml),
  [exchange admission](../../dtcc_core/model/exchange.py),
  [strict CityJSON admission](../../dtcc_core/io/cityjson/admission.py),
  [terrain builder](../../dtcc_core/builder/geometry_builders/terrain.py),
  [road vector importer](../../dtcc_core/io/roadnetwork.py).

The conceptual association tables and published XSD must not be conflated. For
example, conceptual tables 392/414/426 display `1..*` on several transport roles,
whereas the corresponding XSD property elements use `minOccurs="0"`. We record
the encoding's optional repeated roles below; these are not automatically proposed
as mandatory DTCC properties. An eventual GML adapter needs an explicit conformance
target rather than mechanically copying one table into LinkML.

The inspected CityJSON 2.0.2 machine schema admits MultiLineString, MultiSurface
and CompositeSurface for transportation geometry. It is broader than some
prose rules: semantic-object extra properties are unrestricted, while §3.3 describes
scalar values and the transport example includes an array. This discrepancy needs
an explicit adapter policy; passing the JSON Schema alone is not proof of every
prose constraint. DTCC's current strict semantic metadata subset is scalar-only.

## Relief property crosswalk

| CityGML fact | CityJSON representation | Existing DTCC carrier | Decision / missing contract |
| --- | --- | --- | --- |
| `ReliefFeature.reliefComponent`, one or more; aggregate `lod`, integer 0–3 | No matching aggregate wrapper in the bounded TIN feature form | Terrain/Object with child Objects | Keep a single terrain tile simple. Introduce component Objects only when they need independent IDs, extent, LoD or source lineage. Do not equate every Terrain with this aggregate. |
| `AbstractReliefComponent.lod`, required integer 0–3 | Geometry `lod` string | GeometryRepresentation.lod | Preserve exact source LoD; do not deduce resolution or accuracy from it. A GML component LoD and representation LoD need an explicit adapter rule. |
| Component `extent`, optional 2D polygon surface | Feature `geographicalExtent` is an axis-aligned numeric box | Surface representation versus derived Bounds | The polygon domain is not interchangeable with Bounds. Preserve it explicitly if admitted; never replace it with a bounding box. |
| `TINRelief.tin`, one triangulated surface | `TINRelief` + `CompositeSurface`; text also permits nontriangular polygons | Mesh for triangles; MultiSurface for polygons | Strict conversion is absent. A triangular subset can use Mesh with no new numerical class. Nontriangular input needs explicit polygon preservation or rejection, not hidden triangulation. |
| `RasterRelief.grid`, one grid coverage; XSD uses RectifiedGridCoverage | No standard RasterRelief feature in this CityJSON subset | Raster(data, georef, nodata, crs), optionally on Terrain | Native arrays and affine survive. Missing elevation-band meaning, unit, vertical reference, cell/point sampling convention and complete coverage metadata. No automatic CityJSON raster-to-TIN conversion. |
| `MassPointRelief.reliefPoints` and `pointCloud`, each optional | No corresponding standard relief feature | PointCloud with numerical classifications and fields | A generic cloud is not certified bare earth. Preserve source point classification and define the terrain-selection method before classifying it as relief. External point-cloud references need a separate adapter. |
| `BreaklineRelief.breaklines` and `ridgeOrValleyLines`, each optional | No corresponding standard relief feature | MultiLineString representations with distinct roles | Keep the two roles distinct; a line collection alone does not identify a breakline or whether it is a ridge/valley. Strict external mapping is absent. |

The new building `ElevationMeasurement` establishes how a scalar carries a vertical
reference. It does not add a unit, band definition or sampling model to Raster, nor
does it establish the interpretation of every geometry Z coordinate. Raster elevation
semantics should be designed around a real DEM workflow, not thousands of per-cell
measurement dictionaries.

Native raster admission preserves dtype/shape, six affine coefficients, nodata and
CRS. It does not interpret pixel values as elevations. The current terrain raster
builder selects classifications 2/9 when class 2 exists, otherwise uses all points;
sets nodata to zero and fills holes; and does not assign source CRS in that function.
These are concrete issues to address before certifying a DEM-building workflow:
zero can be valid terrain elevation, interpolation changes provenance, and a
requested ground-only operation can currently include non-ground points. They are
recorded here, not silently changed as part of this semantic mapping milestone.

The permissive CityJSON reader's `get_terrain_mesh` returns only the first of
multiple TIN objects. The permissive exporter has a terrain path, but neither is
the strict preservation contract. A future strict path must keep every terrain
feature and representation or fail explicitly. Adding `TINRelief` to YAML alone
cannot resolve these conversion issues. No implicit `GroundSurface` labels should
be manufactured for relief triangles from the building surface vocabulary.

## Transportation property crosswalk

| CityGML fact | CityJSON difference | DTCC decision / gap |
| --- | --- | --- |
| Road, Railway, Track, Square, Waterway specialize AbstractTransportationSpace | Road, Railway, Waterway, TransportSquare; Section/Intersection/Track are flattened rather than corresponding feature types | Use generic Objects with schema-declared semantic types for physical features. Map Square ↔ TransportSquare explicitly at the boundary. RoadNetwork stays a separate computational graph. |
| Optional `class`; repeated `function` and `usage` on Road and other main features | Feature attributes allow arbitrary JSON values | Reuse qualified codes and distinct intended/actual use when declaring these schema types. Do not infer function from a native RoadType enum or coerce arbitrary supplier values automatically. |
| Optional repeated Section and Intersection links on main transport features in XSD | Sections can be modeled through attributes and feature hierarchy | DTCC children express single ownership; relations reference IDs. Shared intersections require relations, not duplicated children. No current strict hierarchy mapping. |
| Optional repeated TrafficSpace, AuxiliaryTrafficSpace, Hole and Marking roles in XSD | Flattened representation differs from full spatial-object decomposition | Generic Objects can retain independent IDs when needed. Surface regions serve surface-only source facts; do not pretend a region is a whole traffic volume. |
| TrafficSpace and AuxiliaryTrafficSpace require `granularity`: `lane` or `way` | No automatic equivalent supplied by the flattened form | Add a schema enum only with a workflow that knows this fact. Neither graph edge count nor LoD determines lane/way granularity. |
| `trafficDirection`: optional `forwards`, `backwards`, `both`, relative to linear geometry | Supplier attributes require explicit interpretation | Proposed DTCC `traffic_direction` must identify the directed representation. Never interpret edge endpoint order as legal travel direction. |
| TrafficSpace `predecessor` and `successor`, optional repeated references in XSD | No full traffic-space reference structure supplied by default | Object relations can carry declared feature IDs. Native graph vertex indices refer to vertices, not TrafficSpace IDs, and are not an equivalent relationship. |
| TrafficArea / AuxiliaryTrafficArea are thematic ground surfaces; optional class/material and repeated function/usage | TrafficArea/AuxiliaryTrafficArea semantic primitives; Marking/Hole become TransportationMarking/TransportationHole | SemanticRegion is the right native carrier for surface membership. Add vocabulary and scalar naming rules at the adapter boundary. Full feature identity/provenance may require Objects instead. |
| `surfaceMaterial`, optional gml:CodeType on thematic areas | Examples and semantic-attribute rules differ on arrays | Proposed `surface_material` needs a declared scalar/qualified-code policy. Current strict semantic adapter accepts scalar metadata only; nested QualifiedCode cannot simply be copied into a standard CityJSON Semantic Object. |
| TrafficSpace clearanceSpace links, optional repeated, to ClearanceSpace objects | Example `clearanceSpace` attribute is a scalar | A scalar clearance is not the full geometry of a clearance space. Keep a qualified distance separate from a spatial Object; no inference between them. |
| Occupancy records on transportation/traffic spaces | Arbitrary attributes do not establish the full record semantics | Temporal population/vehicle occupancy remains outside current schema. Do not derive it from VehicleCollection or a snapshot count. |

The current strict adapter admits none of the transport feature types, no line
collections, and only building semantic surface names. `Road` with an unfamiliar
URI can be retained as generic native data, but that does not mean its transport
meaning was validated or that strict CityJSON can export it. The vocabulary,
record rules and adapter mappings are three separate work items.

## RoadNetwork audit and ordinary access

[RoadNetwork](../../dtcc_core/model/object/roadnetwork.py) exposes `vertices`,
`edges`, `length` and optional `multilinestrings`. Native `.dtcc` preserves the
graph arrays and checks endpoint indices and finite nonnegative edge-aligned
lengths. No graph geometry, road width, traffic surface or legal direction is
created by the semantic schema.

The vector importer currently builds one edge per source line from its endpoints,
optionally rounds coordinates (two decimals by default), retains full line geometry
when requested, and stores source properties as lists in attributes. Inspection
shows `id_field` is accepted but does not create independently identified edge
Objects; list alignment and stable source-ID uniqueness are not enforced by the
generic Object schema. Graph vertices come from a set, so their index ordering is
not a stable external identity. Geometric crossings are not automatically junctions;
grade separation and source direction are not resolved. Length comes from the
source line in its coordinate units, not necessarily metres in a geographic CRS.

These are boundaries to settle for a routing workflow, not reasons to replace
the efficient arrays with an IFC-like object graph. Keep graph access direct.
If a user needs a particular physical road's polygon, store an explicitly identified
Road Object with a footprint/surface representation and a declared connection to
the graph where supplied. Do not create that polygon from the graph silently.

## Recommended next implementation slice

Start with **strict triangular CityJSON TINRelief import/export**:

1. Declare TINRelief as an explicit semantic specialization compatible with native
   Terrain; retain generic Terrain without asserting it is always a TIN or DTM.
2. Map triangular CompositeSurface to existing Mesh, preserving each feature ID,
   attributes, triangle order, coordinates, CRS and exact representation LoD.
   Preserve multiple terrain objects and representations. Reject polygon holes,
   nontriangles and unhandled semantics in this first subset.
3. Use ordinary `terrain.get_geometry(id=...)` for direct NumPy vertices/faces and
   exercise default `.dtcc` and canonical package round trips plus strict CityJSON.
   Match numerical comparison to CityJSON quantization; keep native arrays exact.
4. Verify a mixed building/terrain example and a meaningful rejected input without
   dropping data or overwriting a valid artifact. Adding raster conversion, new
   geometry base classes or a general GML runtime is outside this slice.

Then choose one explicit transport workflow, preferably a source Road surface with
TrafficArea regions and qualified feature classifications, before extending graph
semantics. A later DEM workflow must address source filtering, nodata, CRS, vertical
units and sampling together. Schema-only names cannot substitute for these facts.

This ordering preserves generic, fast numerical carriers and makes each semantic
extension earn its place through a concrete user workflow. It does not commit DTCC
to the entire CityGML hierarchy.
