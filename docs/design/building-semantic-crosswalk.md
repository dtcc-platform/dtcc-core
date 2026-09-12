# Building semantics: CityGML, CityJSON and DTCC

Current implementation: standard **0.9.0** retains snake_case DTCC
properties and explicit CityJSON mappings. See [the naming contract](model-naming.md).
The crosswalk and behavioral evidence below retain their original **0.2.0**
spellings and observed behavior at the time of the audit. The subsequent
[height milestone](building-height.md) resolves scalar height access and separates
modelling estimates. The [basic building milestone](building-attributes.md) adds
attributes, sibling schema classes and the six missing surface declarations.
[Strict CityJSON schema admission](cityjson-schema-io.md) is also implemented.
The first [qualified height/code records](qualified-values.md) are implemented.
Signed [elevations and an explicit 3DBAG attribute mapping](3dbag-attribute-mapping.md)
are implemented. The [terrain and transportation crosswalk](terrain-transportation-crosswalk.md)
continues the property audit. Other qualified records and adapters remain open.

Status: reviewed 11 September 2026 against DTCC Model schema **0.2.0** and wire
version **6**. This is the building/part, inherited construction, and boundary/opening
crosswalk. It is not a complete CityGML audit or a new executable schema.
Recommendations below are proposals for subsequent implementation; current behavior
is stated separately. No native API or published schema version changes in this step.

The existing generic model and Protobuf value representation can carry most of the
missing building metadata. The main work is defining meaning, units and access,
then declaring the relevant rules in a new schema version. Generating CityGML's
Python class hierarchy would not solve those problems.

## Sources and how to read the comparison

Primary sources inspected:

- [CityGML 3.0 Conceptual Model, OGC 20-010](https://docs.ogc.org/is/20-010/20-010.html):
  Core feature/space properties, Construction, Building and their referenced data
  types. The full HTML was downloaded and the relevant sections read locally.
- OGC's encoding schemas, each declaring version 3.0.0:
  [Building](https://schemas.opengis.net/citygml/building/3.0/building.xsd),
  [Construction](https://schemas.opengis.net/citygml/construction/3.0/construction.xsd),
  [Core](https://schemas.opengis.net/citygml/3.0/core.xsd).
  These cross-check encoding cardinalities; they do not replace the conceptual model.
- [GML basic types](https://schemas.opengis.net/gml/3.2.1/basicTypes.xsd), particularly
  CodeType, MeasureType and MeasureOrNilReasonListType. The served file declares
  version 3.2.2 despite its 3.2.1 directory.
- [CityJSON 2.0.2](https://www.cityjson.org/specs/2.0.2/), especially §§2.1, 2.3,
  3.3, 4, 5 and 10. Its [implementation notes](https://www.cityjson.org/citygml/v30/)
  provide context; use the versioned specification and OGC definitions for precise
  property claims. For example, the notes mention measuredHeight while the reviewed
  CityGML 3.0 Construction definition has qualified height records.

`0..1` means optional scalar, `0..*` optional collection, and `1` required.
Source cardinality and the proposed DTCC rule are different columns deliberately.
Code lists are not universal enums: a code and its code-list context can matter.
Optional means absence is allowed; it does not mean zero, empty text or JSON null
has the same meaning as absence.

One source discrepancy must not become a DTCC requirement: conceptual catalogue
tables 668 and 674 print `1..*` for associations such as address and buildingPart,
whereas the published Building XSD declares `minOccurs="0"`, `maxOccurs="unbounded"`.
Similar differences occur for boundary and fillingSurface. The relationship table
below reports this explicitly. We do not infer that every building must have an
address, part, room or window. Resolving formal CityGML conformance would require
reconciling the conceptual UML/catalogue and encoding; this milestone makes no such
conformance claim.

Current implementation authorities are the
[standard schema](../../schemas/archive/model/0.2.0/schema.yaml),
[native Building API](../../dtcc_core/model/object/building.py),
[semantic projection](../../dtcc_core/model/profiles.py),
[LinkML backend](../../dtcc_core/model/_profile_backend.py),
[strict CityJSON admission](../../dtcc_core/io/cityjson/admission.py), and
[surface mapping](../../dtcc_core/io/cityjson/semantics.py).

## Current validation boundary

The induced Building slots in schema 0.2.0 are exactly `id`, `parent`, `name`,
`measuredHeight` and `roofType`; BuildingPart inherits those semantic rules.
Most other building attributes are preserved by the native typed-value codec but
are not projected into LinkML evaluation. Their JSON-like structure and finiteness
are checked; their domain meaning, units, date validity and cardinality are not.

CityJSON leaves feature attributes open to arbitrary JSON values. Consequently,
a CityJSON file can be structurally valid while disagreeing with DTCC's selected
building schema. Current strict CityJSON import/export invokes native admission,
not standard LinkML evaluation. Canonical `.dtcc` save/load applies both by default.
This remains an adapter-adoption gap, not evidence that arbitrary attributes have
already been semantically validated. [CityJSON §2.1](https://www.cityjson.org/specs/2.0.2/).

## Identity and descriptive properties

Source definitions: [AbstractFeature](https://docs.ogc.org/is/20-010/20-010.html#AbstractFeature-section),
[AbstractCityObject](https://docs.ogc.org/is/20-010/20-010.html#AbstractCityObject-section),
[ExternalReference](https://docs.ogc.org/is/20-010/20-010.html#ExternalReference-section).
All recommendations concern DTCC's own contract, not a universal source conversion.

| Source concept / type / cardinality | Current DTCC storage and checks | Recommended DTCC decision |
| --- | --- | --- |
| featureID: ID, `1`; CityJSON object dictionary key | `Object.id`; native identity/reference checks; `.dtcc` preserves it | Keep one model-local identity authority. Do not require UUID syntax or mistake an ID for a globally persistent real-world identifier. |
| identifier: ScopedName, `0..1` | Can be carried as generic metadata; no scoped-identifier rule | Keep stable external identity distinct from local `id`. Define namespace/scope before declaring a standard field; do not generate a new identity silently. |
| name: GenericName, `0..*` | `attributes['name']` is optional scalar string under the DTCC schema; a list fails | Keep a convenient scalar display name as an intentional DTCC simplification. Preserve additional names/language/scope separately when a workflow needs them; do not claim scalar equivalence to all source names. |
| description: CharacterString, `0..1` | Generic attribute; no Building type rule | Add optional string `description` in YAML. It does not require a new Python class or wire field. |
| Generic attributes | Nested typed values survive; generic validity only for undeclared keys | Keep open metadata. Adding a declared slot in a later schema can make previously accepted data invalid; document that change. |
| externalReference: records with required targetResource URI and optional informationSystem/relationType URIs | Structured attributes can preserve the record; `Object.relations` targets must resolve inside the model | Use structured metadata for external references. Never put external URLs into model-local ID relations and then weaken dangling-reference checks. |

CityJSON metadata.identifier identifies the **dataset**, not every building. Its
metadata.referenceDate likewise must not be interpreted as each building's
construction date. Dataset Context and feature attributes have different owners.
[CityJSON §5](https://www.cityjson.org/specs/2.0.2/).

## Building and construction properties

The building-specific rows follow
[AbstractBuilding, table 669](https://docs.ogc.org/is/20-010/20-010.html#AbstractBuilding-section).
Inherited construction rows follow
[AbstractConstruction, table 538](https://docs.ogc.org/is/20-010/20-010.html#AbstractConstruction-section).
The types/cardinalities were also checked against the Building and Construction
XSDs linked above. CityJSON does not impose these types on arbitrary attributes.

| Source property / type / cardinality / unit | Current DTCC schema 0.2.0 | Recommendation / intentional limitation |
| --- | --- | --- |
| class: BuildingClassValue, `0..1`, code | Generic metadata only | Optional code string in the initial simple subset; preserve code-list context when supplied. No invented universal enum or automatic numeric-code conversion. |
| function: BuildingFunctionValue, `0..*`, codes | Generic; scalar and list can both survive | Declare a list of strings for intended purposes. Source scalar-to-list normalization belongs in an explicit adapter policy, not implicit core coercion. |
| usage: BuildingUsageValue, `0..*`, codes | Generic | Declare a list of strings for actual uses, independently of function. Never derive one from the other. |
| roofType: RoofTypeValue, `0..1`, code | Optional string, no universal code list | Keep current simple rule; document that code-space-qualified values require a separate representation decision. |
| storeysAboveGround: Integer, `0..1`, count | Generic; even a negative value is currently accepted | Declare optional integer ≥0 as a DTCC count rule, excluding booleans. Preserve absence. |
| storeysBelowGround: Integer, `0..1`, count | Generic | Same count rule, independently optional. Do not infer missing below-ground count as zero. |
| storeyHeightsAboveGround: MeasureOrNilReasonList, `0..1`, unit-bearing list | Generic nested/scalar values; no list/units/order rule | Defer until a unit-bearing ordered-list contract exists. Retain missing-entry reasons; do not replace them with zero. |
| storeyHeightsBelowGround: MeasureOrNilReasonList, `0..1`, unit-bearing list | Generic | Preserve order from nearest ground outward, separately from the above-ground list. Do not automatically sum lists into building height or require a matching count without a chosen profile rule. |
| conditionOfConstruction: enumeration, `0..1` | Generic | A later declared slot can use the actual versioned source enum. Keep it distinct from arbitrary condition/quality scores. |
| dateOfConstruction: Date, `0..1` | Generic string; calendar syntax not checked | Define an optional ISO calendar date. Preserve year-only source precision separately; never invent January 1. |
| dateOfDemolition: Date, `0..1` | Generic | Same date policy. Any ordering rule should apply only to comparable declared dates; future/planned facts need their own context. |
| constructionEvent: ConstructionEvent, `0..*` | Generic nested metadata | Deferred structured records; source has event code and date, plus optional description. Not one free-text replacement for all lifecycle facts. |
| height: Height, `0..*`, qualified Length | Generic metadata; no qualified-height schema | Preserve measurement qualifiers and units as described below. Do not feed this record collection into the current scalar `Building.height` accessor. |
| elevation: Elevation, `0..*` | Generic metadata; no elevation record rule | Defer qualified position/reference handling. A height difference and an elevation relative to a vertical reference are different facts. |
| occupancy: Occupancy, `0..*` | Generic metadata | Defer qualified occupancy records; do not reduce them to a single unqualified population integer. |
| measuredHeight: example/source attribute spelling in CityJSON; not the reviewed CityGML 3.0 Height record | Optional numeric scalar ≥0, **metres** by DTCC declaration | Retain as a useful deliberately restricted DTCC measurement. It cannot express several height definitions or their capture status. Never claim automatic CityGML 3.0 Height equivalence. |
| yearOfConstruction / yearOfDemolition: source conventions, distinct from the reviewed Date properties | Generic metadata | Retain original precision and names until explicitly mapped; do not silently rename them to dates. |

Units are part of meaning. GML MeasureType and MeasureOrNilReasonListType require
`uom`; CodeType can carry `codeSpace`. DTCC's present measuredHeight rule instead
fixes the unit in its schema. A number imported from CityJSON does not prove its
author intended metres. Conversion from another known unit must be explicit;
unknown units must not be guessed from the coordinate system.
[GML basic types](https://schemas.opengis.net/gml/3.2.1/basicTypes.xsd).

### Height deserves the first implementation decision

CityGML 3.0 [Height, table 619](https://docs.ogc.org/is/20-010/20-010.html#Height-section)
requires four facts per record: highReference, lowReference, status and a Length
value. The [Construction XSD](https://schemas.opengis.net/citygml/construction/3.0/construction.xsd)
confirms those required members and a collection of height records on constructions.

DTCC currently has three different sources of an apparent height:

| Access | Observed meaning |
| --- | --- |
| `building.attributes['measuredHeight']` | Declared optional nonnegative scalar in metres; checked at canonical persistence |
| `building.attributes['height']` | Open metadata, used by existing builders; not a declared Building measurement in schema 0.2.0 |
| `building.height` | Returns the `height` attribute if present; otherwise computes `bounds.zmax - bounds.zmin`; ignores measuredHeight |

A probe saved measuredHeight=12.5 with geometry spanning Z=0..3. After canonical
load, `building.height` returned **3.0**. Adding the separate height attribute with
value 8 made it return **8.0**. All are possible today without a semantic conflict
diagnostic. The getter's metres claim is also not established merely by subtracting
local Z coordinates: CRS units and transforms need to support that interpretation.

Recommended direction: one declared scalar measurement authority, with simple
Python access; geometry-based estimation must be an explicit operation. Do not
silently choose a roof/ground reference pair from several qualified measurements.
Do not write a computed extent back as a supposedly measured value. Resolve the
exact `building.height` behavior while auditing its existing builder producers and
consumers, before changing that public property. Missing measurements should remain
missing, not become zero or an undocumented fallback.

A future qualified-measurement record can use native typed attributes, keeping
`value`, `uom`, references and status together without extending Protobuf. The
standard semantic backend currently treats class-valued slots as non-inlined ID
references and explicitly rejects inlined class slots. Therefore declaring such
records is **not yet a YAML-only change**: it needs a small generic backend
extension, followed by focused nested-record validation. It does not need a new
native class for every CityGML data type or a parallel validator.

## Other inherited facts relevant to buildings

Definitions:
[lifespan](https://docs.ogc.org/is/20-010/20-010.html#AbstractFeatureWithLifespan-section),
[space](https://docs.ogc.org/is/20-010/20-010.html#AbstractSpace-section),
[city-object properties](https://docs.ogc.org/is/20-010/20-010.html#AbstractCityObject-section).

| Source property / cardinality | Current DTCC handling | Decision |
| --- | --- | --- |
| creationDate / terminationDate: each DateTime `0..1` | Generic feature attributes; Dataset Context is separate | Preserve model/database lifecycle separately from the physical object's lifespan. Defer temporal rules until that distinction is explicit. |
| validFrom / validTo: each DateTime `0..1` | Generic feature attributes | Preserve physical validity separately from construction events and simulation time. Timezone/interval policy remains to be designed. |
| relativeToTerrain / relativeToWater: each enumeration `0..1` | Generic attributes | Do not infer these classifications from a cached bounding box. Add only with an explicit use case and declared enum. |
| spaceType: enumeration `0..1` | Generic | Retain as metadata pending a workflow; do not create a Python space hierarchy solely for it. |
| area / volume: qualified records `0..*` | Generic attributes, or computed geometry values outside this schema | Area needs m² or a declared alternative; volume needs m³. Preserve the kind of area/volume; footprint area, floor area and envelope volume are not interchangeable. |

## Objects, containment and surface relationships

Building and BuildingPart are **siblings** under AbstractBuilding in CityGML.
DTCC's semantic schema currently makes BuildingPart a subtype of Building to reuse
rules; native BuildingPart instead subclasses Object and does not inherit all of
Building's convenience methods. This is a DTCC modeling shortcut, not source
ontology equivalence. A shared abstract **schema** class can express shared facts
without requiring an additional native Python class.

| Source relationship / cardinality | Current DTCC mapping and validation | Decision |
| --- | --- | --- |
| Building → buildingPart: CM table 674 `1..*`; XSD `0..*` | Building children retain BuildingPart identity; single-owner acyclic containment checked. Native roots may be standalone parts. Strict CityJSON requires a part's one parent and reciprocal links. | Keep optional parts and allow standalone native results. Apply external completeness requirements at the adapter boundary. |
| Part within another Part | Native/schema/strict adapter permit it. Reviewed CityGML 3.0 BuildingPart has no buildingPart property of its own. | Record as a DTCC hierarchy choice; a future GML adapter must address it explicitly rather than claiming a direct encoding. |
| Building/Part → address: CM table 668 `1..*`; XSD `0..*`; CityJSON optional address array | Strict adapter rejects top-level `address`; generic attributes can retain address-like data but have no standard ownership/location rule | Defer full addresses. Preserve structure, multiplicity and optional location when implementing; moving it into attributes without an explicit reverse mapping is not faithful CityJSON support. |
| Building/Part → installation / furniture / constructive element / room / subdivision: CM table 668 `1..*`; XSD `0..*` | No declared type/relationship rules or strict external mapping for these types | Prefer generic Objects plus schema classes when needed. Never flatten them into a parent and lose identity. |
| Building/Part → boundary: conceptual construction table `1..*`; Core XSD `0..*` | SemanticRegion membership selects existing surfaces/faces; regions are optional and local to a representation | Keep this efficient representation for the implemented surface slice. It does not represent every independently identified, shared, multi-LoD CityGML boundary feature. |
| Construction surface → fillingSurface: CM table 540 `1..*`; XSD `0..*` | Optional region.parent gives one local host; index bounds and cycles checked. Known Window/Door hosts must be WallSurface or RoofSurface under the DTCC schema. | Keep host distinct from object containment. Missing host is not inferred from a polygon hole; current host type restriction is DTCC's bounded subset. |
| Window/Door physical feature → WindowSurface/DoorSurface | DTCC Window/Door schema names denote **surfaces**, following CityJSON spelling | Do not use those URIs for physical frames, glazing or door elements. Add separate semantic classes if physical elements become required. |
| generalizesTo / relatedTo: optional repeated encoding relations | Named `Object.relations` preserve local target IDs and validate target existence; these meanings are undeclared | A typed local relation can be schema-driven. Qualified relation records and links to a region are not currently expressible as ordinary object-ID relations. |

References for this table are the
[Building catalogue](https://docs.ogc.org/is/20-010/20-010.html#AbstractBuilding-section),
[construction surfaces](https://docs.ogc.org/is/20-010/20-010.html#AbstractConstructionSurface-section),
[physical Window](https://docs.ogc.org/is/20-010/20-010.html#Window-section),
[WindowSurface](https://docs.ogc.org/is/20-010/20-010.html#WindowSurface-section),
the encoding XSDs above, and [CityJSON §§2.3, 3.3](https://www.cityjson.org/specs/2.0.2/).

The strict adapter recognizes eleven building surface spellings, while the standard
schema declares only RoofSurface, WallSurface, GroundSurface, Window and Door.
ClosureSurface, OuterCeilingSurface, OuterFloorSurface, InteriorWallSurface,
CeilingSurface and FloorSurface currently receive generic region checks. Importing
one successfully does not mean all its surface-specific meaning was validated.
These six classes are a bounded schema-vocabulary gap; introducing Python subclasses
is unnecessary.

## Geometry, footprint access and intentional limits

The existing representation store is the right authority: geometry plus explicit
representation ID, LoD and optional role. Geometry arrays stay out of LinkML.
This is consistent with separating a feature's meaning from its representations,
without claiming every source geometry relationship is already mapped.

| Topic | Current DTCC behavior / limit | Decision |
| --- | --- | --- |
| Multiple representations at one LoD | Stored separately; ambiguous singular selectors fail | Keep direct ID/role selectors. Do not assume one LoD implies one geometry. |
| Footprint | `get_geometry(id=...)` retrieves stored data. `building.footprint()` defaults to LoD0, creates a polygon-derived Surface and can retain only the largest disconnected polygon. | Distinguish exact stored footprint access from derived extraction. Do not equate every LoD0 representation with the authoritative footprint or discard components when exact exchange is intended. |
| Strict Building/Part CityJSON geometry | MultiSurface and Solid supported; CompositeSurface/CompositeSolid and other native geometry require explicit mappings | Keep unsupported encodings explicit. Strict import does not establish geometric watertightness. |
| Representation IDs / roles / region IDs | Preserved canonically; strict CityJSON export restricts IDs to its generated cityjson-N convention, rejects roles and explicit region IDs | A browser-readable `.dtcc` building is not automatically exportable as strict CityJSON. Define extension mappings before widening the adapter. |
| CRS and transforms | Native matrices/SRS survive; strict CityJSON dequantizes coordinates and exports only an agreed global frame/CRS | Dequantization is not reprojection. The source integer grid is not the native transform authority. Vertical reference and unit interpretation remain explicit concerns. |
| Point clouds, terrain-intersection curves and implicit geometry | Native PointCloud/LineString carriers exist; strict building mappings and instancing are absent | Do not infer complete CityGML support from available numerical carriers. Defer until required by a consumer. |
| Appearance, shared boundary identity, interior subdivisions and temporal simulation axes | No complete canonical/external contracts for this building slice | Deliberately outside this milestone. Preserve the inventory gaps and avoid speculative classes. |

CityJSON's representation/geometry rules are in
[§§2.3, 3 and 4](https://www.cityjson.org/specs/2.0.2/);
inherited CityGML geometry roles are in
[AbstractSpace](https://docs.ogc.org/is/20-010/20-010.html#AbstractSpace-section) and
[AbstractPhysicalSpace](https://docs.ogc.org/is/20-010/20-010.html#AbstractPhysicalSpace-section).

## Proposed implementation sequence

1. **Resolve height authority and Python access.** Audit the existing height
   producers/consumers, decide the declared measurement versus explicit geometry
   estimate behavior, and exercise ordinary building construction and canonical
   save/load. Do this before adding a second competing height abstraction.
2. **Publish the next local schema version for the small exterior-building slice.**
   Add optional description, class, function, usage and storey counts with the
   decisions above; complete the six missing surface declarations and use a shared
   abstract building schema class if it clarifies inheritance. New declarations
   constrain previously open metadata. Keep 0.2.0 unchanged and document the
   validation change. No production dependency or `.proto` change is needed for
   these scalar/list rules. Validate boolean-versus-integer behavior explicitly.
3. **Apply the selected standard schema at the strict CityJSON boundary.** Preserve
   explicit bypass semantics and faithful unsupported-content failures. First test
   a source with declared units/attribute types; do not label all conforming
   CityJSON as conforming to DTCC's more specific semantic contract.
4. **Add qualified records only for a concrete workflow.** Measurements, code-space
   records, external references or addresses justify one generic inlined-record
   capability in the schema backend, not bespoke validators per feature. Date and
   nil-reason policies belong in the relevant versioned contract.

Permanent vocabulary URIs still need a project-owned namespace before publication.
The current `https://example.org/dtcc/` URIs are experimental; the public schema ID
and schema version are different identifiers. A namespace decision must cover
surface-versus-physical Window/Door names and future stable concept identity, not
just replace a string in one YAML file.

## Observed behavioral evidence

The local public-API probe completed with these results:

- measuredHeight=12.5 survived `.dtcc`, while the current height accessor returned
  a geometry span of 3.0, then 8.0 when a height attribute was added.
- Negative storey count, malformed date text, scalar function and a structured
  height attribute survived canonical persistence as undeclared metadata.
- List-valued name and structured measuredHeight were rejected by their declared
  scalar slots, demonstrating the boundary between preservation and validation.
- Negative measuredHeight passed strict CityJSON import/export but failed the
  subsequent default `.dtcc` save with an attribute-path diagnostic.
- CityJSON address was explicitly rejected by the strict adapter; an undeclared
  ClosureSurface URI passed canonical generic-region admission.

Evidence: `/private/tmp/dtcc-building-crosswalk-evidence/report.json`; the standalone
probe is `/private/tmp/dtcc-building-crosswalk-probe.py`. These are design probes,
not regression tests that enshrine behavior we intend to improve. Use the public
`io.save_model(value, path)` for roots without a type-specific `.save` convenience.
No runtime/schema behavior was changed. The
[implementation plan](../../.agent/plans/2026-09-11-building-crosswalk.md) records
completion and source evidence; the [native inventory](model-inventory.md) retains
the remaining thematic and numerical work.
