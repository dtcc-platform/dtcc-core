# DTCC buildings profile: standards mapping and first workflow

This milestone document records historical implementation evidence. For current
wire support and schema coverage, see [the native inventory](model-inventory.md)
and [standard I/O contract](standard-schema-io.md). Old wire readers and legacy
serialization paths described below have since been removed.

Status: first local subset, 10 September 2026. This implements the buildings
milestone of issue #85; it does not certify full CityGML or CityJSON conformance.
[Core Design](../../DESIGN.md) remains the authority. The user's controlling
constraints are a simple generic implementation, efficient numerical exchange,
direct Python access and independently editable versioned semantic schemas.

## Read and edit the schema

The current entry point is
[`schemas/profiles/buildings/0.2.0/schema.yaml`](../../schemas/profiles/buildings/0.2.0/schema.yaml).
The unchanged [0.1.0 schema](../../schemas/profiles/buildings/0.1.0/schema.yaml)
records the first subset described below; [0.2.0 adds opening surfaces and hosts](building-openings.md).
It contains all its own classes, slots and primitive type definitions; it has no
imports and needs no network resolution. The `example.org` namespace explicitly
marks this as an experimental vocabulary. Choose a durable namespace before a
public schema release; no external ontology equivalence is asserted by these URIs.

Keep released versions immutable. A changed vocabulary or constraint set gets a
new versioned file, retaining the profile ID. A validation operation selects its
schema file explicitly and reuses that loaded view. Two versions can coexist.
There is no process-global profile replacement or automatic remote schema fetch.
Standalone native data and the canonical codec remain usable without LinkML.

Adding a semantic subtype or attribute normally edits only this YAML. A new native
numerical representation, reference scope or axis still needs an implementation
and possibly a wire-version change. Schema editability cannot manufacture missing
geometry operations or guarantee their correctness. This distinction prevents
semantic schema changes from turning into a generated Python class hierarchy.

The first profile declares City, Building, BuildingPart and roof/wall/ground
regions; 0.2.0 also declares window/door filling surfaces. Height and classification are optional; the restrictive sensor and
low-rise rules in the earlier experiment are not universal building requirements.
The `measuredHeight` property uses metres by explicit profile convention. An
external numeric value with unknown units must not be assumed to have that meaning.
Additional source attributes are retained and allowed by this profile's open
validation mode. Domain code lists and richer quantity representations need
explicit source mappings before tighter validation is introduced.

## What we take from the standards

References are the [CityGML 3.0 conceptual model](https://docs.ogc.org/is/20-010/20-010.html)
and [CityJSON 2.0.2](https://www.cityjson.org/specs/2.0.2/). CityJSON describes its
[implementation differences](https://www.cityjson.org/citygml/v30/) separately.
CityGML provides the conceptual definitions; CityJSON supplies the concrete
interchange fixture. The profile is a DTCC application of selected concepts.

| Concept | Standard distinction | Native decision and current coverage |
| --- | --- | --- |
| Model aggregate | CityJSON holds features in its document's CityObjects map. | DTCC City aggregates native children. Source feature IDs survive. Its generated aggregate UUID is not a CityJSON feature ID. |
| Building / BuildingPart | Whole building versus a physical or functional subdivision. | Existing concrete classes, IDs, attributes and recursive children. Schema inheritance does not force a change to Python inheritance. BuildingPart now offers the same `building_parts` convenience as Building. |
| Boundary role | Roof, wall and ground boundaries describe meaning, independently of geometric primitives. | Generic SemanticRegion with a semantic URI. No RoofSurface Python class or per-triangle semantic object. |
| Assignment | Several primitives can refer to one semantic entity; a primitive may be unclassified. | A region stores one NumPy index array. Unselected elements remain unclassified. Multiple assignments to one element are rejected in this subset. |
| Openings | Semantic entities may have local parent/child relationships. | Implemented by the [openings milestone](building-openings.md): one optional native region parent index, with inverse children derived. |
| Geometry | MultiSurface and Solid have different structure; solids include shells. | Preserve MultiSurface polygon rings/holes and Solid exterior/interior shells. Composite geometries remain unsupported. |
| LoD | One object may have several geometries, including several at one LoD. | Preserve exact string LoDs and repeated representations using local IDs; selection raises on ambiguity. See the representation contract. |
| Properties | Source attributes, code lists and measures need explicit interpretation. | Typed attributes remain native. Optional `measuredHeight` and `roofType` are declared; no universal height/usage requirement or invented classification enumeration. |
| Coordinates | CityJSON dequantizes integer vertices with scale and translation. | Convert to global float64 coordinates on import. Preserve CRS in native transforms. Strict export requires identity affines and one CRS, then applies requested quantization. |
| Surface identity | Source semantic properties belong to the shared entity, not every polygon. | Preserve scalar properties once per region. Native region IDs are geometry-local; exporting a native region ID requires a future explicit CityJSON mapping. |
| Context | Document metadata differs from intrinsic object facts. | This strict input subset accepts referenceSystem and checked derived extents; explicit extent recomputation retains discrepant summaries in Dataset Context. Other source metadata fails explicitly pending Dataset Context mapping. Canonical packages preserve attached root Context. |

The strict boundary validates representability and intrinsic structure. It does
not claim to implement every normative CityJSON geometry test. Planarity,
self-intersection and the full thematic vocabulary require appropriate geometry
and profile checks. The profile validator is separate from this admission check.
The fixture and exported example also pass the published CityJSON 2.0.2 JSON
Schema from the [official schema distribution](https://3d.bk.tudelft.nl/schemas/cityjson/2.0.2/).
An object may have no geometry; an encoded MultiSurface must have nonempty
boundaries. Empty native MultiSurfaces remain valid in canonical exchange.

## Generic native representation

`SemanticRegion(semantic_type, indices, id=None, attributes={}, parent=None)` is a small value
attached through `Geometry.regions`. On MultiSurface, indices address its Surface
list; on Solid they address its flat Surface list; on Mesh they address triangular faces. The region owns metadata, not copied
vertices. One roof can cover two polygons or many triangles. Optional region IDs
are unique within their geometry; they are not added to the Object reference graph.
The optional parent indexes the same geometry's region list; see the
[reference scope and editing rules](building-openings.md).

The generic codec checks URI spelling, typed attributes, index dtype/range,
uniqueness and non-overlap. It does not know the roof vocabulary. Region order and
index array dtype are preserved in canonical exchange. Surface arrays, normals,
holes and each nested geometry's own transform/fields are also preserved.

A topology change must explicitly preserve/reindex regions or report that it
cannot. MultiSurface.merge offsets copied element and parent indices. The subsequent
[real-building milestone](real-buildings.md) adds region-preserving triangulation
through `MultiSurface.mesh()` and batch meshing, using the existing per-surface
backends. Cleaning, welding/snapping and field interpolation remain unsupported
on that path. The raw C++ geometry conversion boundary still rejects regions;
the meshing wrapper transfers membership after geometry conversion. Other
geometry algorithms still require an explicit
region-preservation audit before their outputs can be treated as faithful semantic
transformations.
Geometry-only copies intentionally produce derived geometry without region facts.
Raw NumPy/list mutations remain the caller's responsibility; no observer framework
or automatic topology tracking has been introduced.

## Python acceptance path

From the repository root in the DTCC environment:

```python
from dtcc_core import io

city = io.load_city("sandbox/model_profiles/fixtures/buildings.city.json", strict=True)
building = city.buildings[0]
print(building.id)                              # building-1
print(building.attributes["name"])              # Courtyard building
polygon = building.footprint().to_polygon()     # 96 m², with its courtyard hole
part = building.building_parts[0]
print(part.id)                                  # part-1
roof = part.lod2.regions_of("https://example.org/dtcc/RoofSurface")[0]
roof_polygons = [part.lod2.surfaces[i] for i in roof.indices]
print(roof.attributes["solar-potential"])        # 42.5, source attribute

city.save("/tmp/buildings.dtcc")
restored = io.load_model("/tmp/buildings.dtcc")
restored.save("/tmp/buildings.city.json", strict=True)
```

The fixture deliberately labels its LoD0 ground geometry. The existing footprint
operation projects the chosen geometry to XY; its default selects LoD0. LoD0
alone does not prove a ground-boundary interpretation, and the existing footprint
method can select the largest component of a disconnected result. Resolving
multi-component/roof-print footprint semantics is a separate explicit API decision.

`strict=True` opts into this documented CityJSON subset on the existing public
load/save path. It never catches a building error and returns a successful partial
city. Ordinary legacy import/export remains available during migration. Both modes
share geometry conversion and the building exporter traversal; strict admission
adds the explicit representation boundary.

The CityJSON round trip preserves source feature IDs, attributes, containment,
polygon rings/holes, surface grouping, unclassified elements and CRS. Canonical
round trips are exact for these native facts. CityJSON export normalizes vertex
indexing/deduplication, recomputes checked extent summaries, preserves exact LoD
strings and quantizes coordinates with the requested scale (default 0.001).
It may write explicit all-null assignments for previously absent semantics and
generates a fresh native aggregate City UUID on reimport. A schema/profile identity
has no plain CityJSON mapping in this subset: export fails if one is attached.
The example separately demonstrates its preservation through canonical packages.

Unsupported strict inputs include other feature types, composite geometries,
appearances, geometry templates, extensions, addresses, additional document
metadata and physical filling-element objects. These are recorded coverage gaps, not
redefinitions of valid CityJSON. Required future platform data must expand coverage.
Strict export also rejects unrepresented native fields, stored polygon normals,
profile facts and named references rather than dropping them.

## Versioning and validation

The semantic profile starts at 0.1.0. Canonical **model wire version 2** adds Surface,
MultiSurface and regions. [Wire version 3](geometry-representations.md) adds ordered
representations and Solid. [Wire version 4](building-openings.md) adds region
parents. The [mixed-city milestone](model-semantic-coverage.md) adds native Tree
and Landuse state in v5; current readers accept versions 1–5 and writers emit 5. Profile 0.2.0
adds opening surface types independently of this wire version.
A payload claiming version 1 cannot contain the new representation. The v1 Point
fixture remains readable, and package readers check manifest model versions against
the actual decoded artifact. The canonical package manifest remains v3 because its
existing format already records each artifact's model version.

LinkML stays in the isolated development environment. Its generic adapter has
moved to `sandbox/model_profiles/linkml_profile.py`, shared by both experiments.
Native geometry arrays do not enter its semantic projection. Validation selects
the schema explicitly; ordinary property access and encoding/decoding do not load
or invoke LinkML. This is not yet a production profile-loader API or dependency
adoption. Native files can preserve a semantic type unknown to an older profile;
that older profile reports it as unsupported when asked to validate it.

Reproduce the two-process workflow:

```bash
../venv/bin/python sandbox/model_profiles/buildings_example.py /tmp/dtcc-buildings-example.json
PYSTOW_HOME=/tmp/dtcc-profile-pystow /tmp/dtcc-model-profile-85-env/bin/python sandbox/model_profiles/evaluate_buildings.py /tmp/dtcc-buildings-example.json
../venv/bin/python sandbox/model_profiles/benchmark_buildings.py --output /tmp/dtcc-buildings-benchmark.json
```

The evaluator copies only the YAML to an unrelated directory, validates the
restored records, then creates version 0.1.1 with SolarRoofSurface and required
`efficiency`. The same native SemanticRegion and wire shape survive. Version 0.1.0
rejects the new type; 0.1.1 accepts it and rejects a missing required efficiency.
Both versions coexist. The seven validity checks include optional and invalid
provided height. The original 18-case experiment remains separate regression evidence.

## Measurement and next coverage

Initial local measurement: Python 3.12.12, macOS 26.6.2 arm64, 500 repeated synthetic
buildings and 4,000 polygons. Median of three warmed runs with collection between
runs; model construction, disk I/O and LinkML excluded. With regions: 1,776,445 bytes,
154.8 ms encode and 234.3 ms decode (10.95 and 7.23 MiB/s). Without regions on the
same geometry: 1,616,945 bytes, 138.8 ms and 210.6 ms. This is a small repeatable
baseline, not a real-city scale/peak-memory or browser/C++ performance claim.

The [representation milestone](geometry-representations.md) now exercises the
complete real 3DBAG tile with fractional/repeated LoDs and explicit shells.
Opening relations are now covered by the [next milestone](building-openings.md).
Remaining decisions include further source metadata mapping
and region-aware topology changes beyond triangulation.
The [real-building evaluation](real-buildings.md) retains the initial data evidence. Benchmark actual geometry before changing
the packed array representation or adding caches.
