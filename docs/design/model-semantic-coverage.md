> Historical mixed-city slice, retained as implementation evidence. For the current
> authoritative coverage matrix, see [exterior-city profile](exterior-city-profile.md).
> Version numbers, unsupported features and acceptance claims below describe that
> earlier milestone and do not define current behavior.

# DTCC semantic coverage and the mixed-city slice

This milestone document records historical implementation evidence. For current
wire support and schema coverage, see [the native inventory](model-inventory.md)
and [standard I/O contract](standard-schema-io.md). Old wire readers and legacy
serialization paths described below have since been removed.
The current [terrain/transportation crosswalk](terrain-transportation-crosswalk.md)
supersedes the historical terrain/transport rows: their numerical carriers are now
canonical, while strict standards conversion still needs implementation.
The subsequent [strict triangular TINRelief slice](strict-tin-relief.md) now
implements that one terrain conversion; other relief forms and transport remain open.

Status: implemented bounded mixed-city slice, 11 September 2026. Authority:
[Core Design](../../DESIGN.md) and the [model contract](model-contract.md).
The current city schema is the single self-contained file
[city/0.1.0/schema.yaml](../../sandbox/model_profiles/profiles/city/0.1.0/schema.yaml).
The buildings 0.1.0 and 0.2.0 files remain usable and unchanged.

## Coverage map

The table distinguishes a semantic vocabulary gap from missing native numerical
facts and external-format mappings. It is an inventory, not a commitment to
implement the entire standards hierarchy. The current acceptance slice comprises
buildings, land use, individual plants, explicit native trees and city furniture.

| Area | Standards concept | Existing DTCC authority | Implemented coverage and remaining gap |
| --- | --- | --- | --- |
| Aggregate | CityJSON document containing CityObjects | City and Object.children | City aggregate and typed children round-trip. BuildingPart containment is supported; other part hierarchies need explicit rules. |
| Buildings | Building, BuildingPart | Existing typed classes, ordered geometry representations | Buildings profile, exact LoDs, polygon/Solid geometry, openings and canonical exchange work. Rooms/installations, addresses and other metadata remain gaps. |
| Boundary/opening surfaces | Roof/wall/ground and filling surfaces | SemanticRegion membership and local parent index | Buildings 0.2.0 rules also appear in the city profile. Other thematic surfaces require schema definitions; geometry operations need their own semantic mapping. |
| Land use | LandUse | Landuse, its geometry and native landuses list | New canonical typed-state support and strict MultiSurface mapping. Source classification attributes remain separate from native enum codes; neither park nor grass is inferred from the other. |
| Individual vegetation | SolitaryVegetationObject | Generic Object for unqualified source plants; Tree for explicit native trees | New profile types; strict source plants preserve attributes and supported geometry. Native Tree position/height/radius round-trip canonically. No automatic source-to-Tree specialization or fabricated measurements. |
| Vegetation cover | PlantCover | Geometry can be carried by Object; no dedicated implementation evaluated | Vocabulary/profile rules and adapter mapping remain. A land-use polygon does not establish a measured vegetation cover. |
| City furniture | CityFurniture | Generic Object with semantic URI, attributes and geometry | Profile and strict mapping added. Bench function needs no Python subclass. Schema-only ChargingBench extension is demonstrated. Templates, appearance and richer geometry remain outside strict coverage. |
| Terrain | Relief / TINRelief | Terrain with Mesh/Raster representations | Existing legacy paths remain; Terrain, Raster and CompositeSurface are not in the canonical/strict subset. Requires numerical and adapter work, beyond adding a semantic name. |
| Transport | Road, Railway, Waterway, TransportSquare | RoadNetwork graph vertices/edges/lengths | Graph state and line geometry are not yet canonical. A routing graph is not automatically a road's physical extent; define a concrete workflow before mapping these concepts. |
| Water | WaterBody and water boundaries | Generic Object, geometry, and numerical fields | No dedicated profile/strict mapping certified. A native WATER land-use code does not supply a water volume or water-level interpretation. |
| Sensors/simulation | Observations and computational values extend this city vocabulary | Object relations, Field and geometry | Generic canonical relationships/fields exist; earlier optional sensor experiment remains evidence. General temporal/scenario axes and further numerical types are separate gaps. |

Primary references: [CityJSON 2.0.2 feature and geometry definitions](https://www.cityjson.org/specs/2.0.2/)
and [CityGML 3.0 conceptual model](https://docs.ogc.org/is/20-010/20-010.html).
CityJSON permits MultiSurface/CompositeSurface for LandUse and any Geometry Object
for CityFurniture and SolitaryVegetationObject. TINRelief uses CompositeSurface.
CityGML individual vegetation includes plants other than trees. These distinctions
motivate the bounded mappings below; DTCC's computational types are not declared
formally equivalent to every related standards concept.

## What changes only the schema

The new city profile includes the buildings 0.2.0 definitions directly, adding
LandUse, SolitaryVegetationObject, Tree and CityFurniture. It has no imports or
remote resolution. Adding ChargingBench with a required nonnegative chargingPower
attribute changes only an extended YAML file; the same Object, geometry and wire
continue to work. No Python class registry or generated runtime is introduced.
Profile identity/version remains explicit and independent of wire version.

Tree is a semantic subtype of SolitaryVegetationObject in this profile. Generic
objects can carry semantic subtypes; a Python class is not chosen from a URI on
canonical reads. Concrete native types come from the wire kind. The DTCC metre
conventions for height, diameter and radius are stated in YAML; the validator
checks values but cannot discover or convert an unknown source unit convention.
Missing source measurements remain absent.

The native LanduseClasses list is an existing computational API, not a universal
land-use classification. Its current names are described by NativeLanduseCode in
the profile. Adding a new native enum member still needs a code/wire-contract
decision; ordinary source classification strings and new semantic subtypes do not.
Self-contained versioned profiles repeat stable definitions intentionally. Edits
create new versions; no runtime sharing or remote imports are required.

## Canonical v5: preserve existing typed state

The prior codec rejected Tree and Landuse. V5 admits them with explicit kind tags
and matching typed state. It does not move their native properties into attributes.

| Native fact | Canonical representation |
| --- | --- |
| Tree.position | Existing numeric Array encoding; dtype, shape and exact finite values retained. XYZ or existing empty forms are supported. |
| Tree.height, Tree.crown_radius | Finite nonnegative doubles, admitted only when conversion is exact. Native zero defaults are preserved as supplied facts. |
| Landuse.landuses | Ordered list of native enum names, preserving empty lists and repetitions; unknown names fail. |
| Other Object facts | Existing geometry, transforms, attributes, children, references and profile identity remain independent and preserved. |

State must match the concrete kind and be present for Tree/Landuse. Old version
declarations cannot carry these new kinds or payloads. Writers emit v5; readers
accept v1–v5. Frozen v1–v4 fixtures remain readable. Unknown/malformed state fails
before file replacement. Canonical package manifests remain v3 and identify the
embedded model version separately.

The existing Tree class behavior is unchanged: position is its separate native
property, and no point representation is automatically added. Tree positions do
not currently contribute to Object's geometry-derived bounds. This milestone
preserves that data without claiming new bounding, reprojection or geometric
interpretation behavior. The source point-object path has ordinary geometry bounds.

The optional semantic projection emits Tree measurements and native land-use code
names from their actual properties. It rejects attribute names colliding with
those projected fields. This is disposable validation input; canonical exchange
can still preserve both independently named native facts and arbitrary attributes.

## Strict CityJSON mapping

| Source feature | Native class | Geometry admitted in this slice |
| --- | --- | --- |
| Building / BuildingPart | Building / BuildingPart | MultiSurface, Solid |
| LandUse | Landuse | MultiSurface |
| SolitaryVegetationObject | Object with the source semantic URI | MultiSurface, Solid, singleton MultiPoint mapped to Point |
| CityFurniture | Object with the source semantic URI | MultiSurface, Solid, singleton MultiPoint mapped to Point |

All selected features are direct City children except BuildingPart. Their source
IDs, attributes, geometry order and LoD labels are preserved. The strict adapter
does not interpret function/class/species strings to manufacture a Python class.
Source LandUse leaves native landuses empty and retains source codes in attributes.
No source missing height is replaced by Tree's default zero.

Geometry collections, multi-point cardinalities other than one, point semantic
regions and non-building surface-region vocabularies remain unsupported. Only
building-part child hierarchies are currently admitted. Native Tree extra state,
nonempty native Landuse code lists, profile identity and new semantic URI subtypes
fail strict export until a faithful external mapping is defined. Users must not
silently downcast or erase these facts to claim a lossless conversion.

The native typed example and source example are deliberately distinct. Each
round-trips through canonical file/package APIs. The source example additionally
round-trips through CityJSON. It would be incorrect to claim that the full native
Tree contract has a faithful plain CityJSON mapping today.

## Python access and reproduction

```python
from dtcc_core import io
from dtcc_core.model import Landuse, Object

city = io.load_city("sandbox/model_profiles/fixtures/mixed-city.city.json", strict=True)
footprint = city.buildings[0].footprint().to_polygon()  # 96 m²
park = city.get_children(Landuse)[0]
park_polygon = park.lod0.surfaces[0].to_polygon()       # 360 m²
objects = {obj.id: obj for obj in city.get_children(Object)}
bench = objects["bench-1"]
point = bench.lod1
plant_species = objects["tree-1"].attributes["species"]

native = io.load_model("/tmp/dtcc-mixed-city/native.dtcc")
tree = native.trees[0]
print(tree.position, tree.height, tree.crown_radius)
```

This uses existing accessors and one ordinary dictionary for ID lookup. It adds
no persistent ID cache or alternative geometry collection. Landuse.surfaces still
refers to its legacy MultiSurface role; explicitly named/LoD representations use
the ordinary representation accessors shown here.

From the repository root, with the established native and isolated tooling environments:

```bash
../venv/bin/python sandbox/model_profiles/mixed_city_example.py /tmp/dtcc-mixed-city
../venv/bin/python sandbox/model_profiles/buildings_example.py /tmp/dtcc-buildings.json
PYSTOW_HOME=/tmp/dtcc-profile-pystow /tmp/dtcc-model-profile-85-env/bin/python sandbox/model_profiles/evaluate_mixed_city.py /tmp/dtcc-mixed-city/profile.json /tmp/dtcc-buildings.json
../venv/bin/python -m pytest tests/io/test_mixed_city.py -q
```

The source fixture's tree and bench use point locations, not detailed physical
geometry. The separately authored native example supplies an explicit Tree's
position, height and radius and a native GRASS code. These deliberate authoring
choices are not an automatic classification or diameter-to-radius conversion.
The optional ChargingBench experiment proves a new furniture subtype without
editing DTCC classes, geometry or wire schema. Evidence for this run is recorded
in the [implementation plan](../../.agent/plans/2026-09-11-mixed-city.md).

Verification: 1,332 affected regression tests passed (1 skipped, 53 deselected),
followed by 65 passing focused checks after the final changes. Both source/export
fixtures passed the recorded official CityJSON 2.0.2 schema. Mixed-city, buildings
and openings profile checks passed. The complete 3DBAG v3 tile preserves its
24,504,507 bytes after restoring only the prior envelope version. No new latency
or memory benchmark is claimed. Runtime artifacts are under
`/private/tmp/dtcc-mixed-city-evidence/`.
