# First exterior-city profile

DTCC Model's current standard contract is schema **0.9.0**, with semantic names
under `https://github.com/dtcc-platform/dtcc-core/schemas/model#`. This milestone
covers a useful exterior urban model on the existing generic native carriers.
It does not claim full CityGML conformance or complete CityJSON support.

The authority for declared semantic rules is
[`schema.yaml`](../../dtcc_core/schemas/model/0.9.0/schema.yaml). Numerical arrays,
indices, native types and persistence integrity remain Core responsibilities.
The CityJSON adapter separately checks whether the selected external format can
preserve each representation and value. These checks remain active with
`validate_schema=False`.

## Current coverage

This table describes **strict CityJSON 2.0 exchange**, not the broader inventory of
native storage types. Polygon surfaces preserve exterior rings and holes; Solid
preserves exterior and cavity shell partitions. The implementation does not certify
polygon planarity, aggregate connectedness, watertightness or manifold validity.

| Feature type | Native carrier | Strict representations | Semantic boundary |
|---|---|---|---|
| Building, BuildingPart | Building, BuildingPart | MultiSurface, CompositeSurface, Solid | Existing building surface/opening regions and selected attributes; building-part containment |
| LandUse | Landuse | MultiSurface, CompositeSurface | Generic source attributes plus existing selected rules; native `landuses` enum state has no automatic CityJSON mapping |
| SolitaryVegetationObject | Object | Singleton MultiPoint, MultiSurface, CompositeSurface, Solid | Individual plant metadata; procedural native Tree state has no strict CityJSON mapping |
| CityFurniture | Object | Singleton MultiPoint, MultiSurface, CompositeSurface, Solid | Identifiable furniture and existing selected metadata rules |
| TINRelief | Terrain | Triangular CompositeSurface | Mesh triangles and exact LoD; no relief regions, marker/normal/field conversion or nontriangular facets |
| Road, Railway, TransportSquare | Object | MultiLineString, MultiSurface, CompositeSurface | TrafficArea, AuxiliaryTrafficArea, TransportationMarking, TransportationHole regions |
| Waterway | Object | MultiLineString, MultiSurface, CompositeSurface | Physical transportation feature; no regions in this bounded mapping |
| WaterBody | Object | MultiLineString, MultiSurface, CompositeSurface, Solid | WaterSurface, WaterGroundSurface, WaterClosureSurface regions |
| PlantCover | Object | MultiSurface, CompositeSurface, Solid | Optional average_height and common classification metadata; no surface-region vocabulary |

There are **12 supported external feature types**. Coverage within each type remains
bounded by this table. Except BuildingPart, admitted features are top-level; the
adapter does not reinterpret new theme relationships as containment. A City is the
native aggregate corresponding to the CityJSON document, not an extra CityObject.

A native terrain raster also has an explicit elevation interpretation workflow:
unit, vertical reference, representation identity and cell-centre sampling remain
associated with the Terrain in native files and canonical packages. Rasters remain
rasters. There is no implicit raster-to-CityJSON conversion. See the
[runnable DEM example](../../sandbox/model_profiles/dem_example.py) for classified
ground selection, negative elevations, CRS and nodata handling.

## Geometry and Python access

The new features use `Object`, `MultiLineString`, `MultiSurface`, `Solid` and
`SemanticRegion`; no new native feature hierarchy or Protobuf geometry type is
needed. Native `.dtcc` remains an ordinary `DTCC.ModelFile` Protobuf message using
wire version 6.

An imported polygon CompositeSurface uses a MultiSurface carrier and
`representation.role == 'cityjson:CompositeSurface'`. This explicit role preserves
the source aggregate declaration through native exchange and restores its exact
CityJSON spelling on export. A polygon MultiSurface has no role. The existing
triangular TIN mapping uses Mesh and no role. Other representation roles remain
unsupported by strict CityJSON. Preserving the declaration does not prove the
CompositeSurface's topology.

Strict representation IDs remain `cityjson-0`, `cityjson-1`, and so on, in source
geometry order. LoDs retain their exact strings. Both MultiSurface and
CompositeSurface can occur at the same LoD, so selecting by representation ID is
appropriate when the LoD is ambiguous.

```python
from dtcc_core import io

city = io.load_city("exterior.city.json", strict=True)
features = {obj.id: obj for group in city.children.values() for obj in group}
road = features["road"]

line_vertices = road.get_geometry(lod="0").linestrings[0].vertices
surface = road.get_geometry(id="cityjson-1")
polygon_vertices = surface.surfaces[0].vertices
material = surface.regions[0].attributes["surface_material"]
bench_id = features["bench-1"].id

city.save("exterior.dtcc")
city.save("exterior.city.json", strict=True)
```

The arrays are ordinary NumPy arrays. RoadNetwork retains its separate native
vertices/edges/lengths routing-graph meaning; no routing topology is inferred from
these physical road features. Waterway similarly denotes transportation, while
WaterBody denotes the water feature itself.

## Metadata and semantic regions

The six new feature classes use optional scalar `class` and list-valued `function`
and `usage`. Codes may be plain strings or QualifiedCode records containing a
value and explicit code_space. The adapter preserves values without scalar/list
coercion, code conversion or remote code-list resolution. Existing LandUse and
CityFurniture scalar conventions remain intact.

TrafficArea and AuxiliaryTrafficArea use scalar string `class`, `function`, `usage`
and `surface_material` in this deliberately restricted region contract.
WaterSurface's optional `water_level` is a string classification, **not a numeric
elevation**. TransportationMarking and TransportationHole preserve scalar metadata
without inventing further domain rules. Regions keep their assignments and source parent/children grouping, normalized to
native parent indices. This preserves declared grouping without inferring a physical
space hierarchy; same-theme host rules are schema-checked. Unsupported cross-theme
owner vocabulary fails even with schema bypass.

Boundary aliases are explicit:

| CityJSON | DTCC | Where |
|---|---|---|
| `surfaceMaterial` | `surface_material` | TrafficArea, AuxiliaryTrafficArea |
| `waterLevel` | `water_level` | WaterSurface |
| `averageHeight` | `average_height` | PlantCover |

PlantCover average_height is a nonnegative height in metres under this DTCC
profile. It is not calculated from bounds. Other existing building/plant aliases
remain documented in the [model contract](model-contract.md). Ambiguous duplicate
or reserved spellings fail; arbitrary source metadata is not globally renamed.

CityJSON 2.0.2 section 3.3 permits scalar semantic-object attributes, while a
transportation example elsewhere uses an array for surfaceMaterial. This adapter
follows the explicit scalar rule and rejects array/object/null region values
without flattening them. Feature-level attributes may retain structured values.
Section 3.3 also omits Waterway from the transportation semantic-surface owner list;
this milestone conservatively excludes Waterway regions. These are documented
adapter boundaries, not claims that the wider standard has no other interpretations.
See [CityJSON 2.0.2](https://www.cityjson.org/specs/2.0.2/).

The CityGML concepts inform the crosswalk rather than generate a parallel Python
object hierarchy. See the [CityGML 3.0 conceptual model](https://docs.ogc.org/is/20-010/20-010.html)
and the [terrain/transportation comparison](terrain-transportation-crosswalk.md).

## Reproducible acceptance workflow

Run from the repository with its existing environment:

```sh
PYSTOW_HOME=/private/tmp/dtcc-profile-pystow ../venv/bin/python -m sandbox.model_profiles.exterior_city_example /private/tmp/dtcc-exterior-city-example
PYSTOW_HOME=/private/tmp/dtcc-profile-pystow ../venv/bin/python -m sandbox.model_profiles.dem_example /private/tmp/dtcc-exterior-dem-example
```

The [exterior example](../../sandbox/model_profiles/exterior_city_example.py)
combines the existing building/part, park, individual plant, bench and TIN fixtures
with all six new feature types. It checks direct Python access, exact native and
canonical package payloads, package DatasetContext, IDs/attributes, representation
types/roles/LoDs, region assignments and geometry after strict CityJSON exchange.
Coordinates permit only the declared half-grid-cell quantization tolerance. An
invalid road region is rejected under schema bypass without changing an existing
file. Negative average vegetation height is rejected at the default native and
strict CityJSON boundaries while preserving both existing files; the example also
exercises explicit semantic bypass in memory. Each example writes a concise
`report.json` alongside its artifacts.

The source is synthetic tutorial data, not a validation corpus for every supported
external producer. DatasetContext belongs to the package manifest; standalone
native payloads do not include it. CityJSON does not carry DTCC schema selection;
strict import selects the current default. CityObjects map order has no external
meaning and is not a native child-order preservation claim.

## Finite exclusions

This milestone defers interiors and installations; bridges, tunnels and other
constructions; CityObjectGroup and additional relationship mappings; appearances,
textures and geometry instancing; MultiSolid and CompositeSolid; transport Solid,
lane/routing topology; nontriangular or non-TIN CityJSON relief; raster conversion;
full addresses, temporal/versioning concepts and arbitrary supplier mappings.
There is no CityGML XML adapter or full CityGML/CityJSON conformance claim.

The native model still permits unfamiliar semantic URIs as generic data, with
native integrity checks. Supporting an additional strict external mapping or
claiming new domain rules requires an explicit future scope decision; the wider
CityGML inventory is not an implicit backlog for this milestone.

## Recorded verification (12 September 2026)

The model/I/O and affected terrain, semantic meshing and reprojection regression
run passed **717 tests, with 1 skipped**. Independent subagent reviews found no
remaining blocker in the bounded adapter and DEM workflows. Three focused
metadata tests passed after the final schema-configuration guard. The final wheel
was installed without new dependencies, and both examples passed from outside the
source repository; its installed schema and backend match the source exactly. The synthetic source
and export also pass the cached official CityJSON 2.0.2 JSON schema; this is
structural evidence, not geometric or complete standards certification.

The unchanged, pinned real 3DBAG tile was rerun under schema 0.9.0: 1,110 buildings,
5,550 qualified elevation records, exact native/package state and strict CityJSON
round trips with maximum coordinate error 0.000375 m. Missing vertical reference
was rejected without replacing an existing file. See the reproducible
[3DBAG example](../../sandbox/model_profiles/three_d_bag_example.py).

Warm in-memory medians on this machine, measured without concurrent heavy jobs:

| Input | Measurement | Seconds |
|---|---|---:|
| Repeated synthetic exterior fixture: 1,300 features, 1,700 representations | Semantic validation | 0.063 |
| Same exterior fixture | Encode with validation / bypass | 0.166 / 0.102 |
| Same exterior fixture | Decode with validation / bypass | 0.231 / 0.170 |
| Real enriched 3DBAG tile | Semantic validation | 0.687 |
| Real enriched 3DBAG tile | Validated encode / decode | 1.977 / 3.064 |
| 250,000 ground points, 251,001 raster cells | Raw grid backend / hardened Raster builder | 0.0054 / 0.0151 |

The mixed fixture repeats the same source coordinates and exercises metadata and
representation costs, not spatial variety. Timings exclude disk I/O and cold
schema initialization; small differences across runs are not regression claims.
The DEM measurement includes the extra support pass needed to distinguish valid
zero elevations from missing cells. The real native payload is 27,354,026 bytes;
the longer semantic namespace adds about 0.50 MB versus the prior 0.7 example,
without changing the numerical arrays or Protobuf wire version.
