# Strict triangular TINRelief exchange

Implemented 12 September 2026, using [standard schema 0.8.0](../../schemas/archive/model/0.8.0/schema.yaml)
and unchanged Protobuf wire v6. This implements the bounded terrain slice from
the [terrain/transportation crosswalk](terrain-transportation-crosswalk.md).

## Public workflow

```python
from dtcc_core import io
from dtcc_core.model import Terrain

city = io.load_city("mixed.city.json", strict=True)
terrain = next(t for t in city.get_children(Terrain) if t.id == "terrain-west")
mesh = terrain.get_geometry(lod="2.2")
print(mesh.vertices, mesh.faces)  # NumPy coordinate and triangle-index arrays

city.save("mixed.dtcc")
restored = io.load_city("mixed.dtcc")
restored.save("mixed-out.city.json", strict=True)
```

Each source TINRelief becomes the existing native Terrain with explicit semantic
URI `https://example.org/dtcc/TINRelief`. Each triangular CompositeSurface becomes
an existing Mesh, with no Surface object per triangle. Multiple terrain features
and all their ordered representations survive. Their attributes stay unchanged.
Representation IDs follow the existing `cityjson-0`, `cityjson-1`, … convention;
LoD strings such as `1.2` and `2.2` stay exact. Ambiguous geometry selection retains
the existing error behavior, so use an ID when several representations share a LoD.

Schema 0.8.0 declares TINRelief as a semantic specialization of Terrain. A generic
Terrain remains generic native data. Strict export requires the explicit TINRelief
URI and supported Mesh representations, rather than classifying every terrain or
mesh automatically. No new Python classes, wire fields or dependencies are added.
The schema describes the classification; coordinate buffers and connectivity
remain in native admission and external representation checks, not LinkML.

## Preservation boundary

| Fact | Strict mapping |
| --- | --- |
| Feature identity, attributes, semantic classification | Preserved for every terrain feature; no first-feature selection or merge |
| Representation order, LoD and CRS | Preserved using existing representation and transform authorities |
| Triangle order and winding | Preserved; no triangulation, cleaning, snapping or reorientation |
| Shared source indices within a representation | Compacted to a local Mesh vertex table in first-use order on import |
| Coordinates | Source transform is applied once; native exchange is exact; CityJSON export uses the selected quantization scale |
| Terrain surface semantics | No labels invented; source `semantics`, including empty assignments, is rejected in this first subset |
| Dataset Context | Preserved by canonical package export when supplied; standalone native ModelFile does not carry root context |

CityJSON's global vertex numbering is not a persistent native identifier. Unused
document vertices are excluded from each imported representation; the writer shares
equal quantized coordinates through its existing global indexer. Comparing a
CityJSON round trip therefore compares ordered triangle coordinates, rather than
raw local vertex numbering or duplicate-vertex identity. Native `.dtcc` preserves
the actual arrays, dtypes and connectivity exactly.

The public save/load schema defaults and `validate_schema=False` work unchanged.
Bypass skips semantic evaluation only: it cannot bypass source indices, supported
geometry, representability or quantization checks. Source extent discrepancies
retain the existing explicit `extent_policy="recompute"` mechanism.

## Rejected inputs and exports

TIN geometry must be a nonempty CompositeSurface whose polygons each contain one
ring of exactly three valid integer indices, referring to three distinct positions.
Nontriangular polygons, holes, invalid/repeated indices and unsupported semantic
assignments fail. A malformed later terrain feature fails the entire import; it
does not return a partial city.

Strict export rejects nonempty Mesh markers, normals, fields or semantic regions,
unused mesh vertices, unsupported representation IDs/roles, nested/non-global
transforms and conflicting CRSs. These facts remain usable and serializable in
native DTCC; an external mapping is required before CityJSON can preserve them.
If quantization merges vertices of a triangle, export fails before opening the
destination. The existing file remains intact on these validation failures.

The subset checks representation integrity, not a complete geometric validity
certificate. It does not certify connected/manifold terrain, triangle area,
bare-earth classification, vertical datum correctness, accuracy or suitability for
simulation. No raster-to-TIN conversion, DEM builder correction, arbitrary polygonal
CompositeSurface support, GML XML adapter or transportation mapping is included.
The older permissive terrain reader/exporter remains outside this strict contract.

## Reproduction and evidence

From the repository, run:

```sh
PYSTOW_HOME=/private/tmp/dtcc-profile-pystow ../venv/bin/python \
  -m sandbox.model_profiles.tin_relief_example /private/tmp/dtcc-tin-relief
```

The synthetic example contains a building and part, two terrain features, three
terrain representations and five triangles. It checks every representation, direct
NumPy access, exact native/package exchange, package context, strict CityJSON and
rejection of unmapped markers without file replacement. It is a reproducible
carrier example, not survey data. Focused tests additionally cover invalid later
features, unsupported mesh state and collapsing quantization.

The [implementation plan](../../.agent/plans/2026-09-12-strict-tin-relief.md) records
observed checks and limitations. Previous schema 0.7.0 is archived unchanged;
new roots select 0.8.0. No automatic schema migration or legacy wire reader is added.

Observed verification: 29 focused checks passed, followed by 690 model/I/O checks
passed and one skipped. The example report is
`/private/tmp/dtcc-tin-relief/report.json`. A separate synthetic 20,000-triangle,
10,201-vertex grid measured warm median import/export of 0.416/0.070 s over three
runs on the development Mac. Those timings include default schema evaluation but
exclude first-use initialization, JSON parsing/stringification and file I/O (the
input/output is an in-memory CityJSON dictionary). Samples are in
`/private/tmp/dtcc-tin-timing.json`; they are workload evidence, not a scaling guarantee.
