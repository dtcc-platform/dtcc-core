# Delft flagship development model

`flagship.dtcc` is a rich, reproducible development fixture using the standard
DTCC Model schema and ordinary Protobuf exchange. It combines 64 real buildings
with an explicitly synthetic laboratory beside the selected neighbourhood.
It is useful for developing readers, viewers, selection tools, field inspection,
semantic queries and save/load workflows. It is not survey or simulation evidence.

For the ordinary Python preview, use `dtcc.load_city(...).plot()`. The built-in
preview also supports `city.plot(representation='dem')` and
`city.plot(field='velocity')`; see the native model preview guide in
`docs/model-preview.md`. The generated `preview.png` below remains the custom
two-panel illustration from the flagship generator.

## Files

Generated in `data/flagship/` (ignored by Git):

- `flagship.dtcc`: the complete native model, with source attribution and
  generation details embedded in `city.attributes['flagship']`.
- `flagship.dtccpkg`: the same model plus Dataset Context and package integrity.
- `source-neighbourhood.city.json`: the original selected CityJSON features with
  a compacted vertex table. This is the real-data input subset, not an export of
  the enriched native model.
- `inventory.json`: counts, geometry types, LoDs, region vocabulary, field
  locations/dtypes/shapes/units, checksum and verification results.
- `preview.png`: an overview derivative and pavilion detail. It selects one
  representation per object and simplifies polygon holes/interior shells for
  display. Overview heights are exaggerated ×2; pavilion detail has true
  proportions and triangulated openings. Complete geometry remains in the file.
- `README.md`: a copy of this guide to accompany the artifacts.

## Contents

The real buildings retain source IDs, BuildingPart containment, footprint rings,
LoDs 0/1.2/1.3/2.2, Solid shell topology, roof/wall/ground regions and original
supplier attributes. `io.load_3dbag` adds qualified roof codes and NAP elevation
records without inventing measured building heights. One source BuildingPart also
has a derived boundary mesh with normals, markers, true triangle areas and a
clearly synthetic irradiance field. Its original representations stay intact.

The synthetic area contains a pavilion (LoDs 0, 1 and 3), two Solid shells,
interior surfaces, wall holes with Window/Door regions and host relationships,
qualified height metadata, a park parcel, plant cover, twelve procedural trees,
a pond, road and routing graph, railway, waterway, square, bench, sensors and
vehicles. The bench demonstrates a nonidentity geometry-local transform.
A solar-panel Surface demonstrates an unfamiliar semantic URI admitted as a
generic Object: its custom classification is preserved but is not a declared
class in the standard schema.

Numerical examples include a terrain triangle mesh, classified point samples,
a DEM with deliberately missing ground samples, a 2D Grid, 3D VolumeGrid and
1,152 tetrahedra. Scalar/vector fields cover vertex, edge, face, cell, sample
and geometry associations, float32/float64 arrays and explicit units. Sensor/vehicle collections
and named ID relationships demonstrate containment versus cross references.
All sensor readings and flow/temperature/irradiance fields are authored examples,
not observations or solver outputs. Each synthetic object says so in its attributes.

Coordinates use EPSG:7415 (Amersfoort / RD New + NAP height), in metres. The
synthetic landscape is authored in that coordinate frame; it is not fitted to
real terrain. The bench's shape carries its own local-to-CRS affine. Generic
aggregate bounds must not be interpreted as transformed world-frame bounds.

## Reproduce

Use the normal Core development environment and the checkout root. No extra
dependencies are needed. Download the openly licensed source once:

```sh
mkdir -p data/flagship-source
curl --fail --location \
  https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/9-284-556.city.json \
  --output data/flagship-source/3dbag.city.json
python scripts/generate_flagship_model.py
```

After the source is downloaded, no command-line arguments are needed. From the
`scripts` directory, run `python generate_flagship_model.py`. Both invocations
read `data/flagship-source/3dbag.city.json` and write `data/flagship/` under the
checkout root, regardless of the working directory. An optional source argument
and `--output` override these defaults; explicit relative paths use the working
directory.

The generator checks SHA-256
`2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba`
before creating outputs. It uses the 64 nearest building footprint-envelope
centres to `NL.IMBAG.Pand.0503100000000030`, breaking ties by ID. The source
release is unspecified; do not relabel this cached sample as the latest 3DBAG.
Generation itself makes no network requests. Re-running intentionally replaces
the output files. Native bytes are deterministic in the same environment;
mesher/dependency changes may change derived triangulations. Inventory timing
and package creation metadata are not byte-reproducibility promises.

## Python access

```python
from dtcc_core import io

city = io.load_model('data/flagship/flagship.dtcc')  # validation is on

def walk(obj):
    yield obj
    for group in obj.children.values():
        for child in group:
            yield from walk(child)

features = {obj.id: obj for obj in walk(city)}
building = features['NL.IMBAG.Pand.0503100000000030']
footprint = building.get_geometry(lod='0').surfaces[0]
polygon = footprint.vertices            # NumPy array; also inspect footprint.holes
solid = building.building_parts[0].get_geometry(id='cityjson-2')
# In this fixture cityjson-2 is the source LoD-2.2 Solid. The additional
# flagship_boundary_mesh also has LoD 2.2, so use an explicit representation ID.

pavilion = features['synthetic-pavilion']
detailed = pavilion.get_geometry(id='detailed_shells')
window = next(r for r in detailed.regions if r.semantic_type.endswith('#Window'))
window_polygon = detailed.surfaces[window.indices[0]].vertices
host_wall = detailed.regions[window.parent]

bench_id = features['synthetic-bench-01'].id
terrain = features['synthetic-terrain']
dem = terrain.get_geometry(id='dem')
elevations, georeference = dem.data, dem.georef

tetra = features['synthetic-flow-domain'].get_geometry(id='tetrahedra')
velocity = next(f for f in tetra.fields if f.name == 'velocity')
vectors, unit, association = velocity.values, velocity.unit, velocity.association
points, temperatures = features['synthetic-sensors'].to_arrays('air_temperature')
source = city.attributes['flagship']  # retained in the standalone .dtcc
```

Use explicit representation IDs when LoD alone is ambiguous. Geometry coordinates,
field values and metadata remain ordinary Python/NumPy data. Dataset Context is
additionally available after `dtcc_core.datasets.load_model_package(...)` on the
companion package. The complete enriched model has no lossless CityJSON mapping;
strict export must reject unsupported native additions rather than discard them.

## Verification and attribution

The generator uses default schema validation for native save/load and package
exchange. It checks exact native re-encoding, exact package model/context,
preservation of all selected source facts, opening access, and rejection of an
invalid storey-count edit without replacing the valid file. The development test
also checks field dimensions and dangling-reference rejection. This does not
certify watertightness or geometric validity of the original 3DBAG tile.

Source geometry and source attributes: **© 3DBAG by tudelft3d and 3DGI**, licensed
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
See the [required attribution and terms](https://docs.3dbag.nl/en/copyright/) and
the [CityJSON sample catalogue](https://www.cityjson.org/datasets/).
Changes are the spatial subset, DTCC mapping/enrichment, derived mesh and synthetic
laboratory described above. Preserve this credit, the license link and the change
description when redistributing the file or its visual derivatives.
