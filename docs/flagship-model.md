# Delft flagship development model

The flagship is a **2 × 2 km** district of Delft containing **10,344 real buildings** and
**138 canal/water features**, with synthetic terrain and numerical fields. It includes
the Nieuwe Kerk, Oude Kerk, the TU Delft EWI tower, dense historic streets, and
the original 64-building reference neighbourhood. Detached demonstration
architecture (pavilion, artificial park, road, bench, etc.) is no longer included.

## Inspect locally

From the Core checkout, using its normal development environment:

```sh
python scripts/inspect_flagship_model.py
```

The dashboard loads `data/flagship/flagship.dtcc` and displays a building-height
map with real water outlines and landmark labels, district wind vectors/speed,
tracer concentration, and a vertical temperature section. A slider selects four
stored snapshots, 60 seconds apart. Colour scales are fixed between snapshots;
grey field samples mean masked solids/ground or missing data, not zero values.
The dashed rectangle marks the fine reference patch.

Additional views:

```sh
python scripts/inspect_flagship_model.py --view map
python scripts/inspect_flagship_model.py --view city
python scripts/inspect_flagship_model.py --view fine --frame 2
python scripts/inspect_flagship_model.py --view tetra
python scripts/inspect_flagship_model.py --view solar
python scripts/inspect_flagship_model.py --view wind --frame 2
python scripts/inspect_flagship_model.py --view dem
python scripts/inspect_flagship_model.py --view streamlines
python scripts/inspect_flagship_model.py --save /tmp/flagship.png
```

`city` is Core's rotatable native `Model.plot()` preview, with numerical domains
excluded and an explicit display budget large enough for this district. Native
polygon previews outline faces with holes; the fine reference buildings also
have triangulated previews. `fine` shows the densely sampled reference patch.
`tetra` intersects the stored tetrahedra with a vertical plane and displays their
connectivity, interpolated vertex pressure and cell tracer. It does not build a
replacement mesh or recompute the fields. Its vertical scale is exaggerated for
legibility. These are Matplotlib previews, without Twin's future volume renderer.

All views read stored geometry/arrays and work offline. `--save` renders a PNG
headlessly; interactive use requires a desktop Matplotlib backend. Cold native
loading and plotting thousands of buildings can take appreciably longer than the
small former fixture.

## Domain and real sources

- Horizontal bounds in RD New metres: **E 83,570–85,570; N 445,800–447,800**.
- CRS: **EPSG:7415**, RD New + NAP height. All main-scene geometry uses world
  coordinates and identity affines. Vectors follow projected east, north, up.
- Include complete building envelopes contained in the square. Buildings that
  straddle its boundary are omitted; solids are never cut into open shells.
- Preserve the original 64 buildings from the reference CityJSON sample, including
  their parts, all source LoDs, attributes, surface semantics and geometry. These
  take precedence over newer buildings with the same BAG IDs.
- Fill the remaining district from **15 3DBAG v20250903 tiles**, retaining footprints, the finest available source LoD and all attributes.
  Coarser alternate LoDs outside the reference patch are omitted to fit the
  native 256 MiB limit; detailed roof geometry is not simplified. This is an explicitly mixed-vintage
  fixture, not a claim that all buildings were observed on one date.
- Use BGT/PDOK `waterdeel` outlines requested in EPSG:28992, the horizontal part
  of EPSG:7415. Retain current records only (no end-registration or termination
  date), clip outlines to the square, and preserve source attributes and holes.
  The water plane at −0.5 m NAP and bed at −0.6 m are **synthetic elevations**.

The source manifest is [`scripts/flagship-sources.json`](../scripts/flagship-sources.json).
It lists the source URLs and local filenames. `building_source_files` on the
City identifies each building's input. Supplier title/contact/version metadata,
which the strict Core CityJSON profile does not accept at its root, is retained
in `source_metadata` on the City. Building geometry/attributes still use the
ordinary strict `load_3dbag` admission and its qualified elevation mapping.
The versioned 2025 tiles use the supported `b3_h_dak_*` convention; a general mapping
for later attribute conventions is not introduced by this example.

3DBAG credit: **© 3DBAG by tudelft3d and 3DGI**, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
See [attribution](https://docs.3dbag.nl/en/copyright/) and
[delivery documentation](https://docs.3dbag.nl/en/delivery/webservices/).
Water outlines: [BGT / PDOK](https://api.pdok.nl/lv/bgt/ogc/v1/collections/waterdeel).
The selected church/tower locations are labelled using public PDOK address
coordinates; reconstructed roof elevations need not match architectural heights.

## Fields and meshes

A shared analytic sampler supplies wind shear, local building wakes/deflection,
a moving vortex, vertical motion, changing wind direction, an advecting tracer
pulse, pressure, and air temperature with illustrative cooling over water.
Spatial indices restrict obstacle/wake work to nearby buildings. Footprint prisms
approximate obstacles, with a smooth finite wake support. These fields are
**synthetic illustrations**, not a fluid solution: conservation and no-slip
conditions are not claimed. Related grids, sensors and slices sample the same
function. Terrain interpolates building-base elevations and pins footprint
interiors to their bases; it is not a surveyed terrain model.

| Representation | Standard | Stress |
|---|---:|---:|
| District terrain/surface sampling | 8 m | 8 m |
| Maximum volume spacing (X, Y, Z) | 20, 20, 10 m | 20, 20, 8 m |
| Volume grid cells (X × Y × Z) | 100 × 100 × 14 | 100 × 100 × 18 |
| Tetrahedra | 840,000 | 1,080,000 |
| Longest/shortest grid-cell dimension | ≈2.02 | ≈2.60 |
| Reference-patch field samples | ≤2 m | ≤2 m |
| Reference boundary mesh triangle target | 2 m | 1 m |
| Stored snapshot times | 0, 60, 120, 180 s | same |

Both volume representations share the same vertical intervals: 14 in standard
and 18 in stress. There are 15 or 19 vertex planes respectively. Actual vertical
spacing is about 9.90 m or 7.70 m over the 138.576 m height of this domain. The
reported axis ratio describes the Cartesian cells, not a tetrahedron quality
metric. Neither the domain height nor the terrain/surface sampling is changed
when volume spacing is overridden.

Both profiles include:

- Terrain TIN and complete coarse overview TIN, classified point samples, and a
  DEM with a deliberate 3 × 3 no-data corner.
- A 2D grid with scalar runoff proxy and two-component downhill direction.
- A 3D structured grid with velocity, speed, temperature, pressure, tracer,
  boolean air mask and integer flow-zone codes.
- A **tetrahedral mesh covering the entire district volume**. Each Cartesian
  box is divided into six tetrahedra along a consistent body diagonal, sharing
  vertices and faces with neighbours. This is direct deterministic connectivity
  construction in the generator, not a call to TetGen or the urban volume mesher.
- Vertex velocity, pressure and validity on the tetrahedra; cell tracer,
  cell validity and flow-zone markers; a six-component symmetric velocity
  outer product (`xx, yy, zz, xy, xz, yz`, in m²/s²). This algebraic tensor example
  is not a physical stress measurement.
- Four snapshot Objects with district samples at 2 m and 15 m above terrain,
  fine samples at 2 m above terrain over the reference patch, a vertical section,
  steady streamlines, and eight sensor stations at stable locations.
- On the 64 reference buildings: source-derived boundary meshes with normals,
  region markers, triangle area, solar irradiance, surface temperature and
  approximate footprint-prism shadow masks. Other buildings keep their real
  source geometry without expensive dense surface subdivision.
- Water-surface temperature illustrating surface-associated scalar values.

The tetrahedral mesh **does not conform to terrain or building walls/roofs**.
It includes cells crossing or lying inside obstacles. Validity is evaluated at
sample locations: `air_mask` at vertices/grid cells and `cell_air_mask` at tetra
centres. Physical fields are NaN at invalid samples. A centre mask is not a
cut-cell geometry guarantee. A mesh suitable for CFD with fitted boundaries
would require a separate meshing workflow.

Flow-zone codes are 0 solid/ground, 1 open air and 2 wake. Times/component labels
use ordinary Object metadata and relations; the native format has no dedicated
time/tensor axes. Full 3D data is stored at t=0; subsequent snapshots store slices,
streamlines and sensors to control size. Streamlines are steady integrations at
one snapshot, not time-dependent particle trajectories.

## Reproduce and inspect artifacts

```sh
python scripts/download_flagship_sources.py
python scripts/generate_flagship_model.py
python scripts/generate_flagship_model.py --detail stress --output data/flagship-stress
# Explicit maximum X, Y, Z spacing in metres (overrides the selected profile):
python scripts/generate_flagship_model.py --volume-spacing 20 20 10 --output data/flagship-custom
```

`--volume-spacing DX DY DZ` controls **both** `air_grid` and `tetrahedra`. Values
must be finite and positive. Each axis uses `ceil(domain_length / requested_spacing)`
cells; actual spacing is reduced slightly to fit the unchanged domain exactly.
The generator prints the dimensions, actual spacing and cell axis ratio, and
records them in `flagship_sampling` and `inventory.json` (XYZ order). The six
positively oriented tetrahedra per grid box remain consistent across neighbours.
For near-cubic cells, choose equal spacings; this grows the file much faster than
refining Z alone. The existing 256 MiB native limit still applies. Requests whose
connectivity/coordinates alone exceed it fail before volume allocation; other
oversized complete models are rejected by the canonical writer.

Default paths are checkout-relative; explicit paths are working-directory-relative.
The downloader is the only network step. It downloads missing files and reuses
existing files in `data/flagship-source`. Delete a cached file to download it again.
The generator reads these files offline and validates the source data. Downloads
do not need to match historical checksums: the BGT endpoint is live and includes
response metadata that changes between requests. Water outlines may also change
over time, so retain the same cached inputs when reproducing a particular model.

Both profiles are generated and round-trip checked against the native 256 MiB
limit. Exact sizes and generation/save/load timings are in each inventory.
Rendering all detailed roofs/walls in the Matplotlib `city` view is substantially
slower than the map or field views; start with the dashboard.

Generated outputs are `flagship.dtcc`, `flagship.dtccpkg`, `inventory.json`,
`preview.png`, and a copy of this guide. Inventory records exact building/water
counts, array sizes/types, field ranges/missing values, sampling, source
preservation, round-trip results, file size and generation/save/load timings.
The package also retains Dataset Context; the bare native file contains the
City's embedded provenance but not that separate context.

The generator checks exact native and canonical package round trips, verifies
all selected source building facts after removing added meshes, checks exact 2 km XY
bounds, and verifies a rejected invalid storey edit cannot replace the output.
The focused regression checks positive tetrahedral volumes and their summed
box volume, field associations/masks/time consistency, sensors, fine sampling,
tetrahedral inspection, dashboard interaction and failed-save preservation.
These checks do not certify every supplier solid geometrically.

The flagship deliberately prioritizes a coherent real city. The removed
pavilion's LoD3 openings/interiors, display-area furniture/transport and standalone
affine specimen are no longer coverage claims for this model. Native schema
coverage tests remain separate from the flagship's visual scene.
