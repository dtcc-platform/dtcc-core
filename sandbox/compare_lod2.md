# Höganäs LOD2 comparison

Run `python sandbox/compare_lod2.py` from an installed dtcc-core environment with
its compatible dtcc-mesher build and Shapely 2.1 or newer. Edit the **CAPS constants
at the top of the script**, especially `DATA_DIR`, `OUTPUT_DIR`, and `BOUNDS`.
LM footprints and LiDAR are downloaded through the DTCC service and its normal
cache; tile lookup requires network access even when the LAZ file is cached.

The default 300 × 300 m block covers residential roofs and a cluster of manually
drawn building parts. Bounds are in the supplied data's **EPSG:3008** coordinates,
not LM's EPSG:3006. Buildings crossing the box boundary are excluded, rather than
cut into partial roofs. The common terrain raster must cover the box plus
`TERRAIN_MARGIN`. The 1 km LAS tile contains approximately 76.6 million points;
increasing the area can make dense-cloud fitting substantially slower.

## Outputs

The default destination is `sandbox/output/compare_lod2/` (ignored by Git).
Open the selected files together in ParaView, click **Apply**, and toggle visibility.
The default three methods produce seven files:

| File | Geometry |
|---|---|
| `reference.vtu` | Supplied 3D roof/wall and consolidated building models |
| `catalogue_lm.vtu` | Default catalogue reconstruction from LM points |
| `catalogue_local.vtu` | Default catalogue reconstruction from local Hexagon points |
| `ransac_lm.vtu` | RANSAC reconstruction from LM points |
| `ransac_local.vtu` | RANSAC reconstruction from local Hexagon points |
| `roofer_lm.vtu` | Standalone Roofer reconstruction from LM points |
| `roofer_local.vtu` | Standalone Roofer reconstruction from local Hexagon points |

All files use absolute EPSG:3008 XY and RH2000 Z in metres, with **no recentering or
vertical offset**. These are building skins for inspection, with open bases,
not volume meshes. The common terrain is used during meshing but excluded from
the VTUs. Catalogue and RANSAC call the normal isolated LOD2 surface pipeline.
Roofer runs as a separate executable; no reconstruction algorithm is reimplemented
in this script or added to the dtcc-core API.

Generated meshes have `building` and `is_roof` cell fields. The reference has
`reference_id` and `source_kind` (0: original roof; 1: original wall;
2: consolidated model). IDs are local to each inventory; they do **not** assert a
one-to-one correspondence between LM and municipal building records.

`comparison.json` records configuration, input counts, reference metadata,
reconstruction/fallback counts, and roof error metrics. Its `status` is `complete`
only after all selected variants and comparison finish; use that manifest to distinguish
current outputs from files left by an interrupted run. `*_audit.json` files contain
the pipeline's detailed fits, rejection reasons and meshing diagnostics.
`reference_sources.json` records source attributes and any omitted faces.
`roof_errors.csv` gives errors per LM footprint and method/cloud combination,
including the official LM `objektidentitet` alongside the pipeline model ID.

## Roofer setup

Roofer is included by default. The script needs its standalone CLI executable,
in addition to the installed dtcc-core Python environment; a Python Roofer
package is not needed. Put `roofer` on `PATH` or set its absolute path at the
top of the script:

```python
METHODS = ("catalogue", "ransac", "roofer")
ROOFER_EXECUTABLE = "/absolute/path/to/roofer"  # Or "roofer" when on PATH.
```

`METHODS = ("roofer",)` runs only Roofer with both clouds. Metrics always cover
only the selected methods in that run; the manifest identifies current outputs.
The default runs catalogue, RANSAC, and Roofer with both point clouds. Roofer is
a requirement of this comparison script's default run, not a dtcc-core package
dependency.
The adapter is tested with standalone Roofer **1.0.0 (bb2a85a)**. The executable
is checked before downloading data or overwriting the manifest.

Install the tested version using Roofer's
[official installer](https://innovation.3dbag.nl/roofer/getting_started.html).
On macOS or Linux, run:

```bash
curl -fsSL https://raw.githubusercontent.com/3DBAG/roofer/refs/heads/main/distribution/install.sh | ROOFER_VERSION=1.0.0 sh
export PATH="$HOME/.local/bin:$PATH"
roofer --version
```

The installer chooses the binary for your machine, installs its files under
`~/.local/share/roofer`, and creates the launcher `~/.local/bin/roofer`.
Add the `export PATH` line to `~/.zshrc` (zsh) or `~/.bashrc` (bash) to keep it
for future terminals. Alternatively, set `ROOFER_EXECUTABLE = "~/.local/bin/roofer"`
in the script; then no PATH change is necessary. On Windows, follow the official
instructions to download and extract the release and configure its executable.

From the dtcc-core repository root, run:

```bash
.venv/bin/python sandbox/compare_lod2.py
```

If Roofer cannot be found, the script prints these setup instructions and stops
before downloading data or overwriting outputs. The data directory must exist,
and network access is required for the LM tile lookup. The existing local Python
environment already contains the comparison's Python dependencies, including
Shapely 2.1.

Each building receives its original LM outline and the same extracted roof
observations offered to DTCC. Roofer's footprint simplification and terrain
clipping are disabled. **Roofer is not transferred to DTCC's cleaned outlines**;
this is a comparison of standalone pipelines, including their boundary handling.
Its input LAS labels observations as class 6 solely for its interface; this is not
new evidence that vegetation or facade points have been classified correctly.
LAS coordinates and Roofer's CityJSON output use millimetre precision. The density
ceiling is raised to 10,000 points/m², with the other reconstruction settings left
at Roofer defaults. Common raster ground heights are supplied explicitly.

`ROOFER_TIMEOUT_SECONDS` limits each process (default 90 seconds). Roofer's own
flat-extrusion time limit is two thirds of this, allowing it to finish writing a
fallback. Its `rf_extrusion_mode` and other attributes are recorded in the audit.
A `reconstructed` outcome means Roofer's standard mode, which can still have a
flat roof. Skipped, failed, timed-out and untriangulatable buildings are recorded
as missing outputs, with **no DTCC substitution**. Invalid geometry rejects the
whole building rather than silently omitting a roof face. Zero successful
buildings produces an empty VTU and zero coverage, not a misleading zero error.

Raw CityJSONSeq, per-building LAS/GPKG inputs, commands and logs remain under
`roofer_lm_raw/` and `roofer_local_raw/`. Each invocation uses new subdirectories
so old successful geometry cannot be mistaken for output from a failed run.
The VTU importer applies CityJSON transforms, includes roof/wall surfaces from
all LOD2.2 building parts, and excludes ground caps to match DTCC's open skins.
Roofer's mesh is exported directly after polygon triangulation; it does not use
the DTCC wall-refinement budget. Triangle counts therefore also reflect meshing
and footprint-handling differences, not just roof reconstruction quality.

### Validation on the default block

On the default 82-building block, Roofer 1.0.0 produced 79 building skins from
LM observations (2,134 triangles) and all 82 from local observations (7,085
triangles). The LM run had one building without observations and two without
returned roof geometry. One standard-mode LM reconstruction (cell field
`building=25`) contains roof vertices below ground, reaching Z = −6.498 m;
the export preserves this defect for inspection. Successful execution is not a
guarantee of a physically valid roof.

On 16,611 common one-metre grid samples, mean absolute roof-height errors were:

| Method | LM cloud | Local cloud |
|---|---:|---:|
| Catalogue | 2.510 m | 2.296 m |
| RANSAC | 2.547 m | 2.173 m |
| Roofer | 1.937 m | 0.291 m |

These measurements include each pipeline's fallback and boundary behavior and
the reference/inventory differences described below. They are observations for
this block, not a general ranking. The generated manifest and CSV contain the
full metrics and per-building errors.

The RANSAC values include centered plane fitting and supported plane-subset
assembly. On the same footprints and observations, the before/after results are:

| RANSAC metric | LM before → after | Local before → after |
|---|---:|---:|
| LOD1 fallback buildings | 39 → 26 | 56 → 31 |
| Retained non-flat fits | 36 → 51 | 11 → 38 |
| Mean absolute roof error | 2.707 → 2.547 m | 2.400 → 2.173 m |
| Building-skin triangles | 2,044 → 2,262 | 1,923 → 2,342 |

The same 16,611 reference samples were used before and after. In the dense-cloud
run, projected roof area with slope above 5° grew from 0.53% to 24.97%. Of 64 LM
footprints with reference samples, 20 improved in mean error, 6 worsened, and 38
changed by no more than 1 cm. The default block therefore shows substantially
fewer flat fallbacks but only a modest aggregate error reduction; complex roofs
remain unsupported, and Roofer's dense-cloud result is still more accurate.

The saved validation run retains its baseline in `before_ransac_improvements/`
and paired error tables, mesh checks, and matched visual examples in
`ransac_final/`, under the output directory. These are generated local artifacts,
not bundled input data. The normal script reruns all configured methods with
the current reconstructor.

## Coordinates and provenance

The supplied `läsmig.txt`, shapefile PRJs and LAS header identify SWEREF99 13 30
(EPSG:3008) and RH2000 (EPSG:5613). The raster has the same horizontal projection.
The script checks these definitions. The Esri raster WKT rounds the ellipsoid
flattening; its normalized projection agrees with EPSG:3008.

LM downloads use SWEREF99 TM (EPSG:3006). Only their XY coordinates are transformed.
The cached LM LAZ header in this test has no CRS VLR; the datum follows the
[Lantmäteriet laser-data specification](https://www.lantmateriet.se/globalassets/geodata/geodataprodukter/hojddata/pb_laserdata_nedladdning_skog.pdf),
which specifies SWEREF99 TM and RH2000. A ground-point comparison against the
supplied raster provides an independent alignment check, not a fitted vertical
correction. Absolute roof heights are retained.

The reference is **mixed-quality reference data**, not uniformly hand-surveyed
truth. Source attributes identify image-matched Terratec roofs, CityEngine
models from building-permit drawings, and manually drawn consolidated models.
The supplied notes identify the dense cloud as Hexagon LiDAR from 2025. Some
models elsewhere in the dataset were drawn from that LiDAR, so they are not
independent validation of it. The seven consolidated parts in the default block
are labelled as drawn in ArcGIS Pro from permit drawings.

The notes explicitly say the reference walls extend 1–2 m below ground. Do not
place these models on terrain or add ground elevation to their absolute Z.
Features marked `Aktualitet=0` are excluded. Original features sharing `BDG_ID`
are grouped; consolidated features retain their separate identities. Original
and consolidated geometry can overlap; the height comparison uses the upper
surface at each sample, without guessing which source should replace the other.

Reference Multipatch polygons, including vertical walls, holes and nested
collections, are triangulated in their dominant projection using original XYZ
vertices. Exact zero-area retraced edges are removed. Other invalid faces are
reported and omitted rather than silently reshaped. Nonplanar faces retain their
vertices; their interiors depend on triangulation. In the default block two
self-crossing wall faces are omitted, and 19 faces deviate from planarity by more
than 1 cm (maximum approximately 9.5 cm).

## Fair comparison and limitations

All variants start with the same complete LM footprints of at least 15 m² and
use the same supplied terrain raster. Catalogue and RANSAC share cleaning
settings; the Roofer adapter preserves the original outlines as above. Point extraction removes ground/water through the standard
pipeline; LAS noise classes 7/18 and locally flagged withheld points are excluded.
Class 12 points and points flagged as overlapping remain. Each method receives the full extracted observations
from its cloud, with the standard statistical outlier removal. The reference
geometry is never used for fitting or selecting roof models.

Errors are vertical differences on a fixed 1 m XY grid inside the LM footprints,
using the upper reference surface and generated roof triangles. MAE, RMSE and
signed bias use **the identical samples covered by the reference and all selected
outputs**. Adding a method with missing roofs can shrink that shared sample set;
compare the recomputed metrics together, not against numbers from an earlier run. Counts also show reference coverage and pairwise generated coverage;
missing support is not assigned zero error. Roofs that fall back to LOD1 remain
in the comparison and are counted in the audits. These area-weighted errors
include shape, height, footprint and acquisition-date differences; they do not
measure wall quality or prove the reference is correct. Inspect the VTUs and
per-building CSV alongside aggregate metrics, especially where inventories differ.
