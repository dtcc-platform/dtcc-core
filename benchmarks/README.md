# City cleaning and meshing benchmarks

Use one runner from the dtcc-core checkout:

```sh
.venv/bin/python benchmarks/bench quick
.venv/bin/python benchmarks/bench survey
.venv/bin/python benchmarks/bench sweep
```

| Mode | What it tests | Default scope |
| --- | --- | --- |
| `quick` | Basic operation across different cities | One central 500 m tile in each of 10 cities |
| `survey` | Geographic robustness at fixed default parameters | 100 tiles per city, 1,000 tiles total |
| `sweep` | Quality and cost when one parameter changes | Central tiles, baseline plus parameter variations |

Defaults are **flat and surface city meshes**: quick has 20 tasks, survey has
2,000. Surface and volume meshes default to **LOD1**. Survey does not contain
parameter sweeps. Start with `--city lund` for a bounded first run.

## LOD2 and roof evidence

Select `--lod lod2` for surface or volume meshes. Catalogue is the default roof
method; `--lod2-method ransac` selects the alternative. Flat and terrain-only
meshes have no building LOD selection.

```sh
# All ten central tiles, catalogue roofs, with VTU surface meshes.
.venv/bin/python benchmarks/bench quick --dataset city_surface_mesh --lod lod2 --save-artifacts --output /tmp/lod2-catalogue

# The identical cleaned buildings, terrain and individual roof observations.
.venv/bin/python benchmarks/bench quick --dataset city_surface_mesh --lod lod2 --lod2-method ransac --phase meshing --input /tmp/lod2-catalogue --save-artifacts --output /tmp/lod2-ransac

.venv/bin/python benchmarks/bench compare /tmp/lod2-catalogue /tmp/lod2-ransac
```

Cleaning retains each source building's original footprint. Roof-point
extraction uses those original polygons for both isolated buildings and parts
of merged groups; the cleaned outlines define the meshing domain. Each part
keeps its own points, ground level and height. A merged parent cloud does not
replace its parts' evidence. Empty source clouds are recorded explicitly.

`metrics.mesh_input.roof_evidence` records source IDs, point counts, ground
levels and SHA-256 fingerprints of the original polygons and point arrays.
The prepared `city.dtcc` retains these observations. Saved LOD2 replay verifies
the fingerprints before fitting; altered or incomplete evidence fails clearly.
Prepared artifacts without this provenance are not reused for LOD2:
the saved cleaning is reused and observations are acquired again. Meshing
metrics include retained nonflat fits and block fallbacks. Builder exceptions
retain detailed fitting and fallback reasons in `stage_audit`; worker timeouts
retain completed input phases.

## Automated checks

Ordinary `pytest` runs offline acceptance tests through
`benchmark_datasets.run_dataset`, the same task entry point used by `bench`.
They exercise actual cleaning, roof extraction, catalogue and RANSAC meshing,
and exact saved-evidence replay on a small hip/gable group and isolated shed.
Only input acquisition is replaced with a deterministic fixture. CI also runs
the volume version and the LOD2 boundary/tetrahedron tests after installing
the volume extra and requiring TetGen to be available.

The real ten-city suite needs provider/cache access and can take many minutes
(each task has a 420-second timeout), so it is opt-in rather than part of every
ordinary test run:

```sh
DTCC_RUN_QUICK_BENCHMARKS=1 .venv/bin/python -m pytest tests/common/test_quick_benchmarks_live.py -v
```

This wrapper calls `bench quick` for catalogue and then RANSAC with saved-input
replay. It requires a nondegenerate mesh with retained nonflat roofs for every
tile and identical roof evidence between methods. All task logs and metrics
remain in pytest's temporary directory. The manually dispatched
`quick-lod2-benchmarks` GitHub Actions workflow runs the same wrapper and uploads
the inputs, logs, metrics and meshes, including when it fails.

## Two phases

By default each task prepares raw footprints, cleans them, then meshes the
cleaned result. Each phase has its own quality measurements, timing and report.
The high-level dataset APIs still clean automatically.

```sh
# Footprints only: no LiDAR, terrain or height estimation.
.venv/bin/python benchmarks/bench quick --city lund --phase cleaning --output /tmp/lund-clean

# Mesh exactly those saved cleaned outlines; do not run cleaning again.
.venv/bin/python benchmarks/bench quick --city lund --phase meshing --input /tmp/lund-clean --output /tmp/lund-mesh

# Full cleaning + flat meshing over one city grid.
.venv/bin/python benchmarks/bench survey --city stockholm --dataset city_flat_mesh
```

`--phase cleaning` selects `city_footprints` and produces one task per spatial
case/scenario. Omit `--dataset` for this phase. Meshing-only requires `--input`;
there is no implicit fallback to downloading new footprints or cleaning again.
The selected tasks must have matching bounds and cleaning parameters in the
saved run. Missing or ambiguous inputs fail clearly before execution.

The cleaning phase saves the raw city, the cleaned city returned by
`City.clean_footprints` (with its cleaning record) and the cleaned footprints.
Flat meshing needs only the requested domain bounds and the cleaned city.
Surface/volume meshing also needs terrain and building heights. These are
prepared separately for the **saved cleaned buildings**, preserving their IDs
and order, and saved as a city artifact. A later meshing run reuses that artifact
when the dataset and preparation parameters match. Otherwise it acquires the
additional LiDAR inputs and records that fact. Input preparation has its own
status and time; its failure does not erase completed cleaning results.

The cleaning handoff is the cleaned city. Its cleaning record (scale and a
fingerprint of building IDs and footprints) lets the mesh builders validate and
consume it without cleaning again; setting heights does not invalidate it. Height attachment, domain clipping,
ground-region construction and meshing remain downstream operations. Existing
mesh datasets compose the same builder stages automatically.

Each dataset task currently owns its own input and cleaning result. A combined
flat/surface run therefore cleans separately for those tasks. To compare mesh
outputs from exactly one shared cleaning result, run cleaning first and use
`--phase meshing --input` as above.

## Selection and inspection

```sh
.venv/bin/python benchmarks/bench list
.venv/bin/python benchmarks/bench list cities
.venv/bin/python benchmarks/bench list scenarios
.venv/bin/python benchmarks/bench survey --city lund --dry-run
.venv/bin/python benchmarks/bench sweep --city lund --dataset city_surface_mesh --scenario max_mesh_size_2
```

- `--city NAME`: restrict to one city; default is all 10.
- `--dataset NAME`: select an output; repeat to select several. Available:
  `city_footprints`, `city_flat_mesh`, `city_surface_mesh`, `city_volume_mesh`,
  `terrain_surface_mesh`. The latter is a separate terrain-only dataset check;
  its timing includes input preparation and it has no footprint-cleaning phase.
- `--scenario NAME`: select settings within `sweep`; repeat if needed. Cleaning
  sweeps omit mesh/raster parameters; flat sweeps omit raster parameters; routine
  volume sweeps omit mesh sizes below 5 m.
- `--dry-run`: print exact bounds, parameters and task counts without execution.
- `--show-output`: stream worker logs; logs are always saved.
- `--save-artifacts`: also save meshes. Raw/cleaned inputs, metrics and prepared
  terrain/height inputs are always saved so meshing can be replayed.
- `--output PATH`: choose a **new** run directory. Default:
  `benchmarks/runs/<timestamp>_<mode>`. Existing directories are never overwritten.

Cities: Lund, Stockholm, Gothenburg, Malmo, Uppsala, Linkoping, Orebro, Vasteras,
Helsingborg and Norrkoping. The central case is tile 056. The city catalog is the
single source for all grids, bounds and parameters.

## Results and quality

Each run contains `manifest.json`, `environment.json`, `results.json`,
`summary.md`, and a directory under `tasks/` for each case/dataset/scenario:

```text
task.json
stdout.log
stderr.log
cleaning/
    raw.dtcc              # Original source buildings and IDs
    cleaned.dtcc          # Cleaned city: blocks with the original buildings as parts
    footprints.geojson    # Regions handed to meshing (a block per roof level) with building index, ID, roof
    metrics.json
meshing/
    city.dtcc             # Cleaned buildings with prepared terrain/heights
    input.json            # Preparation metrics and LOD2 source-evidence fingerprints
    metrics.json
artifacts/                # Optional exported mesh
```

The footprint artifact explicitly records **EPSG:3006, metre coordinates** in
its metadata; do not interpret it as WGS84 GeoJSON. Loading checks its version,
CRS, bounds and cleaning parameters, that the cleaned city's cleaning record
is valid, and that the saved footprints match the cleaned buildings. Mesher
input validation remains active. Results from before phase separation do not
contain this handoff and cannot be used for meshing-only runs, and artifacts of
another format version are rejected; rerun their cleaning phase.

| Phase | Measurements |
| --- | --- |
| Cleaning | Input/output/vertex/hole counts, invalid input/output counts, unrepresented sources, added/removed area, maximum boundary displacement, overlap, admissibility, minimum feature size, fidelity, budget and achieved budget, stage contract, cleaning time |
| Meshing | Vertex/element/region counts; min/mean/max element quality, aspect ratio, radius ratio, edge ratio and skewness; 1st-percentile element quality; counts below quality 0.02 and of degenerate cells; stage contracts; meshing time |

Area and displacement compare coverage unions, so overlapping raw polygons do
not inflate area. Invalid raw polygons use GEOS `make_valid` for these
measurements; the original raw geometry and invalid-input count are retained.
Unrepresented sources and changes in holes/area are observations, not automatic
failures: declared scale rules intentionally remove some geometry.

Cleaning metrics include independent
[contract checks](../dtcc_core/builder/cleaning/contract.py). The columns
`Admissible`, `mfs m`, `Fidelity`, `Budget m` and `Achieved m` describe the final
cleaned geometry **before area selection**, in the terms of the
[contract](../docs/design/footprint-cleaning-contract.md): membership in
`A_delta` under the union topology profile, the minimum feature size of the
essential boundary graph (capped at delta, so an admissible result shows delta),
the fidelity status, the fidelity budget used (scale / 2 by default, or the
retry budget 3 scale / 4 when a merge group needed it) and the achieved budget,
the smallest epsilon the result would also satisfy.
`Selected out` counts regions excluded afterward by `min_building_area`.

JSON separates `before_selection_contract`, `selection`, and the selected
handoff's `geometric_contract`. The usual counts and added/removed areas still
compare raw inputs with that selected handoff. Selection reports excluded area
and cleaned-building indices explicitly; it never changes the cleaning
reference. Dataset and benchmark `min_building_area` parameters apply after
cleaning, as in the mesh builders.

A `borderline` fidelity result is inconclusive, not a pass. By default, residual
separation defects alone produce a prominent warning and allow meshing to be
attempted; topology, fidelity, source attribution and mesher-input safety remain
mandatory. `Admissible` remains `no` and the geometric contract remains `fail`;
the pipeline reports `warning`, not geometric conformance. Smaller elements and
greater mesh cost are possible. Warnings survive selection and saved-stage
replay. Older reports are not relabeled and show `-` for missing observations.

Replayed cleaning is marked `reused`, without a new cleaning time or claimed
new cleaning measurements. Consult its source run for those measurements.
Meshing input preparation records `bounds`, `saved`, or `provider/cache`.
Completed input artifacts and phase metrics survive worker timeouts or crashes,
so another method can still use the same prepared observations.

Geometry/pipeline failures and quality-evaluation failures fail the task.
Stage quality warnings and degenerate final cells are explicit warnings; data
coverage/cache misses retain their existing warning classification. Warnings
remain nonfatal to the command. **This is a measurement campaign, not yet a
new numerical quality acceptance policy.** Review quality changes as well as
success counts; a zero exit code does not mean every mesh is acceptable.

## Review, compare and rerun

Plot saved footprints without downloading data or rerunning either phase:

```sh
.venv/bin/python benchmarks/bench plot /tmp/lund-clean --city lund
.venv/bin/python benchmarks/bench plot /tmp/lund-mesh --city lund --dataset city_surface_mesh --show
```

Each selected task produces a PNG under `<run>/plots/`: raw footprints, cleaned
footprints, and an overlay showing added coverage in green and removed coverage
in red. All panels have matching extents. The difference panel uses the same
`make_valid` interpretation of invalid raw polygons as the area measurements.
`--show` also opens each figure interactively with linked pan/zoom; close it to
advance. By default the command only saves images, including on headless systems.
Use `--task <exact-id>` to inspect one survey tile/scenario, or `--output <dir>`
to choose the image directory. Repeating the command replaces those generated
PNGs. Tasks with no saved cleaning result are skipped; an empty selection fails
clearly. A combined flat/surface run has one image per dataset task.

```sh
.venv/bin/python benchmarks/bench report /tmp/lund-mesh
.venv/bin/python benchmarks/bench compare /tmp/before /tmp/after
.venv/bin/python benchmarks/bench rerun /tmp/lund-mesh --failed
.venv/bin/python benchmarks/bench rerun /tmp/lund-mesh --task 'quick:city_flat_mesh:city_grid:lund:056:baseline'
```

`compare` reports phase timing and quality deltas, not just overall runtime.
Use matching input geometry and settings for meaningful before/after comparison.
For LOD2, compare exported roof/wall/terrain quality against both saved production
baselines, separately from the current minimal geometry:

```sh
.venv/bin/python -m benchmarks.lod2_comparison NEW_RUN \
  --production benchmarks/runs/lod2-principled-comparison-20261018-v7 \
  --original-catalogue benchmarks/runs/lod2-central-catalogue-20261018 \
  --original-ransac benchmarks/runs/lod2-central-ransac-20261018
```

This writes `NEW_RUN/comparison/report.md`, explicit historical acceptance
targets, quality CSV, per-source geometry attribution, and matched junction
neighborhood measurements. A minimal-to-refined gain cannot clear a historical
regression. Source evidence and group triangle/point budgets are independent
checks.

`rerun` selects failed/incomplete tasks by default; `--task` selects an exact ID,
including a warning or successful case. Every rerun writes a new run directory.
A meshing-only rerun uses the saved handoff; a full/cleaning rerun reacquires raw
footprints from the provider/cache as the original command did.

The old `run` verb and `smoke`, `regression`, `grid`, `stress`, `triage` suites
have been removed. Use `quick`, `survey --city ...`, `sweep`, and `rerun`.
Historical reports in this directory document the old commands and results;
they are not today's acceptance baseline.
