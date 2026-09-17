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
2,000. Survey does not contain parameter sweeps. These are development tools,
not normal CI tests. Start with `--city lund` for a bounded first run.

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

Flat meshing needs only the requested domain bounds and cleaned footprints.
Surface/volume meshing also needs terrain and building heights. These are
prepared separately from the **saved source buildings**, preserving source
indices, and saved as a city artifact. A later meshing run reuses that artifact
when the dataset and preparation parameters match. Otherwise it acquires the
additional LiDAR inputs and records that fact. Input preparation has its own
status and time; its failure does not erase completed cleaning results.

The cleaning handoff is the final mesher-ready footprint stage, including its
existing normalization, not an earlier intermediate. Mesh builders validate
and consume it without footprint repair. Height attachment, domain clipping,
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
    footprints.geojson    # Cleaned polygons and many-to-many source mapping
    metrics.json
meshing/
    city.dtcc             # Source buildings with prepared terrain/heights
    metrics.json
artifacts/                # Optional exported mesh
```

The footprint artifact explicitly records **EPSG:3006, metre coordinates** in
its metadata; do not interpret it as WGS84 GeoJSON. Loading checks its version,
CRS, bounds, cleaning parameters, polygon geometry and source indices. Mesher
input validation remains active. Results from before phase separation do not
contain this handoff and cannot be used for meshing-only runs.

| Phase | Measurements |
| --- | --- |
| Cleaning | Input/output/vertex/hole counts, invalid input/output counts, unrepresented sources, added/removed area, maximum boundary displacement, overlap, minimum clearance, short edges, stage contract, cleaning time |
| Meshing | Vertex/element/region counts; min/mean/max element quality, aspect ratio, radius ratio, edge ratio and skewness; 1st-percentile element quality; counts below quality 0.02 and of degenerate cells; stage contracts; meshing time |

Area and displacement compare coverage unions, so overlapping raw polygons do
not inflate area. Invalid raw polygons use GEOS `make_valid` for these
measurements; the original raw geometry and invalid-input count are retained.
Unrepresented sources and changes in holes/area are observations, not automatic
failures: declared scale rules intentionally remove some geometry.

Replayed cleaning is marked `reused`, without a new cleaning time or claimed
new cleaning measurements. Consult its source run for those measurements.
Meshing input preparation records `bounds`, `saved`, or `provider/cache`.

Geometry/pipeline failures and quality-evaluation failures fail the task.
Stage quality warnings and degenerate final cells are explicit warnings; data
coverage/cache misses retain their existing warning classification. Warnings
remain nonfatal to the command. **This is a measurement campaign, not yet a
new numerical quality acceptance policy.** Review quality changes as well as
success counts; a zero exit code does not mean every mesh is acceptable.

## Review, compare and rerun

```sh
.venv/bin/python benchmarks/bench report /tmp/lund-mesh
.venv/bin/python benchmarks/bench compare /tmp/before /tmp/after
.venv/bin/python benchmarks/bench rerun /tmp/lund-mesh --failed
.venv/bin/python benchmarks/bench rerun /tmp/lund-mesh --task 'quick:city_flat_mesh:city_grid:lund:056:baseline'
```

`compare` reports phase timing and quality deltas, not just overall runtime.
Use matching input geometry and settings for meaningful before/after comparison.
`rerun` selects failed/incomplete tasks by default; `--task` selects an exact ID,
including a warning or successful case. Every rerun writes a new run directory.
A meshing-only rerun uses the saved handoff; a full/cleaning rerun reacquires raw
footprints from the provider/cache as the original command did.

The old `run` verb and `smoke`, `regression`, `grid`, `stress`, `triage` suites
have been removed. Use `quick`, `survey --city ...`, `sweep`, and `rerun`.
Historical reports in this directory document the old commands and results;
they are not today's acceptance baseline.
