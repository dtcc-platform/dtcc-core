## Benchmarks

Long-running benchmark tools live in the top-level `benchmarks/` directory and
are not part of normal CI. The active benchmark suite is the manifest-driven
runner:

- `cd benchmarks`
- `./bench list`
- `./bench list suites`
- `./bench list datasets`
- `./bench list cities`
- `./bench run smoke --dry-run`
- `./bench run smoke`
- `./bench run smoke --show-output`
- `./bench run smoke --save-artifacts`
- `./bench run grid --city stockholm`
- `./bench run grid --city stockholm --dataset city_surface_mesh`
- `./bench run sweep --city lund`
- `./bench run survey --dry-run`
- `./bench run survey`
- `./bench report runs/<run-id>`

The runner uses dataset names:

- `city_footprints`
- `terrain_surface_mesh`
- `city_flat_mesh`
- `city_surface_mesh`
- `city_volume_mesh`

The benchmark suites are:

- `smoke`: small live sanity check across all datasets
- `regression`: center grid tile across all cities and datasets
- `sweep`: one-axis parameter sweeps for `city_surface_mesh`
- `grid`: full 10x10 grid survey for one explicit city
- `stress`: fine-raster and dataset-specific fine-mesh center grid tiles across
  all cities
- `survey`: multi-hour flat/surface robustness and parameter-envelope survey

A benchmark task is one spatial case run through one dataset and one
scenario. For example, `./bench run grid --city stockholm` is `100` spatial
cases times `5` datasets times `1` scenario, for `500` tasks.

The city catalog is the only location catalog. Every city has a center point,
and every 10x10, 500 m grid is computed from that center. Tile `056` is the
standard center tile for city-grid cases.

Benchmark scenarios use normalized parameter names, such as
`raster_cell_size_0.5`, `max_mesh_size_1`, and `bbox_size_m_500`. The runner
maps those normalized names to dataset-specific arguments where needed.
`./bench list scenarios` shows the default parameter set and marks `baseline`
as the default scenario.

The stress suite intentionally uses dataset-specific mesh-size limits. Surface
mesh stress includes `max_mesh_size_1` and `max_mesh_size_2`. Volume mesh
stress uses `max_mesh_size_5` as its smallest mesh-size scenario; larger 3D
volume meshes below that size are treated as a separate capacity/performance
challenge rather than part of the routine benchmark gate.

The survey suite is the default "big but bounded" status benchmark for this
pipeline. It runs:

- full 10x10 grid baseline for `city_flat_mesh` and `city_surface_mesh` across
  all benchmark cities;
- center-grid flat/surface sweeps for `max_mesh_size`, `raster_cell_size`,
  `min_building_detail`, and `min_building_area`;
- city-center bbox-size sweeps for flat/surface.

By default this is 2400 tasks. Use `--city <name>` to reduce it to 240 tasks
for one city, or `--dataset city_surface_mesh` / `--dataset city_flat_mesh` to
run only one of the two mesh families.

The active suite no longer compares Spade or Triangle as benchmark dimensions.
It uses the default DTCC meshing path exposed through the dataset APIs.

Run outputs are written under:

- `runs/<run-id>/manifest.json`
- `runs/<run-id>/environment.json`
- `runs/<run-id>/results.json`
- `runs/<run-id>/summary.md`

Worker stdout and stderr are saved under `runs/<run-id>/tasks/<task-id>/`.
Use `--show-output` to also echo worker output in the terminal while the run is
active.
Use `--save-artifacts` to write successful dataset outputs under
`runs/<run-id>/tasks/<task-id>/artifacts/`; the artifact paths are recorded in
`results.json` and `summary.md`.

Every completed run prints a status summary and a full per-task result table.
The same summary is saved to `runs/<run-id>/summary.md`. Completed runs also
record `started_at`, `finished_at`, and `elapsed_seconds` in `manifest.json`
and show total wall-clock time near the top of `summary.md`.

Some data availability or local data-read misses are reported as warnings
instead of hard failures. For example, `lidar_coverage` means the requested
bounds are outside available lidar coverage, and `lidar_cache` means the runner
could not read a lidar cache/tile payload. These tasks remain visible in
`results.json` and the summary table, but they do not count as conditioning or
meshing failures and do not make the benchmark command fail by themselves.
