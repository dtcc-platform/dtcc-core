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
- `./bench run grid --city stockholm`
- `./bench run grid --city stockholm --dataset city_surface_mesh`
- `./bench run sweep --city lund`
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
- `stress`: fine-raster and fine-mesh center grid tiles across all cities

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

Every completed run prints a status summary and a full per-task result table.
The same summary is saved to `runs/<run-id>/summary.md`.
