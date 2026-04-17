## Benchmarks

Long-running benchmark tools live in the top-level `benchmarks/` directory and
are not part of normal CI. The main Stockholm meshing benchmarks are:

- `python benchmarks/bench_footprints.py`
- `python benchmarks/bench_mesh_2d.py`
- `python benchmarks/bench_mesh_3d.py`

They share a common core CLI:

- `--cases 15 16 54` to run one or more specific Stockholm tiles
- `--output-dir ...` to override the default output root
- `--delay ...` to throttle multi-case runs
- `--no-plots` to skip PNG generation

Default outputs are written under:

- `benchmarks/output_footprints/`
- `benchmarks/output_mesh_2d/`
- `benchmarks/output_mesh_3d/`
