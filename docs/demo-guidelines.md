# Demo guidelines

Demos teach the public Python API. Each has one clear purpose that does not
repeat another demo; prefer improving an existing example over adding one.

- Run from the repository root with `python demos/<name>.py` in the documented
  environment, without arguments or an execution switch.
- Use a short docstring, editable settings near the top, and linear top-level
  statements. No `main()`, entry-point guard, helper framework or mode loops.
- Aim for at most 40 lines. Show public API calls directly; reusable setup,
  processing and plotting algorithms belong in the library.
- Name dataset examples `dataset_<name>.py`: request a model, call `.info()`,
  then plot or view the relevant result. Direct builder or solver examples
  should teach a distinct workflow. Use domain names for datasets.
- Prefer small, offline inputs where possible. State environments, live-data
  requirements, units and limitations in the [demo catalog](../demos/README.md).
  Missing required dependencies or data must fail clearly, without fake results.
- Show a useful visualization without a flag. Save only useful artifacts under
  `output/`; directory creation belongs in the library. Supporting APIs must
  allow safe repeated runs without manual path edits or deleting valuable cases.
- Keep CLIs, deployment, detached launchers, benchmarks and reproduction tooling
  outside `demos/`. Numerical acceptance and export round-trip checks belong in
  tests. Keep headless setup, mocking and temporary output redirection there too.

Verify a new or revised demo's purpose, argument-free execution, visualization
and repeated-run behavior in its required environment. Use real numerical or
geometric execution for the claims it demonstrates.

## Core specifics

Use `import dtcc_core as dtcc`; workflow examples may also import Shapely's
`box`. Spatial mesh and city examples show `.view()` as well as `.plot()`.
Comparisons use separate explicit calls, each starting from the original input.
The inventory and static checks live in `tests/demos/`; execution instructions
are in [CONTRIBUTING.md](../CONTRIBUTING.md#demos).
