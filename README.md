# DTCC Core

DTCC Core provides the core functionality for DTCC Platform, including data
modeling, data wrangling, data generation, and data input/output.

This project is part of the
[Digital Twin Platform (DTCC Platform)](https://github.com/dtcc-platform/)
developed at the
[Digital Twin Cities Centre](https://dtcc.chalmers.se/)
supported by Sweden’s Innovation Agency Vinnova under Grant No. 2019-421 00041.

## Installation

DTCC Core is installed from source with [uv](https://docs.astral.sh/uv/). It is
not installed from PyPI. The package includes a C++ extension, so a C++17
compiler is required. CMake and Ninja are downloaded during the build if they are
not already installed.

```bash
git clone https://github.com/dtcc-platform/dtcc-core.git
cd dtcc-core
uv sync
```

`uv sync` creates `.venv`, installs the exact dependency versions recorded in
`uv.lock` together with the development tools, and builds DTCC Core in editable
mode. Run commands in that environment with `uv run`, for example
`uv run python my_script.py`.

To depend on DTCC Core from another uv project, add it from Git or from a local
checkout:

```bash
uv add "dtcc-core @ git+https://github.com/dtcc-platform/dtcc-core.git@develop"
uv add --editable ../dtcc-core
```

## Development

| Task | Command |
|---|---|
| Set up or update the environment | `uv sync` |
| Run the tests | `cd tests && uv run pytest` |
| Build a source distribution and wheel into `dist/` | `uv build` |
| Add a dependency | `uv add <package>` |
| Upgrade one locked dependency | `uv lock --upgrade-package <package>` |

Commit `uv.lock` together with any change to dependencies in `pyproject.toml`.

Editing Python files needs no reinstall. After editing C++ sources or CMake files,
run `uv sync` (or any `uv run` command) and uv rebuilds the extension. Build trees
are kept in `build/` so rebuilds are incremental.

### Using a local dtcc-mesher checkout

By default `dtcc-mesher` is installed from the Git commit pinned in
`pyproject.toml`. To test DTCC Core against a local checkout instead, install it
into the environment and pass `--no-sync` so uv does not restore the pinned
version:

```bash
uv pip install -e ../dtcc-mesher
cd tests && uv run --no-sync pytest
```

Run `uv sync` to return to the pinned version.

## Load a city

```python
import dtcc_core as dtcc

city = dtcc.load_city("path/to/city.dtcc")
```

All existing `load_*` functions are available directly: `load_model`, `load_city`,
`load_cityjson`, `load_3dbag`, `load_mesh`, `load_volume_mesh`, `load_mesh_as_city`,
`load_pointcloud`, `load_pointcloud_directory`, `load_raster`, `load_footprints`,
`load_landuse` and `load_roadnetwork`. These lazy aliases use the same loaders as
`dtcc.io`, including default native schema validation. Model classes live under
`dtcc.model`; other I/O operations remain under `dtcc.io`.

```python
print(city)                                 # compact representation
city.info()                                 # detailed tables and dataset context
city.plot()                                  # quick 3D Matplotlib inspection
ax = city.plot(field="velocity", show=False)  # vector magnitudes at their samples
```

See [native model previews](docs/model-preview.md) for representation selection,
supported geometry and preview limitations. Full visualization belongs in DTCC Twin.
See [model display and inspection](docs/model-display.md) for `repr`, `str`,
`.info()` and dataset parameter help.

## Documentation

The DTCC data contract lives in [`dtcc_core/schemas/`](dtcc_core/schemas/):
[`dtcc.yaml`](dtcc_core/schemas/dtcc.yaml) defines the LinkML semantic schema, and
[`dtcc.proto`](dtcc_core/schemas/dtcc.proto) defines the Protobuf wire format.
Their versions are independent; both specifications have stable file paths.

This project is documented as part of the
[DTCC Platform Documentation](https://platform.dtcc.chalmers.se/).

The guiding architecture and design principles are described in
[DESIGN.md](DESIGN.md).

## Authors (in order of appearance)

* [Anders Logg](http://anders.logg.org)
* [Vasilis Naserentin](https://www.chalmers.se/en/Staff/Pages/vasnas.aspx)
* [Dag Wästerberg](https://chalmersindustriteknik.se/sv/medarbetare/dag-wastberg/)
* [Orfeas Eleutheriou](http://orfeasel.com/)
* [Anton Olsson](mailto:anton.j.olsson@bredband.net)
* [Anton Annlöv](mailto:annlova@student.chalmers.se)
* [George Spaias](mailto:gspaiasa@ece.auth.gr)

## License

This project is licensed under the
[MIT license](https://opensource.org/licenses/MIT).

Copyrights are held by the individual authors as listed at the top of
each source file.

## Community guidelines

Comments, contributions, and questions are welcome. Please engage with
us through Issues, Pull Requests, and Discussions on our GitHub page.

## Local Function-Call Check

CI checks that every public Python function exported via `__all__` has at least
one executed body line in test coverage. Importing a function, evaluating its
defaults, or applying its decorators does not count. You can run the same check locally:

- Install the project and test tools:
  - `uv sync`

- Run tests with coverage to produce `tests/coverage.json` (run from the `tests` directory to match CI):
  - `cd tests`
  - `uv run pytest --maxfail=1 --disable-warnings --cov=dtcc_core --cov-report=term-missing --cov-report=json:coverage.json`

- From the project root, run the function-call checker:
  - `cd ..`
  - `uv run python scripts/check_public_api_calls.py --package dtcc_core --coverage-file tests/coverage.json`

Exit status `0` means each discovered function has body-execution evidence. A
non-zero exit prints missed functions and their source locations. Line coverage
cannot distinguish imports from calls for a function written entirely on its
`def` line; give it a separate body line to make it verifiable. Docstring-only
functions and native functions without Python source also cannot be verified.
This check does not measure native C++ coverage, all model methods, or call counts.

## Demos and Examples

We ship several demos to illustrate how to interact with the API. They are separated into two main categories:

1. **Dataset Demos**: small scripts that fetch one dataset and print a summary or
   preview (e.g., `demos/weather.py`, `demos/building_footprints.py`).
2. **Workflow Demos**: end-to-end pipelines that assemble multiple builder
   components to produce meshes, point clouds, and terrain (e.g.,
   `demos/build_city.py`, `demos/build_city_volume_mesh.py`).

The full list, with the preview each demo produces, is in
[docs/datasets/demo-catalog.md](docs/datasets/demo-catalog.md).

**Running Demos:**
All demos can be run directly as Python scripts from the repository root:
```bash
python demos/build_city.py
```
Most workflow demos download footprints and point clouds from the DTCC data
service and run meshing, so they need network access and are not run by normal
CI. `tests/demos/` static-checks every demo (syntax, imports, viewer guards); to
actually execute the workflow demos as a test, opt in with:
```bash
DTCC_RUN_DEMOS=1 pytest tests/demos/test_demo_execution.py
```

**Viewing Output:**
Some workflow demos generate meshes and point clouds. To interactively view the results, append the `--view` flag to the command:
```bash
python demos/build_city_flat_mesh.py --view
```
*Note: The 3D viewer will block execution until the window is closed.*

## Installation Notes

* **Native build and package contents**:
  `uv sync` builds the private C++ extension in an isolated build environment;
  `pybind11` is required there only. Installed wheels contain the extension and
  dependency notices, while source archives retain the C++ and vendor inputs.
  Third-party notices are listed in [licenses/README.md](licenses/README.md).
  The [native code map](docs/native-code.md) lists the retained kernels and their
  Python callers.
  On macOS, `cmake.define.DTCC_USE_HOMEBREW_LLVM=ON` selects Homebrew Clang before
  CMake configures the compiler; use a fresh build directory when changing compilers.

* **Surface meshing backends**:
  Earcut is the built-in fast triangulation used for the lightweight (no-refinement) meshing path. Quality-controlled meshing uses the external `dtcc_mesher` package when installed (preferred by the `auto` mesher). Support for the Triangle backend is optional and disabled by default to keep the standard installation minimal.

* **Enabling Triangle**:
  A Triangle implementation header is bundled in `dtcc_core/cpp/external/triangle`, which is the default `DTCC_TRIANGLE_DIR`. To supply a compatible implementation header, pass `cmake.define.DTCC_TRIANGLE_DIR=/path/to/triangle/prefix`; discovery checks that directory and its `include`, `include/triangle`, and `triangle` subdirectories. This backend compiles the implementation in `triangle.h`; a declarations-only header and a separate library are not supported. Build settings are passed to the CMake build as config settings.

  For the development environment:

  ```
  uv sync --reinstall-package dtcc-core \
    --config-settings-package dtcc-core:cmake.define.DTCC_USE_TRIANGLE=ON
  ```

  Repeat these settings when rebuilding. CMake options can persist in the incremental build directory; explicitly pass `cmake.define.DTCC_USE_TRIANGLE=OFF` when switching back to the default backend configuration.

  For a wheel:

  ```
  uv build --wheel -C cmake.define.DTCC_USE_TRIANGLE=ON
  ```

  If these options are omitted, the build will proceed without Triangle and will use earcut (and `dtcc_mesher`, if installed) depending on configuration and availability.

* **Volume meshing with TetGen**:
  Install the optional volume-meshing dependencies from the lockfile:

  ```
  uv sync --extra volume
  uv run --extra volume python demos/build_city_volume_mesh.py
  ```

  The [`dtcc-tetgen-wrapper`](https://github.com/dtcc-platform/dtcc-tetgen-wrapper)
  dependency is pinned to a Git commit. Its build downloads a fixed upstream
  TetGen source archive, verifies its SHA-256 checksum, and compiles and installs
  the native Python module. A fresh build requires internet access and a C++
  compiler; no separate wrapper clone or vendoring command is needed. TetGen's
  source version and checksum belong to the wrapper, while `uv.lock` records
  the wrapper revision used by dtcc-core.

  Include `--extra volume` when syncing this environment; a plain `uv sync`
  omits the extra and removes the wrapper. The demo above also downloads city
  data and writes its output under `demos/output/`.

  TetGen is the only 3D meshing backend. A missing wrapper raises an installation
  error before geometry preparation or dataset downloads; TetGen errors propagate.
  Base installation still supports non-volume operations.

  The legacy columnar mesher and its options `smoother_max_iterations`,
  `smoothing_relative_tolerance`, `aspect_ratio_threshold`, and `debug_step` have
  been removed from `build_city_volume_mesh()` and `City.build_volume_mesh()`.
  Pass remaining options after `tetgen_switch_overrides` by keyword. Footprint
  `smoothing` and the intermediate 2D `mesher` option remain supported.

  TetGen and its wrapper use the AGPL license. The wrapper includes the relevant
  license and third-party notices; installing it as an optional dependency does
  not change those terms. See the wrapper documentation and
  [TetGen licensing information](https://www.wias-berlin.de/software/tetgen/FAQ-license.jsp?lang=0)
  before distributing applications containing it or deploying a service using it.
  
## Import Time

`import dtcc_core` is kept deliberately cheap, because it runs before any
useful work in every script, test, and CLI invocation that touches the package.

Measure it with:

```
python scripts/measure_import_time.py
```

or `make measure-import-time`. The script reports the wall-clock cost, the
cumulative cost per top-level package, and the most expensive individual
imports. Pass `--module dtcc_core.io` to measure a subpackage, and use
`python -X importtime -c "import dtcc_core"` for the raw trace.

Two rules keep the import from drifting upward:

* **Heavy third-party imports go inside the function that uses them**, not at
  module scope. This applies to `scipy`, `rasterio`, `rasterstats`, `fiona`,
  `laspy`, `geopandas`, and `skimage`. A single module-scope `from
  scipy.spatial.transform import Rotation` once accounted for over half the
  cost of importing the package, for one call in one function.
  `tests/test_import_time.py` fails if any of them reach `sys.modules` on a
  bare import.

* **A module using `@register_model_method` must be imported eagerly.** Those
  decorators attach methods to model classes as an import side effect, so
  deferring such a module silently removes public API — `Raster.slope_aspect()`
  and eleven siblings went missing this way. `tests/test_model_methods.py`
  pins the expected method set per class. Note also that such a module cannot
  use `from __future__ import annotations`: `register_model_method` resolves
  the first parameter's annotation with `issubclass()`, which needs a real
  class rather than a string. Quote individual annotations instead.

### Makefile shortcuts

If you have `make` available, you can use these shortcuts from the repository root:

- `make install` — set up the environment with `uv sync`.
- `make test` — run the test suite.
- `make coverage` — run tests with coverage and write `tests/coverage.json`.
- `make check-public-api` — check that all public API functions were executed.
- `make verify-public-api` — run coverage and then the public API check.
- `make measure-import-time` — report the cost of `import dtcc_core`.
