# Contributing to DTCC Core

Questions, bug reports and pull requests are welcome. Include a small reproducer
when reporting a bug and describe the user-visible behavior of a proposed change.

Name development branches `dev/<short-description>`, unless there is a specific
reason to use another convention.

## Development environment

Follow the [installation instructions](README.md#installation), then run commands
from the repository root unless noted otherwise. `uv sync` creates `.venv`,
installs the locked dependencies and development tools, and builds DTCC Core in
editable mode.

Python edits are available without reinstalling. After changing C++ sources or
CMake files, run `uv sync` to rebuild the extension. Incremental build directories
are kept under `build/`. See the [installation guide](docs/installation.md) for
compiler settings and optional backends.

| Task | Command |
| --- | --- |
| Set up or update the environment | `uv sync` |
| Run the tests | `cd tests && uv run pytest` |
| Build a source distribution and wheel | `uv build` |
| Add a dependency | `uv add <package>` |
| Upgrade a locked dependency | `uv lock --upgrade-package <package>` |

Commit `uv.lock` together with dependency changes in `pyproject.toml`.
Build artifacts are written to `dist/`.

### Using a local dtcc-mesher checkout

To test against a sibling checkout instead of the pinned dependency:

```bash
uv pip install -e ../dtcc-mesher
cd tests
uv run --no-sync pytest
```

Use `--no-sync` while testing the local installation. Run `uv sync` from the Core
repository root to restore the pinned version.

## Tests and coverage

Run tests from `tests/` to match CI. The public-function check uses coverage to
verify that every exported Python function has at least one executed body line.
To run the check locally:

```bash
cd tests
uv run pytest --maxfail=1 --disable-warnings --cov=dtcc_core \
  --cov-report=term-missing --cov-report=json:coverage.json
cd ..
uv run python scripts/check_public_api_calls.py --strict \
  --package dtcc_core --coverage-file tests/coverage.json
```

`--strict` also fails when a package submodule cannot be imported. The checker
reports missed functions and their source locations. Imports, decorators and
default arguments do not count as body execution. Functions need a separate
executable body line to be verifiable; docstring-only functions and native
functions without Python source cannot be checked this way. This check does not
measure native C++ coverage, all model methods or call counts.

### Running demos

The [demo catalog](docs/datasets/demo-catalog.md) lists dataset examples and
complete workflows. Run individual demos from the repository root:

```bash
uv run python demos/dataset_city.py
```

Workflow demos that support `--view` open a viewer until its window is closed.
Many demos download data and run meshing. Normal demo tests perform static
checks; opt in to workflow execution with:

```bash
DTCC_RUN_DEMOS=1 uv run --extra volume pytest tests/demos/test_demo_execution.py
```

These execution tests need network access and the
[volume-meshing dependencies](docs/installation.md#volume-meshing-with-tetgen).

## Import conventions

Keep heavyweight third-party imports inside the functions that use them. This
includes `scipy`, `rasterio`, `rasterstats`, `fiona`, `laspy`, `geopandas` and
`skimage`. `tests/test_import_time.py` checks that a bare Core import does not
load these packages.

Modules using `@register_model_method` must be imported eagerly so their methods
are attached to model classes. `tests/test_model_methods.py` checks the expected
method set. These modules cannot use `from __future__ import annotations`:
registration passes the first parameter's annotation to `issubclass()`, which
requires a class. Quote individual annotations instead.

Measure import cost with:

```bash
uv run python scripts/measure_import_time.py
```

The report includes wall-clock time and expensive imports. Add
`--module dtcc_core.io` to measure a subpackage. For the raw trace, use
`uv run python -X importtime -c "import dtcc_core"`.

## Makefile shortcuts

Run these from the repository root:

| Command | Action |
| --- | --- |
| `make install` | Set up the environment |
| `make test` | Run the test suite |
| `make coverage` | Write test coverage to `tests/coverage.json` |
| `make check-public-api-strict` | Check public functions using existing coverage |
| `make verify-public-api-strict` | Run coverage and the strict public-function check |
| `make measure-import-time PYTHON="uv run python"` | Measure import cost |

See [DESIGN.md](DESIGN.md) for architecture and
[the native code map](docs/native-code.md) for C++ kernels and their Python callers.

## Progress phase timing

Use existing `ProgressTracker.phase()` scopes for observable dataset stages.
Callers can collect wall-clock timings with the context-local
[`phase_timing` hook](docs/progress-timing.md), including in silent mode.
