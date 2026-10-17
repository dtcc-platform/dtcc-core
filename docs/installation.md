# Installation guide

For the standard source checkout, follow the [README](../README.md#installation).
DTCC Core requires Python 3.11 or later and a C++17 compiler. The isolated build
environment supplies Python build dependencies, including pybind11, CMake and
Ninja as needed.

Run the commands below from the Core repository root unless another project is
specified.

## Using Core in another project

From an existing uv project, install Core from Git:

```bash
uv add "dtcc-core @ git+https://github.com/dtcc-platform/dtcc-core.git@develop"
```

Alternatively, use an editable local checkout:

```bash
uv add --editable ../dtcc-core
```

## Volume meshing with TetGen

Install the optional volume-meshing dependencies and run a volume-mesh example:

```bash
uv sync --extra volume
uv run --extra volume python demos/dataset_city_volume_mesh.py
```

The demo downloads city data and writes results under `demos/output/`.
Include `--extra volume` when syncing or running commands in this environment;
a plain `uv sync` omits the extra and removes the wrapper.

Volume meshing uses
[`dtcc-tetgen-wrapper`](https://github.com/dtcc-platform/dtcc-tetgen-wrapper).
Its build downloads and verifies the upstream TetGen sources, so a fresh build
needs internet access and a C++ compiler. Core pins the wrapper dependency in
`pyproject.toml` and `uv.lock`.

TetGen and its wrapper use the AGPL license. See the wrapper's documentation and
[TetGen licensing information](https://www.wias-berlin.de/software/tetgen/FAQ-license.jsp?lang=0)
for their terms. The base Core installation supports operations that do not
require volume meshing.

## Surface meshing

The default installation includes `dtcc-mesher` for quality-controlled surface
meshing. Earcut supplies the lightweight triangulation path without refinement.
Triangle is an optional backend, disabled by default.

### Enabling Triangle

For the development environment:

```bash
uv sync --reinstall-package dtcc-core \
  --config-settings-package dtcc-core:cmake.define.DTCC_USE_TRIANGLE=ON
```

For a wheel:

```bash
uv build --wheel -C cmake.define.DTCC_USE_TRIANGLE=ON
```

Repeat these settings when rebuilding. CMake options can persist in the
incremental build directory; explicitly pass
`cmake.define.DTCC_USE_TRIANGLE=OFF` when switching back to the default
configuration.

A Triangle implementation header is bundled in
`dtcc_core/cpp/external/triangle`, the default `DTCC_TRIANGLE_DIR`. To use another
compatible implementation, also pass
`cmake.define.DTCC_TRIANGLE_DIR=/path/to/triangle/prefix`. Discovery checks that
directory and its `include`, `include/triangle` and `triangle` subdirectories.
The backend compiles the implementation in `triangle.h`; a declarations-only
header and a separate library are not supported.

Triangle has its own license terms. See the
[third-party notices](../licenses/README.md).

## Compiler configuration

On macOS, select Homebrew Clang with:

```bash
uv sync --reinstall-package dtcc-core \
  --config-settings-package dtcc-core:cmake.define.DTCC_USE_HOMEBREW_LLVM=ON
```

Install Homebrew LLVM first and use a fresh build directory when changing
compilers. The option selects the compiler before CMake configures the project.

## Package contents

Wheels include the compiled C++ extension and dependency notices. Source archives
also include the C++ and vendor inputs. See the
[third-party notices](../licenses/README.md) and
[native code map](native-code.md) for details.

For editable development, local dependency overrides and testing, see
[CONTRIBUTING.md](../CONTRIBUTING.md).


## Compact catalogue meshing

Catalogue surfaces with the resolved `dtcc_mesher` backend now use the bounded
compact path by default. Flat/LOD1 and ordinary Mesher behavior remain unchanged.
Read [the public contract](compact_mesher_interface.md) for zero-point controls,
backend restrictions, source/group budgets and triangle-growth ceilings.

The required native `improve` operation and bounded `segment_trials` result are
currently uncommitted in the sibling Mesher checkout. Its base commit is
`0bba17bbd0d5de4b3a65b0451c4e5db303e541ac`; that commit alone is **not** a
compatible release. The portable Core pin remains
`7ada8a89c98dbbedc34213387ee598bd92826f10`. No compatible revision is fabricated
and no absolute local dependency is placed in pyproject.toml or the lockfile.
A normal clean install at that pin cannot yet run compact catalogue meshing;
the early capability check explains the required build instead of silently
using an unrecorded alternative. The remaining release gate is to commit/release
the validated Mesher changes and update the portable pin to that real revision.

For this local integration, from dtcc-core with both checkouts present:

```sh
uv pip install --python .venv/bin/python --no-deps ../dtcc-mesher
.venv/bin/python -c 'from dtcc_core.builder.meshing.compact_surface import require_compact_mesher; require_compact_mesher()'
```

This builds native Mesher and Python bindings together. Source hashes, actual
loaded binary SHA-256, base revision, compiler/test commands and fixed evaluation
inputs are recorded in `.agent/work/compact-integration/validation.json` and its
source snapshot. A C consumer must rebuild against the matching unreleased
improvement-result header because the bounded trial array extends its struct.
The Core native extension is unchanged by this integration. The small native
constraint-recovery fix skips non-flippable crossed edges before adding recovery
points; its generic holed-polygon regression ensures a recoverable zero-refinement
input keeps all 11 points. Existing ordinary Mesher suites remain required.

Use `.venv/bin/python` for this unreleased local build. `uv sync` or an automatic
sync by `uv run` can restore the portable pinned dependency; rebuild the local
Mesher afterward, or use `uv run --no-sync` while validating this checkout.
