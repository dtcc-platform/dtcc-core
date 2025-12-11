# DTCC Core

DTCC Core provides the core functionality for DTCC Platform, including data
modeling, data wrangling, data generation, and data input/output.

This project is part of the
[Digital Twin Platform (DTCC Platform)](https://github.com/dtcc-platform/)
developed at the
[Digital Twin Cities Centre](https://dtcc.chalmers.se/)
supported by Sweden’s Innovation Agency Vinnova under Grant No. 2019-421 00041.

## Documentation

This project is documented as part of the
[DTCC Platform Documentation](https://platform.dtcc.chalmers.se/).




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

CI enforces that every public API function (exported via `__all__`) in `dtcc_core` is executed at least once by the test suite. You can run the same check locally:

- Create a virtual environment and install the project and test tools:
  - `python -m venv .venv && source .venv/bin/activate`  (Windows: `python -m venv .venv && .venv\\Scripts\\activate`)
  - `pip install -e .`
  - `pip install pytest pytest-cov`

- Run tests with coverage to produce `tests/coverage.json` (run from the `tests` directory to match CI):
  - `cd tests`
  - `pytest --maxfail=1 --disable-warnings --cov=dtcc_core --cov-report=term-missing --cov-report=json:coverage.json`

- From the project root, run the function-call checker:
  - `cd ..`
  - `python scripts/check_public_api_calls.py --package dtcc_core --coverage-file tests/coverage.json`

Exit status `0` means all public functions were exercised by tests. A non‑zero exit prints the list of missed functions with their source locations so you can add or adjust tests.

## Installation Notes

* **Surface meshing backends**:
  SPADE is used by default via [`dtcc-pyspade-native`](https://github.com/dtcc-platform/dtcc-pyspade-native). Earcut is provided as a fast alternative for cases where a lightweight triangulation method is preferred. Support for the Triangle backend is optional and disabled by default to keep the standard installation minimal.

* **Enabling Triangle**:
  If you wish to build with Triangle support, ensure that the Triangle library is available on your system (a header-only setup is sufficient) and install `dtcc-core` with:

  ```
  pip install . \
    --config-settings=cmake.define.DTCC_USE_TRIANGLE=ON \
    --config-settings=cmake.define.DTCC_TRIANGLE_DIR=/path/to/triangle/prefix
  ```

  If these options are omitted, the build will proceed without Triangle and will use SPADE or earcut depending on configuration and availability.

* **Volume meshing with TetGen**:
  TetGen can be used for tetrahedral meshing through the minimal wrapper provided in the [`dtcc-tetgen-wrapper`](https://github.com/dtcc-platform/dtcc-tetgen-wrapper) repository:

  ```
  git clone https://github.com/dtcc-platform/dtcc-tetgen-wrapper.git
  cd dtcc-tetgen-wrapper
  pip install .
  ```
  
### Makefile shortcuts

If you have `make` available, you can use these shortcuts from the repository root:

- `make install` — install the package and test tooling.
- `make test` — run the test suite.
- `make coverage` — run tests with coverage and write `tests/coverage.json`.
- `make check-public-api` — check that all public API functions were executed.
- `make verify-public-api` — run coverage and then the public API check.
