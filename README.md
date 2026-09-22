# DTCC Core

DTCC Core is a Python library for working with digital twins of cities. It
provides data models, access to geospatial datasets, tools for building terrain
and city meshes, and support for importing, inspecting and exporting city models.

DTCC Core is part of the [DTCC Platform](https://github.com/dtcc-platform/).

## Installation

Install from source using [uv](https://docs.astral.sh/uv/). You need Python 3.11
or later and a C++17 compiler.

```bash
git clone https://github.com/dtcc-platform/dtcc-core.git
cd dtcc-core
uv sync
```

Run Python scripts in the installed environment with `uv run python`.
See the [installation guide](docs/installation.md) for optional volume meshing,
compiler settings and using DTCC Core in another project.

## Quickstart

Load and inspect the bundled city model of The Hague. From the repository root,
start Python:

```bash
uv run python
```

Then run:

```python
import dtcc_core as dtcc

city = dtcc.load_city("demos/data/DenHaag_01.city.json.zip")
city.info()
city.plot()
```

The example uses local data. `city.plot()` opens an interactive 3D preview and
requires a graphical environment.

## Examples and documentation

- [Demo catalog](docs/datasets/demo-catalog.md): datasets and city-building workflows.
- [Model inspection](docs/model-display.md) and [plotting](docs/model-preview.md).
- [Reprojection](docs/reprojection.md): changing coordinates while preserving fields.
- [Data schemas](dtcc_core/schemas/): model definitions and exchange formats.
- [Architecture and design](DESIGN.md).

## Contributing

Questions, bug reports and contributions are welcome through the repository's
[issues](https://github.com/dtcc-platform/dtcc-core/issues) and
[pull requests](https://github.com/dtcc-platform/dtcc-core/pulls).
See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, tests and coding
conventions.

## Credits

Developed at the [Digital Twin Cities Centre](https://dtcc.chalmers.se/),
supported by Sweden’s Innovation Agency Vinnova under Grant No. 2019-421 00041.

Authors (in order of appearance):

* [Anders Logg](http://anders.logg.org)
* [Vasilis Naserentin](https://www.chalmers.se/en/Staff/Pages/vasnas.aspx)
* [Dag Wästerberg](https://chalmersindustriteknik.se/sv/medarbetare/dag-wastberg/)
* [Orfeas Eleutheriou](http://orfeasel.com/)
* [Anton Olsson](mailto:anton.j.olsson@bredband.net)
* [Anton Annlöv](mailto:annlova@student.chalmers.se)
* [George Spaias](mailto:gspaiasa@ece.auth.gr)

## License

DTCC Core is licensed under the [MIT license](LICENSE). Copyrights are held by
the individual authors as listed in the source files. See the
[third-party notices](licenses/README.md) and
[optional dependency licenses](docs/installation.md#volume-meshing-with-tetgen)
for bundled code and optional meshing dependencies.
