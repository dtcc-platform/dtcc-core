# DTCC Core

> **Note:** New core to match plans of Oct 2024

This project is part of the
[Digital Twin Platform (DTCC Platform)](https://github.com/dtcc-platform/)
developed at the
[Digital Twin Cities Centre](https://dtcc.chalmers.se/)
supported by Sweden’s Innovation Agency Vinnova under Grant No. 2019-421 00041.

# 1. DTCC Model

DTCC Model defines the common data model and data formats for DTCC Platform.
Additionally, DTCC Model provides utilities for working with the data model
and data formats.

# 2. DTCC IO

DTCC IO provides input/output (IO) for DTCC Platform.

# 3. DTCC Builder

DTCC Builder provides functionality for building city models for DTCC Platform.

## Citing

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04928/status.svg)](https://doi.org/10.21105/joss.04928)

# 4. DTCC Common

DTCC Common provides common utilities for DTCC Platform.


## Documentation

This project is documented as part of the
[DTCC Platform Documentation](https://platform.dtcc.chalmers.se/).


## Authors (in order of appearance)

* [Anders Logg](http://anders.logg.org) [1](#1-dtcc-model), [2](#2-dtcc-io), [3](#3-dtcc-builder), [4](#4-dtcc-common)
* [Vasilis Naserentin](https://www.chalmers.se/en/Staff/Pages/vasnas.aspx) [1](#1-dtcc-model), [2](#2-dtcc-io), [3](#3-dtcc-builder)
* [Dag Wästerberg](https://chalmersindustriteknik.se/sv/medarbetare/dag-wastberg/) [1](#1-dtcc-model), [2](#2-dtcc-io), [3](#3-dtcc-builder)
* [Orfeas Eleutheriou](http://orfeasel.com/) [3](#3-dtcc-builder)
* [Anton Olsson](mailto:anton.j.olsson@bredband.net) [3](#3-dtcc-builder)
* [Anton Annlöv](mailto:annlova@student.chalmers.se) [3](#3-dtcc-builder)
* [George Spaias](mailto:gspaiasa@ece.auth.gr) [3](#3-dtcc-builder)

## License

This project is licensed under the
[MIT license](https://opensource.org/licenses/MIT).

Copyright is held by the individual authors as listed at the top of
each source file.


## Initial sketch of structure
```
dtcc-core/
│
├── dtcc_core/
│   ├── __init__.py
│   ├── dtcc_io/
│   │   ├── __init__.py
│   ├── dtcc_common/
│   │   ├── __init__.py
│   └── dtcc_builder/
│       ├── __init__.py
│
├── tests/
│
├── docs/
│
├── examples/
│
├── pyproject.toml
├── README.md
└── LICENSE
```
