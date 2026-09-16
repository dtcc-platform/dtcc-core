# Native dependency notices

These files accompany the compiled `_dtcc_builder` extension. DTCC Core's own
license is the repository-root `LICENSE`; it does not replace the licenses of
the dependencies below.

- **Eigen 3.4.0:** `Eigen-MPL-2.0.txt` is the Mozilla Public License 2.0 text
  from https://www.mozilla.org/media/MPL/2.0/index.txt. `Eigen-notices.txt`
  retains copyright/license comment blocks from the vendored files, grouped by
  source path. The corresponding source is retained in
  `dtcc_core/cpp/external/Eigen` in the Git checkout and source archive.
- **nanoflann and its vector adaptor:** the BSD notices are copied from
  `dtcc_core/cpp/external/nanoflann.hpp` and
  `dtcc_core/cpp/include/KDTreeVectorOfVectorsAdaptor.h` respectively.
- **Mapbox earcut:** `earcut-ISC.txt` comes from
  https://github.com/mapbox/earcut.hpp/blob/master/LICENSE. The vendored
  implementation is `dtcc_core/cpp/external/earcut.hpp`.
- **pybind11:** `pybind11-BSD.txt` comes from
  https://github.com/pybind/pybind11/blob/v3.0.1/LICENSE and covers the 3.0-series
  build dependency. It is used to compile the extension, not imported at runtime.
- **Triangle (optional, disabled by default):** `Triangle-notice.txt` copies
  the complete opening notice from the bundled implementation header. Its terms
  are distinct from DTCC Core's MIT license. The source archive and Git checkout
  retain `dtcc_core/cpp/external/triangle/triangle.h`. An alternate implementation
  supplied by a builder must retain its own applicable notices and source terms.

TetGen is built and installed by the separate `dtcc-tetgen-wrapper` dependency;
its code is not part of DTCC Core's extension or source archive.

Update the corresponding notices when updating native dependencies. Wheels
include this directory in their distribution license metadata; source archives
also retain the native sources and their original notices.
