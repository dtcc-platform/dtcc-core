# Nits

Nit findings from merged pull requests that were accepted but left unfixed.
The pull request that fixes one deletes its line. Keep entries sorted by path.

Entry format; the symbol and PR number are optional:

    - `path/to/file.py` (`symbol`): what is wrong; why it was left open. (#PR)

- `dtcc_core/model/values/dataset.py` (`DatasetValue.__eq__`): `numpy_array == DatasetValue(...)` returns NumPy's elementwise array, not one bool; `__array_ufunc__ = None` would fix it but also disables NumPy arithmetic on DatasetValue.
- `dtcc_core/plotting/face_map.py` (`_layers`): The sampling note can overlap the scale bar, mention absent context and omit that field faces were sampled. This affects annotation clarity. (#197)
- `dtcc_core/plotting/face_map.py` (`_limits`): An invalid range with one explicit limit does not identify the automatic bound in the error message. The invalid range is rejected. (#197)
