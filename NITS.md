# Nits

Nit findings from merged pull requests that were accepted but left unfixed.
The pull request that fixes one deletes its line. Keep entries sorted by path.

Entry format; the symbol and PR number are optional:

    - `path/to/file.py` (`symbol`): what is wrong; why it was left open. (#PR)

- `dtcc_core/model/values/dataset.py` (`DatasetValue.__eq__`): `numpy_array == DatasetValue(...)` returns NumPy's elementwise array, not one bool; `__array_ufunc__ = None` would fix it but also disables NumPy arithmetic on DatasetValue.
