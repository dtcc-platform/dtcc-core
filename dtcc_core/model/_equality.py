"""Array-aware equality shared by DTCC model classes.

The policy is documented in docs/design/model-contract.md ("Python equality"):
exact values, identical array shapes and dtypes, and missing values (NaN, NaT)
equal in matching positions.
"""

from collections.abc import Mapping
from dataclasses import fields

import numpy as np


def values_equal(a, b) -> bool:
    """Return one bool for two field values under the model equality policy."""
    a, b = _unwrap(a), _unwrap(b)
    if a is b:
        return True
    if _is_wrapper(a) or _is_wrapper(b):
        # Only a self-wrapping DatasetValue is left wrapped; identity decides.
        return False
    if isinstance(a, np.void) or isinstance(b, np.void):
        return (
            isinstance(a, np.void)
            and isinstance(b, np.void)
            and _arrays_equal(np.asarray(a), np.asarray(b))
        )
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return _arrays_equal(a, b)
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        # A list never equals a tuple; a namedtuple equals a plain tuple.
        same_kind = (isinstance(a, list) and isinstance(b, list)) or (
            isinstance(a, tuple) and isinstance(b, tuple)
        )
        return (
            same_kind
            and len(a) == len(b)
            and all(values_equal(x, y) for x, y in zip(a, b, strict=True))
        )
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        return a.keys() == b.keys() and all(values_equal(a[k], b[k]) for k in a)
    if _is_nat(a) and _is_nat(b):
        return type(a) is type(b)
    if _is_inexact(a) and _is_inexact(b):
        return _parts_equal(a.real, b.real) and _parts_equal(a.imag, b.imag)
    return bool(a == b)


def _arrays_equal(a, b) -> bool:
    if type(a) is not type(b) or a.shape != b.shape or a.dtype != b.dtype:
        return False
    if a.dtype.names is not None:
        # A masked record's field can be a plain scalar, so compare as values.
        return all(values_equal(a[name], b[name]) for name in a.dtype.names)
    if isinstance(a, np.ma.MaskedArray):
        mask = np.ma.getmaskarray(a)
        if not np.array_equal(mask, np.ma.getmaskarray(b)):
            return False
        return _arrays_equal(a.data[~mask], b.data[~mask])
    if a.dtype.hasobject:
        return all(values_equal(x, y) for x, y in zip(a.flat, b.flat, strict=True))
    if a.dtype.kind == "c":
        return _arrays_equal(a.real, b.real) and _arrays_equal(a.imag, b.imag)
    if a.dtype.kind in "fmM":
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


def _unwrap(value):
    """Innermost value of containers that compare as their contents (DatasetValue).

    Stops at a wrapper already seen and returns it, so a self-wrapping value
    cannot loop.
    """
    seen = set()
    while (unwrap := getattr(type(value), "_eq_unwrap", None)) is not None:
        if id(value) in seen:
            break
        seen.add(id(value))
        value = unwrap(value)
    return value


def _is_wrapper(value) -> bool:
    return getattr(type(value), "_eq_unwrap", None) is not None


def _is_nat(value) -> bool:
    return isinstance(value, (np.datetime64, np.timedelta64)) and bool(np.isnat(value))


def _is_inexact(value) -> bool:
    return isinstance(value, (float, complex, np.floating, np.complexfloating))


def _parts_equal(a, b) -> bool:
    """Compare real numbers exactly, with NaN equal to NaN."""
    return bool(a == b or (np.isnan(a) and np.isnan(b)))


def model_equal(self, other):
    """Compare two instances of the same concrete dataclass field by field.

    Fields declared with ``compare=False`` are skipped. Each value is read
    through ``_eq_field_value`` when the class defines it.
    """
    if other is self:
        return True
    if type(other) is not type(self):
        return NotImplemented
    read = getattr(type(self), "_eq_field_value", None)
    for f in fields(self):
        if not f.compare:
            continue
        if read is None:
            mine, theirs = getattr(self, f.name), getattr(other, f.name)
        else:
            mine, theirs = read(self, f.name), read(other, f.name)
        if not values_equal(mine, theirs):
            return False
    return True
