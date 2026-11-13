# Copyright(C) 2026 Dag Wästberg
# Licensed under the MIT License

import datetime as _dt
import re
from dataclasses import dataclass, field

import numpy as np

from ..model import Model
from .field import Field

_NS = np.dtype("datetime64[ns]")
_TD_NS = np.dtype("timedelta64[ns]")
_MIN_NS = np.iinfo(np.int64).min + 1
_MAX_NS = np.iinfo(np.int64).max
_TIME_RANGE_ERROR = (
    "Time is outside the datetime64[ns] range (1677-09-21 to 2262-04-11)."
)
_SUB_NS_UNITS = {"ps": 1000, "fs": 10**6, "as": 10**9}
_ASSOCIATIONS = ("vertex", "edge", "face", "cell", "sample", "geometry")
_METHODS = ("linear", "previous", "nearest", "pchip")
_AGGREGATES = ("mean", "min", "max", "sum", "first", "last", "count")
_DUPLICATE_POLICIES = ("raise", "keep_old", "keep_new")
_DURATION = re.compile(r"^\s*(\d+)?\s*(ns|us|ms|s|min|h|d|w)\s*$")
_DURATION_UNITS = {
    "ns": "ns",
    "us": "us",
    "ms": "ms",
    "s": "s",
    "min": "m",
    "h": "h",
    "d": "D",
    "w": "W",
}


def _is_pandas(value):
    return type(value).__module__.split(".")[0] == "pandas"


def _checked_times(value):
    """Convert a datetime array to nanoseconds without wrapping finite times."""
    array = np.asarray(value)
    if array.dtype == _NS:
        return array
    unit, count = np.datetime_data(array.dtype)
    if unit == "generic":
        return array.astype(_NS)
    if unit in ("Y", "M"):
        dtype = f"datetime64[{unit}]"
        lower = int(np.datetime64(_MIN_NS, "ns").astype(dtype).astype(np.int64)) + 1
        upper = int(np.datetime64(_MAX_NS, "ns").astype(dtype).astype(np.int64))
        lower, upper = -(-lower // count), upper // count
    else:
        divisor = _SUB_NS_UNITS.get(unit, 1)
        multiplier = (
            count
            if divisor > 1
            else count * int(np.timedelta64(1, unit).astype(_TD_NS).astype(np.int64))
        )
        lower = -(-_MIN_NS * divisor // multiplier)
        upper = _MAX_NS * divisor // multiplier
    ticks = array.astype(np.int64)
    outside = (ticks < max(lower, _MIN_NS - 1)) | (ticks > min(upper, _MAX_NS))
    if (outside & ~np.isnat(array)).any():
        raise ValueError(_TIME_RANGE_ERROR)
    if unit in _SUB_NS_UNITS and count != 1:
        # Python integers avoid an overflowing intermediate product for scaled
        # subnanosecond dtypes whose final nanosecond count is representable.
        nanos = ticks.astype(object) * count // _SUB_NS_UNITS[unit]
        return np.where(np.isnat(array), _MIN_NS - 1, nanos).astype(np.int64).view(_NS)
    return array.astype(_NS)


def _datetime_time(value, extra_ns=0):
    """Normalize a Python date or datetime with an extra nanosecond count."""
    if isinstance(value, _dt.datetime) and value.tzinfo is not None:
        try:
            value = value.astimezone(_dt.UTC).replace(tzinfo=None)
        except OverflowError as exc:
            raise ValueError(_TIME_RANGE_ERROR) from exc
    nanos = int(np.datetime64(value, "us").astype(np.int64)) * 1000 + extra_ns
    if not _MIN_NS <= nanos <= _MAX_NS:
        raise ValueError(_TIME_RANGE_ERROR)
    return np.datetime64(nanos, "ns")


def _iso_time(value):
    """Parse ISO time, retaining nanoseconds in the clock and UTC offset."""
    try:
        parsed = _dt.datetime.fromisoformat(value)
    except ValueError:
        stamp = np.datetime64(value)
        # Parsing directly into an inferred nanosecond dtype can itself wrap.
        if not np.isnat(stamp) and stamp.astype("datetime64[D]") != np.datetime64(
            value, "D"
        ):
            raise ValueError(_TIME_RANGE_ERROR) from None
        return _checked_times(stamp)[()]
    zone_start = (
        max(value.rfind("+", 10), value.rfind("-", 10))
        if parsed.tzinfo is not None
        else -1
    )
    extra_ns = 0
    parts, cursor = [], 0
    for match in re.finditer(r"[.,](\d+)(?=$|[Zz+-])", value):
        if match.start() <= 10:
            continue
        remainder = int(match.group(1)[:9].ljust(9, "0"))
        if zone_start != -1 and match.start() > zone_start:
            remainder *= -1 if value[zone_start] == "+" else 1
        extra_ns += remainder
        parts.append(value[cursor : match.start()])
        cursor = match.end()
    if parts:
        parsed = _dt.datetime.fromisoformat("".join(parts) + value[cursor:])
    return _datetime_time(parsed, extra_ns)


def _scalar_time(value):
    """Convert one time value to a naive UTC ``datetime64[ns]``."""
    if isinstance(value, np.datetime64):
        return _checked_times(value)[()]
    if _is_pandas(value):
        if getattr(value, "tzinfo", None) is not None:
            value = value.tz_convert("UTC").tz_localize(None)
        return _checked_times(value.to_datetime64())[()]
    if isinstance(value, str):
        return _iso_time(value)
    if isinstance(value, _dt.date):
        return _datetime_time(value)
    raise TypeError(f"Cannot interpret {type(value).__name__} as a time.")


def _to_times(value):
    """Return ``(times, is_scalar)`` with times a 1-D ``datetime64[ns]`` array.

    Accepts NumPy datetime64 of any unit, ``datetime``, ISO-8601 strings and
    pandas timestamps or indexes. Time-zone-aware input is converted to UTC;
    naive input is taken as UTC. Finite times outside the nanosecond storage
    range and arrays with more than one dimension raise ValueError.
    """
    if isinstance(value, np.datetime64):
        return np.array([_scalar_time(value)]), True
    if _is_pandas(value) and hasattr(value, "__len__"):
        if value.dtype.kind != "M":
            return _to_times(value.to_numpy(dtype=object))
        if hasattr(value, "dt"):
            if value.dt.tz is not None:
                value = value.dt.tz_convert("UTC").dt.tz_localize(None)
        elif value.tz is not None:
            value = value.tz_convert("UTC").tz_localize(None)
        return _checked_times(value.to_numpy()), False
    if isinstance(value, np.ndarray):
        if value.ndim > 1:
            raise ValueError(
                f"Times must be one time or a 1-D sequence; got shape {value.shape}."
            )
        if value.dtype.kind == "M":
            return np.atleast_1d(_checked_times(value)), value.ndim == 0
        if value.ndim == 0:
            return np.array([_scalar_time(value.item())]), True
        return np.array([_scalar_time(v) for v in value.ravel()], dtype=_NS), False
    if isinstance(value, (list, tuple)):
        return np.array([_scalar_time(v) for v in value], dtype=_NS), False
    return np.array([_scalar_time(value)]), True


def _numpy_duration(value, what):
    """Convert a fixed-length duration without overflowing nanosecond storage."""
    unit, count = np.datetime_data(value.dtype)
    if unit == "generic":
        raise ValueError(f"{what} needs a unit, such as np.timedelta64(1, 'h').")
    if unit in ("Y", "M"):
        raise ValueError(
            f"{what} must be a fixed-length duration; months and years vary in length."
        )
    if np.isnat(value):
        return value.astype(_TD_NS)
    nanos = int(value.astype(np.int64)) * count
    if unit in _SUB_NS_UNITS:
        nanos //= _SUB_NS_UNITS[unit]
    else:
        nanos *= int(np.timedelta64(1, unit).astype(_TD_NS).astype(np.int64))
    if not _MIN_NS <= nanos <= _MAX_NS:
        raise ValueError(
            f"{what} is outside the timedelta64[ns] range (about 292 years)."
        )
    return np.timedelta64(nanos, "ns")


def _to_duration(value, what):
    """Convert a positive fixed-length duration to ``timedelta64[ns]``.

    Strings take a whole count and a lowercase unit: ``ns``, ``us``, ``ms``,
    ``s``, ``min``, ``h``, ``d`` or ``w``. Months and years are rejected
    because their length varies, and so is a ``np.timedelta64`` without a
    unit. NaT and durations longer than ``timedelta64[ns]`` can hold (about
    292 years) raise ValueError.
    """
    try:
        if isinstance(value, np.timedelta64):
            result = value
        elif _is_pandas(value):
            import pandas as pd

            if value is pd.NaT:
                raise ValueError(f"{what} must not be NaT.")
            if not hasattr(value, "to_timedelta64"):
                raise TypeError(
                    f"Cannot interpret {type(value).__name__} as a duration."
                )
            result = value.to_timedelta64()
        elif isinstance(value, _dt.timedelta):
            result = np.timedelta64(value)
        elif isinstance(value, str):
            match = _DURATION.match(value)
            if match is None:
                raise ValueError(f"Cannot interpret {value!r} as a duration.")
            count = int(match.group(1) or 1)
            unit = _DURATION_UNITS[match.group(2)]
            result = np.timedelta64(count, unit)
        else:
            raise TypeError(f"Cannot interpret {type(value).__name__} as a duration.")
        result = _numpy_duration(result, what)
    except OverflowError as exc:
        raise ValueError(
            f"{what} is outside the timedelta64[ns] range (about 292 years)."
        ) from exc
    if np.isnat(result):
        raise ValueError(f"{what} must not be NaT.")
    if result <= np.timedelta64(0, "ns"):
        raise ValueError(f"{what} must be a positive duration.")
    return result


def _span(start, stop):
    """Exact ``stop - start`` in nanoseconds as uint64, for ``stop >= start``.

    ``start`` and ``stop`` are int64 nanosecond times. Their int64 difference
    overflows for spans longer than about 292 years; the difference of their
    uint64 views wraps back to the exact value, which is below 2**64.
    """
    with np.errstate(over="ignore"):
        return np.asarray(stop).astype(np.uint64) - np.asarray(start).astype(np.uint64)


def _interval(nanos):
    """A ``timedelta64`` for a nanosecond count of any size.

    Counts that ``timedelta64[ns]`` cannot hold are given in microseconds,
    truncated.
    """
    if nanos <= _MAX_NS:
        return np.timedelta64(nanos, "ns")
    return np.timedelta64(nanos // 1000, "us")


def _strictly_increasing(times):
    return len(times) < 2 or bool((times[1:] > times[:-1]).all())


def _checked_cast(values, dtype, what):
    """Cast ``values`` to ``dtype``, refusing casts that lose integers or overflow.

    Kinds follow NumPy ``same_kind`` casting, except that any integer may go
    into any integer dtype. Integers outside the range of an integer
    ``dtype``, integers a float ``dtype`` would round, and finite values that
    would overflow to infinity raise ValueError. Floats cast to a narrower
    float round as NumPy rounds them.
    """
    if values.dtype == dtype:
        return values
    integers = dtype.kind in "iu" and values.dtype.kind in "biu"
    if not integers and not np.can_cast(values.dtype, dtype, casting="same_kind"):
        raise TypeError(f"Cannot cast {what} {values.dtype} values to {dtype}.")
    if dtype.kind in "iu" and values.size:
        info = np.iinfo(dtype)
        if values.min() < info.min or values.max() > info.max:
            raise ValueError(
                f"{what.capitalize()} values do not fit in {dtype} "
                f"(range {info.min} to {info.max})."
            )
    with np.errstate(over="ignore"):
        cast = values.astype(dtype)
    if dtype.kind == "f" and (np.isinf(cast) & ~np.isinf(values)).any():
        raise ValueError(f"{what.capitalize()} values overflow {dtype}.")
    if dtype.kind == "f" and values.dtype.kind in "iu":
        # Integers below 2**(nmant + 1) are exact; compare the others as
        # Python integers, which never round.
        mask = np.abs(cast) >= 2.0 ** (np.finfo(dtype).nmant + 1)
        for original, converted in zip(
            values[mask].tolist(), cast[mask].tolist(), strict=True
        ):
            if int(converted) != original:
                raise ValueError(
                    f"{what.capitalize()} values hold the integer {original}, "
                    f"which {dtype} cannot represent exactly."
                )
    return cast


def _format_duration(duration):
    """Readable text for a ``timedelta64`` in nanoseconds or coarser units."""
    if np.datetime_data(duration.dtype)[0] == "ns":
        ns = int(duration.astype(np.int64))
        if ns % 1000:
            return f"{ns} ns"
    return str(duration.astype("timedelta64[us]").item())


def _interpolation_gap(method, fill, max_gap):
    """Check the interpolation options and return ``max_gap`` as a duration."""
    if method not in _METHODS:
        raise ValueError(
            f"method must be one of {', '.join(_METHODS)}; got {method!r}."
        )
    if fill not in (None, "edge"):
        raise ValueError(f"fill must be None or 'edge'; got {fill!r}.")
    return None if max_gap is None else _to_duration(max_gap, "max_gap")


def _pchip_slopes(t, y, knots):
    """Shape-preserving (Fritsch-Carlson) slopes of ``y`` at the given knots.

    ``t`` holds int64 nanosecond times, shape ``(T,)`` with ``T >= 2``, and
    ``y`` has time on axis 0. Only the samples around each knot are read, so
    the cost grows with ``len(knots)``, not ``T``. A NaN sample makes the
    slopes that depend on it NaN.
    """
    shape = (-1,) + (1,) * (y.ndim - 1)
    if len(t) == 2:
        slope = (y[1].astype(np.float64) - y[0]) / (_span(t[0], t[1]) / 1e9)
        return np.broadcast_to(slope, (len(knots),) + y.shape[1:])
    # Each knot uses the two secants around it; the end knots use the first
    # or last two secants.
    first = np.clip(knots - 1, 0, len(t) - 3)
    y0, y1, y2 = (np.asarray(y[first + i], dtype=np.float64) for i in range(3))
    h0 = (_span(t[first], t[first + 1]) / 1e9).reshape(shape)
    h1 = (_span(t[first + 1], t[first + 2]) / 1e9).reshape(shape)
    m0, m1 = (y1 - y0) / h0, (y2 - y1) / h1
    w1, w2 = 2 * h1 + h0, h1 + 2 * h0
    with np.errstate(divide="ignore", invalid="ignore"):
        harmonic = (w1 + w2) / (w1 / m0 + w2 / m1)
    flat = (np.sign(m0) != np.sign(m1)) | (m0 == 0) | (m1 == 0)
    slopes = np.where(flat, 0.0, harmonic)

    def edge(ha, hb, ma, mb):
        d = ((2 * ha + hb) * ma - ha * mb) / (ha + hb)
        d = np.where(np.sign(d) != np.sign(ma), 0.0, d)
        steep = (np.sign(ma) != np.sign(mb)) & (np.abs(d) > 3 * np.abs(ma))
        return np.where(steep, 3 * ma, d)

    at_start = (knots == 0).reshape(shape)
    at_end = (knots == len(t) - 1).reshape(shape)
    slopes = np.where(at_start, edge(h0, h1, m0, m1), slopes)
    slopes = np.where(at_end, edge(h1, h0, m1, m0), slopes)
    slopes[np.isnan(m0) | np.isnan(m1)] = np.nan
    return slopes


def _aggregate(values, starts, how):
    """Aggregate ``values`` over the time groups beginning at ``starts``.

    NaN samples are skipped; a group with no valid sample gives NaN, or 0 for
    ``count``.
    """
    values = np.asarray(values, dtype=np.float64)
    valid = ~np.isnan(values)
    count = np.add.reduceat(valid.astype(np.int64), starts, axis=0)
    if how == "count":
        return count
    empty = count == 0
    if how in ("mean", "sum"):
        total = np.add.reduceat(np.where(valid, values, 0.0), starts, axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            result = total / count if how == "mean" else total
    elif how in ("min", "max"):
        ufunc = np.fmin if how == "min" else np.fmax
        result = ufunc.reduceat(values, starts, axis=0)
    else:
        steps = np.arange(len(values)).reshape((-1,) + (1,) * (values.ndim - 1))
        padded = np.concatenate([values, np.full((1,) + values.shape[1:], np.nan)])
        if how == "first":
            index = np.minimum.reduceat(
                np.where(valid, steps, len(values)), starts, axis=0
            )
        else:
            index = np.maximum.reduceat(np.where(valid, steps, -1), starts, axis=0)
            index[index < 0] = len(values)
        result = np.take_along_axis(padded, index, axis=0)
    result[empty] = np.nan
    return result


@dataclass(repr=False, eq=False, init=False)
class TimeSeries(Model):
    """Represents values (scalar or vector) on a set of elements over time.

    A time series is a field with a time axis: at each of ``T`` strictly
    increasing UTC timestamps it holds one value per element, for ``N``
    elements with ``dim`` components each. Historic observations and
    predictions use the same type and are told apart by ``role``.

    Times are stored as ``datetime64[ns]`` in UTC. Every method that accepts
    times also accepts ``datetime`` objects, ISO-8601 strings and pandas
    timestamps; time-zone-aware input is converted to UTC and naive input is
    taken as UTC. ISO strings retain nanosecond precision. Times outside the
    finite ``datetime64[ns]`` range raise ``ValueError``. ``NaN`` in ``values``
    marks a missing value.

    ``times`` and ``values`` are views into storage that grows geometrically,
    so appending a step is amortized O(1). The ``times`` view is read-only so
    the ordering cannot be broken in place; the ``values`` view is writable.
    A view taken before an append may stop reflecting the series once the
    append reallocates storage. The series copies arrays it is given, so it
    never shares memory with the caller. Assigning a new array to either
    attribute revalidates and copies it. Equality is identity.

    Attributes
    ----------
    name: str
        Name of the time series.
    unit: str
        Unit of measurement of the values.
    description: str
        Description of the time series.
    times: np.ndarray
        Strictly increasing ``datetime64[ns]`` UTC timestamps, shape ``(T,)``.
    values: np.ndarray
        Real numeric values, shape ``(T, N)`` for ``dim == 1`` and
        ``(T, N, dim)`` otherwise. A 1-D array of shape ``(T,)`` is accepted
        for ``dim == 1`` and stored as ``(T, 1)``.
    dim: int
        Number of components per element.
    association: str or None
        Element type the values live on: vertex, edge, face, cell, sample, or
        geometry, with the same meaning as ``Field.association``.
    role: str
        Free label such as ``"observed"`` or ``"predicted"``; empty means
        unspecified.
    """

    # times and values are dataclass fields whose defaults are the properties
    # defined below, so dataclasses.fields() and replace() include the data.
    name: str = ""
    unit: str = ""
    description: str = ""
    times: np.ndarray
    values: np.ndarray
    dim: int = 1
    association: str | None = field(default=None, kw_only=True)
    role: str = field(default="", kw_only=True)

    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def __init__(
        self,
        name="",
        unit="",
        description="",
        times=None,
        values=None,
        dim=1,
        *,
        association=None,
        role="",
    ):
        self.name = name
        self.unit = unit
        self.description = description
        self.dim = dim
        self.association = association
        self.role = role
        if times is None:
            if values is not None:
                raise ValueError("TimeSeries values were given without times.")
            times = np.empty(0, dtype=_NS)
        else:
            times, _ = _to_times(times)
        if values is None:
            if len(times):
                raise ValueError("TimeSeries times were given without values.")
            values = np.empty((0, 0) if dim == 1 else (0, 0, dim))
        self._set_data(times, values)

    # Storage

    def __setattr__(self, name, value):
        """Validate ``dim`` and ``association`` on every assignment."""
        if name == "dim":
            if (
                isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer))
                or value < 1
            ):
                raise ValueError("TimeSeries.dim must be a positive integer.")
            storage = self.__dict__.get("_values")
            if storage is not None and value != self.dim:
                if self._length or storage.shape[1]:
                    raise ValueError(
                        "TimeSeries.dim cannot change once the series has "
                        "elements; build a new TimeSeries."
                    )
                shape = (0, 0) if value == 1 else (0, 0, value)
                self.__dict__["_values"] = np.empty(shape, storage.dtype)
        elif name == "association":
            if value is not None and value not in _ASSOCIATIONS:
                raise ValueError(
                    "TimeSeries.association must be one of "
                    f"{', '.join(_ASSOCIATIONS)} or None; got {value!r}."
                )
        super().__setattr__(name, value)

    def _canonical_values(self, values, steps):
        values = np.asarray(values)
        if values.dtype.kind not in "biuf":
            raise ValueError("TimeSeries.values must be a real numeric array.")
        if self.dim == 1 and values.ndim == 1:
            values = values.reshape(-1, 1)
        expected_ndim = 2 if self.dim == 1 else 3
        if values.ndim != expected_ndim or (
            self.dim > 1 and values.shape[2] != self.dim
        ):
            raise ValueError(
                "TimeSeries.values must have shape (T, N) for dim=1 or "
                f"(T, N, dim) otherwise; got shape {values.shape} with "
                f"dim={self.dim}."
            )
        if values.shape[0] != steps:
            raise ValueError(
                f"TimeSeries.values has {values.shape[0]} steps but times has {steps}."
            )
        return values

    def _check_data(self, times, values):
        """Check the time and value invariants; return the canonical values."""
        if np.isnat(times).any():
            raise ValueError("TimeSeries.times must not contain NaT.")
        if not _strictly_increasing(times):
            raise ValueError("TimeSeries.times must be strictly increasing.")
        return self._canonical_values(values, len(times))

    def _set_data(self, times, values, copy=True):
        # copy=False adopts the arrays as storage; the caller must own them.
        values = self._check_data(times, values)
        if copy:
            times, values = np.array(times, dtype=_NS), np.array(values)
        self._times = times
        self._values = values
        self._length = len(times)

    def __copy__(self):
        """Return a copy with its own storage, so appends never interfere."""
        result = object.__new__(type(self))
        result.__dict__.update(self.__dict__)
        result._times = self._times[: self._length].copy()
        result._values = self._values[: self._length].copy()
        return result

    def __getstate__(self):
        """Return the state for deepcopy and pickle without spare capacity."""
        state = dict(self.__dict__)
        state["_times"] = self._times[: self._length]
        state["_values"] = self._values[: self._length]
        return state

    @property
    def times(self) -> np.ndarray:
        """Timestamps as a read-only ``datetime64[ns]`` view, shape ``(T,)``."""
        view = self._times[: self._length]
        view.flags.writeable = False
        return view

    @times.setter
    def times(self, times):
        times, _ = _to_times(times)
        if len(times) != self._length:
            raise ValueError(
                f"New times has {len(times)} steps but the series has "
                f"{self._length}; build a new TimeSeries to change the length."
            )
        self._set_data(times, self.values)

    @property
    def values(self) -> np.ndarray:
        """Values as a writable view, shape ``(T, N)`` or ``(T, N, dim)``."""
        return self._values[: self._length]

    @values.setter
    def values(self, values):
        self._set_data(self.times, values)

    def _reserve(self, extra):
        needed = self._length + extra
        capacity = len(self._times)
        if needed <= capacity and self._values.shape[0] >= needed:
            return
        capacity = max(needed, 2 * capacity, 16)
        times = np.empty(capacity, dtype=_NS)
        times[: self._length] = self._times[: self._length]
        values = np.empty((capacity,) + self._values.shape[1:], self._values.dtype)
        values[: self._length] = self._values[: self._length]
        self._times, self._values = times, values

    # Inspection

    def __len__(self):
        return self._length

    @property
    def start(self):
        """First timestamp as ``np.datetime64``, or None for an empty series."""
        return self._times[0] if self._length else None

    @property
    def end(self):
        """Last timestamp as ``np.datetime64``, or None for an empty series."""
        return self._times[self._length - 1] if self._length else None

    @property
    def num_elements(self) -> int:
        """Number of elements ``N`` each step holds a value for."""
        return self._values.shape[1]

    def is_regular(self):
        """Return ``(True, interval)`` when steps are evenly spaced.

        ``interval`` is a ``timedelta64[ns]``, or a ``timedelta64[us]``
        truncated to whole microseconds when it is longer than
        ``timedelta64[ns]`` can hold (about 292 years). Returns
        ``(False, None)`` for uneven spacing and for fewer than two steps.
        """
        if self._length < 2:
            return False, None
        t = self.times.view(np.int64)
        steps = _span(t[:-1], t[1:])
        if np.all(steps == steps[0]):
            return True, _interval(int(steps[0]))
        return False, None

    def _summary_items(self):
        return [
            ("name", self.name),
            ("unit", self.unit),
            ("shape", self.values.shape),
            ("dtype", str(self.values.dtype)),
            ("association", self.association),
            ("role", self.role),
        ]

    def _info_sections(self):
        sections = super()._info_sections()
        regular, interval = self.is_regular()
        sections[0][2].extend(
            [
                ("Dimension", self.dim),
                ("Description", self.description),
                ("Steps", self._length),
                ("Start", str(self.start) if self._length else "-"),
                ("End", str(self.end) if self._length else "-"),
                ("Interval", _format_duration(interval) if regular else "irregular"),
            ]
        )
        return sections

    # Append and merge

    def _field_values(self, value):
        """Return a Field's values after checking it matches this series."""
        if value.dim != self.dim:
            raise ValueError(
                f"Field dim {value.dim} does not match series dim {self.dim}."
            )
        if value.unit != self.unit:
            raise ValueError(
                f"Field unit {value.unit!r} does not match series unit {self.unit!r}."
            )
        if value.association != self.association:
            raise ValueError(
                f"Field association {value.association!r} does not match "
                f"series association {self.association!r}."
            )
        return value.values

    def _step_values(self, values, count, single):
        """Normalize appended values to shape ``(count, N[, dim])``."""
        if isinstance(values, Field):
            values = self._field_values(values)
        values = np.asarray(values)
        if values.dtype.kind not in "biuf":
            raise ValueError("Appended values must be real numeric.")
        if single:
            if self.dim == 1:
                if values.ndim > 2 or (values.ndim == 2 and values.shape[1] != 1):
                    raise ValueError(
                        "One step of a dim=1 series takes values of shape (N,); "
                        f"got shape {values.shape}."
                    )
                values = values.reshape(1, -1)
            else:
                values = values[np.newaxis]
        values = self._canonical_values(values, count)
        if self._length or self.num_elements:
            if values.shape[1] != self.num_elements:
                raise ValueError(
                    f"Appended values have {values.shape[1]} elements; the series "
                    f"has {self.num_elements}."
                )
        return _checked_cast(values, self._values.dtype, "appended")

    def append(self, times, values):
        """Append one step or a batch of steps after the end of the series.

        Parameters
        ----------
        times : time or sequence of times
            One time, or ``k`` strictly increasing times, all later than
            ``end``.
        values : array_like or Field
            For one time: the per-step values, shape ``(N,)`` or ``(N, dim)``
            (a scalar when ``N == dim == 1``), or a ``Field`` with matching
            dim, unit and association. For ``k`` times: values with leading
            axis ``k``.

        An empty batch with empty values does nothing.

        Raises
        ------
        ValueError
            If a time is not after ``end``, the batch is not strictly
            increasing, shapes do not match, or a value does not fit the
            series dtype. The series is unchanged.
        TypeError
            If the values are of a kind the series dtype cannot hold under
            ``same_kind`` casting, such as floats for an integer series. The
            series is unchanged.
        """
        new_times, single = _to_times(times)
        if not len(new_times):
            if isinstance(values, Field) or np.size(values):
                raise ValueError("An empty batch of times takes empty values.")
            return
        if np.isnat(new_times[0] if single else new_times).any():
            raise ValueError("Appended times must not contain NaT.")
        if not _strictly_increasing(new_times):
            raise ValueError("Appended times must be strictly increasing.")
        if self._length and new_times[0] <= self._times[self._length - 1]:
            raise ValueError(
                f"Appended time {new_times[0]} is not after the series end "
                f"{self.end}; use merge() for out-of-order data."
            )
        new_values = self._step_values(values, len(new_times), single)
        if not self._length and new_values.shape[1] != self.num_elements:
            self._values = np.empty((0,) + new_values.shape[1:], self._values.dtype)
        self._reserve(len(new_times))
        stop = self._length + len(new_times)
        self._times[self._length : stop] = new_times
        self._values[self._length : stop] = new_values
        self._length = stop

    def merge(self, other, on_duplicate="raise"):
        """Merge another series into this one, in time order.

        ``other`` must have the same dim, association, unit and element
        count; ``name`` and ``role`` may differ and this series keeps its own.

        Parameters
        ----------
        other : TimeSeries
            Series to merge in. It is not modified.
        on_duplicate : {"raise", "keep_old", "keep_new"}
            What to do with timestamps present in both series. With
            ``"raise"`` a duplicate raises ``ValueError`` and nothing changes.
        """
        if not isinstance(other, TimeSeries):
            raise TypeError(f"Cannot merge {type(other).__name__} into a TimeSeries.")
        if on_duplicate not in _DUPLICATE_POLICIES:
            raise ValueError(
                f"on_duplicate must be one of {', '.join(_DUPLICATE_POLICIES)}; "
                f"got {on_duplicate!r}."
            )
        for key in ("unit", "dim", "association"):
            if getattr(other, key) != getattr(self, key):
                raise ValueError(
                    f"Cannot merge: {key} {getattr(other, key)!r} does not match "
                    f"{getattr(self, key)!r}."
                )
        if not len(other):
            return
        incoming = other.values
        if (self._length or self.num_elements) and (
            incoming.shape[1] != self.num_elements
        ):
            raise ValueError(
                f"Cannot merge: {incoming.shape[1]} elements does not match "
                f"{self.num_elements}."
            )
        incoming = _checked_cast(incoming, self._values.dtype, "merged")
        if not self._length:
            self._set_data(other.times, incoming)
            return
        times = np.concatenate([self.times, other.times])
        values = np.concatenate([self.values, incoming])
        order = np.argsort(times, kind="stable")
        times, values = times[order], values[order]
        duplicate = times[1:] == times[:-1]
        if duplicate.any():
            if on_duplicate == "raise":
                raise ValueError(
                    f"Cannot merge: {int(duplicate.sum())} duplicate timestamps; "
                    "pass on_duplicate='keep_old' or 'keep_new'."
                )
            keep = np.ones(len(times), dtype=bool)
            # Stable sort puts this series' sample first in each duplicate pair.
            if on_duplicate == "keep_old":
                keep[1:][duplicate] = False
            else:
                keep[:-1][duplicate] = False
            times, values = times[keep], values[keep]
        self._set_data(times, values)

    # Access

    def _index_of(self, key):
        if isinstance(key, (int, np.integer)) and not isinstance(key, (bool, np.bool_)):
            index = int(key)
            if index < 0:
                index += self._length
            if not 0 <= index < self._length:
                raise IndexError(
                    f"Step {key} is out of range for {self._length} steps."
                )
            return index
        time, _ = _to_times(key)
        index = int(np.searchsorted(self.times, time[0]))
        if index == self._length or self._times[index] != time[0]:
            raise KeyError(f"No step at {time[0]}; use at() to interpolate.")
        return index

    def field_at(self, key) -> Field:
        """Return one step as a ``Field`` whose values are a zero-copy view.

        Parameters
        ----------
        key : int or time
            Step index (negative counts from the end) or an exact timestamp.

        Raises
        ------
        IndexError
            If an index is out of range.
        KeyError
            If no step has exactly this time.
        """
        index = self._index_of(key)
        values = self._values[index]
        return Field(
            name=self.name,
            unit=self.unit,
            description=self.description,
            values=values,
            dim=self.dim,
            association=self.association,
        )

    def _like(self, times, values, copy=True):
        result = TimeSeries(
            name=self.name,
            unit=self.unit,
            description=self.description,
            dim=self.dim,
            association=self.association,
            role=self.role,
        )
        result._set_data(times, values, copy=copy)
        return result

    def between(self, start=None, end=None) -> "TimeSeries":
        """Return a copy of the steps in the half-open range ``[start, end)``.

        ``None`` leaves that side open; a NaT bound raises ValueError.
        """

        def bound(value):
            time = _to_times(value)[0][0]
            if np.isnat(time):
                raise ValueError("between() bounds must not be NaT; use None.")
            return time

        times = self.times
        lo = 0 if start is None else np.searchsorted(times, bound(start))
        hi = len(times) if end is None else np.searchsorted(times, bound(end))
        hi = max(lo, hi)
        return self._like(times[lo:hi].copy(), self.values[lo:hi].copy())

    # Interpolation

    def at(self, times, method="linear", fill=None, max_gap=None) -> np.ndarray:
        """Interpolate values at arbitrary times.

        Query times that match a step exactly return that step's values for
        every method. NaT query times give NaN.

        Parameters
        ----------
        times : time or sequence of times
            Query times, in any order.
        method : {"linear", "previous", "nearest", "pchip"}
            ``"previous"`` holds the last value (states, setpoints,
            counters); ``"nearest"`` takes the closer step, the earlier one
            on a tie; ``"pchip"`` is a shape-preserving cubic without
            overshoot (Fritsch-Carlson, as in SciPy's ``PchipInterpolator``).
            Angles such as wind direction are interpolated without
            wrap-around.
        fill : None or "edge"
            Outside ``[start, end]``: ``None`` gives NaN, ``"edge"`` holds the
            first or last step.
        max_gap : duration, optional
            Queries inside an interval between steps longer than this give
            NaN.

        Returns
        -------
        np.ndarray
            Float values of shape ``(Q, N)`` or ``(Q, N, dim)`` for ``Q``
            query times, or the per-step shape for a single time. A NaN
            sample makes every query whose interpolation uses it NaN.
        """
        gap = _interpolation_gap(method, fill, max_gap)
        query, single = _to_times(times)
        result = self._interpolate(query, method, fill, gap)
        return result[0] if single else result

    def align_to(
        self, target, method="linear", fill=None, max_gap=None
    ) -> "TimeSeries":
        """Return this series interpolated onto the time steps of ``target``.

        The result has the times of ``target`` and the metadata of this
        series, with float values as returned by ``at`` with the same
        options. Steps of ``target`` outside ``[start, end]`` are NaN unless
        ``fill="edge"``, so the result always has ``len(target)`` steps; use
        ``between`` to crop it to the overlap.

        Target times that match steps of this series are copied, not
        interpolated. When the target times are a contiguous run of this
        series' steps, as for two series on the same clock, the values are
        copied as one block without searching.

        Parameters
        ----------
        target : TimeSeries or times
            A series whose times to align to, or strictly increasing times.
        method, fill, max_gap
            As for ``at``.
        """
        gap = _interpolation_gap(method, fill, max_gap)
        if isinstance(target, TimeSeries):
            times = target.times.copy()
        else:
            times = np.array(_to_times(target)[0])
            if np.isnat(times).any():
                raise ValueError("align_to() target times must not contain NaT.")
            if not _strictly_increasing(times):
                raise ValueError("align_to() target times must be strictly increasing.")
        run = self._matching_run(times)
        if run is None:
            values = self._interpolate(times, method, fill, gap)
        else:
            values = self.values[run].astype(np.float64)
        return self._like(times, values, copy=False)

    def _matching_run(self, times):
        """Return the slice of steps whose times equal ``times``, or None."""
        t = self.times.view(np.int64)
        q = times.view(np.int64)
        if not len(q):
            return slice(0, 0)
        if len(q) > len(t):
            return None
        lo = int(np.searchsorted(t, q[0]))
        hi = lo + len(q)
        if hi > len(t) or not np.array_equal(t[lo:hi], q):
            return None
        return slice(lo, hi)

    def _interpolate(self, query, method, fill, gap):
        # Only the rows each query needs are read and converted to float, so
        # a call costs O(Q log T), not O(T).
        t = self.times.view(np.int64)
        q = query.view(np.int64)
        values = self.values
        n = len(t)
        out = np.full((len(q),) + values.shape[1:], np.nan)
        if n == 0:
            return out
        missing = np.isnat(query)
        right = np.searchsorted(t, q, side="right")
        exact = (right > 0) & (t[np.maximum(right - 1, 0)] == q)
        inside = (q > t[0]) & (q < t[-1]) & ~exact
        if inside.any():
            j = right[inside] - 1
            lo, hi = t[j], t[j + 1]
            span = _span(lo, hi)
            w = _span(lo, q[inside]).astype(np.float64) / span.astype(np.float64)
            w = w.reshape((-1,) + (1,) * (values.ndim - 1))
            left = np.asarray(values[j], dtype=np.float64)
            right_values = np.asarray(values[j + 1], dtype=np.float64)
            if method == "linear":
                inner = left * (1.0 - w) + right_values * w
            elif method == "previous":
                inner = left
            elif method == "nearest":
                inner = np.where(w <= 0.5, left, right_values)
            else:
                h = (span / 1e9).reshape(w.shape)
                w2, w3 = w * w, w * w * w
                inner = (
                    (2 * w3 - 3 * w2 + 1) * left
                    + (w3 - 2 * w2 + w) * h * _pchip_slopes(t, values, j)
                    + (3 * w2 - 2 * w3) * right_values
                    + (w3 - w2) * h * _pchip_slopes(t, values, j + 1)
                )
            if gap is not None:
                inner[span > np.uint64(gap.astype(np.int64))] = np.nan
            out[inside] = inner
        out[exact] = values[right[exact] - 1]
        if fill == "edge":
            out[(q < t[0]) & ~missing] = values[0]
            out[q > t[-1]] = values[-1]
        return out

    # Resampling and aggregation

    def resample(self, freq, how="mean") -> "TimeSeries":
        """Aggregate into regular bins ``[t, t + freq)``.

        Bins are aligned to the UTC epoch and labelled by their left edge.
        NaN samples are skipped. A bin with no valid sample gives NaN, or 0
        for ``count``. Output values are float, except int64 for ``count``.

        The output has one row for every bin from the first sample to the
        last, empty or not, so its size is ``(end - start) / freq`` rather
        than the number of samples. Sparse data over a long span with a small
        ``freq`` can allocate a very large array.

        Raises ``ValueError`` when the first bin would start before the
        ``datetime64[ns]`` range, which only happens for data within ``freq``
        of 1677-09-21.

        Parameters
        ----------
        freq : duration
            Bin width, such as ``"1h"``, ``"15min"``, ``np.timedelta64`` or
            ``datetime.timedelta``.
        how : {"mean", "min", "max", "sum", "first", "last", "count"}
            Aggregate per bin.
        """
        if how not in _AGGREGATES:
            raise ValueError(
                f"how must be one of {', '.join(_AGGREGATES)}; got {how!r}."
            )
        width = _to_duration(freq, "freq").astype(np.int64)
        if not self._length:
            dtype = np.int64 if how == "count" else np.float64
            return self._like(
                np.empty(0, dtype=_NS), self.values.astype(dtype), copy=False
            )
        bins = self.times.view(np.int64) // width
        if int(bins[0]) * int(width) < _MIN_NS:
            raise ValueError(
                "The first resample bin starts before 1677-09-21, outside the "
                "datetime64[ns] range."
            )
        starts = np.flatnonzero(np.r_[True, bins[1:] != bins[:-1]])
        occupied = bins[starts] - bins[0]
        aggregated = _aggregate(self.values, starts, how)
        count = int(bins[-1]) - int(bins[0]) + 1
        fill = 0 if how == "count" else np.nan
        values = np.full((count,) + aggregated.shape[1:], fill, aggregated.dtype)
        values[occupied] = aggregated
        times = ((bins[0] + np.arange(count)) * width).view(_NS)
        return self._like(times, values, copy=False)

    def _reduce(self, how):
        if not self._length:
            return np.full(self.values.shape[1:], np.nan)
        return _aggregate(self.values, np.array([0]), how)[0]

    def mean(self) -> np.ndarray:
        """Mean over time per element, skipping NaN."""
        return self._reduce("mean")

    def min(self) -> np.ndarray:
        """Minimum over time per element, skipping NaN."""
        return self._reduce("min")

    def max(self) -> np.ndarray:
        """Maximum over time per element, skipping NaN."""
        return self._reduce("max")

    def sum(self) -> np.ndarray:
        """Sum over time per element, skipping NaN; NaN if no valid sample."""
        return self._reduce("sum")

    # pandas interop

    def to_pandas(self):
        """Return a pandas object indexed by a UTC ``DatetimeIndex``.

        A series with one scalar element gives a ``Series`` named after the
        series. Otherwise the result is a ``DataFrame`` with one column per
        element for ``dim == 1``, or ``(element, component)`` MultiIndex
        columns for ``dim > 1``.
        """
        import pandas as pd

        index = pd.DatetimeIndex(self.times).tz_localize("UTC")
        values = self.values
        if self.dim == 1 and self.num_elements == 1:
            return pd.Series(
                values[:, 0], index=index, name=self.name or None, copy=True
            )
        if self.dim == 1:
            return pd.DataFrame(values, index=index, copy=True)
        columns = pd.MultiIndex.from_product(
            [range(self.num_elements), range(self.dim)],
            names=["element", "component"],
        )
        flat = values.reshape(len(values), self.num_elements * self.dim)
        return pd.DataFrame(flat, index=index, columns=columns, copy=True)

    @classmethod
    def from_pandas(
        cls,
        obj,
        name=None,
        unit="",
        description="",
        *,
        association=None,
        role="",
    ) -> "TimeSeries":
        """Build a TimeSeries from the pandas layout ``to_pandas`` produces.

        A ``Series`` gives one scalar element and supplies the name when
        ``name`` is None. A ``DataFrame`` with flat columns gives one element
        per column, in column order; ``(element, component)`` MultiIndex
        columns give ``dim`` components per element, ordered by label
        whatever the column order. MultiIndex columns must hold every
        ``(element, component)`` pair exactly once.
        """
        dim = 1
        if name is None:
            name = obj.name if obj.ndim == 1 and obj.name is not None else ""
        if obj.ndim == 2 and obj.columns.nlevels > 1:
            obj, elements, dim = cls._element_component_frame(obj)
            shape = (len(obj), elements) if dim == 1 else (len(obj), elements, dim)
            values = obj.to_numpy().reshape(shape)
        else:
            values = obj.to_numpy()
        return cls(
            name=str(name),
            unit=unit,
            description=description,
            times=obj.index,
            values=values,
            dim=dim,
            association=association,
            role=role,
        )

    @staticmethod
    def _element_component_frame(frame):
        """Return ``frame`` with columns in (element, component) label order."""
        import pandas as pd

        columns = frame.columns
        if len(columns):
            columns = columns.remove_unused_levels()
        if columns.nlevels != 2:
            raise ValueError(
                "DataFrame columns must be flat or an (element, component) "
                f"MultiIndex; got {columns.nlevels} levels."
            )
        elements, components = columns.levels
        expected = pd.MultiIndex.from_product([elements, components])
        if (
            not columns.is_unique
            or len(columns) != len(expected)
            or not expected.isin(columns).all()
        ):
            raise ValueError(
                "MultiIndex columns must hold every (element, component) pair "
                "exactly once."
            )
        frame = frame.set_axis(columns, axis=1)
        return frame.loc[:, expected], len(elements), len(components)
