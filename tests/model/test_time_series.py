"""Behaviour of the TimeSeries value container."""

import copy
import dataclasses
import datetime as dt
import time

import numpy as np
import pytest

from dtcc_core.model import Field, TimeSeries

T0 = np.datetime64("2024-01-01T00:00", "ns")
H = np.timedelta64(3600 * 10**9, "ns")


def hours(*offsets):
    return T0 + np.array(offsets) * H


def scalar(values, offsets=None, **kwargs):
    values = np.asarray(values)
    if offsets is None:
        offsets = range(len(values))
    return TimeSeries(
        name="temperature",
        unit="C",
        times=hours(*offsets),
        values=values,
        association="sample",
        **kwargs,
    )


# Construction and validation


def test_empty_series():
    ts = TimeSeries()
    assert len(ts) == 0
    assert ts.times.dtype == np.dtype("datetime64[ns]")
    assert ts.values.shape == (0, 0)
    assert ts.start is None and ts.end is None
    assert TimeSeries(dim=3).values.shape == (0, 0, 3)


def test_scalar_shorthand_is_normalized():
    ts = scalar([1.0, 2.0, 3.0])
    assert ts.values.shape == (3, 1)
    assert ts.num_elements == 1
    assert ts.start == T0 and ts.end == T0 + 2 * H


def test_vector_series_shape():
    ts = TimeSeries(times=hours(0, 1), values=np.zeros((2, 4, 3)), dim=3)
    assert ts.values.shape == (2, 4, 3)
    assert ts.num_elements == 4


@pytest.mark.parametrize(
    "times,values,dim,match",
    [
        (hours(0, 1), np.zeros(3), 1, "times"),
        (hours(1, 0), np.zeros(2), 1, "strictly increasing"),
        (hours(0, 0), np.zeros(2), 1, "strictly increasing"),
        (hours(0, 1), np.zeros((2, 4)), 3, "shape"),
        (hours(0, 1), np.zeros((2, 4, 1)), 1, "shape"),
        (hours(0, 1), np.array(["a", "b"]), 1, "real numeric"),
        (hours(0, 1), np.zeros(2), 0, "positive integer"),
        (hours(0, 1), np.zeros(2), True, "positive integer"),
    ],
)
def test_constructor_rejects_invalid_data(times, values, dim, match):
    with pytest.raises(ValueError, match=match):
        TimeSeries(times=times, values=values, dim=dim)


def test_constructor_rejects_nat_and_bad_metadata():
    with pytest.raises(ValueError, match="NaT"):
        TimeSeries(
            times=np.array(["2024-01-01", "NaT"], dtype="datetime64[ns]"),
            values=np.zeros(2),
        )
    with pytest.raises(ValueError, match="association"):
        TimeSeries(association="pixel")
    with pytest.raises(ValueError, match="values"):
        TimeSeries(values=np.zeros(2))


def test_time_input_forms_are_equivalent():
    expected = np.array(["2024-01-01T00:00", "2024-01-01T01:00"], dtype="M8[ns]")
    utc_plus_1 = dt.timezone(dt.timedelta(hours=1))
    forms = [
        expected.astype("M8[s]"),
        [dt.datetime(2024, 1, 1, 0), dt.datetime(2024, 1, 1, 1)],
        [dt.datetime(2024, 1, 1, 1, tzinfo=utc_plus_1), "2024-01-01T01:00:00Z"],
        ["2024-01-01T00:00", "2024-01-01T02:00+01:00"],
    ]
    for form in forms:
        ts = TimeSeries(times=form, values=np.zeros(2))
        np.testing.assert_array_equal(ts.times, expected)


@pytest.mark.parametrize("suffix", ["", "Z", "+01:00", "-01:30"])
def test_iso_time_preserves_nanoseconds_across_public_methods(suffix):
    stamps = [f"2024-01-01T01:30:00.12345678{i}{suffix}" for i in range(3)]
    offset = {"": 0, "Z": 0, "+01:00": 60, "-01:30": -90}[suffix]
    expected = np.array(
        [s.removesuffix(suffix) if suffix else s for s in stamps], dtype="M8[ns]"
    ) - np.timedelta64(offset, "m")
    ts = TimeSeries(times=stamps[:2], values=[1.0, 2.0])
    ts.append(stamps[2], 3.0)
    np.testing.assert_array_equal(ts.times, expected)
    assert ts.field_at(stamps[1]).values[0] == 2.0
    np.testing.assert_array_equal(ts.at(stamps)[:, 0], [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(ts.between(stamps[1], stamps[2]).values, [[2.0]])
    ts.times = stamps
    np.testing.assert_array_equal(ts.times, expected)


@pytest.mark.parametrize(
    "stamp,expected",
    [
        ("20240101T010203.123456789", "2024-01-01T01:02:03.123456789"),
        ("2024-W01-1T01:02:03.123456789", "2024-01-01T01:02:03.123456789"),
        ("2024-01-01T01:02:03,123456789", "2024-01-01T01:02:03.123456789"),
        (
            "2024-01-01T00:00:00.123456789+01:00:00.000000001",
            "2023-12-31T23:00:00.123456788",
        ),
        (
            "2024-01-01T00:00:00.123456789-01:00:00.000000001",
            "2024-01-01T01:00:00.123456790",
        ),
        ("2024-01-01T00:00:00.123456789+00:00:00.123456789", "2024-01-01T00:00:00"),
        ("2024", "2024-01-01"),
        ("2024-01", "2024-01-01"),
    ],
)
def test_iso_time_forms_preserve_fractional_seconds(stamp, expected):
    assert TimeSeries(times=stamp, values=[1.0]).start == np.datetime64(expected, "ns")


@pytest.mark.parametrize(
    "stamp", ["2024-01-01T00:00:00.123.456", "2024-01-01T00:00:00+01:00:00.123.456"]
)
def test_malformed_iso_fractions_are_rejected(stamp):
    with pytest.raises(ValueError):
        TimeSeries(times=stamp, values=[1.0])


@pytest.mark.parametrize("form", ["NaT", np.datetime64("NaT")])
def test_scalar_nat_is_rejected(form):
    with pytest.raises(ValueError, match="NaT"):
        TimeSeries(times=form, values=[1.0])


@pytest.mark.parametrize("stamp", ["1500-01-01", "2500-01-01"])
@pytest.mark.parametrize("form", ["str", "datetime", "date", "numpy", "array"])
def test_out_of_range_time_is_rejected(stamp, form):
    forms = {
        "str": stamp,
        "datetime": dt.datetime.fromisoformat(stamp),
        "date": dt.date.fromisoformat(stamp),
        "numpy": np.datetime64(stamp, "s"),
        "array": np.array([stamp], dtype="M8[s]"),
    }
    with pytest.raises(ValueError, match="range"):
        TimeSeries(times=forms[form], values=[1.0])


@pytest.mark.parametrize("unit", ["Y", "M", "W", "D", "h", "m", "s", "ms", "us"])
def test_coarse_numpy_time_range_is_checked_before_conversion(unit):
    for stamp in ("1500-01-01", "2500-01-01"):
        with pytest.raises(ValueError, match="range"):
            TimeSeries(times=np.array([stamp], dtype=f"M8[{unit}]"), values=[1.0])
    expected = np.array(["1900-01-01", "2024-01-01"], dtype=f"M8[{unit}]")
    ts = TimeSeries(times=expected, values=[1.0, 2.0])
    np.testing.assert_array_equal(ts.times, expected.astype("M8[ns]"))


@pytest.mark.parametrize(
    "unit,ticks,expected",
    [
        ("ps", -1001, -2),
        ("fs", 1000001, 1),
        ("as", 1000000001, 1),
        ("1234ps", 10**18, 1234 * 10**15),
        ("1234ps", -(10**18), -1234 * 10**15),
        ("2ns", 10**9, 2 * 10**9),
        ("2Y", 1, 63072000000000000),
    ],
)
def test_numpy_time_units_convert_without_intermediate_overflow(unit, ticks, expected):
    array = np.array([ticks], dtype=np.int64).view(f"M8[{unit}]")
    for form in (array, array[0], array.reshape(())):
        assert TimeSeries(times=form, values=[1.0]).start == np.datetime64(
            expected, "ns"
        )
    with pytest.raises(ValueError, match="NaT"):
        TimeSeries(times=np.array(["NaT"], dtype=f"M8[{unit}]"), values=[1.0])


@pytest.mark.parametrize(
    "unit,first,last,before,after",
    [
        ("Y", "1678", "2262", "1677", "2263"),
        ("M", "1677-10", "2262-04", "1677-09", "2262-05"),
        (
            "us",
            "1677-09-21T00:12:43.145225",
            "2262-04-11T23:47:16.854775",
            "1677-09-21T00:12:43.145224",
            "2262-04-11T23:47:16.854776",
        ),
    ],
)
def test_coarse_numpy_time_range_endpoints(unit, first, last, before, after):
    valid = np.array([first, last], dtype=f"M8[{unit}]")
    np.testing.assert_array_equal(
        TimeSeries(times=valid, values=[1.0, 2.0]).times, valid.astype("M8[ns]")
    )
    for stamp in (before, after):
        with pytest.raises(ValueError, match="range"):
            TimeSeries(times=np.datetime64(stamp, unit), values=[1.0])


@pytest.mark.parametrize(
    "stamp",
    ["1677-09-21T00:12:43.145224193", "2262-04-11T23:47:16.854775807"],
)
def test_nanosecond_range_endpoints_are_preserved(stamp):
    for form in (stamp, stamp + "Z", np.datetime64(stamp, "ns")):
        assert TimeSeries(times=form, values=[1.0]).start == np.datetime64(stamp, "ns")


@pytest.mark.parametrize(
    "stamp",
    [
        "1677-09-21T00:12:43.145224192",
        "2262-04-11T23:47:16.854775808",
        "1677-09-21T00:12:43.145224193+00:01",
        "2262-04-11T23:47:16.854775807-00:01",
        "10000-01-01T00:00:00.123456789",
        "0000-01-01T00:00:00.123456789",
    ],
)
def test_times_outside_nanosecond_endpoints_are_rejected(stamp):
    with pytest.raises(ValueError, match="range"):
        TimeSeries(times=stamp, values=[1.0])


@pytest.mark.parametrize(
    "stamp,expected",
    [
        ("1677-09-21T00:11:43.145224193-00:01", "1677-09-21T00:12:43.145224193"),
        ("2262-04-11T23:48:16.854775807+00:01", "2262-04-11T23:47:16.854775807"),
    ],
)
def test_range_is_checked_after_timezone_normalization(stamp, expected):
    assert TimeSeries(times=stamp, values=[1.0]).start == np.datetime64(expected, "ns")


def test_out_of_range_time_is_rejected_by_all_time_apis_without_mutation():
    ts = scalar([1.0, 2.0])
    original_times, original_values = ts.times.copy(), ts.values.copy()
    with pytest.raises(ValueError, match="range"):
        ts.times = ["2500-01-01", "2500-01-02"]
    with pytest.raises(ValueError, match="range"):
        ts.append(["2024-01-02", "2500-01-01"], [3.0, 4.0])
    with pytest.raises(ValueError, match="range"):
        ts.field_at("2500-01-01")
    with pytest.raises(ValueError, match="range"):
        ts.at(["2024-01-01", "2500-01-01"])
    with pytest.raises(ValueError, match="range"):
        ts.between(end="2500-01-01")
    np.testing.assert_array_equal(ts.times, original_times)
    np.testing.assert_array_equal(ts.values, original_values)


@pytest.mark.parametrize("stamp", ["1500-01-01", "2500-01-01"])
def test_pandas_out_of_range_times_are_rejected(stamp):
    pd = pytest.importorskip("pandas")
    value = pd.Timestamp(stamp)
    index = pd.DatetimeIndex(np.array([stamp], dtype="M8[s]"))
    for form in (value, index, pd.Series(index), index.tz_localize("UTC")):
        with pytest.raises(ValueError, match="range"):
            TimeSeries(times=form, values=[1.0])
    with pytest.raises(ValueError, match="range"):
        TimeSeries.from_pandas(pd.Series([1.0], index=index))


def test_pandas_nanosecond_times_are_preserved():
    pd = pytest.importorskip("pandas")
    stamp = "2024-01-01T01:00:00.123456789+01:00"
    index = pd.DatetimeIndex([stamp])
    expected = np.datetime64("2024-01-01T00:00:00.123456789", "ns")
    for form in (pd.Timestamp(stamp), index, pd.Series(index)):
        assert TimeSeries(times=form, values=[1.0]).start == expected
    assert TimeSeries.from_pandas(pd.Series([1.0], index=index)).start == expected


def test_pandas_time_input():
    pd = pytest.importorskip("pandas")
    index = pd.date_range("2024-01-01 01:00", periods=2, freq="h", tz="Europe/Paris")
    ts = TimeSeries(times=index, values=np.zeros(2))
    assert ts.start == T0
    ts.append(pd.Timestamp("2024-01-01 02:00", tz="UTC"), 1.0)
    assert ts.end == T0 + 2 * H


def test_pandas_series_time_input():
    pd = pytest.importorskip("pandas")
    naive = pd.Series(pd.date_range("2024-01-01", periods=2, freq="h"))
    aware = naive.dt.tz_localize("Europe/Paris")
    np.testing.assert_array_equal(
        TimeSeries(times=naive, values=[1, 2]).times, hours(0, 1)
    )
    np.testing.assert_array_equal(
        TimeSeries(times=aware, values=[1, 2]).times, hours(-1, 0)
    )


def test_pandas_string_time_input():
    pd = pytest.importorskip("pandas")
    strings = ["2024-01-01T00:00", "2024-01-01T02:00+01:00"]
    for form in (pd.Series(strings), pd.Index(strings)):
        ts = TimeSeries(times=form, values=[1.0, 2.0])
        np.testing.assert_array_equal(ts.times, hours(0, 1))


def test_times_view_is_read_only_and_values_view_is_writable():
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError):
        ts.times[0] = T0 - H
    ts.values[0, 0] = 10.0
    assert ts.values[0, 0] == 10.0


def test_assigning_arrays_revalidates():
    ts = scalar([1.0, 2.0])
    ts.values = np.array([5.0, 6.0])
    assert ts.values.shape == (2, 1)
    with pytest.raises(ValueError, match="steps"):
        ts.values = np.zeros(3)
    ts.times = hours(10, 11)
    assert ts.start == T0 + 10 * H
    with pytest.raises(ValueError, match="strictly increasing"):
        ts.times = hours(11, 10)
    with pytest.raises(ValueError, match="times"):
        ts.times = hours(1, 2, 3)


def test_is_regular():
    assert scalar([1.0, 2.0, 3.0]).is_regular() == (True, H)
    assert scalar([1.0, 2.0, 3.0], offsets=[0, 1, 3]).is_regular() == (False, None)
    assert scalar([1.0]).is_regular() == (False, None)


def test_repr_and_info():
    ts = scalar([1.0, 2.0], role="observed")
    text = repr(ts)
    assert "TimeSeries" in text and "temperature" in text
    labels = [row[0] for row in ts._info_sections()[0][2]]
    assert "Steps" in labels and "Interval" in labels


def test_copy_is_independent():
    ts = scalar([1.0, 2.0])
    other = ts.copy()
    other.append(T0 + 5 * H, 3.0)
    other.values[0, 0] = 99.0
    assert len(ts) == 2 and ts.values[0, 0] == 1.0


def test_equality_is_identity_and_series_are_hashable():
    a, b = scalar([1.0, 2.0]), scalar([1.0, 2.0])
    assert a == a and a != b
    assert [a, b].index(b) == 1
    assert len({a, b}) == 2


def test_constructor_copies_its_input():
    arr = np.array([1.0, 2.0])
    ts = TimeSeries(times=hours(0, 1), values=arr)
    arr[0] = 100.0
    assert ts.values[0, 0] == 1.0
    readonly = np.array([1.0, 2.0])
    readonly.flags.writeable = False
    ts = TimeSeries(times=hours(0, 1), values=readonly)
    ts.values[0, 0] = 5.0
    assert ts.values[0, 0] == 5.0


def test_shallow_copy_has_its_own_storage():
    a = scalar([1.0, 2.0, 3.0])
    a.append(T0 + 3 * H, 4.0)  # leaves spare capacity
    b = copy.copy(a)
    b.append(T0 + 4 * H, 40.0)
    a.append(T0 + 5 * H, 50.0)
    assert b.end == T0 + 4 * H and b.values[-1, 0] == 40.0
    assert a.end == T0 + 5 * H and len(a) == len(b) == 5


def test_dataclasses_replace_keeps_the_data():
    ts = scalar([1.0, 2.0], role="observed")
    renamed = dataclasses.replace(ts, name="renamed")
    assert renamed.name == "renamed" and renamed.role == "observed"
    np.testing.assert_array_equal(renamed.values, ts.values)
    np.testing.assert_array_equal(renamed.times, ts.times)


def test_info_shows_a_readable_interval():
    rows = dict(scalar([1.0, 2.0])._info_sections()[0][2])
    assert rows["Interval"] == "1:00:00"


# Append


def test_append_single_steps():
    ts = TimeSeries(name="t", association="sample")
    ts.append(T0, 1.0)
    ts.append(dt.datetime(2024, 1, 1, 1), np.array([2.0]))
    ts.append("2024-01-01T02:00", 3)
    np.testing.assert_array_equal(ts.values[:, 0], [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(ts.times, hours(0, 1, 2))


def test_append_batch():
    ts = scalar([1.0])
    ts.append(hours(1, 2, 3), np.array([2.0, 3.0, 4.0]))
    assert len(ts) == 4
    vec = TimeSeries(times=hours(0), values=np.zeros((1, 2, 3)), dim=3)
    vec.append(hours(1, 2), np.ones((2, 2, 3)))
    assert vec.values.shape == (3, 2, 3)


def test_append_field():
    ts = TimeSeries(
        name="t",
        unit="C",
        times=hours(0),
        values=np.zeros((1, 3)),
        association="vertex",
    )
    ts.append(T0 + H, Field(unit="C", values=np.ones(3), association="vertex"))
    np.testing.assert_array_equal(ts.values[1], np.ones(3))
    with pytest.raises(ValueError, match="unit"):
        ts.append(T0 + 2 * H, Field(unit="K", values=np.ones(3), association="vertex"))
    with pytest.raises(ValueError, match="association"):
        ts.append(T0 + 2 * H, Field(unit="C", values=np.ones(3), association="face"))


@pytest.mark.parametrize(
    "times,values,error",
    [
        (T0, 9.0, ValueError),  # not after end
        (T0 - H, 9.0, ValueError),
        (hours(3, 2), np.array([9.0, 9.0]), ValueError),  # batch unsorted
        (T0 + 5 * H, np.array([9.0, 9.0]), ValueError),  # wrong element count
        (hours(3, 4), np.array([9.0]), ValueError),  # batch length mismatch
    ],
)
def test_rejected_append_leaves_series_unchanged(times, values, error):
    ts = scalar([1.0, 2.0])
    with pytest.raises(error):
        ts.append(times, values)
    assert len(ts) == 2
    np.testing.assert_array_equal(ts.values[:, 0], [1.0, 2.0])


def test_single_step_rejects_batch_shaped_values():
    ts = scalar(np.zeros((1, 2)))
    with pytest.raises(ValueError, match="shape"):
        ts.append(T0 + H, np.zeros((2, 2)))
    assert len(ts) == 1


def test_append_rejects_unsafe_cast():
    ts = TimeSeries(times=hours(0), values=np.array([1], dtype=np.int32))
    with pytest.raises(TypeError, match="cast"):
        ts.append(T0 + H, 1.5)
    ts.append(T0 + H, np.int64(2))
    assert ts.values.dtype == np.int32


@pytest.mark.parametrize(
    "dtype,value",
    [(np.int8, 300), (np.uint8, -1), (np.int32, 2**40), (np.float16, 1e6)],
)
def test_append_and_merge_reject_values_that_do_not_fit(dtype, value):
    ts = TimeSeries(times=hours(0), values=np.array([1], dtype=dtype))
    with pytest.raises(ValueError, match="fit|overflow"):
        ts.append(T0 + H, value)
    other = TimeSeries(times=hours(5), values=np.array([value]))
    with pytest.raises(ValueError, match="fit|overflow"):
        ts.merge(other)
    assert len(ts) == 1
    ts.append(T0 + H, 2)
    assert ts.values.dtype == dtype


def test_first_append_sets_element_count():
    ts = TimeSeries(dim=2)
    ts.append(T0, np.ones((4, 2)))
    assert ts.values.shape == (1, 4, 2)
    with pytest.raises(ValueError):
        ts.append(T0 + H, np.ones((3, 2)))


def test_append_does_not_grow_quadratically():
    def per_append(n):
        ts = TimeSeries()
        times = T0 + np.arange(n) * H
        start = time.perf_counter()
        for t in times:
            ts.append(t, 1.0)
        return (time.perf_counter() - start) / n

    per_append(1_000)  # warm up
    small = min(per_append(10_000) for _ in range(3))
    large = min(per_append(100_000) for _ in range(3))
    assert large < 3 * small


def test_empty_batch_append_is_a_no_op():
    ts = scalar([1.0, 2.0])
    ts.append(np.array([], dtype="datetime64[ns]"), np.empty((0, 1)))
    ts.append([], [])
    assert len(ts) == 2
    with pytest.raises(ValueError, match="empty"):
        ts.append([], [1.0])


# Merge


def test_merge_interleaves():
    ts = scalar([1.0, 3.0], offsets=[0, 2])
    ts.merge(scalar([2.0, 4.0], offsets=[1, 3]))
    np.testing.assert_array_equal(ts.values[:, 0], [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(ts.times, hours(0, 1, 2, 3))


@pytest.mark.parametrize(
    "policy,expected", [("keep_old", [1.0, 2.0, 30.0]), ("keep_new", [1.0, 20.0, 30.0])]
)
def test_merge_duplicate_policies(policy, expected):
    ts = scalar([1.0, 2.0])
    ts.merge(scalar([20.0, 30.0], offsets=[1, 2]), on_duplicate=policy)
    np.testing.assert_array_equal(ts.values[:, 0], expected)


def test_merge_raise_leaves_series_unchanged():
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="duplicate"):
        ts.merge(scalar([20.0], offsets=[1]))
    np.testing.assert_array_equal(ts.values[:, 0], [1.0, 2.0])
    with pytest.raises(ValueError, match="on_duplicate"):
        ts.merge(scalar([20.0], offsets=[5]), on_duplicate="newest")


def test_merge_rejects_incompatible_series():
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="unit"):
        ts.merge(TimeSeries(unit="K", times=hours(5), values=np.zeros(1)))
    with pytest.raises(ValueError, match="element"):
        ts.merge(scalar(np.zeros((1, 2)), offsets=[5]))
    with pytest.raises(TypeError):
        ts.merge(np.zeros(3))


def test_merge_into_empty_and_ignores_role():
    ts = TimeSeries(name="t", unit="C", association="sample", role="observed")
    ts.merge(scalar([1.0, 2.0], role="predicted"))
    assert len(ts) == 2 and ts.role == "observed"


# Access


def test_field_at_is_zero_copy_view():
    ts = TimeSeries(
        name="t",
        unit="C",
        description="d",
        times=hours(0, 1),
        values=np.zeros((2, 3)),
        association="vertex",
    )
    f = ts.field_at(1)
    assert isinstance(f, Field)
    assert (f.name, f.unit, f.description, f.dim, f.association) == (
        "t",
        "C",
        "d",
        1,
        "vertex",
    )
    f.values[0] = 7.0
    assert ts.values[1, 0] == 7.0
    assert ts.field_at(T0 + H).values[0] == 7.0
    assert ts.field_at(-1).values[0] == 7.0


def test_field_at_errors():
    ts = scalar([1.0, 2.0])
    with pytest.raises(KeyError):
        ts.field_at(T0 + H / 2)
    with pytest.raises(IndexError):
        ts.field_at(2)


def test_between_is_half_open_copy():
    ts = scalar([0.0, 1.0, 2.0, 3.0])
    part = ts.between(T0 + H, T0 + 3 * H)
    np.testing.assert_array_equal(part.values[:, 0], [1.0, 2.0])
    assert (part.name, part.unit, part.association) == ("temperature", "C", "sample")
    part.values[0, 0] = 99.0
    assert ts.values[1, 0] == 1.0
    assert len(ts.between(None, T0 + H)) == 1
    assert len(ts.between(T0 + 2 * H, None)) == 2


# Interpolation


def test_at_linear_and_exact_hits():
    ts = scalar([0.0, 10.0, np.nan, 30.0])
    out = ts.at(hours(0, 0.5, 1, 1.5, 3))
    np.testing.assert_array_equal(out[:, 0], [0.0, 5.0, 10.0, np.nan, 30.0])
    assert out.shape == (5, 1)
    assert ts.at(T0 + H / 2).shape == (1,)


def test_at_methods():
    ts = scalar([0.0, 10.0, 20.0], offsets=[0, 1, 3])
    q = hours(0.25, 0.75, 2)
    np.testing.assert_array_equal(ts.at(q, method="previous")[:, 0], [0.0, 0.0, 10.0])
    np.testing.assert_array_equal(ts.at(q, method="nearest")[:, 0], [0.0, 10.0, 10.0])
    pchip = ts.at(q, method="pchip")[:, 0]
    assert np.all(np.diff(pchip) > 0) and 0.0 < pchip[0] < 10.0 < pchip[2] < 20.0
    with pytest.raises(ValueError, match="method"):
        ts.at(q, method="cubic")


def test_pchip_nan_stays_local_to_its_column():
    values = np.column_stack([np.arange(6.0), np.arange(6.0)])
    values[3, 1] = np.nan
    ts = TimeSeries(times=hours(0, 1, 2, 3, 4, 5), values=values)
    out = ts.at(hours(0.5, 2.5), method="pchip")
    np.testing.assert_allclose(out[:, 0], [0.5, 2.5])
    assert np.isnan(out[1, 1])


def test_pchip_matches_scipy():
    from scipy.interpolate import PchipInterpolator

    offsets = np.array([0, 1, 2.5, 3, 5, 6])
    values = np.column_stack([[0, 1, 4, 4, 2, 3], [5, 3, 2, 1, 1, 0.5]]).astype(float)
    ts = TimeSeries(times=hours(*offsets), values=values)
    query = np.linspace(0, 6, 61)
    expected = PchipInterpolator(offsets, values, axis=0)(query)
    np.testing.assert_allclose(ts.at(hours(*query), method="pchip"), expected)
    two = TimeSeries(times=hours(0, 2), values=[0.0, 4.0])
    np.testing.assert_allclose(two.at(hours(0.5), method="pchip")[:, 0], [1.0])


def test_pchip_on_integer_and_three_step_series():
    from scipy.interpolate import PchipInterpolator

    ts = TimeSeries(times=hours(0, 1, 3), values=np.array([0, 4, 5], dtype=np.int16))
    query = np.array([0.5, 2.0, 2.9])
    expected = PchipInterpolator([0, 1, 3], [0.0, 4.0, 5.0])(query)
    np.testing.assert_allclose(ts.at(hours(*query), method="pchip")[:, 0], expected)


def test_at_fill_and_max_gap():
    ts = scalar([0.0, 10.0, 20.0], offsets=[0, 1, 5])
    outside = hours(-1, 6)
    assert np.isnan(ts.at(outside)).all()
    np.testing.assert_array_equal(ts.at(outside, fill="edge")[:, 0], [0.0, 20.0])
    gapped = ts.at(hours(0.5, 3), max_gap=np.timedelta64(2, "h"))[:, 0]
    np.testing.assert_array_equal(gapped, [5.0, np.nan])
    assert np.isnan(ts.at(hours(3), max_gap="2h", method="previous")).all()
    with pytest.raises(ValueError, match="fill"):
        ts.at(outside, fill="zero")


def test_at_integer_series_returns_float_and_vectors():
    ts = TimeSeries(times=hours(0, 2), values=np.array([[0, 10]], dtype=int).T)
    assert ts.at(hours(1)).dtype == np.float64
    vec = TimeSeries(
        times=hours(0, 2),
        values=np.stack([np.zeros((2, 3)), np.full((2, 3), 2.0)]),
        dim=3,
    )
    np.testing.assert_array_equal(vec.at(hours(1)), np.ones((1, 2, 3)))


def test_at_on_empty_and_single_step():
    empty = TimeSeries(dim=2)
    assert empty.at(hours(0, 1)).shape == (2, 0, 2)
    one = scalar([4.0])
    np.testing.assert_array_equal(one.at(hours(0, 1))[:, 0], [4.0, np.nan])
    np.testing.assert_array_equal(one.at(hours(1), fill="edge")[:, 0], [4.0])


# Resampling and aggregation


def test_resample_hows():
    ts = scalar([1.0, np.nan, 3.0, 5.0], offsets=[0, 0.5, 1.25, 3.5])
    expect = {
        "mean": [1.0, 3.0, np.nan, 5.0],
        "min": [1.0, 3.0, np.nan, 5.0],
        "max": [1.0, 3.0, np.nan, 5.0],
        "sum": [1.0, 3.0, np.nan, 5.0],
        "first": [1.0, 3.0, np.nan, 5.0],
        "last": [1.0, 3.0, np.nan, 5.0],
        "count": [1, 1, 0, 1],
    }
    for how, values in expect.items():
        out = ts.resample("1h", how=how)
        np.testing.assert_array_equal(out.times, hours(0, 1, 2, 3))
        np.testing.assert_array_equal(out.values[:, 0], values)
    assert ts.resample("1h", how="count").values.dtype == np.int64
    with pytest.raises(ValueError, match="how"):
        ts.resample("1h", how="median")


def test_resample_groups_and_aligns_to_epoch():
    ts = scalar([1.0, 2.0, 4.0, 8.0], offsets=[0.5, 1, 2.5, 3])
    out = ts.resample(dt.timedelta(hours=2), how="sum")
    np.testing.assert_array_equal(out.times, hours(0, 2))
    np.testing.assert_array_equal(out.values[:, 0], [3.0, 12.0])
    first = ts.resample("2h", how="first").values[:, 0]
    last = ts.resample("2h", how="last").values[:, 0]
    np.testing.assert_array_equal(first, [1.0, 4.0])
    np.testing.assert_array_equal(last, [2.0, 8.0])
    assert len(TimeSeries().resample("1h")) == 0
    with pytest.raises(ValueError, match="positive"):
        ts.resample("0h")


@pytest.mark.parametrize("freq", ["MS", "1H", "D", "1.5h", "1M"])
def test_duration_strings_reject_other_units_and_case(freq):
    with pytest.raises(ValueError, match="duration"):
        scalar([1.0, 2.0]).resample(freq)


@pytest.mark.parametrize("unit", ["M", "Y"])
def test_calendar_timedelta_is_rejected(unit):
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="fixed-length"):
        ts.resample(np.timedelta64(1, unit))
    with pytest.raises(ValueError, match="fixed-length"):
        ts.at(hours(0.5), max_gap=np.timedelta64(1, unit))


def test_reductions():
    ts = TimeSeries(
        times=hours(0, 1, 2),
        values=np.array([[1.0, np.nan], [3.0, np.nan], [5.0, 1.0]]),
    )
    np.testing.assert_array_equal(ts.mean(), [3.0, 1.0])
    np.testing.assert_array_equal(ts.min(), [1.0, 1.0])
    np.testing.assert_array_equal(ts.max(), [5.0, 1.0])
    np.testing.assert_array_equal(ts.sum(), [9.0, 1.0])
    assert np.isnan(TimeSeries(times=hours(0), values=[np.nan]).mean()).all()


# pandas


@pytest.mark.parametrize(
    "values,dim",
    [
        (np.array([1.0, 2.0, 3.0]), 1),
        (np.arange(6, dtype=np.int32).reshape(3, 2), 1),
        (np.arange(18.0).reshape(3, 2, 3), 3),
        (np.arange(9.0).reshape(3, 1, 3), 3),
    ],
)
def test_pandas_round_trip(values, dim):
    pytest.importorskip("pandas")
    ts = TimeSeries(name="x", times=hours(0, 1, 2), values=values, dim=dim)
    obj = ts.to_pandas()
    assert str(obj.index.tz) == "UTC"
    back = TimeSeries.from_pandas(obj, unit="C", association="sample")
    np.testing.assert_array_equal(back.times, ts.times)
    np.testing.assert_array_equal(back.values, ts.values)
    assert back.values.dtype == ts.values.dtype
    assert (back.dim, back.unit, back.association) == (dim, "C", "sample")


def test_from_pandas_orders_columns_by_label():
    pytest.importorskip("pandas")
    ts = TimeSeries(times=hours(0, 1), values=np.arange(12.0).reshape(2, 2, 3), dim=3)
    shuffled = ts.to_pandas().sort_index(axis=1, level=1)
    np.testing.assert_array_equal(TimeSeries.from_pandas(shuffled).values, ts.values)


def test_from_pandas_component_subset_and_incomplete_columns():
    pytest.importorskip("pandas")
    ts = TimeSeries(times=hours(0, 1), values=np.arange(12.0).reshape(2, 2, 3), dim=3)
    frame = ts.to_pandas()
    subset = TimeSeries.from_pandas(frame.loc[:, (slice(None), [0, 1])])
    assert subset.dim == 2
    np.testing.assert_array_equal(subset.values, ts.values[:, :, :2])
    with pytest.raises(ValueError, match="exactly once"):
        TimeSeries.from_pandas(frame.iloc[:, :-1])


@pytest.mark.parametrize("elements", [1, 3])
@pytest.mark.parametrize("steps", [0, 2])
def test_from_pandas_single_component_has_scalar_shape(elements, steps):
    pytest.importorskip("pandas")
    values = np.arange(steps * elements * 3, dtype=np.int32).reshape(steps, elements, 3)
    ts = TimeSeries(times=hours(*range(steps)), values=values, dim=3)
    frame = ts.to_pandas().loc[:, (slice(None), [1])].iloc[:, ::-1]
    back = TimeSeries.from_pandas(
        frame, unit="m/s", association="vertex", role="observed"
    )
    assert back.dim == 1 and back.values.shape == (steps, elements)
    assert (back.unit, back.association, back.role) == ("m/s", "vertex", "observed")
    assert back.values.dtype == values.dtype
    np.testing.assert_array_equal(back.values, values[:, :, 1])
    np.testing.assert_array_equal(back.times, ts.times)


def test_pandas_round_trip_of_empty_vector_series():
    pytest.importorskip("pandas")
    ts = TimeSeries(dim=3)
    back = TimeSeries.from_pandas(ts.to_pandas())
    assert back.values.shape == (0, 0, 3) and back.dim == 3


def test_to_pandas_scalar_is_named_series():
    pd = pytest.importorskip("pandas")
    obj = scalar([1.0, 2.0]).to_pandas()
    assert isinstance(obj, pd.Series) and obj.name == "temperature"
    assert TimeSeries.from_pandas(obj).name == "temperature"


def test_serialization_is_not_supported_yet():
    with pytest.raises(NotImplementedError):
        scalar([1.0]).to_proto()


# Range and validation edge cases

LONG_SPAN = ["1700-01-01", "2200-01-01"]
LONG_DAYS = (dt.date(2200, 1, 1) - dt.date(1700, 1, 1)).days


@pytest.mark.parametrize("method", ["linear", "pchip"])
def test_interpolation_across_spans_longer_than_int64_nanoseconds(method):
    ts = TimeSeries(times=LONG_SPAN, values=[0.0, 10.0])
    mid = "1950-01-01"
    expected = 10.0 * (dt.date(1950, 1, 1) - dt.date(1700, 1, 1)).days / LONG_DAYS
    np.testing.assert_allclose(ts.at(mid, method=method), [expected])
    assert np.isnan(ts.at(mid, method=method, max_gap="1w")).all()


def test_pchip_across_spans_longer_than_int64_nanoseconds():
    ts = TimeSeries(times=["1680-01-01", "2000-01-01", "2250-01-01"], values=[0, 5, 10])
    value = ts.at("1840-01-01", method="pchip")[0]
    assert 0.0 < value < 5.0


def test_is_regular_and_info_on_spans_longer_than_int64_nanoseconds():
    ts = TimeSeries(times=LONG_SPAN, values=[0.0, 10.0])
    regular, interval = ts.is_regular()
    assert regular and interval == np.timedelta64(LONG_DAYS, "D")
    rows = dict(ts._info_sections()[0][2])
    assert rows["Interval"] == str(dt.timedelta(days=LONG_DAYS))
    uneven = TimeSeries(
        times=["1700-01-01", "1701-01-01", "2200-01-01"], values=[0, 1, 2]
    )
    assert uneven.is_regular() == (False, None)
    assert dict(uneven._info_sections()[0][2])["Interval"] == "irregular"


@pytest.mark.parametrize("fill", [None, "edge"])
@pytest.mark.parametrize("method", ["linear", "previous", "nearest", "pchip"])
def test_nat_query_times_give_nan(fill, method):
    ts = scalar([1.0, 2.0, 3.0])
    out = ts.at([T0, np.datetime64("NaT", "ns"), "NaT"], method=method, fill=fill)[:, 0]
    np.testing.assert_array_equal(out, [1.0, np.nan, np.nan])
    assert np.isnan(ts.at("NaT", method=method, fill=fill)).all()


def test_between_rejects_nat_bounds():
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="NaT"):
        ts.between("NaT", T0 + H)
    with pytest.raises(ValueError, match="NaT"):
        ts.between(T0, np.datetime64("NaT", "ns"))


@pytest.mark.parametrize(
    "times",
    [
        hours(0, 1, 2).reshape(-1, 1),
        np.array([["2024-01-01", "2024-01-02"]], dtype=object),
    ],
)
def test_multidimensional_times_are_rejected(times):
    with pytest.raises(ValueError, match="1-D"):
        TimeSeries(times=times, values=np.zeros(times.size))
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="1-D"):
        ts.at(times)
    with pytest.raises(ValueError, match="1-D"):
        ts.append(times, np.zeros(times.size))
    with pytest.raises(ValueError, match="1-D"):
        ts.times = times
    assert len(ts) == 2


def test_nat_durations_are_rejected():
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="max_gap must not be NaT"):
        ts.at(hours(0.5), max_gap=np.timedelta64("NaT", "s"))
    with pytest.raises(ValueError, match="freq must not be NaT"):
        ts.resample(np.timedelta64("NaT", "h"))


def test_pandas_nat_duration_is_rejected():
    pd = pytest.importorskip("pandas")
    with pytest.raises(ValueError, match="freq must not be NaT"):
        scalar([1.0, 2.0]).resample(pd.NaT)
    out = scalar([1.0, 2.0]).resample(pd.Timedelta("2h"))
    np.testing.assert_array_equal(out.values[:, 0], [1.5])


@pytest.mark.parametrize(
    "duration",
    [
        "100000w",
        "9" * 30 + "s",
        np.timedelta64(10**17, "s"),
        dt.timedelta(days=999_999),
    ],
)
def test_durations_outside_the_nanosecond_range_raise_value_error(duration):
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="range"):
        ts.resample(duration)
    with pytest.raises(ValueError, match="range"):
        ts.at(hours(0.5), max_gap=duration)


@pytest.mark.parametrize(
    "dtype,value",
    [(np.float64, 2**53 + 1), (np.float32, 2**24 + 1), (np.float64, -(2**63) + 1)],
)
def test_integers_that_would_round_are_rejected(dtype, value):
    ts = TimeSeries(times=hours(0), values=np.array([1.0], dtype=dtype))
    with pytest.raises(ValueError, match="exactly"):
        ts.append(T0 + H, np.array([value]))
    other = TimeSeries(times=hours(5), values=np.array([value]))
    with pytest.raises(ValueError, match="exactly"):
        ts.merge(other)
    assert len(ts) == 1
    ts.append(T0 + H, np.array([2**53 + 2 if dtype == np.float64 else 2**24 + 2]))
    assert len(ts) == 2


def test_metadata_assignment_is_validated():
    ts = scalar([1.0, 2.0])
    with pytest.raises(ValueError, match="dim"):
        ts.dim = 3
    with pytest.raises(ValueError, match="dim"):
        ts.copy(dim=4)
    with pytest.raises(ValueError, match="positive integer"):
        ts.dim = 0
    with pytest.raises(ValueError, match="association"):
        ts.association = "bogus"
    assert ts.dim == 1 and ts.association == "sample"
    ts.dim = 1
    ts.association = "vertex"
    assert ts.copy(association=None).association is None


def test_dim_of_an_empty_series_can_change():
    ts = TimeSeries()
    ts.dim = 3
    assert ts.values.shape == (0, 0, 3)
    ts.append(T0, np.ones((2, 3)))
    assert ts.values.shape == (1, 2, 3)
    shaped = TimeSeries(times=hours(), values=np.zeros((0, 2)))
    with pytest.raises(ValueError, match="dim"):
        shaped.dim = 3


def test_resample_bins_before_the_time_range_raise():
    first = "1677-09-21T00:12:43.145224193"
    ts = TimeSeries(times=[first, "1677-09-22"], values=[1.0, 2.0])
    with pytest.raises(ValueError, match="range"):
        ts.resample("1d")
    near = TimeSeries(times=["1677-09-21T00:13", "1677-09-21T00:15"], values=[1, 2])
    out = near.resample("1min")
    assert out.times[0] == np.datetime64("1677-09-21T00:13", "ns")
    np.testing.assert_array_equal(out.values[:, 0], [1.0, np.nan, 2.0])
