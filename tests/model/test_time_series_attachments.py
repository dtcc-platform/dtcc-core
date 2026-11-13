"""Time series attached to geometries: display."""

import numpy as np

from dtcc_core.model import Mesh, Object, Point, TimeSeries


def series(elements=1, association="sample", name="temperature"):
    return TimeSeries(
        name,
        "C",
        "hourly mean",
        times=["2024-01-01T00:00", "2024-01-01T01:00"],
        values=np.zeros((2, elements)),
        association=association,
        role="observed",
    )


# Display


def test_point_with_time_series_is_not_shown_as_a_complete_literal():
    point = Point(x=1.0, time_series=[series()])
    assert "num_time_series=1" in repr(point)
    assert repr(point) != repr(Point(x=1.0))
    assert dict(point._summary_items())["num_time_series"] == 1


def test_info_lists_attached_time_series():
    mesh = Mesh(vertices=np.zeros((3, 3)), time_series=[series(3, "vertex")])
    sections = {title: rows for title, _, rows in mesh._info_sections()}
    (row,) = sections["Time series"]
    ts = mesh.time_series[0]
    assert row == (
        "temperature",
        "C",
        "vertex",
        "observed",
        (2, 3),
        ts.start,
        ts.end,
        np.dtype("float64"),
        "hourly mean",
    )
    assert "Time series" not in {s[0] for s in Mesh()._info_sections()}


def test_verbose_tree_lists_attached_time_series(capsys):
    obj = Object(id="station")
    obj.add_geometry(Point(time_series=[series()]), "location")
    obj.tree(verbose=True)
    out = capsys.readouterr().out
    assert "time series 'temperature' (C, observed): hourly mean" in out
