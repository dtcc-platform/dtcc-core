"""Time series attached to geometries: display and exporters without a mapping."""

import numpy as np
import pytest

from dtcc_core.model import (
    DeSO,
    FieldSlice,
    GeometryType,
    LineString,
    Mesh,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    RoadNetwork,
    StreamlineCollection,
    Surface,
    TimeSeries,
    VolumeMesh,
)
from dtcc_core.model.object import SensorCollection

SQUARE = [(320000, 6399000), (320100, 6399000), (320100, 6399100), (320000, 6399100)]


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


# Exporters that write fields refuse time series instead of dropping them


def _sensors():
    collection = SensorCollection()
    station = Object(attributes={"station_id": "s1"})
    station.add_geometry(
        Point(x=320000.0, y=6399000.0, time_series=[series()]), "location"
    )
    collection.add_station(station)
    collection.transform.srs = "EPSG:3006"
    return collection


def _slice():
    return FieldSlice(
        points=np.array([[320000.0, 6399000.0, 10.0]]),
        crs="EPSG:3006",
        time_series=[series()],
    )


def _streamlines():
    line = LineString(
        vertices=np.array([[320000.0, 6399000.0, 1.0], [320010.0, 6399020.0, 2.0]]),
        time_series=[series(2, "vertex")],
    )
    return StreamlineCollection(lines=[line], crs="EPSG:3006")


def _road():
    road = RoadNetwork()
    road.vertices = np.array([SQUARE[0], SQUARE[1]], dtype=float)
    road.edges = np.array([[0, 1]])
    road.transform.srs = "EPSG:3006"
    line = MultiLineString(
        linestrings=[LineString(vertices=road.vertices)],
        time_series=[series(1, "sample")],
    )
    road.add_geometry(line, GeometryType.MULTILINESTRING)
    return road


def _deso():
    deso = DeSO()
    deso.transform.srs = "EPSG:3006"
    area = Object(attributes={"desokod": "1480A0001"})
    area.id = "1480A0001"
    surface = Surface(vertices=np.array([[x, y, 0.0] for x, y in SQUARE]))
    area.add_geometry(
        MultiSurface(surfaces=[surface], time_series=[series(1, "geometry")]),
        GeometryType.LOD0,
    )
    deso.add_child(area)
    return deso


@pytest.mark.parametrize(
    "make", [_sensors, _slice, _streamlines, _road, _deso], ids=lambda f: f.__name__
)
def test_geojson_exports_refuse_time_series(make):
    model = make()
    with pytest.raises(NotImplementedError, match="Time series .* GeoJSON"):
        model.to_geojson()


def _mesh():
    return Mesh(
        vertices=np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        faces=np.array([[0, 1, 2]]),
        time_series=[series(3, "vertex")],
    )


def _volume_mesh():
    return VolumeMesh(
        vertices=np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        cells=np.array([[0, 1, 2, 3]]),
        time_series=[series(1, "cell")],
    )


@pytest.mark.parametrize(
    "make, suffix",
    [(_mesh, ".vtk"), (_mesh, ".vtu"), (_volume_mesh, ".vtu"), (_volume_mesh, ".xdmf")],
)
def test_mesh_files_refuse_time_series(tmp_path, make, suffix):
    mesh = make()
    path = tmp_path / f"mesh{suffix}"
    with pytest.raises(NotImplementedError, match="Time series"):
        mesh.save(path)
    assert not path.exists()
    mesh.time_series = []
    mesh.save(path)
    assert path.exists()
