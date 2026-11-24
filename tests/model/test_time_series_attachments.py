"""Attached time series in display, ordinary exports and canonical storage."""

from copy import deepcopy

import h5py
import meshio
import numpy as np
import pygltflib
import pytest

from dtcc_core.io import load_model
from dtcc_core.model import (
    DeSO,
    Field,
    FieldSlice,
    GeometryType,
    LineString,
    Mesh,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    PointCloud,
    RoadNetwork,
    SemanticRegion,
    StreamlineCollection,
    Surface,
    TimeSeries,
    VolumeMesh,
    exchange,
)
from dtcc_core.model.object import SensorCollection, VehicleCollection

SQUARE = [(320000, 6399000), (320100, 6399000), (320100, 6399100), (320000, 6399100)]


def series(elements=1, association="sample", name="temperature"):
    return TimeSeries(
        name,
        "C",
        "hourly mean",
        times=["2024-01-01T00:00", "2024-01-01T01:00"],
        values=np.arange(2 * elements, dtype=float).reshape(2, elements),
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


# Ordinary exports keep their documented subset and leave attachments unchanged


def _sensors():
    collection = SensorCollection()
    station = Object(attributes={"station_id": "s1"})
    point = Point(x=320000.0, y=6399000.0, time_series=[series()])
    station.add_geometry(point, "location")
    collection.add_station(station)
    collection.transform.srs = "EPSG:3006"
    return collection, point


def _vehicles():
    collection = VehicleCollection()
    vehicle = Object(attributes={"vehicle_id": "v1"})
    point = Point(x=320000.0, y=6399000.0, time_series=[series()])
    vehicle.add_geometry(point, "location")
    collection.add_vehicle(vehicle)
    collection.transform.srs = "EPSG:3006"
    return collection, point


def _slice():
    model = FieldSlice(
        points=np.array([[320000.0, 6399000.0, 10.0]]),
        crs="EPSG:3006",
        time_series=[series()],
    )
    return model, model


def _streamlines():
    line = LineString(
        vertices=np.array([[320000.0, 6399000.0, 1.0], [320010.0, 6399020.0, 2.0]]),
        time_series=[series(2, "vertex")],
    )
    return StreamlineCollection(lines=[line], crs="EPSG:3006"), line


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
    return road, line


def _deso():
    deso = DeSO()
    deso.transform.srs = "EPSG:3006"
    area = Object(attributes={"desokod": "1480A0001"})
    area.id = "1480A0001"
    surface = Surface(vertices=np.array([[x, y, 0.0] for x, y in SQUARE]))
    geometry = MultiSurface(surfaces=[surface], time_series=[series(1, "geometry")])
    area.add_geometry(geometry, GeometryType.LOD0)
    deso.add_child(area)
    return deso, geometry


def _map_mesh():
    mesh = _mesh()
    mesh.transform.srs = "EPSG:3006"
    return mesh, mesh


@pytest.mark.parametrize(
    "make",
    [_sensors, _vehicles, _slice, _streamlines, _road, _deso, _map_mesh],
    ids=lambda f: f.__name__,
)
def test_geojson_exports_keep_their_subset_with_time_series(make):
    model, owner = make()
    baseline, baseline_owner = deepcopy((model, owner))
    baseline_owner.time_series = []
    before = exchange.dumps(model)
    attached = owner.time_series[0]

    payload = model.to_geojson()

    assert payload == baseline.to_geojson()
    assert payload["type"] == "FeatureCollection"
    assert len(payload["features"]) == 1
    assert owner.time_series[0] is attached
    assert exchange.dumps(model) == before


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


@pytest.mark.parametrize("attached_input", ["receiver", "incoming", "both"])
def test_point_cloud_merge_rejects_series_without_mutating_either_input(attached_input):
    first = PointCloud(points=np.array([[0.0, 0.0, 0.0]]))
    second = PointCloud(points=np.array([[1.0, 0.0, 0.0]]))
    if attached_input in ("receiver", "both"):
        first.time_series = [series()]
    if attached_input in ("incoming", "both"):
        second.time_series = [series()]
    before = exchange.dumps(first), exchange.dumps(second)

    with pytest.raises(NotImplementedError, match="time series"):
        first.merge(second)
    assert (exchange.dumps(first), exchange.dumps(second)) == before

    if second.time_series:
        empty = PointCloud()
        with pytest.raises(NotImplementedError, match="time series"):
            empty.merge(second)
        assert len(empty) == 0
    else:
        assert first.merge(PointCloud()) is first
        assert exchange.dumps(first) == before[0]


@pytest.mark.parametrize("attached_input", ["receiver", "incoming"])
def test_multisurface_merge_rejects_collection_series_without_mutation(attached_input):
    first = MultiSurface(surfaces=[Surface(vertices=_mesh().vertices)])
    second = MultiSurface(surfaces=[Surface(vertices=_mesh().vertices + [0, 0, 1])])
    owner = first if attached_input == "receiver" else second
    owner.time_series = [series(1, "face")]
    before = exchange.dumps(first), exchange.dumps(second)

    with pytest.raises(NotImplementedError, match="time series"):
        first.merge(second)
    assert (exchange.dumps(first), exchange.dumps(second)) == before


def test_multisurface_merge_preserves_series_on_the_incoming_surface():
    surface = Surface(vertices=_mesh().vertices, time_series=[series(1, "geometry")])
    result = MultiSurface().merge(MultiSurface(surfaces=[surface]))
    assert exchange.dumps(result.surfaces[0]) == exchange.dumps(surface)
    assert len(exchange.loads(exchange.dumps(result)).surfaces[0].time_series) == 1


def test_mesh_processing_rejects_series_before_losing_it():
    from dtcc_core.builder.meshing.meshing import disjoint_meshes, merge_meshes
    from dtcc_core.builder.model_conversion import (
        mesh_to_builder_mesh,
        volume_mesh_to_builder_volume_mesh,
    )

    mesh = _mesh()
    before = exchange.dumps(mesh)
    for operation in (
        lambda: merge_meshes([mesh]),
        lambda: mesh.snap_vertices(0.01),
        lambda: disjoint_meshes(mesh),
        lambda: mesh_to_builder_mesh(mesh),
        lambda: volume_mesh_to_builder_volume_mesh(_volume_mesh()),
    ):
        with pytest.raises(NotImplementedError, match="time series"):
            operation()
    assert exchange.dumps(mesh) == before


@pytest.mark.parametrize("with_regions", [False, True])
@pytest.mark.parametrize("attachment", ["collection", "surface"])
def test_collection_meshing_rejects_series_before_cleaning(with_regions, attachment):
    from dtcc_core.builder.meshing.meshing import mesh_multisurface, mesh_multisurfaces

    surface = Surface(vertices=_mesh().vertices)
    collection = MultiSurface(surfaces=[surface])
    if with_regions:
        collection.regions = [SemanticRegion("urn:roof", indices=np.array([0]))]
    owner = collection if attachment == "collection" else surface
    owner.time_series = [series(1, "geometry")]
    before = exchange.dumps(collection)

    with pytest.raises(NotImplementedError, match="time series"):
        mesh_multisurface(collection, clean=True)
    with pytest.raises(NotImplementedError, match="time series"):
        mesh_multisurfaces([collection], clean=True)
    if owner is surface:
        with pytest.raises(NotImplementedError, match="time series"):
            surface.mesh(clean=True)
    assert exchange.dumps(collection) == before


@pytest.mark.parametrize(
    "make, suffix",
    [
        *[
            (_mesh, suffix)
            for suffix in (".obj", ".ply", ".stl", ".vtk", ".vtu", ".xdmf")
        ],
        (_volume_mesh, ".vtu"),
        (_volume_mesh, ".xdmf"),
    ],
)
@pytest.mark.parametrize("center_on_origin", [False, True])
def test_mesh_exports_write_geometry_with_time_series(
    tmp_path, make, suffix, center_on_origin
):
    mesh = make()
    path = tmp_path / f"mesh{suffix}"
    before = exchange.dumps(mesh)

    mesh.save(path, center_on_origin=center_on_origin)

    kind = "triangle" if isinstance(mesh, Mesh) else "tetra"
    if suffix == ".xdmf":
        with h5py.File(path.with_suffix(".h5"), "r") as saved:
            saved_vertices = np.asarray(saved["Mesh/mesh/geometry"])
            saved_cells = np.asarray(saved["Mesh/mesh/topology"])
        xdmf = path.read_text()
        assert f"{path.stem}.h5:/Mesh/mesh/geometry" in xdmf
        assert f"{path.stem}.h5:/Mesh/mesh/topology" in xdmf
    else:
        saved = meshio.read(path)
        saved_vertices = saved.points
        saved_cells = saved.cells_dict[kind]
    vertices = mesh.vertices.copy()
    if center_on_origin:
        vertices[:, :2] -= (
            vertices[:, :2].min(axis=0) + vertices[:, :2].max(axis=0)
        ) / 2
    np.testing.assert_allclose(saved_vertices, vertices)
    np.testing.assert_array_equal(
        saved_cells, mesh.faces if kind == "triangle" else mesh.cells
    )
    assert exchange.dumps(mesh) == before


@pytest.mark.parametrize("suffix", [".gltf", ".gltf2", ".glb"])
def test_gltf_exports_write_geometry_with_time_series(tmp_path, suffix):
    mesh = _mesh()
    path = tmp_path / f"mesh{suffix}"
    before = exchange.dumps(mesh)

    mesh.save(path)

    saved = pygltflib.GLTF2.load(path)
    primitive = saved.meshes[0].primitives[0]
    assert saved.accessors[primitive.attributes.POSITION].count == len(mesh.vertices)
    assert saved.accessors[primitive.indices].count == mesh.faces.size
    assert exchange.dumps(mesh) == before


@pytest.mark.parametrize("make", [_mesh, _volume_mesh])
def test_native_mesh_files_preserve_time_series_and_fields(tmp_path, make):
    mesh = make()
    mesh.add_field(
        Field(
            name="static_temperature",
            unit="C",
            values=np.arange(len(mesh.vertices), dtype=float),
            association="vertex",
        )
    )
    before = exchange.dumps(mesh)
    path = tmp_path / "mesh.dtcc"

    mesh.save(path)
    restored = load_model(path, expected_type=type(mesh))

    assert path.exists()
    assert len(restored.time_series) == 1
    assert len(restored.fields) == 1
    assert exchange.dumps(restored) == before
    assert exchange.dumps(mesh) == before


@pytest.mark.parametrize("suffix", [".dae", ".fbx"])
def test_assimp_export_reports_unimplemented_with_time_series(
    tmp_path, monkeypatch, suffix
):
    from dtcc_core.io import meshes

    monkeypatch.setitem(meshes._save_formats[Mesh], suffix, meshes._save_assimp_mesh)
    path = tmp_path / f"mesh{suffix}"
    mesh = _mesh()
    before = exchange.dumps(mesh)

    with pytest.raises(RuntimeError, match="Not implemented"):
        mesh.save(path)

    assert not path.exists()
    assert exchange.dumps(mesh) == before
