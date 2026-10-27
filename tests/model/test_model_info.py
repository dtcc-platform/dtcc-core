"""Public compact-display and explicit-inspection contracts."""

import inspect
import re
from io import StringIO

import numpy as np
import pytest

from dtcc_core import model
from dtcc_core.model import (
    Bounds,
    Building,
    City,
    DatasetCollection,
    Field,
    GeometryRepresentation,
    Mesh,
    Object,
    Point,
    PointCloud,
    Raster,
    SemanticRegion,
    SensorCollection,
    VehicleCollection,
)
from dtcc_core.model.model import Model


def test_all_concrete_models_share_display_contract(capsys):
    classes = {
        getattr(model, name)
        for name in model.__all__
        if inspect.isclass(getattr(model, name))
    }
    for cls in sorted(classes, key=lambda cls: cls.__name__):
        if not issubclass(cls, Model) or inspect.isabstract(cls):
            continue
        obj = cls()
        text = repr(obj)
        assert str(obj) == text
        assert "\n" not in text
        assert len(text) < 250
        assert text.lstrip("<").startswith(f"{cls.__name__}(")
        report = obj.info(print=False)
        assert report.startswith(cls.__name__ + "\n")
        assert "│" in report
        assert "\x1b" not in report
    assert capsys.readouterr().out == ""


def test_complete_small_values_and_omitted_state_markers():
    for obj in (Point(x=1.25, y=2, z=np.float64(3)), Bounds(0, 1, 2, 3, 4, 5)):
        text = repr(obj)
        assert not text.startswith("<")
        restored = eval(text, {"Point": Point, "Bounds": Bounds})
        assert restored._summary_items() == obj._summary_items()
    point = Point(x=1, fields=[Field(name="temperature", values=np.array([21.0]))])
    assert repr(point).startswith("<Point(")
    assert "num_fields=1" in repr(point)
    point = Point()
    point.transform.set_translation(1, 0, 0)
    assert repr(point).startswith("<Point(")
    point = Point()
    point.transform.srs = "EPSG:3006"
    assert repr(point).startswith("<Point(")
    for obj in (Point(x=np.nan), Bounds(xmax=np.inf), Point(x=10**120)):
        assert repr(obj).startswith("<")
    obj = Bounds()
    obj.dataset_context = object()
    assert repr(obj).startswith("<Bounds(")


def test_repr_does_not_scan_arrays_or_expand_children(monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("Compact representation must not inspect geometry or children")

    city = City(id="example")
    city.children[Building] = [Building(id="child")] * 10000
    pc = PointCloud(points=np.zeros((100000, 3)))
    original_bounds = pc._bounds
    monkeypatch.setattr(PointCloud, "calculate_bounds", unexpected)
    monkeypatch.setattr(Building, "__repr__", unexpected)
    assert "num_children=10000" in repr(city)
    assert len(repr(city)) < 150
    assert repr(pc) == "<PointCloud(num_points=100000, num_fields=0)>"
    assert pc._bounds is original_bounds
    sensors = SensorCollection()
    monkeypatch.setattr(SensorCollection, "stations", unexpected)
    assert repr(sensors) == "<SensorCollection(num_stations=0)>"


def test_info_print_return_and_file_are_identical(capsys):
    mesh = Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        faces=np.array([[0, 1, 2]]),
        fields=[
            Field(
                name="temperature",
                unit="°C",
                association="vertex",
                values=np.array([20.0, 21.0, 22.0]),
            )
        ],
    )
    text = mesh.info(print=False)
    assert capsys.readouterr().out == ""
    assert re.search(r"│Vertices\s*│3\s*│", text)
    assert re.search(r"│Faces\s*│1\s*│", text)
    assert "temperature" in text and "°C" in text and "vertex" in text
    assert mesh.info() is None
    assert capsys.readouterr().out == text + "\n"
    stream = StringIO()
    assert mesh.print_info(file=stream) is None
    assert stream.getvalue() == text + "\n"
    assert capsys.readouterr().out == ""


def test_info_treats_attributes_as_literal_text_and_limits_rows():
    obj = Object(id="literal", attributes={"[bold]": "[red]literal[/red]"})
    text = obj.info(print=False)
    assert "[bold]" in text and "[red]literal[/red]" in text
    obj.attributes = {f"attribute_{i}": i for i in range(100)}
    text = obj.info(print=False)
    assert "attribute_19" in text
    assert "attribute_20" not in text
    assert "80 more row(s)" in text
    samples = DatasetCollection(items=list(range(100))).info(print=False)
    assert "20 shown, 80 omitted" in samples


def test_sensor_details_and_vehicle_guidance_live_in_info():
    sensors = SensorCollection(attributes={"phenomenon": "temperature"})
    station = Object(id="station", attributes={"station_name": "Central"})
    station.add_geometry(
        Point(
            x=1,
            y=2,
            fields=[Field(name="temperature", unit="°C", values=np.array([21.0]))],
        ),
        id="location",
    )
    sensors.add_station(station)
    text = sensors.info(print=False)
    assert "Measurement statistics: temperature" in text
    assert "21.00" in text and "Central" in text and "°C" in text
    assert "Central" not in repr(sensors)
    vehicles = VehicleCollection(
        attributes={
            "partial_result": True,
            "configuration_help": ["Set [API_KEY] to enable the provider."],
            "upstream_errors": [
                {"failure_class": "credentials", "message": "Missing key"}
            ],
        }
    )
    text = vehicles.info(print=False)
    assert "Set [API_KEY]" in text and "Missing key" in text
    assert "configuration_help" not in repr(vehicles)


def test_records_do_not_expand_geometry_or_region_indices():
    record = GeometryRepresentation(Mesh(vertices=np.zeros((10000, 3))), lod="2.2")
    region = SemanticRegion("urn:roof", indices=np.arange(10000))
    assert str(record) == repr(record)
    assert str(region) == repr(region)
    assert "geometry_type='Mesh'" in repr(record)
    assert "num_elements=10000" in repr(region)
    assert len(repr(record)) < 150 and len(repr(region)) < 150


def test_tree_short_and_long_modes_cover_mixed_attachments(capsys):
    root = Object(id="root", attributes={"note": "first\nsecond"})
    child = Object(id="child")
    dem = Raster(data=np.zeros((2, 3)))
    extent = Bounds(xmax=3, ymax=2)
    sample = Point(
        fields=[Field(name="temperature", unit="°C", description="Air temperature")]
    )
    child.add_geometry(dem, id="dem", role="elevation")
    child.add_geometry(extent, id="extent")
    child.add_geometry(sample, id="sample")
    root.add_child(child)
    last = Object(id="last")
    root.add_child(last)

    assert root.tree() is None
    assert capsys.readouterr().out.splitlines() == [
        repr(root),
        f"├── {child!r}",
        f"│   ├── 'dem': {dem!r}",
        f"│   ├── 'extent': {extent!r}",
        f"│   └── 'sample': {sample!r}",
        f"└── {last!r}",
    ]

    assert root.tree(verbose=True) is None
    text = capsys.readouterr().out
    assert "├── attribute 'note': 'first\\nsecond'" in text
    assert "│   ├── 'dem' (role='elevation'): <Raster(" in text
    assert "│       └── field 'temperature' (°C): Air temperature" in text
    assert text.endswith(f"└── {last!r}\n")

    root.tree(max_depth=1)
    assert capsys.readouterr().out.splitlines() == [
        repr(root),
        f"├── {child!r} ...",
        f"└── {last!r}",
    ]
    root.tree(max_depth=0)
    assert capsys.readouterr().out == repr(root) + " ...\n"
    dem.tree()
    assert capsys.readouterr().out == repr(dem) + "\n"


def test_tree_rejects_invalid_depth_without_printing(capsys):
    with pytest.raises(ValueError, match="nonnegative"):
        Object().tree(max_depth=-1)
    with pytest.raises(TypeError, match="integer"):
        Object().tree(max_depth=True)
    assert capsys.readouterr().out == ""
