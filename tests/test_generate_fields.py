"""Exercise the synthetic fields CLI and its native/ParaView agreement."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import meshio
import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import Bounds, City, VolumeGrid, VolumeMesh


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "generate_fields.py"


def test_cli_defaults_use_twice_the_equal_scale_vertical_count(monkeypatch):
    spec = importlib.util.spec_from_file_location("generate_fields", SCRIPT)
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    bounds = Bounds(xmin=83570., ymin=445800., zmin=-13.65049923706055,
                    xmax=85570., ymax=447800., zmax=263.5015007629395)
    monkeypatch.setattr(script, "flagship_domain", lambda _: (bounds, "EPSG:7415"))
    calls = []
    monkeypatch.setattr(script, "generate", lambda *args: calls.append(args))
    monkeypatch.setattr(sys, "argv", [str(SCRIPT)])
    script.main()
    assert calls[0][:3] == (128, 128, 36)
    assert calls[0][2] == 2*int(np.ceil(128*bounds.depth/bounds.width))


def test_fields_cli_round_trip_and_time_series(tmp_path):
    # Reuse an already-expanded flagship domain exactly, without doubling it again.
    bounds = dict(xmin=83570., ymin=445800., zmin=-17.,
                  xmax=85570., ymax=447800., zmax=283.)
    city = City(id="flagship-test")
    city.transform.srs = "EPSG:7415"
    domain = VolumeGrid(width=1, height=1, depth=1)
    domain.bounds = Bounds(**bounds)
    domain.transform.srs = city.transform.srs
    city.add_geometry(domain, id="air_grid")
    flagship = tmp_path / "flagship.dtcc"
    io.save_model(city, flagship)
    output = tmp_path / "fields"
    subprocess.run([sys.executable, str(SCRIPT), "-nx", "3", "-ny", "4", "-nz", "5",
                    "-nt", "5", "-T", "8", "--output-dir", str(output),
                    "--flagship-model", str(flagship)],
                   cwd=SCRIPT.parent, check=True, capture_output=True, text=True)
    series = json.loads((output / "fields-series.json").read_text())
    frames = series["frames"]
    assert series["bounds"] == bounds and series["crs"] == "EPSG:7415"
    assert [frame["time_seconds"] for frame in frames] == [0., 2., 4., 6., 8.]
    assert len(list(output.glob("*.dtcc"))) == 10
    first_fields = {}
    for name, native_type, cell_type in (("tet", VolumeMesh, "tetra"),
                                         ("grid", VolumeGrid, "hexahedron")):
        datasets = ET.parse(output / f"fields_{name}.pvd").findall("./Collection/DataSet")
        assert [float(item.get("timestep")) for item in datasets] == [0., 2., 4., 6., 8.]
        assert [item.get("file") for item in datasets] == [frame[name]["vtu"] for frame in frames]
        snapshots = []
        for frame in frames:
            native = io.load_model(output / frame[name]["dtcc"])
            reference = meshio.read(output / frame[name]["vtu"])
            assert type(native) is native_type
            assert native.transform.srs == "EPSG:7415"
            assert {key: getattr(native.bounds, key) for key in bounds} == bounds
            points = native.vertices if name == "tet" else native.coordinates()
            np.testing.assert_array_equal(reference.points, points)
            assert [field.name for field in native.fields] == ["velocity", "pressure"]
            assert set(reference.point_data) == {"velocity", "pressure"}
            for field in native.fields:
                assert field.association == "vertex" and field.values.dtype == np.float32
                assert field.values.shape == ((120, 3) if field.name == "velocity" else (120,))
                assert np.isfinite(field.values).all()
                np.testing.assert_array_equal(reference.point_data[field.name], field.values)
            assert len(reference.cells_dict[cell_type]) == (360 if name == "tet" else 60)
            snapshots.append({field.name: field.values for field in native.fields})
            if name == "tet":
                np.testing.assert_array_equal(reference.cells_dict["tetra"], native.cells)
            else:
                assert (native.width, native.height, native.depth) == (3, 4, 5)
        for field_name in ("velocity", "pressure"):
            np.testing.assert_allclose(snapshots[0][field_name], snapshots[-1][field_name], atol=1e-6)
            assert not np.allclose(snapshots[0][field_name], snapshots[1][field_name])
        first_fields[name] = snapshots[0]
    for field_name in ("velocity", "pressure"):
        np.testing.assert_array_equal(first_fields["tet"][field_name], first_fields["grid"][field_name])
    assert np.all(np.ptp(first_fields["tet"]["velocity"], axis=0) > 1)
    assert first_fields["tet"]["pressure"].min() < 0 < first_fields["tet"]["pressure"].max()

    mesh = io.load_model(output / frames[0]["tet"]["dtcc"])
    v = mesh.vertices[mesh.cells]
    volumes = np.linalg.det(v[:, 1:] - v[:, :1]) / 6
    assert np.all(volumes > 0)
    assert volumes.sum() == pytest.approx(2000*2000*300)
    # Every interior triangle is shared; unmatched faces must lie on the box.
    faces = np.concatenate([mesh.cells[:, face] for face in ((0, 1, 2), (0, 1, 3),
                                                            (0, 2, 3), (1, 2, 3))])
    unique, counts = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
    assert np.all((counts == 1) | (counts == 2))
    boundary = mesh.vertices[unique[counts == 1]]
    on_wall = (np.all(np.isclose(boundary, [bounds[key] for key in ("xmin", "ymin", "zmin")], rtol=0), axis=1)
               | np.all(np.isclose(boundary, [bounds[key] for key in ("xmax", "ymax", "zmax")], rtol=0), axis=1))
    assert np.all(np.any(on_wall, axis=1))


@pytest.mark.parametrize("arguments, message", [
    (["-nx", "0"], "positive integer"),
    (["-nt", "-1"], "positive integer"),
    (["-T", "nan"], "finite and positive"),
    (["-nx", "1000000"], "Protobuf limit"),
])
def test_invalid_cli_input_writes_nothing(tmp_path, monkeypatch, capsys, arguments, message):
    spec = importlib.util.spec_from_file_location("generate_fields", SCRIPT)
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), *arguments, "--output-dir", str(tmp_path / "absent")])
    with pytest.raises(SystemExit) as error:
        script.main()
    assert error.value.code == 2
    assert message in capsys.readouterr().err
    assert not (tmp_path / "absent").exists()


def test_missing_flagship_writes_nothing(tmp_path, monkeypatch, capsys):
    spec = importlib.util.spec_from_file_location("generate_fields", SCRIPT)
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    missing = tmp_path / "missing.dtcc"
    output = tmp_path / "absent"
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--flagship-model", str(missing),
                                      "--output-dir", str(output)])
    with pytest.raises(SystemExit) as error:
        script.main()
    assert error.value.code == 1
    # OSError reports the filename via repr(), which escapes Windows backslashes.
    assert repr(str(missing)) in capsys.readouterr().err
    assert not output.exists()
