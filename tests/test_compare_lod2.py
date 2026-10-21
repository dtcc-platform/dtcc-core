"""Geometry and comparison contracts for the standalone Höganäs experiment."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("shapely", minversion="2.1")
from pyproj import CRS
from shapely.geometry import Polygon

from dtcc_core.model import Building, GeometryType, Surface

spec = importlib.util.spec_from_file_location(
    "compare_lod2", Path(__file__).parents[1] / "sandbox" / "compare_lod2.py"
)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_reference_vertical_wall_preserves_absolute_xyz():
    ring = np.array(
        [
            [91000.0, 6231000.0, 2.0],
            [91010.0, 6231000.0, 2.0],
            [91010.0, 6231000.0, 8.0],
            [91000.0, 6231000.0, 7.0],
        ]
    )
    vertices, faces, deviation = comparison.triangulate_reference(
        [np.vstack([ring, ring[0]])]
    )
    np.testing.assert_array_equal(vertices, ring)
    assert len(faces) == 2
    assert deviation == 0
    triangles = vertices[faces]
    assert (
        np.linalg.norm(
            np.cross(
                triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
            ),
            axis=1,
        ).sum()
        / 2
        == 55
    )


def test_reference_hole_not_filled_and_nonplanar_vertices_preserved():
    rings = [
        np.array(
            [[0.0, 0.0, 10.0], [4.0, 0.0, 11.0], [4.0, 4.0, 11.0], [0.0, 4.0, 10.0]]
        ),
        np.array(
            [[1.0, 1.0, 10.25], [1.0, 3.0, 10.25], [3.0, 3.0, 10.75], [3.0, 1.0, 10.75]]
        ),
    ]
    v, f, _ = comparison.triangulate_reference(rings)
    z = comparison.top_heights(v, f, np.array([[0.5, 0.5], [2.0, 2.0]]))
    assert z[0] == pytest.approx(10.125)
    assert np.isnan(z[1])
    rings[0][2, 2] += 0.1
    v, _, deviation = comparison.triangulate_reference(rings)
    np.testing.assert_array_equal(v, np.vstack(rings))
    assert deviation > 0


def test_top_envelope_ignores_winding_and_buried_base():
    v = np.array([[0.0, 0.0, 8.0], [2.0, 0.0, 10.0], [0.0, 2.0, 8.0]])
    v = np.vstack([v, v - [0, 0, 5]])
    z = comparison.top_heights(
        v, np.array([[2, 1, 0], [3, 4, 5]]), np.array([[0.5, 0.5], [3.0, 3.0]])
    )
    assert z[0] == pytest.approx(8.5)
    assert np.isnan(z[1])


def test_crs_rejects_tm_instead_of_local_projection():
    comparison.check_crs(CRS(3008), "local")
    with pytest.raises(ValueError, match="Unexpected CRS"):
        comparison.check_crs(CRS(3006), "wrong")
    with pytest.raises(ValueError, match="Missing CRS"):
        comparison.check_crs(None, "missing")


def test_metrics_use_identical_samples_and_report_missing_coverage(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(comparison, "OUTPUT_DIR", tmp_path)
    b = Building(id="lm-1")
    s = Surface()
    s.from_polygon(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]), 0)
    b.add_geometry(s, GeometryType.LOD0)
    v = np.array(
        [[0.0, 0.0, 10.0], [2.0, 0.0, 10.0], [2.0, 2.0, 10.0], [0.0, 2.0, 10.0]]
    )
    f = np.array([[0, 1, 2], [0, 2, 3]])
    result = comparison.compare(
        [b],
        (v, f),
        {"full": (v + [0, 0, 2], f), "partial": (v + [0, 0, 1], f[:1])},
        (0, 0, 2, 2),
    )
    assert result["lm_footprint_samples"] == 4
    assert result["common_samples"] == 3
    assert (
        result["variants"]["full"]["samples"]
        == result["variants"]["partial"]["samples"]
        == 3
    )
    assert result["variants"]["full"]["rmse_m"] == 2
    assert result["variants"]["partial"]["rmse_m"] == 1
    assert result["variants"]["partial"]["generated_samples"] == 3
    assert (tmp_path / "roof_errors.csv").exists()


def test_reference_retraced_spur_does_not_drop_roof_area():
    ring = np.array(
        [
            [0.0, 0.0, 10.0],
            [4.0, 0.0, 10.0],
            [4.0, 4.0, 10.0],
            [2.0, 4.0, 10.0],
            [2.0, 4.01, 10.0],
            [2.0, 4.0, 10.0],
            [0.0, 4.0, 10.0],
        ]
    )
    v, f, _ = comparison.triangulate_reference([ring])
    assert sum(Polygon(t[:, :2]).area for t in v[f]) == pytest.approx(16)
    assert np.all(v[:, 2] == 10)


def test_ground_alignment_keeps_vertical_bias_visible():
    from affine import Affine

    from dtcc_core.model import PointCloud, Raster

    cloud = PointCloud(
        points=np.array([[0.5, 1.5, 3.0], [1.5, 0.5, 3.0], [5.0, 5.0, 0.0]]),
        classification=np.array([2, 2, 2]),
    )
    terrain = Raster(data=np.full((2, 2), 3.1), georef=Affine(1, 0, 0, 0, -1, 2))
    result = comparison.ground_alignment(cloud, terrain)
    assert result["samples"] == 2
    assert result["median_raster_minus_cloud_m"] == pytest.approx(0.1)
    np.testing.assert_array_equal(cloud.points[:, 2], [3, 3, 0])


def test_roofer_cityjson_transform_semantics_and_building_parts(tmp_path):
    import json

    header = dict(
        type="CityJSON",
        transform=dict(
            scale=[0.001, 0.001, 0.001], translate=[91000.0, 6231000.0, 2.0]
        ),
        vertices=[],
        CityObjects={},
    )
    geometry = dict(
        type="Solid",
        lod="2.2",
        boundaries=[[[[0, 1, 2, 3]], [[4, 5, 6, 7]]]],
        semantics=dict(
            surfaces=[dict(type="GroundSurface"), dict(type="RoofSurface")],
            values=[[0, 1]],
        ),
    )
    feature = dict(
        type="CityJSONFeature",
        vertices=[
            [0, 0, 0],
            [1000, 0, 0],
            [1000, 1000, 0],
            [0, 1000, 0],
            [0, 0, 5000],
            [1000, 0, 5000],
            [1000, 1000, 5000],
            [0, 1000, 5000],
        ],
        CityObjects={
            "b": dict(
                type="Building",
                attributes=dict(rf_success=True, rf_extrusion_mode="standard"),
                children=["b-0"],
            ),
            "b-0": dict(type="BuildingPart", parents=["b"], geometry=[geometry]),
        },
    )
    (tmp_path / "model.city.jsonl").write_text(
        json.dumps(header) + "\n" + json.dumps(feature) + "\n"
    )
    (v, f, roof), attrs = comparison.read_roofer(tmp_path, "b")
    assert len(f) == 2 and roof.all()
    np.testing.assert_allclose(v.min(0), [91000, 6231000, 7])
    np.testing.assert_allclose(v.max(0), [91001, 6231001, 7])
    assert attrs["rf_extrusion_mode"] == "standard"
    with pytest.raises(ValueError, match="Unexpected Roofer building ID"):
        comparison.read_roofer(tmp_path, "other")


def test_roofer_no_geometry_is_missing_not_flat_roof(tmp_path):
    import json

    feature = dict(
        CityObjects={
            "b": dict(type="Building", attributes=dict(rf_extrusion_mode="skip"))
        },
        vertices=[],
    )
    (tmp_path / "empty.jsonl").write_text(json.dumps(feature))
    mesh, attrs = comparison.read_roofer(tmp_path, "b")
    assert mesh is None and attrs["rf_extrusion_mode"] == "skip"
    z = comparison.top_heights(
        np.empty((0, 3)), np.empty((0, 3), dtype=int), np.zeros((2, 2))
    )
    assert np.isnan(z).all()


def test_roofer_unavailable_fails_before_output_changes(tmp_path, monkeypatch):
    monkeypatch.setattr(comparison, "METHODS", ("roofer",))
    monkeypatch.setattr(
        comparison, "ROOFER_EXECUTABLE", str(tmp_path / "no-such-roofer")
    )
    monkeypatch.setattr(comparison, "OUTPUT_DIR", tmp_path / "output")
    with pytest.raises(SystemExit, match="ROOFER_EXECUTABLE") as exc:
        comparison.main()
    message = str(exc.value)
    assert "ROOFER_VERSION=1.0.0 sh" in message
    assert 'export PATH="$HOME/.local/bin:$PATH"' in message
    assert "roofer --version" in message
    assert not (tmp_path / "output").exists()


def test_roofer_process_failure_is_audited_without_dtcc_fallback(tmp_path, monkeypatch):
    import json
    import subprocess

    from dtcc_core.model import City, PointCloud

    monkeypatch.setattr(comparison, "OUTPUT_DIR", tmp_path)
    b = Building(id="b")
    s = Surface()
    s.from_polygon(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]), 8.0)
    b.add_geometry(s, GeometryType.LOD0)
    b.add_geometry(
        PointCloud(
            points=np.array([[0.5, 0.5, 8.0], [1.0, 1.0, 8.0], [1.5, 1.5, 8.0]])
        ),
        GeometryType.POINT_CLOUD,
    )
    b.attributes["ground_height"] = 2.0
    city = City()
    city.add_buildings([b], remove_outside_terrain=False)

    def fail(command, **kwargs):
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(comparison.subprocess, "run", fail)
    v, f, summary = comparison.reconstruct_roofer(
        city, "roofer_test", dict(executable="/fake/roofer", version="test")
    )
    assert v.shape == f.shape == (0, 3)
    assert summary["outcomes"] == {"failed": 1}
    audit = json.loads((tmp_path / "roofer_test_audit.json").read_text())
    assert audit["buildings"][0]["source_id"] == "b"
    assert (tmp_path / "roofer_test.vtu").exists()


def test_default_run_dispatches_all_six_variants(tmp_path, monkeypatch):
    """Keep Roofer in the real script orchestration without live data/CLI needs."""
    from dtcc_core.model import City, PointCloud

    assert comparison.METHODS == ("catalogue", "ransac", "roofer")
    monkeypatch.setattr(comparison, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(comparison, "BOUNDS", (0.0, 0.0, 2.0, 2.0))
    monkeypatch.setattr(comparison, "MIN_BUILDING_AREA", 1.0)
    building = Building(id="b")
    surface = Surface()
    surface.from_polygon(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]), 0.0)
    building.add_geometry(surface, GeometryType.LOD0)
    cloud = PointCloud(points=np.array([[1.0, 1.0, 8.0]]), classification=np.array([6]))
    building.add_geometry(cloud, GeometryType.POINT_CLOUD)
    city = City()
    city.add_buildings([building], remove_outside_terrain=False)
    v = np.array([[0.0, 0.0, 8.0], [2.0, 0.0, 8.0], [0.0, 2.0, 8.0]])
    f = np.array([[0, 1, 2]])
    calls = []
    tool = dict(executable="/fake/roofer", version="test")
    monkeypatch.setattr(comparison, "resolve_roofer", lambda: tool)
    monkeypatch.setattr(comparison, "load_terrain", lambda bounds: None)
    monkeypatch.setattr(
        comparison,
        "load_reference",
        lambda bounds: (v, f, dict(buildings=1, omitted_faces=[])),
    )
    monkeypatch.setattr(
        comparison.dtcc.io.data, "download_footprints", lambda bounds: [building]
    )
    monkeypatch.setattr(
        comparison.dtcc.io.data, "download_pointcloud", lambda bounds: cloud
    )
    monkeypatch.setattr(comparison, "reproject_surface", lambda surface, *args: surface)
    monkeypatch.setattr(comparison, "reproject_array", lambda points, *args: points)
    monkeypatch.setattr(comparison, "local_pointcloud", lambda bounds: (cloud, {}))
    monkeypatch.setattr(comparison, "ground_alignment", lambda *args: {})
    monkeypatch.setattr(comparison, "prepare_city", lambda *args: city)

    def reconstruct(city, method, name):
        calls.append((method, name))
        return v, f, dict(triangles=1)

    def roofer(city, name, selected_tool):
        assert selected_tool == tool
        return reconstruct(city, "roofer", name)

    monkeypatch.setattr(comparison, "reconstruct", reconstruct)
    monkeypatch.setattr(comparison, "reconstruct_roofer", roofer)
    comparison.main()
    assert calls == [
        (method, f"{method}_{cloud}")
        for cloud in ("lm", "local")
        for method in ("catalogue", "ransac", "roofer")
    ]
    import json

    manifest = json.loads((tmp_path / "comparison.json").read_text())
    assert manifest["status"] == "complete"
    assert set(manifest["variants"]) == {name for _, name in calls}


def test_roofer_success_preserves_coordinates_and_open_skin(tmp_path, monkeypatch):
    import json
    import subprocess

    import meshio

    from dtcc_core.model import City, PointCloud

    monkeypatch.setattr(comparison, "OUTPUT_DIR", tmp_path)
    b = Building(id="b")
    surface = Surface()
    surface.from_polygon(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]), 8.0)
    b.add_geometry(surface, GeometryType.LOD0)
    b.add_geometry(
        PointCloud(points=np.array([[0.5, 0.5, 8.0], [1.0, 1.0, 8.0]])),
        GeometryType.POINT_CLOUD,
    )
    b.attributes["ground_height"] = 2.0
    city = City()
    city.add_buildings([b], remove_outside_terrain=False)
    feature = dict(
        vertices=[[0, 0, 8], [2, 0, 8], [2, 2, 8], [0, 2, 8], [0, 0, 2], [2, 0, 2]],
        CityObjects={
            "b": dict(
                type="Building",
                attributes=dict(rf_extrusion_mode="standard"),
                geometry=[
                    dict(
                        type="MultiSurface",
                        lod="2.2",
                        boundaries=[[[0, 1, 2, 3]], [[4, 5, 1, 0]]],
                        semantics=dict(
                            surfaces=[
                                dict(type="RoofSurface"),
                                dict(type="WallSurface"),
                            ],
                            values=[0, 1],
                        ),
                    )
                ],
            )
        },
    )

    def run(command, **kwargs):
        assert "--no-simplify" in command and "--no-clip-terrain" in command
        assert (
            kwargs["check"] and kwargs["timeout"] == comparison.ROOFER_TIMEOUT_SECONDS
        )
        destination = Path(command[-1])
        destination.mkdir()
        (destination / "model.jsonl").write_text(json.dumps(feature))
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(comparison.subprocess, "run", run)
    v, roof, summary = comparison.reconstruct_roofer(
        city, "roofer_test", dict(executable="/fake/roofer")
    )
    assert summary["outcomes"] == {"reconstructed": 1}
    assert len(roof) == 2 and np.all(v[roof, 2] == 8.0)
    saved = meshio.read(tmp_path / "roofer_test.vtu")
    np.testing.assert_array_equal(saved.points, v)
    assert len(saved.cells_dict["triangle"]) == 4  # Roof and wall, with no base.
    np.testing.assert_array_equal(
        saved.cell_data_dict["is_roof"]["triangle"], [1, 1, 0, 0]
    )
