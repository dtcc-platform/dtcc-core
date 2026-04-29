from __future__ import annotations

import numpy as np
from affine import Affine

from dtcc_core.builder.geometry_builders import terrain as terrain_module
from dtcc_core.model import Mesh, Raster


class _FakeCppMesh:
    def __init__(self, mesh: Mesh):
        self._mesh = mesh

    def from_cpp(self) -> Mesh:
        return self._mesh


def _make_raster() -> Raster:
    raster = Raster()
    raster.data = np.zeros((3, 3), dtype=np.float64)
    raster.georef = Affine.translation(0.0, 3.0) * Affine.scale(1.0, -1.0)
    raster.nodata = np.nan
    return raster


def test_build_terrain_surface_mesh_defaults_to_dtcc_mesher(monkeypatch):
    raster = _make_raster()
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [3.0, 3.0, 0.0],
                [0.0, 3.0, 0.0],
            ],
            dtype=np.float64,
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64),
        markers=np.array([-2, -2], dtype=np.int64),
    )
    terrain_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [3.0, 3.0, 0.0],
                [0.0, 3.0, 0.0],
            ],
            dtype=np.float64,
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64),
        markers=np.array([-2, -2], dtype=np.int64),
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(terrain_module, "resolve_2d_mesher", lambda mesher=None: "dtcc_mesher")
    monkeypatch.setattr(terrain_module, "raster_to_builder_gridfield", lambda raster: "grid")
    def fake_build_ground_mesh(**kwargs):
        captured["ground_request"] = kwargs
        return ground_mesh

    monkeypatch.setattr(
        terrain_module,
        "build_city_flat_mesh_from_coverage",
        fake_build_ground_mesh,
    )
    monkeypatch.setattr(terrain_module, "mesh_to_builder_mesh", lambda mesh: mesh)
    monkeypatch.setattr(
        terrain_module._dtcc_builder,
        "build_terrain_surface_mesh_from_ground_mesh",
        lambda builder_mesh, grid, smoothing: (
            captured.update(
                {
                    "builder_mesh": builder_mesh,
                    "grid": grid,
                    "smoothing": smoothing,
                }
            )
            or _FakeCppMesh(terrain_mesh)
        ),
    )
    monkeypatch.setattr(
        terrain_module._dtcc_builder,
        "build_terrain_surface_mesh",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError(
                "Explicit builder backend path should not run for the default terrain mesher."
            )
        ),
    )

    result = terrain_module.build_terrain_surface_mesh(
        raster,
        max_mesh_size=2.0,
        min_mesh_angle=20.0,
        smoothing=2,
        report_mesh_quality=False,
    )

    assert result is terrain_mesh
    assert captured["ground_request"]["backend"] == "dtcc_mesher"
    assert captured["ground_request"]["bounds"] == raster.bounds.tuple
    assert captured["builder_mesh"] is ground_mesh
    assert captured["grid"] == "grid"
    assert captured["smoothing"] == 2


def test_build_terrain_surface_mesh_forwards_explicit_mesher(monkeypatch):
    raster = _make_raster()
    captured: dict[str, object] = {}

    def fake_resolve(mesher=None):
        captured["requested_mesher"] = mesher
        return "triangle"

    monkeypatch.setattr(terrain_module, "resolve_2d_mesher", fake_resolve)
    monkeypatch.setattr(terrain_module, "raster_to_builder_gridfield", lambda raster: "grid")
    monkeypatch.setattr(terrain_module, "create_builder_polygon", lambda polygon: polygon)
    monkeypatch.setattr(
        terrain_module._dtcc_builder,
        "build_terrain_surface_mesh",
        lambda *args, **kwargs: Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0],
                    [3.0, 3.0, 0.0],
                    [0.0, 3.0, 0.0],
                ],
                dtype=np.float64,
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64),
            markers=np.array([-2, -2], dtype=np.int64),
        ),
    )
    monkeypatch.setattr(terrain_module, "builder_mesh_to_mesh", lambda mesh: mesh)

    terrain_module.build_terrain_surface_mesh(
        raster,
        mesher="triangle",
        report_mesh_quality=False,
    )

    assert captured["requested_mesher"] == "triangle"
