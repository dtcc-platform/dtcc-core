"""Volume-generation repair forwarding, final reporting and optional-backend contracts."""

from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from dtcc_core.builder import build_city_volume_mesh
from dtcc_core.builder.geometry_builders import meshes as pipeline
from dtcc_core.builder.meshing import tetgen
from dtcc_core.model import City, Mesh

from .test_cleaning_meshing_integration import make_flat_city
from .test_volume_mesh_corners import tag, two_cells


@pytest.fixture
def backend(monkeypatch):
    source = tag(two_cells())

    def generate(**kwargs):
        output = kwargs["return_boundary_faces"]
        return SimpleNamespace(
            points=source.vertices.copy(),
            tets=source.cells.copy(),
            boundary_tri_faces=source.boundary_faces.copy() if output else None,
            boundary_tri_markers=source.boundary_markers.copy() if output else None,
        )

    call = Mock(side_effect=generate)
    monkeypatch.setattr(tetgen, "HAS_TETGEN", True)
    monkeypatch.setattr(
        tetgen, "tetwrap", SimpleNamespace(tetrahedralize=call), raising=False
    )
    monkeypatch.setattr(tetgen, "_tetgen_switch_module", None)
    return source, call


@pytest.mark.parametrize("enabled", [False, True])
@pytest.mark.parametrize("output", [False, True])
def test_low_level_generation_flag(backend, monkeypatch, enabled, output):
    from dtcc_core.builder.meshing import volume_corners

    source, call = backend
    surface = Mesh(
        vertices=source.vertices,
        faces=source.boundary_faces,
        markers=source.boundary_markers,
    )
    plc = SimpleNamespace(
        vertices=surface.vertices,
        shell_faces=surface.faces,
        boundary_facets=[[0, 1, 2]],
        boundary_facet_markers=[-2],
        diagnostics=SimpleNamespace(errors=[], warnings=[]),
    )
    monkeypatch.setattr(
        tetgen.tetgen_utils, "format_tetgen_plc_diagnostics", lambda _: ""
    )
    repair = Mock(wraps=volume_corners.repair_boundary_corners)
    monkeypatch.setattr(volume_corners, "repair_boundary_corners", repair)
    result = tetgen.build_volume_mesh(
        surface,
        prebuilt_plc=plc,
        return_boundary_faces=output,
        repair_boundary_corners=enabled,
        switches_params={"max_volume": 3.5},
    )
    assert result.num_cells == (8 if enabled else 2)
    assert result.find_boundary_corners().num_corners == (0 if enabled else 2)
    assert len(result.boundary_faces) == (6 if output else 0)
    assert result.boundary_corner_repair["detected_corners"] == 2
    assert result.boundary_corner_repair["detected_parent_cells"] == 2
    assert result.boundary_corner_repair["repaired_parent_cells"] == (
        2 if enabled else 0
    )
    assert repair.call_count == int(enabled)
    assert call.call_args.kwargs["switches_params"]["max_volume"] == 3.5
    np.testing.assert_array_equal(source.cells, two_cells().cells)


@pytest.mark.parametrize("enabled", [False, True])
def test_city_forwards_once_and_audits_final_mesh(backend, monkeypatch, enabled):
    from dtcc_core.builder.meshing import volume_corners

    monkeypatch.setattr(pipeline, "_require_tetgen", lambda: None)
    repair = Mock(wraps=volume_corners.repair_boundary_corners)
    monkeypatch.setattr(volume_corners, "repair_boundary_corners", repair)
    forwarded = Mock(wraps=tetgen.build_volume_mesh)
    monkeypatch.setattr(pipeline, "tetgen_build_volume_mesh", forwarded)
    monkeypatch.setattr(pipeline, "get_default_tetgen_switches", lambda: {})
    audit = {}
    mesh = build_city_volume_mesh(
        make_flat_city([]),
        max_mesh_size=20.0,
        domain_height=40.0,
        boundary_face_markers=False,
        repair_boundary_corners=enabled,
        stage_audit=audit,
        report_mesh_quality=True,
        max_volume=12.5,
    )
    assert forwarded.call_args.kwargs["repair_boundary_corners"] is enabled
    assert forwarded.call_args.kwargs["switches_params"]["max_volume"] == 12.5
    assert repair.call_count == int(enabled)
    attempt = audit["attempts"][-1]
    assert attempt["config"]["repair_boundary_corners"] is enabled
    assert (
        attempt["stages"]["volume_mesh"]["num_cells"]
        == mesh.num_cells
        == (8 if enabled else 2)
    )
    assert attempt["stages"]["boundary_corner_repair"] == mesh.boundary_corner_repair
    assert mesh.quality()["num_cells"] == mesh.num_cells
    assert attempt["result"]["status"] == "success"


@pytest.mark.parametrize("enabled", [False, True])
def test_isolated_city_volume_remains_explicitly_unsupported(enabled):
    with pytest.raises(ValueError, match="isolation"):
        build_city_volume_mesh(None, isolated=True, repair_boundary_corners=enabled)


@pytest.mark.parametrize("value", [1, None, "true"])
def test_generation_requires_boolean_before_backend_work(value):
    with pytest.raises(TypeError, match="repair_boundary_corners must be a boolean"):
        tetgen.build_volume_mesh(None, repair_boundary_corners=value)
    with pytest.raises(TypeError, match="repair_boundary_corners must be a boolean"):
        build_city_volume_mesh(None, repair_boundary_corners=value)


def test_city_model_forwards_flag(monkeypatch):
    from dtcc_core import builder

    call = Mock(return_value=two_cells())
    monkeypatch.setattr(builder, "build_city_volume_mesh", call)
    City().build_volume_mesh(repair_boundary_corners=True)
    assert call.call_args.kwargs["repair_boundary_corners"] is True


@pytest.mark.skipif(not tetgen.is_tetgen_available(), reason="TetGen is not installed")
@pytest.mark.parametrize("output", [False, True])
def test_actual_city_generation_removes_corners(output):
    mesh = build_city_volume_mesh(
        make_flat_city([]),
        max_mesh_size=20.0,
        domain_height=40.0,
        boundary_face_markers=output,
        repair_boundary_corners=True,
        report_mesh_quality=False,
    )
    assert mesh.num_cells > 0
    assert mesh.find_boundary_corners().num_corners == 0
    assert bool(len(mesh.boundary_faces)) == output


def test_suitable_generation_is_unchanged_with_flag(backend, monkeypatch):
    source, _ = backend
    suitable = source.repair_boundary_corners()
    source.vertices, source.cells = suitable.vertices, suitable.cells
    surface = Mesh(
        vertices=source.vertices,
        faces=source.boundary_faces,
        markers=source.boundary_markers,
    )
    plc = SimpleNamespace(
        vertices=surface.vertices,
        shell_faces=surface.faces,
        boundary_facets=[[0, 1, 2]],
        boundary_facet_markers=[-2],
        diagnostics=SimpleNamespace(errors=[], warnings=[]),
    )
    monkeypatch.setattr(
        tetgen.tetgen_utils, "format_tetgen_plc_diagnostics", lambda _: ""
    )
    result = tetgen.build_volume_mesh(
        surface, prebuilt_plc=plc, repair_boundary_corners=True
    )
    np.testing.assert_array_equal(result.vertices, source.vertices)
    np.testing.assert_array_equal(result.cells, source.cells)
    assert result.boundary_corner_repair["detected_corners"] == 0
    assert result.boundary_corner_repair["repaired_parent_cells"] == 0
