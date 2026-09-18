"""Region identity and membership survive the public polygon-to-mesh workflow."""

from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.builder.meshing import mesh_multisurfaces
from dtcc_core.builder.meshing.backends import available_2d_meshers
from dtcc_core.model import Field, Solid, Surface, SemanticRegion, exchange


@pytest.fixture
def geometry():
    path = (
        Path(__file__).resolve().parents[2]
        / "sandbox/model_profiles/fixtures/buildings.city.json"
    )
    value = io.load_city(path, strict=True).buildings[0].building_parts[0].lod2
    value.regions[1].id = "shared-roof"
    value.regions[1].attributes["details"] = {"source": ["survey"]}
    value.transform.affine[0, 3] = 12.0
    return value


def areas(mesh):
    triangles = mesh.vertices[mesh.faces]
    return (
        np.linalg.norm(
            np.cross(
                triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
            ),
            axis=1,
        )
        / 2
    )


@pytest.mark.parametrize("mesher", available_2d_meshers())
def test_mesh_preserves_regions_holes_and_unclassified_faces(
    geometry, mesher, tmp_path
):
    before = exchange.dumps(geometry)
    mesh = geometry.mesh(mesher=mesher)
    assert exchange.dumps(geometry) == before
    np.testing.assert_array_equal(mesh.transform.affine, geometry.transform.affine)
    assert mesh.transform.srs == geometry.transform.srs
    assert [r.semantic_type for r in mesh.regions] == [
        r.semantic_type for r in geometry.regions
    ]
    assert mesh.regions[1].id == "shared-roof"
    assert mesh.regions[1].attributes == geometry.regions[1].attributes
    face_areas = areas(mesh)
    # The courtyard hole is excluded; the two roof triangles share one region.
    assert face_areas[mesh.regions[0].indices].sum() == pytest.approx(96)
    assert face_areas[mesh.regions[1].indices].sum() == pytest.approx(100)
    assert len(mesh.regions[1].indices) == 2
    classified = np.concatenate([r.indices for r in mesh.regions])
    assert len(np.unique(classified)) == len(classified)
    assert len(classified) < len(mesh.faces)
    mesh.save(tmp_path / "roof.dtcc")
    assert exchange.dumps(io.load_model(tmp_path / "roof.dtcc")) == exchange.dumps(mesh)
    # Batch entry points have the same behavior; metadata must not alias input.
    batch = mesh_multisurfaces([geometry], mesher=mesher)[0]
    assert exchange.dumps(batch) == exchange.dumps(mesh)
    mesh.regions[1].attributes["details"]["source"].append("derived")
    assert exchange.dumps(geometry) == before


def test_meshing_reports_unmapped_operations_before_mutating_input(geometry):
    before = exchange.dumps(geometry)
    for options in ({"clean": True}, {"weld": True}, {"snap": 0.1}):
        with pytest.raises(NotImplementedError, match="face mapping"):
            geometry.mesh(**options)
    with pytest.raises(NotImplementedError, match="face mapping"):
        mesh_multisurfaces([geometry], clean=True)
    assert exchange.dumps(geometry) == before
    geometry.surfaces[0].transform.affine[0, 3] = 1
    with pytest.raises(NotImplementedError, match="geometry frame"):
        geometry.mesh()
    geometry.surfaces[0].transform.affine[0, 3] = 0
    geometry.fields = [
        Field(
            name="temperature",
            values=np.zeros(len(geometry.surfaces)),
            association="face",
        )
    ]
    with pytest.raises(NotImplementedError, match="fields"):
        geometry.mesh()


def test_invalid_region_membership_fails_before_meshing(geometry):
    geometry.regions[1].indices[0] = len(geometry.surfaces)
    with pytest.raises(ValueError, match="out of range"):
        geometry.mesh()


def test_zero_area_source_ring_fails_before_meshing(geometry):
    geometry.surfaces[1].vertices[2] = geometry.surfaces[1].vertices[0]
    with pytest.raises(ValueError, match="nonzero-area"):
        geometry.mesh()


def test_surface_regions_cannot_be_silently_cleaned_away(geometry):
    surface = geometry.surfaces[0]
    surface.regions = geometry.regions
    with pytest.raises(NotImplementedError, match="MultiSurface"):
        surface.mesh(clean=True)
    with pytest.raises(NotImplementedError, match="MultiSurface"):
        mesh_multisurfaces([geometry], clean=True)


@pytest.fixture
def solid():
    vertices = np.array(
        [
            [0, 0, 0],
            [4, 0, 0],
            [4, 4, 0],
            [0, 4, 0],
            [0, 0, 4],
            [4, 0, 4],
            [4, 4, 4],
            [0, 4, 4],
        ],
        dtype=float,
    )
    faces = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    value = Solid(
        surfaces=[Surface(vertices=vertices[f]) for f in faces]
        + [Surface(vertices=(1 + vertices[f[::-1]] / 2)) for f in faces],
        shells=[np.arange(6), np.arange(6, 12)],
        regions=[
            SemanticRegion(
                "https://github.com/dtcc-platform/dtcc-core/schemas/model#WallSurface",
                indices=np.array([2, 3, 8, 9]),
                id="walls",
                attributes={"source": ["survey"]},
            )
        ],
    )
    value.transform.set_translation(10, 20, 30)
    return value


@pytest.mark.parametrize("mesher", available_2d_meshers())
def test_solid_mesh_includes_cavity_and_preserves_source(solid, mesher, tmp_path):
    before = exchange.dumps(solid)
    mesh = solid.mesh(mesher=mesher)
    assert exchange.dumps(solid) == before
    np.testing.assert_array_equal(mesh.transform.affine, solid.transform.affine)
    assert areas(mesh).sum() == pytest.approx(120)  # 96 exterior + 24 cavity.
    assert areas(mesh)[mesh.regions[0].indices].sum() == pytest.approx(40)
    assert mesh.regions[0].id == "walls"
    assert mesh.regions[0].attributes == solid.regions[0].attributes
    mesh.fields = [
        Field(name="triangle_area", unit="m2", association="face", values=areas(mesh))
    ]
    mesh.save(tmp_path / "boundary.dtcc")
    assert exchange.dumps(io.load_model(tmp_path / "boundary.dtcc")) == exchange.dumps(
        mesh
    )
    mesh.regions[0].attributes["source"].append("derived")
    mesh.vertices[0] += 1
    assert exchange.dumps(solid) == before


def test_solid_mesh_rejects_unmapped_data_and_invalid_rings_without_regions(solid):
    solid.regions = []
    before = exchange.dumps(solid)
    mesh = solid.mesh()
    assert areas(mesh).sum() == pytest.approx(120)
    np.testing.assert_array_equal(mesh.transform.affine, solid.transform.affine)
    for options in ({"clean": True}, {"weld": True}, {"snap": 0.1}):
        with pytest.raises(NotImplementedError, match="face mapping"):
            solid.mesh(**options)
    solid.fields = [Field(name="temperature", association="face", values=np.zeros(12))]
    with pytest.raises(NotImplementedError, match="fields"):
        solid.mesh()
    solid.fields = []
    assert exchange.dumps(solid) == before
    solid.surfaces[0].vertices[:] = 0
    with pytest.raises(ValueError, match="nonzero-area"):
        solid.mesh()


@pytest.mark.skipif(
    "dtcc_mesher" not in available_2d_meshers(), reason="dtcc_mesher not installed"
)
def test_nearly_collinear_leading_vertices_do_not_tilt_mesh_projection():
    # Millimetre quantization near a long straight edge must not turn an almost
    # horizontal roof into a steep plane. This reproduces the 3DBAG area loss.
    vertices = np.array(
        [
            [0, 0, 0],
            [5, 0, 0.001],
            [6, 0.001, 0.001],
            [10, 0, 0],
            [10, 10, 0],
            [0, 10, 0],
        ],
        dtype=float,
    )
    surface = Surface(vertices=vertices + [84000, 446000, 15])
    before = exchange.dumps(surface)
    mesh = surface.copy().mesh(mesher="dtcc_mesher")
    assert areas(mesh).sum() == pytest.approx(99.9975, abs=1e-5)
    assert np.ptp(mesh.vertices[:, 2]) < 0.002
    assert exchange.dumps(surface) == before
    forward = surface.calculate_normal().copy()
    surface.vertices = surface.vertices[::-1].copy()
    np.testing.assert_allclose(surface.calculate_normal(), -forward, atol=1e-12)
    surface.vertices[:] = [84000, 446000, 15]
    with pytest.raises(ValueError, match="nonzero finite area"):
        surface.calculate_normal()
