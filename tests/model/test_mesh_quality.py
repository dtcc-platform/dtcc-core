"""Tests for mesh quality metrics (Mesh.quality() and VolumeMesh.quality())."""

import pytest
import numpy as np
from dtcc_core.model import Mesh, VolumeMesh
from dtcc_core.model.mixins.mesh.quality import (
    triangle_mesh_quality,
    tetrahedron_mesh_quality,
    tri_element_quality,
    tri_aspect_ratio,
    tri_edge_ratio,
    tri_skewness,
    tet_element_quality,
    tet_aspect_ratio,
    tet_edge_ratio,
    tet_skewness,
    format_quality,
    report_quality,
)


# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------

@pytest.fixture
def equilateral_tri_mesh():
    """Single equilateral triangle (side = 2)."""
    h = np.sqrt(3.0)
    vertices = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, h, 0.0]])
    faces = np.array([[0, 1, 2]])
    return Mesh(vertices=vertices, faces=faces)


@pytest.fixture
def right_tri_mesh():
    """Single right triangle (3-4-5)."""
    vertices = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
    faces = np.array([[0, 1, 2]])
    return Mesh(vertices=vertices, faces=faces)


@pytest.fixture
def two_tri_mesh():
    """Two triangles forming a unit square."""
    vertices = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=np.float64
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    return Mesh(vertices=vertices, faces=faces)


@pytest.fixture
def regular_tet_mesh():
    """Single regular tetrahedron."""
    vertices = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ]
    )
    cells = np.array([[0, 1, 2, 3]])
    return VolumeMesh(vertices=vertices, cells=cells)


@pytest.fixture
def flat_tet_mesh():
    """A degenerate (nearly flat) tet — poor quality."""
    vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [5.0, 10.0, 0.0],
            [5.0, 5.0, 0.01],
        ]
    )
    cells = np.array([[0, 1, 2, 3]])
    return VolumeMesh(vertices=vertices, cells=cells)


# -----------------------------------------------------------------------
# Triangle mesh quality tests
# -----------------------------------------------------------------------


class TestTriangleMeshQuality:
    def test_equilateral_element_quality(self, equilateral_tri_mesh):
        """Equilateral triangle should have element quality ≈ 1.0."""
        q = tri_element_quality(
            equilateral_tri_mesh.vertices, equilateral_tri_mesh.faces
        )
        assert pytest.approx(q[0], abs=1e-10) == 1.0

    def test_equilateral_aspect_ratio(self, equilateral_tri_mesh):
        """Equilateral triangle should have aspect ratio ≈ 1.0."""
        ar = tri_aspect_ratio(
            equilateral_tri_mesh.vertices, equilateral_tri_mesh.faces
        )
        assert pytest.approx(ar[0], abs=1e-10) == 1.0

    def test_equilateral_edge_ratio(self, equilateral_tri_mesh):
        """Equilateral triangle should have edge ratio ≈ 1.0."""
        er = tri_edge_ratio(
            equilateral_tri_mesh.vertices, equilateral_tri_mesh.faces
        )
        assert pytest.approx(er[0], abs=1e-10) == 1.0

    def test_equilateral_skewness(self, equilateral_tri_mesh):
        """Equilateral triangle should have skewness ≈ 0.0."""
        sk = tri_skewness(
            equilateral_tri_mesh.vertices, equilateral_tri_mesh.faces
        )
        assert pytest.approx(sk[0], abs=1e-10) == 0.0

    def test_right_triangle_quality_less_than_equilateral(self, right_tri_mesh):
        """Right triangle quality should be < 1.0."""
        q = tri_element_quality(right_tri_mesh.vertices, right_tri_mesh.faces)
        assert 0.0 < q[0] < 1.0

    def test_right_triangle_aspect_ratio_gt_one(self, right_tri_mesh):
        """Non-equilateral triangle should have aspect ratio > 1."""
        ar = tri_aspect_ratio(right_tri_mesh.vertices, right_tri_mesh.faces)
        assert ar[0] > 1.0

    def test_two_triangles_dict_shape(self, two_tri_mesh):
        """quality() should return correct dict structure for 2 triangles."""
        q = two_tri_mesh.quality()
        assert q["num_cells"] == 2
        for key in ("element_quality", "aspect_ratio", "edge_ratio", "skewness"):
            assert "min" in q[key]
            assert "max" in q[key]
            assert "mean" in q[key]

    def test_quality_values_ranges(self, two_tri_mesh):
        """All metrics should be within their valid ranges."""
        q = two_tri_mesh.quality()
        assert 0.0 < q["element_quality"]["min"] <= 1.0
        assert q["aspect_ratio"]["min"] >= 1.0
        assert q["edge_ratio"]["min"] >= 1.0
        assert 0.0 <= q["skewness"]["min"] <= 1.0


# -----------------------------------------------------------------------
# Tetrahedron mesh quality tests
# -----------------------------------------------------------------------


class TestTetrahedronMeshQuality:
    def test_regular_tet_element_quality(self, regular_tet_mesh):
        """Regular tetrahedron should have element quality ≈ 1.0."""
        q = tet_element_quality(regular_tet_mesh.vertices, regular_tet_mesh.cells)
        assert pytest.approx(q[0], abs=1e-6) == 1.0

    def test_regular_tet_aspect_ratio(self, regular_tet_mesh):
        """Regular tetrahedron should have aspect ratio ≈ 1.0."""
        ar = tet_aspect_ratio(regular_tet_mesh.vertices, regular_tet_mesh.cells)
        assert pytest.approx(ar[0], abs=1e-6) == 1.0

    def test_regular_tet_edge_ratio(self, regular_tet_mesh):
        """Regular tetrahedron should have edge ratio ≈ 1.0."""
        er = tet_edge_ratio(regular_tet_mesh.vertices, regular_tet_mesh.cells)
        assert pytest.approx(er[0], abs=1e-10) == 1.0

    def test_regular_tet_skewness(self, regular_tet_mesh):
        """Regular tetrahedron should have skewness ≈ 0.0."""
        sk = tet_skewness(regular_tet_mesh.vertices, regular_tet_mesh.cells)
        assert pytest.approx(sk[0], abs=1e-6) == 0.0

    def test_flat_tet_poor_quality(self, flat_tet_mesh):
        """Nearly-flat tet should have poor quality (< 0.1)."""
        q = tet_element_quality(flat_tet_mesh.vertices, flat_tet_mesh.cells)
        assert q[0] < 0.1

    def test_flat_tet_high_aspect_ratio(self, flat_tet_mesh):
        """Nearly-flat tet should have high aspect ratio."""
        ar = tet_aspect_ratio(flat_tet_mesh.vertices, flat_tet_mesh.cells)
        assert ar[0] > 5.0

    def test_quality_dict_shape(self, regular_tet_mesh):
        """quality() should return correct dict structure."""
        q = regular_tet_mesh.quality()
        assert q["num_cells"] == 1
        for key in ("element_quality", "aspect_ratio", "edge_ratio", "skewness"):
            assert "min" in q[key]
            assert "max" in q[key]
            assert "mean" in q[key]


# -----------------------------------------------------------------------
# Format helper test
# -----------------------------------------------------------------------


class TestFormatQuality:
    def test_format_quality(self, equilateral_tri_mesh):
        q = equilateral_tri_mesh.quality()
        s = format_quality(q)
        assert "Number of cells: 1" in s
        assert "Element quality" in s
        assert "Aspect ratio" in s
        assert "Edge ratio" in s
        assert "Skewness" in s

    def test_report_quality(self, equilateral_tri_mesh):
        q = equilateral_tri_mesh.quality()
        logged = []
        report_quality(q, title="Test mesh", log_fn=logged.append)
        assert any("Test mesh" in line for line in logged)
        assert any("Element quality" in line for line in logged)
        assert any("Number of cells: 1" in line for line in logged)
