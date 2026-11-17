import os
import pathlib
import subprocess
import sys
import tempfile

import meshio
import numpy as np
import pygltflib
import pytest

from dtcc_core import io
from dtcc_core.model import City, Field, Mesh, VolumeMesh


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).parent / ".." / "data" / "meshes"


@pytest.fixture
def mesh_paths(data_dir):
    return {
        "stl": str((data_dir / "cube.stl").resolve()),
        "vtk": str((data_dir / "cube.vtk").resolve()),
        "fbx": str((data_dir / "cube.fbx").resolve()),
        "quad_cube": str((data_dir / "quad_cube.obj").resolve()),
        "quad_tri_cube": str((data_dir / "quad_and_tri_cube.obj").resolve()),
        "nine_cubes": str((data_dir / "9_cubes.obj").resolve()),
    }


@pytest.fixture
def basic_mesh(mesh_paths):
    return io.load_mesh(mesh_paths["stl"])


@pytest.fixture
def volume_mesh_cube():
    vertices = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )
    cells = np.array(
        [
            [0, 1, 2, 4],
            [1, 2, 3, 5],
            [2, 3, 0, 6],
            [3, 0, 1, 7],
            [0, 4, 5, 1],
            [2, 6, 7, 3],
        ]
    )
    return VolumeMesh(vertices=vertices, cells=cells)


# Basic mesh loading tests
def test_load_mesh_stl(mesh_paths):
    mesh = io.load_mesh(mesh_paths["stl"])
    assert len(mesh.vertices) == 24
    assert len(mesh.faces) == 44


def test_load_mesh_vtk(mesh_paths):
    mesh = io.load_mesh(mesh_paths["vtk"])
    assert len(mesh.vertices) == 24
    assert len(mesh.faces) == 44


# @pytest.mark.skip(reason="This test segfaults")
def test_load_mesh_fbx(mesh_paths):
    if not io.meshes.has_assimp():
        pytest.skip("pyassimp not found, skipping test")
    mesh = io.load_mesh(mesh_paths["fbx"])
    assert len(mesh.vertices) == 24
    assert len(mesh.faces) == 44


def test_import_mesh_io_without_assimp_is_quiet():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys\n"
            "sys.modules['pyassimp'] = None\n"
            "from dtcc_core import io\n"
            "assert not io.meshes.has_assimp()\n",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert "assimp" not in (result.stdout + result.stderr).lower()


@pytest.mark.parametrize("extension", [".fbx", ".dae"])
def test_load_assimp_mesh_reports_missing_dependency(monkeypatch, tmp_path, extension):
    monkeypatch.setattr(io.meshes, "HAS_ASSIMP", False)
    with pytest.raises(
        RuntimeError, match="requires pyassimp and the native Assimp library"
    ):
        io.load_mesh(tmp_path / f"mesh{extension}")


# Save and load tests for different formats
@pytest.mark.parametrize("file_extension", [".stl", ".vtk"])
def test_save_load_mesh(basic_mesh, file_extension):
    with tempfile.NamedTemporaryFile(suffix=file_extension, delete=False) as tmp_file:
        path = tmp_file.name

    try:
        io.save_mesh(basic_mesh, path)
        loaded_mesh = io.load_mesh(path)
        assert len(loaded_mesh.vertices) == 24
        assert len(loaded_mesh.faces) == 44
    finally:
        os.unlink(path)


@pytest.mark.parametrize("file_extension", [".glb", ".gltf"])
def test_save_mesh_only(basic_mesh, file_extension):
    with tempfile.NamedTemporaryFile(suffix=file_extension, delete=False) as tmp_file:
        path = pathlib.Path(tmp_file.name)

    bin_path = path.with_suffix(".bin")

    try:
        if file_extension == ".gltf":
            with pytest.warns(UserWarning, match=r"saved to a \.bin file"):
                basic_mesh.save(path)
            assert bin_path.exists()
        else:
            basic_mesh.save(path)
            assert not bin_path.exists()

        assert os.path.exists(path)
    finally:
        if bin_path.exists():
            bin_path.unlink()
        if path.exists():
            path.unlink()


def _read_gltf_arrays(path):
    gltf = pygltflib.GLTF2().load(str(path))
    if path.suffix == ".glb":
        blob = gltf.binary_blob()
    else:
        blob = (path.parent / gltf.buffers[0].uri).read_bytes()
    index_accessor, position_accessor = gltf.accessors
    index_view = gltf.bufferViews[index_accessor.bufferView]
    position_view = gltf.bufferViews[position_accessor.bufferView]
    assert index_accessor.componentType == pygltflib.UNSIGNED_INT
    assert position_accessor.componentType == pygltflib.FLOAT
    assert index_view.byteLength == 4 * index_accessor.count
    assert position_view.byteLength == 12 * position_accessor.count
    indices = np.frombuffer(
        blob, np.uint32, index_accessor.count, index_view.byteOffset
    )
    positions = np.frombuffer(
        blob, np.float32, 3 * position_accessor.count, position_view.byteOffset
    ).reshape(-1, 3)
    assert np.allclose(positions.min(axis=0), position_accessor.min)
    assert np.allclose(positions.max(axis=0), position_accessor.max)
    return indices, positions, np.array(gltf.nodes[0].translation)


@pytest.mark.parametrize("file_extension", [".glb", ".gltf"])
def test_save_gltf_round_trip(basic_mesh, tmp_path, file_extension):
    mesh = Mesh(
        vertices=basic_mesh.vertices * 10.0 + [319370.0, 6397790.0, 12.0],
        faces=basic_mesh.faces,
    )
    path = tmp_path / f"mesh{file_extension}"
    if file_extension == ".gltf":
        with pytest.warns(UserWarning, match=r"saved to a \.bin file"):
            mesh.save(path)
    else:
        mesh.save(path)

    indices, positions, translation = _read_gltf_arrays(path)

    assert np.array_equal(indices.reshape(-1, 3), mesh.faces)
    y_up = positions.astype(np.float64) + translation
    z_up = np.column_stack([y_up[:, 0], -y_up[:, 2], y_up[:, 1]])
    assert np.allclose(z_up, mesh.vertices, rtol=0.0, atol=1e-3)


@pytest.mark.parametrize("bad_index", [-1, 24])
def test_save_gltf_rejects_out_of_range_index(basic_mesh, tmp_path, bad_index):
    faces = basic_mesh.faces.copy()
    faces[0, 0] = bad_index
    mesh = Mesh(vertices=basic_mesh.vertices, faces=faces)
    with pytest.raises(ValueError, match="face indices"):
        mesh.save(tmp_path / "mesh.glb")


def test_save_gltf_rejects_mesh_without_faces(basic_mesh, tmp_path):
    mesh = Mesh(vertices=basic_mesh.vertices, faces=np.zeros((0, 3), dtype=int))
    with pytest.raises(ValueError, match="no faces"):
        mesh.save(tmp_path / "mesh.glb")


# Specialized mesh tests
def test_load_quad_mesh(mesh_paths):
    mesh = io.load_mesh(mesh_paths["quad_cube"])
    assert len(mesh.vertices) == 8
    assert len(mesh.faces) == 12
    assert mesh.faces.shape[1] == 3


def test_load_quad_tri_mesh(mesh_paths):
    mesh = io.load_mesh(mesh_paths["quad_tri_cube"])
    assert len(mesh.vertices) == 8
    assert len(mesh.faces) == 12
    assert mesh.faces.shape[1] == 3


def test_load_disjoint_mesh(mesh_paths):
    mesh = io.load_mesh(mesh_paths["nine_cubes"])
    assert len(mesh.vertices) == 9 * 8
    assert len(mesh.faces) == 9 * 12


@pytest.mark.parametrize("merge_coplanar", [False, True])
def test_load_city_mesh(mesh_paths, merge_coplanar):
    city = io.load_mesh_as_city(
        mesh_paths["nine_cubes"], merge_coplanar_surfaces=merge_coplanar
    )
    assert isinstance(city, City)
    assert len(city.buildings) == 9

    if merge_coplanar:
        assert len(city.buildings[0].lod1.surfaces) == 6


# Volume mesh tests
def test_write_read_volume_mesh(volume_mesh_cube):
    with tempfile.NamedTemporaryFile(suffix=".vtk", delete=False) as tmp_file:
        path = tmp_file.name

    try:
        volume_mesh_cube.save(path)
        assert os.path.exists(path)
        mesh = io.load_volume_mesh(path)
        assert len(mesh.vertices) == 8
        assert len(mesh.cells) == 6
        assert mesh.cells.dtype == np.int64
    finally:
        os.unlink(path)


def test_save_volume_mesh_vtu_includes_point_fields(volume_mesh_cube, tmp_path):
    velocity = np.ones((len(volume_mesh_cube.vertices), 3))
    speed = np.ones((len(volume_mesh_cube.vertices), 1))
    volume_mesh_cube.add_field(
        Field(name="velocity", unit="m/s", values=velocity, dim=3)
    )
    volume_mesh_cube.add_field(Field(name="speed", unit="m/s", values=speed, dim=1))

    path = tmp_path / "volume_mesh.vtu"
    volume_mesh_cube.save(path)

    saved = meshio.read(path)
    assert "velocity" in saved.point_data
    assert "speed" in saved.point_data
    assert saved.point_data["velocity"].shape == (len(volume_mesh_cube.vertices), 3)
    assert saved.point_data["speed"].shape == (len(volume_mesh_cube.vertices),)


def test_save_xdmf_volume_mesh_writes_mesh_and_markers_to_single_file(
    volume_mesh_cube, tmp_path
):
    path = tmp_path / "volume_mesh.xdmf"
    volume_mesh_cube.boundary_faces = np.array([[0, 1, 2], [0, 1, 4]], dtype=int)
    volume_mesh_cube.boundary_markers = np.array([-1, -2], dtype=int)

    volume_mesh_cube.save(path)

    h5_path = tmp_path / "volume_mesh.h5"

    assert path.exists()
    assert h5_path.exists()

    main_xdmf = path.read_text()

    assert 'Grid Name="mesh"' in main_xdmf
    assert 'Grid Name="boundary_markers"' in main_xdmf
    assert 'NumberType="Float" Precision="8"' in main_xdmf
    assert 'NumberType="Int" Precision="8"' in main_xdmf

    loaded = io.load_volume_mesh(path)
    assert len(loaded.vertices) == len(volume_mesh_cube.vertices)
    assert len(loaded.cells) == len(volume_mesh_cube.cells)
    assert np.array_equal(loaded.boundary_faces, volume_mesh_cube.boundary_faces)
    assert np.array_equal(loaded.boundary_markers, volume_mesh_cube.boundary_markers)


def test_save_xdmf_volume_mesh_removes_stale_marker_sidecar(volume_mesh_cube, tmp_path):
    path = tmp_path / "volume_mesh.xdmf"
    marker_path = tmp_path / "volume_mesh_boundary_markers.xdmf"
    marker_path.write_text("stale")

    volume_mesh_cube.save(path)

    assert path.exists()
    assert not marker_path.exists()
