import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg", force=True)

from dtcc_core.model import Field, Mesh, Point, exchange
from dtcc_core.plotting.style import require_matplotlib


@pytest.fixture
def mesh():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]]),
    )
    mesh.vertices[:, :2] += [319845.0, 6397785.0]
    mesh.transform.srs = "EPSG:3006"
    mesh.fields = [
        Field(
            name="temperature",
            unit="K",
            association="face",
            dim=1,
            values=np.array([[0.0], [0.2], [50.0], [50.0]]),
        )
    ]
    yield mesh
    require_matplotlib().close("all")


def plot(mesh, **kwargs):
    return mesh.plot(projection="xy", field="temperature", show=False, **kwargs)


def test_face_mapping_context_perimeter_render_and_immutability(mesh, tmp_path):
    before = exchange.dumps(mesh)
    ax = plot(
        mesh,
        context_mask=np.array([False, False, True, True]),
        coordinate_unit="m",
        scale_bar=True,
        title="Temperature at t = 2.5 s",
    )
    ax.figure.canvas.draw()
    assert ax.name == "rectilinear" and ax.get_aspect() == 1.0
    layers = {c.get_gid(): c for c in ax.collections}
    values = layers["face_field"]
    np.testing.assert_array_equal(values.get_array(), [0.0, 0.2])
    assert values.get_clim() == (0.0, 0.2)
    for path, face in zip(values.get_paths(), mesh.faces[:2], strict=True):
        np.testing.assert_array_equal(path.vertices[:3], mesh.vertices[face, :2])
    assert len(layers["face_context"].get_paths()) == 2
    segments = layers["context_outline"].get_segments()
    assert len(segments) == 4
    diagonal = mesh.vertices[[1, 5], :2]
    assert not any(np.array_equal(segment, diagonal) for segment in segments)
    assert ax.figure.axes[1].get_ylabel() == "temperature [K]"
    assert ax.get_xlabel() == "X [m]"
    assert ax.get_title() == "Temperature at t = 2.5 s"
    path = tmp_path / "map.png"
    ax.figure.savefig(path)
    from PIL import Image

    pixels = np.asarray(Image.open(path).convert("RGB"))
    assert len(np.unique(pixels.reshape(-1, 3), axis=0)) > 100
    assert exchange.dumps(mesh) == before


def test_native_roundtrip_transform_and_supplied_axes(mesh):
    mesh.transform.set_translation(10.0, 20.0, 3.0)
    restored = exchange.loads(exchange.dumps(mesh))
    plt = require_matplotlib()
    _, ax = plt.subplots()
    assert plot(restored, ax=ax, vmin=0.0, vmax=60.0) is ax
    np.testing.assert_array_equal(
        ax.collections[0].get_paths()[0].vertices[:3],
        mesh.vertices[mesh.faces[0], :2] + [10.0, 20.0],
    )
    assert ax.get_xlabel() == "X"
    assert any(t.get_text() == "EPSG:3006" for t in ax.texts)


def test_zero_empty_markers_and_show(mesh, monkeypatch):
    mesh.fields[0].values[:] = 0
    mesh.markers = np.array([], dtype=np.int32)
    shown = []
    monkeypatch.setattr(require_matplotlib(), "show", lambda: shown.append(True))
    ax = mesh.plot(projection="xy", field="temperature")
    assert shown == [True]
    assert any(t.get_text() == "All selected values are zero" for t in ax.texts)
    lo, hi = ax.collections[0].get_clim()
    assert lo < hi and np.isfinite([lo, hi]).all()
    assert plot(mesh, vmin=0.0, vmax=0.5).collections[0].get_clim() == (0.0, 0.5)


@pytest.mark.parametrize(
    "options,match",
    [
        ({"vmax": float("nan")}, "vmax"),
        ({"vmin": True}, "vmin"),
        ({"vmin": 1, "vmax": 0}, "vmin"),
        ({"vmin": 0, "vmax": 0}, "vmin"),
        ({"context_mask": [True]}, "context_mask"),
        ({"face_mask": [1, 1, 0, 0]}, "face_mask"),
        ({"context_mask": [True] * 4}, "at least one"),
        ({"context_mask": [True] * 4, "face_mask": [True] * 4}, "disjoint"),
        ({"scale_bar": True}, "coordinate_unit"),
        ({"coordinate_unit": ""}, "coordinate_unit"),
        ({"max_elements": 0}, "max_elements"),
    ],
)
def test_invalid_options(mesh, options, match):
    with pytest.raises(ValueError, match=match):
        plot(mesh, **options)


def test_missing_vector_wrong_association_and_count(mesh):
    with pytest.raises(ValueError, match="exactly one"):
        mesh.plot(projection="xy", field="missing")
    mesh.fields[0].association = "vertex"
    with pytest.raises(ValueError, match="association"):
        plot(mesh)
    mesh.fields[0].association = "face"
    mesh.fields[0].dim = 2
    with pytest.raises(ValueError, match="scalar"):
        plot(mesh)
    mesh.fields[0].dim = 1
    mesh.fields[0].values = np.array([[1.0]])
    with pytest.raises(ValueError, match="count"):
        plot(mesh)


def test_axes_projection_and_transform_errors(mesh):
    ax = require_matplotlib().figure().add_subplot(projection="3d")
    with pytest.raises(ValueError, match="2D axes"):
        plot(mesh, ax=ax)
    with pytest.raises(ValueError, match="projection"):
        mesh.plot(projection="invalid")
    with pytest.raises(ValueError, match="projection"):
        mesh.plot(vmin=0.0)
    with pytest.raises(ValueError, match="standalone"):
        Point().plot(projection="xy", field="temperature")
    mesh.transform.affine[3, 0] = 1
    with pytest.raises(ValueError, match="final row"):
        plot(mesh)


def test_selection_budget_and_affine_rotation(mesh):
    mesh.transform.affine[:2, :2] = [[0.0, -1.0], [1.0, 0.0]]
    ax = plot(mesh, max_elements=3, context_mask=np.array([False, False, True, True]))
    count = sum(len(c.get_paths()) for c in ax.collections[:2])
    count += len(ax.collections[2].get_segments())
    assert count <= 3
    assert any("sampled" in t.get_text() for t in ax.texts)
    expected = mesh.vertices[mesh.faces[0], :2] @ np.array([[0.0, 1.0], [-1.0, 0.0]])
    np.testing.assert_array_equal(
        ax.collections[0].get_paths()[0].vertices[:3], expected
    )
    selected = plot(mesh, face_mask=np.array([True, False, False, False]))
    assert len(selected.collections[0].get_paths()) == 1
    assert selected.collections[0].get_array().tolist() == [0.0]


def test_context_union_handles_nonconforming_edges(mesh):
    from dtcc_core.plotting.face_map import _perimeter

    xy = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [1.0, 0.0], [1.0, -1.0]])
    faces = np.array([[0, 1, 2], [0, 3, 4]])
    segments = _perimeter(xy, faces)
    assert len(segments) == 5
    assert not any(
        np.array_equal(np.sort(segment, axis=0), [[0.0, 0.0], [1.0, 0.0]])
        for segment in segments
    )
