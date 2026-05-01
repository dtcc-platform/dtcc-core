import numpy as np

from dtcc_core.model import MultiSurface, Surface
from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from dtcc_core.builder.geometry_builders.lod2 import _fit_plane, _ransac_planes


def _surface(coords):
    return Surface(vertices=np.array(coords, dtype=float))


def _closed_box():
    return MultiSurface(
        surfaces=[
            _surface([[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]),
            _surface([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]),
            _surface([[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]),
            _surface([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]]),
            _surface([[1, 1, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1]]),
            _surface([[0, 1, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]]),
        ]
    )


def test_watertight_validator_accepts_closed_box():
    assert is_watertight(_closed_box())


def test_watertight_validator_rejects_open_shell():
    shell = _closed_box()
    shell.surfaces.pop()

    assert not is_watertight(shell)


def test_watertight_validator_rejects_non_manifold_edge():
    shell = _closed_box()
    duplicate_top = shell.surfaces[0].copy(geometry_only=True)
    shell.surfaces.append(duplicate_top)

    assert not is_watertight(shell)


def test_fit_plane_predicts_roof_z_values():
    points = np.array(
        [
            [0, 0, 10],
            [1, 0, 11],
            [0, 1, 12],
            [1, 1, 13],
            [2, 1, 14],
        ],
        dtype=float,
    )

    plane = _fit_plane(points)

    assert np.isclose(plane.z_at(2, 2), 16.0)


def test_ransac_planes_is_deterministic_for_two_planes():
    left = np.array([[x, y, 10 + 0.2 * x] for x in range(5) for y in range(5)], dtype=float)
    right = np.array([[x + 6, y, 12 - 0.2 * x] for x in range(5) for y in range(5)], dtype=float)
    points = np.vstack([left, right])

    first = _ransac_planes(points, seed=7)
    second = _ransac_planes(points, seed=7)

    assert [len(plane.inliers) for plane in first] == [len(plane.inliers) for plane in second]
    assert len(first) == 2
