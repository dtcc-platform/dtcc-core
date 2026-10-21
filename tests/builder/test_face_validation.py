"""Area allowance follows measured boundary coincidence, not an arbitrary epsilon."""

import json
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import box

from dtcc_core.builder.meshing.face_validation import validate_face_geometry
from dtcc_core.builder.meshing.surface_faces import SurfaceFace


def face(polygon, *, origin=None, basis=None, normal=None, source="wall-source"):
    return SurfaceFace(
        0,
        "wall",
        source,
        [0],
        np.empty(0, dtype=int),
        np.zeros(3) if origin is None else np.asarray(origin),
        np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        if basis is None
        else np.asarray(basis),
        np.array([0.0, -1.0, 0.0]) if normal is None else np.asarray(normal),
        polygon,
        [],
        [],
        0.0,
    )


def test_real_wall_coincidence_propagates_length_error_to_area():
    data = json.loads(
        (
            Path(__file__).parents[1] / "data/lod2/wall_boundary_coincidence.json"
        ).read_text()
    )
    domain = face(
        shapely.from_wkt(data["polygon"]),
        origin=data["origin"],
        basis=data["basis"],
        normal=data["normal"],
        source=data["source_key"],
    )
    result = validate_face_geometry(
        domain, np.array(data["vertices"]), np.array(data["triangles"])
    )
    assert 1e-7 < result["coverage_error_m2"] < 1.3e-7
    assert result["boundary_deviation_m"] < 1e-8
    assert result["coverage_error_m2"] < result["coverage_tolerance_m2"]
    assert result["overlap_m2"] == 0.0
    assert len(data["observations"]) == 881


@pytest.mark.parametrize(
    "corruption", ["missing_strip", "duplicate", "wrong_region", "missing_all"]
)
def test_real_wall_defects_are_not_covered_by_numerical_area_allowance(corruption):
    domain = face(box(0, 0, 1, 1e-7))
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1e-7], [0.0, 0.0, 1e-7]]
    )
    triangles = np.array([[0, 1, 2], [0, 2, 3]])
    validate_face_geometry(domain, vertices, triangles)
    # The missing strip is only 2e-8 m2: the previous fixed 1e-7 m2
    # area allowance could miss it, whereas the 10 nm boundary contract cannot.
    if corruption == "missing_strip":
        vertices[:2, 2] = 2e-8
    elif corruption == "duplicate":
        triangles = np.vstack([triangles, triangles])
    elif corruption == "wrong_region":
        vertices[:, 0] += 20.0
    else:
        triangles = np.empty((0, 3), dtype=int)
    with pytest.raises(ValueError):
        validate_face_geometry(domain, vertices, triangles)
