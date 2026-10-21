import numpy as np
import pytest

from dtcc_core.model import MultiSurface, Object, SemanticRegion, Solid, Surface
from dtcc_core.reproject.reproject import reproject_object


def test_reprojection_preserves_attachment_identity_and_descriptors():
    obj = Object()
    surface = Surface(
        vertices=np.array(
            [
                [500000.0, 6500000.0, 100.0],
                [500100.0, 6500000.0, 100.0],
                [500000.0, 6500100.0, 100.0],
            ]
        )
    )
    obj.add_geometry(surface, id="survey", lod="2.2", role="footprint")
    result = reproject_object(obj, "EPSG:3006", "EPSG:4326")
    assert list(result.geometry) == ["survey"]
    geometry = result.get_geometry(lod="2.2", role="footprint")
    assert geometry is result.geometry["survey"].geometry
    assert not np.allclose(geometry.vertices[:, :2], surface.vertices[:, :2])
    assert obj.get_geometry(id="survey") is surface


@pytest.mark.parametrize("shape", [MultiSurface, Solid])
def test_reprojection_rejects_unsupported_semantic_shapes(shape):
    geometry = shape(
        surfaces=[
            Surface(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                )
            )
        ]
    )
    geometry.regions = [
        SemanticRegion(
            semantic_type="https://github.com/dtcc-platform/dtcc-core/schemas/model#RoofSurface",
            indices=np.array([0]),
        )
    ]
    if shape is Solid:
        geometry.shells = [np.array([0])]
    obj = Object()
    obj.add_geometry(geometry, id="survey", lod="2.2")
    with pytest.raises(NotImplementedError, match="reprojection|Reprojection"):
        reproject_object(obj, "EPSG:3006", "EPSG:4326")
    assert obj.get_geometry(id="survey") is geometry
