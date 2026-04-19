import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface
from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType


def _make_building_with_footprint(xmin, ymin, xmax, ymax):
    verts = np.array(
        [[xmin, ymin, 0], [xmax, ymin, 0], [xmax, ymax, 0], [xmin, ymax, 0]],
        dtype=float,
    )
    b = Building()
    b.add_geometry(Surface(vertices=verts), GeometryType.LOD0)
    return b


def test_extracted_points_have_classification():
    """After extraction, per-building PointCloud should carry classification."""
    b = _make_building_with_footprint(0, 0, 10, 10)
    pts = np.array(
        [
            [5, 5, 10],   # inside, class 6
            [5, 5, 0.5],  # inside, class 2 (ground -- filtered by non-ground mask)
            [50, 50, 10],  # outside
            [3, 3, 8],    # inside, class 6
        ],
        dtype=float,
    )
    classification = np.array([6, 2, 6, 6])
    pc = PointCloud(points=pts, classification=classification)

    from dtcc_core.builder.geometry_builders.buildings import extract_roof_points

    extract_roof_points([b], pc, statistical_outlier_remover=False)

    bpc = b.point_cloud
    assert bpc is not None
    assert bpc.classification is not None
    assert len(bpc.classification) == len(bpc.points)


def test_index_alignment_with_missing_footprints():
    """Buildings without footprints should not misalign roof point assignment."""
    b1 = Building()  # no footprint
    b2 = _make_building_with_footprint(0, 0, 10, 10)

    pts = np.array([[5, 5, 10], [3, 3, 8]], dtype=float)
    pc = PointCloud(points=pts)

    from dtcc_core.builder.geometry_builders.buildings import extract_roof_points

    extract_roof_points([b1, b2], pc, statistical_outlier_remover=False)

    # b1 should NOT have a point cloud (no footprint)
    assert b1.point_cloud is None
    # b2 should have the roof points
    assert b2.point_cloud is not None
