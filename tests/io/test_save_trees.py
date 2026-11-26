"""Tests for writing trees to vector files with ``save_trees``."""

import fiona
import numpy as np
import pyproj
import pytest
import shapely.geometry
import shapely.ops

from dtcc_core.io.trees import save_trees
from dtcc_core.model import Tree


def _tree(tree_id="73003757-40eb-4251-9892-43900d17e0c6", radius=2.0, srs=None):
    tree = Tree(
        position=np.array([319900.0, 6397900.0, 10.0]),
        height=8.0,
        crown_radius=radius,
    )
    tree.id = tree_id
    if srs is not None:
        tree.transform.srs = srs
    return tree


def _read(path):
    with fiona.open(path) as src:
        return [
            (shapely.geometry.shape(f["geometry"]), dict(f["properties"])) for f in src
        ]


def test_geojson_circle_radius_is_metres(tmp_path):
    """A WGS84 circle covers the crown radius in metres, not in degrees."""
    path = tmp_path / "circles.geojson"
    save_trees([_tree(radius=2.0)], str(path), as_circles=True)

    [(circle, _)] = _read(path)
    to_sweref = pyproj.Transformer.from_crs(
        "EPSG:4326", "EPSG:3006", always_xy=True
    ).transform
    circle_m = shapely.ops.transform(to_sweref, circle)

    minx, miny, maxx, maxy = circle_m.bounds
    assert maxx - minx == pytest.approx(4.0, abs=0.01)
    assert maxy - miny == pytest.approx(4.0, abs=0.01)
    assert circle_m.centroid.x == pytest.approx(319900.0, abs=0.01)
    assert circle_m.centroid.y == pytest.approx(6397900.0, abs=0.01)


def test_gpkg_circle_radius_is_metres(tmp_path):
    """A circle kept in the projected source CRS covers the crown radius."""
    path = tmp_path / "circles.gpkg"
    save_trees([_tree(radius=2.0)], str(path), as_circles=True)

    [(circle, _)] = _read(path)
    minx, miny, maxx, maxy = circle.bounds
    assert maxx - minx == pytest.approx(4.0, abs=1e-6)
    assert maxy - miny == pytest.approx(4.0, abs=1e-6)


def _wgs84_tree(radius=2.0):
    to_wgs84 = pyproj.Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)
    lon, lat = to_wgs84.transform(319900.0, 6397900.0)
    tree = _tree(radius=radius, srs="EPSG:4326")
    tree.position = np.array([lon, lat, 10.0])
    return tree


def test_circle_from_geographic_source_to_projected_output(tmp_path):
    """Trees in degrees get metre circles when the output CRS is projected."""
    path = tmp_path / "circles.gpkg"
    save_trees([_wgs84_tree(radius=2.0)], str(path), True, output_crs="EPSG:3006")

    [(circle, _)] = _read(path)
    minx, miny, maxx, maxy = circle.bounds
    assert maxx - minx == pytest.approx(4.0, abs=0.01)
    assert maxy - miny == pytest.approx(4.0, abs=0.01)
    assert circle.centroid.x == pytest.approx(319900.0, abs=0.01)
    assert circle.centroid.y == pytest.approx(6397900.0, abs=0.01)


@pytest.mark.parametrize("suffix", [".geojson", ".gpkg"])
def test_circle_with_no_projected_crs_is_rejected(tmp_path, suffix):
    """With degrees in and out, there is no metre CRS to buffer the crown in."""
    with pytest.raises(ValueError, match="projected"):
        save_trees([_wgs84_tree()], str(tmp_path / f"c{suffix}"), as_circles=True)


@pytest.mark.parametrize("suffix", [".geojson", ".gpkg"])
@pytest.mark.parametrize("as_circles", [False, True])
def test_string_ids_are_kept(tmp_path, suffix, as_circles):
    """String tree ids are written unchanged, so distinct ids stay distinct."""
    ids = ["73003757-40eb-4251-9892-43900d17e0c6", "tree-a", "tree-b"]
    path = tmp_path / f"trees{suffix}"
    save_trees([_tree(tree_id=i) for i in ids], str(path), as_circles=as_circles)

    assert [props["id"] for _, props in _read(path)] == ids
