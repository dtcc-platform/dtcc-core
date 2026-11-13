"""Geometry and placement of printable city tile packages."""

import json
from importlib import import_module

import numpy as np
import pytest
import shapely
from pyproj import Transformer
from shapely.geometry import box

from dtcc_core import datasets, io
from dtcc_core.builder.geometry_builders.printable_tiles import (
    _build_printable_tile,
    _validate_solid,
)
from dtcc_core.model import Mesh, Object
from tests.datasets.test_export_matrix import SITE, X0, Y0, _site


def test_layout_projects_the_selection_without_acquisition(no_network):
    bounds = (319720, 6397660, 320220, 6398160)
    preview = datasets.printable_city_tiles.layout(bounds=bounds)
    assert [feature["properties"]["id"] for feature in preview["features"]] == [
        "tile_0_0",
        "tile_1_0",
        "tile_0_1",
        "tile_1_1",
    ]
    assert preview["features"][-1]["properties"]["bounds"][2:] == list(bounds[2:])

    geographic = (11.98, 57.685, 11.983, 57.688)
    preview = datasets.printable_city_tiles.layout(
        bounds=geographic, crs="EPSG:4326", grid_origin=bounds[:2]
    )
    transform = Transformer.from_crs(4326, 3006, always_xy=True)
    drawn = box(*transform.transform_bounds(*geographic))
    tiles = [box(*feature["properties"]["bounds"]) for feature in preview["features"]]
    assert shapely.union_all(tiles).covers(drawn)
    lon, lat = preview["features"][0]["geometry"]["coordinates"][0][0]
    assert transform.transform(lon, lat) == pytest.approx(tiles[0].bounds[:2], abs=1e-7)


@pytest.mark.parametrize("lod", ["lod1", "lod2"])
def test_adjacent_print_stls_match_the_layout_and_reference_pockets(
    lod, monkeypatch, tmp_path, no_network
):
    for target, provider in _site(slope=0.02).items():
        monkeypatch.setattr(target, provider)

    result = datasets.printable_city_tiles(
        bounds=SITE, scale=100, tile_size_mm=200, lod=lod
    )
    package = result.export(tmp_path / "tiles", canonical=True)
    artifacts = package.path / "artifacts"
    layout = json.loads((artifacts / "layout.json").read_text())
    assert layout["base_tile"]["count"] == sum(
        t["base_stack_count"] for t in layout["tiles"]
    )
    assert len(layout["tiles"]) == 9
    assert (layout["width_mm"], layout["height_mm"]) == (600, 600)
    meshes = {}
    for tile in layout["tiles"]:
        mesh = io.load_mesh(artifacts / tile["stl"])
        _validate_solid(mesh)
        assert mesh.vertices.min(axis=0) == pytest.approx((0, 0, 0), abs=1e-6)
        assert mesh.vertices.max(axis=0)[:2] == pytest.approx((200, 200))
        meshes[tile["id"]] = mesh

    first = meshes["tile_0_0"]
    triangles = first.vertices[first.faces]
    footprint = box(0, 0, 200, 200)
    underside = shapely.union_all(
        shapely.polygons(triangles[(triangles[:, :, 2] == 0).all(axis=1), :, :2])
    )
    holes = shapely.get_parts(footprint.difference(underside))
    centers = sorted((round(p.centroid.x, 4), round(p.centroid.y, 4)) for p in holes)
    assert centers == [(8, 8), (8, 192), (192, 8), (192, 192)]
    assert all(p.area == pytest.approx(78.036, abs=1e-3) for p in holes)
    roof_triangles = triangles[(triangles[:, :, 2] == 2.5).all(axis=1), :, :2]
    assert shapely.union_all(shapely.polygons(roof_triangles)).area == pytest.approx(
        4 * 78.036, abs=0.005
    )
    preview = json.loads((artifacts / "preview.geojson").read_text())
    building_tiles = {
        f["properties"]["tile_id"]
        for f in preview["features"]
        if f["properties"]["kind"] == "building"
    }
    assert {"tile_0_0", "tile_1_0"} <= building_tiles
    seam_profiles = []
    for index, side in enumerate(("tile_0_0", "tile_1_0")):
        tile = next(t for t in layout["tiles"] if t["id"] == side)
        assert tile["bounds"][0] == X0 + tile["column"] * 20
        assert tile["bounds"][1] == Y0
        triangles = meshes[side].vertices[meshes[side].faces].astype(float)
        normals = np.cross(
            triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
        )
        upper = triangles[normals[:, 2] > 0].reshape(-1, 3)
        cut = upper[upper[:, 0] == 200 * (1 - index), 1:]
        seam_profiles.append(np.unique(cut + tile["print_origin_mm"][1:], axis=0))
    np.testing.assert_allclose(*seam_profiles, atol=3e-5, rtol=0)


def test_base_changes_pitch_and_retains_the_corner_geometry(monkeypatch, tmp_path):
    for target, provider in _site().items():
        monkeypatch.setattr(target, provider)
    result = datasets.printable_city_tiles(bounds=SITE, tile_size_mm=160, scale=375)
    package = result.export(tmp_path / "resized", canonical=True)
    base = io.load_mesh(package.path / "artifacts" / "base_tile.stl")
    _validate_solid(base)
    assert base.vertices.max(axis=0) == pytest.approx((160, 160, 20))
    ring = base.vertices[np.isclose(base.vertices[:, 2], 20.0)]
    for x in (8, 152):
        for y in (8, 152):
            radius = np.linalg.norm(ring[:, :2] - (x, y), axis=1)
            assert np.count_nonzero(np.isclose(radius, 5.0, atol=2e-5)) == 32


def test_export_stls_writes_the_print_files_directly(monkeypatch, tmp_path, no_network):
    for target, provider in _site().items():
        monkeypatch.setattr(target, provider)
    dataset = datasets.printable_city_tiles
    selection = dict(bounds=SITE, scale=300)
    folder = dataset.export_stls(tmp_path / "print", **selection)
    assert {path.name for path in folder.iterdir()} == {
        "tile_0_0.stl",
        "base_tile.stl",
        "layout.json",
    }
    for name in ("tile_0_0.stl", "base_tile.stl"):
        _validate_solid(io.load_mesh(folder / name))
    layout = json.loads((folder / "layout.json").read_text())
    assert layout["base_tile"]["count"] == 1

    tiles = dataset(**selection)
    assert type(tiles) is Object
    assert isinstance(tiles.get_geometry(id="base_tile"), Mesh)
    [tile] = tiles.get_children(Object)
    assert isinstance(tile.get_geometry(id="print_mesh"), Mesh)
    package = tiles.export(tmp_path / "assembly.dtccpkg", canonical=True)
    restored = datasets.load_model_package(package.path)
    with monkeypatch.context() as offline:
        offline.setattr(type(dataset), "build", None)
        direct = tiles.export_stl(folder=tmp_path / "from-model")
        reexport = restored.export_stl(folder=tmp_path / "from-package")
        for output in (direct, reexport):
            for path in folder.iterdir():
                assert (output / path.name).read_bytes() == path.read_bytes()
    with pytest.raises(ValueError, match="realized print layout"):
        Object().export_stl(folder=tmp_path / "empty")
    with pytest.raises(ValueError, match="either built tiles or Dataset arguments"):
        dataset.export_stls(tmp_path / "ambiguous", tiles=tiles, **selection)

    with pytest.raises(FileExistsError, match="new print output folder"):
        tiles.export_stl(folder=folder)
    assert json.loads((folder / "layout.json").read_text()) == layout


def test_invalid_layout_is_rejected_before_fetch(no_network):
    with pytest.raises(ValueError, match="64 tiles"):
        datasets.printable_city_tiles(bounds=SITE, scale=10)
    with pytest.raises(ValueError, match="finite"):
        datasets.printable_city_tiles.layout(bounds=(X0, Y0, float("inf"), Y0 + 1))
    with pytest.raises(ValueError, match="precise placement"):
        datasets.printable_city_tiles.layout(
            bounds=(319720, 6397660, 320220, 6398160), grid_origin=(1e18, 1e18)
        )


def test_fractional_tile_size_preserves_cut_planes():
    size = 140.12345
    pitch = size / 1000
    vertices = np.array(
        [[-1.0, -1.0, 5.0], [2.0, -1.0, 5.0], [2.0, 2.0, 5.0], [-1.0, 2.0, 5.0]]
    )
    triangles = vertices[[[0, 1, 2], [0, 2, 3]]]
    mesh, _, _ = _build_printable_tile(
        triangles, np.array([-2, -2]), (0, 0, pitch, pitch), size, 1, 0
    )
    assert np.ptp(mesh.vertices, axis=0)[:2] == pytest.approx((size, size))


def test_print_formats_describe_assembly_exports(tmp_path, no_network):
    dataset = datasets.printable_city_tiles
    described = dataset.describe()
    assert described["multi_file_formats"] == ["stl"]
    assert not described["serialization"]["bytes_when_format_is_set"]
    assert not described["serialization"]["format_parameter"]
    with pytest.raises(ValueError, match="export_stls"):
        dataset.export(tmp_path / "print.stl", bounds=SITE)


def test_print_export_matches_tile_ids(monkeypatch, tmp_path, no_network):
    for target, provider in _site().items():
        monkeypatch.setattr(target, provider)
    tiles = datasets.printable_city_tiles(bounds=SITE, scale=150)
    expected = tiles.export_stl(folder=tmp_path / "expected")
    tiles.children[Object].reverse()
    reordered = tiles.export_stl(folder=tmp_path / "reordered")
    assert {p.name: p.read_bytes() for p in expected.iterdir()} == {
        p.name: p.read_bytes() for p in reordered.iterdir()
    }

    layout = tiles.attributes["layout"]
    layout["tiles"][1] = layout["tiles"][0]
    with pytest.raises(ValueError, match="one unique mesh per tile ID"):
        tiles.export_stl(folder=tmp_path / "duplicate")
    assert not (tmp_path / "duplicate").exists()


def test_print_export_preserves_a_destination_created_during_staging(
    monkeypatch, tmp_path, no_network
):
    for target, provider in _site().items():
        monkeypatch.setattr(target, provider)
    tiles = datasets.printable_city_tiles(bounds=SITE, scale=300)
    destination = tmp_path / "print"
    module = import_module("dtcc_core.datasets.printable_city_tiles")
    write_print_files = module._write_print_files

    def stage(assembly, directory):
        write_print_files(assembly, directory)
        destination.mkdir()

    monkeypatch.setattr(module, "_write_print_files", stage)
    with pytest.raises(FileExistsError):
        tiles.export_stl(folder=destination)
    assert list(destination.iterdir()) == []
