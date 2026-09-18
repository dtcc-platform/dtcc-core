"""Strict TIN conversion preserves mixed data and rejects unrepresentable facts."""

import json

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.model import Terrain, Mesh, Field, SemanticRegion, exchange
from sandbox.model_profiles.tin_relief_example import source, verify


def west(city):
    return next(t for t in city.get_children(Terrain) if t.id == "terrain-west")


def test_public_mixed_terrain_roundtrip_and_direct_access(tmp_path):
    data = source()
    path = tmp_path / "source.city.json"
    path.write_text(json.dumps(data))
    city = io.load_city(path, strict=True)
    terrain = west(city)
    mesh = terrain.get_geometry(lod="1.2")
    assert type(terrain) is Terrain and type(mesh) is Mesh
    assert (
        terrain.semantic_type
        == "https://github.com/dtcc-platform/dtcc-core/schemas/model#TINRelief"
    )
    assert terrain.get_geometry(id="cityjson-0") is mesh
    assert mesh.vertices.shape == (4, 3)
    np.testing.assert_array_equal(mesh.faces, [[0, 1, 2], [0, 2, 3]])
    assert list(terrain.geometry) == ["cityjson-0", "cityjson-1"]
    # Every source triangle, in order, after the document's exact dequantization.
    coordinates = (
        np.asarray(data["vertices"]) * data["transform"]["scale"]
        + data["transform"]["translate"]
    )
    for t in city.get_children(Terrain):
        original = data["CityObjects"][t.id]
        assert t.attributes == original["attributes"]
        for record, g in zip(t.geometry.values(), original["geometry"]):
            np.testing.assert_array_equal(
                record.geometry.vertices[record.geometry.faces],
                coordinates[np.array(g["boundaries"])[:, 0]],
            )
            assert record.lod == g["lod"]
    payload = exchange.dumps(city)
    city.save(tmp_path / "mixed.dtcc")
    assert exchange.dumps(io.load_city(tmp_path / "mixed.dtcc")) == payload
    city.dataset_context = DatasetContext(
        identity={"name": "test-tin", "title": "Synthetic TIN"},
        metadata={},
        provenance={},
        presentation={},
        request={"dataset_name": "test-tin"},
    )
    package = city.export(tmp_path / "mixed.dtccpkg", canonical=True)
    assert exchange.dumps(load_model_package(package.path)) == payload
    city.save(tmp_path / "mixed.city.json", strict=True)
    verify(city, io.load_city(tmp_path / "mixed.city.json", strict=True))
    exported = json.loads((tmp_path / "mixed.city.json").read_text())
    assert set(exported["CityObjects"]) == set(data["CityObjects"])
    for t in city.get_children(Terrain):
        assert all(
            g["type"] == "CompositeSurface" and "semantics" not in g
            for g in exported["CityObjects"][t.id]["geometry"]
        )


def test_malformed_or_unsupported_tin_is_never_partially_imported():
    # Corrupt the second terrain: success must not return only the first feature.
    for boundaries in (
        [],
        [[[]]],
        [[[0, 1, 2, 3]]],
        [[[0, 1, 2], [3, 4, 5]]],
        [[[0, 0, 1]]],
        [[[True, 1, 2]]],
        [[[-1, 1, 2]]],
        [[[0, 1, 99999]]],
    ):
        data = source()
        data["CityObjects"]["terrain-east"]["geometry"][0]["boundaries"] = boundaries
        with pytest.raises(ValueError, match="TINRelief"):
            io.load_cityjson(data, strict=True, validate_schema=False)
    for key, value in [
        ("semantics", {"surfaces": [], "values": [None]}),
        ("type", "MultiSurface"),
    ]:
        data = source()
        data["CityObjects"]["terrain-east"]["geometry"][0][key] = value
        with pytest.raises(NotImplementedError, match="TINRelief"):
            io.load_cityjson(data, strict=True)


def test_unmapped_mesh_state_and_generic_terrain_do_not_overwrite(tmp_path):
    city = io.load_cityjson(source(), strict=True)
    path = tmp_path / "terrain.city.json"
    city.save(path, strict=True)
    before = path.read_bytes()
    for field, value, message in [
        ("markers", np.array([1, 2]), "markers"),
        ("normals", np.zeros((2, 3)), "normals"),
        (
            "fields",
            [Field(name="z", unit="m", association="vertex", values=np.arange(4.0))],
            "fields",
        ),
        (
            "regions",
            [
                SemanticRegion(
                    "https://github.com/dtcc-platform/dtcc-core/schemas/model#GroundSurface",
                    np.array([0]),
                )
            ],
            "semantic regions",
        ),
        (
            "vertices",
            np.vstack([west(city).get_geometry(lod="1.2").vertices, [1, 2, 3]]),
            "Unused",
        ),
    ]:
        model = exchange.loads(exchange.dumps(city))
        setattr(west(model).get_geometry(lod="1.2"), field, value)
        with pytest.raises(NotImplementedError, match=message):
            model.save(path, strict=True, validate_schema=False)
        assert path.read_bytes() == before
    west(city).semantic_type = None
    # Generic Terrain remains valid native data, without guessing a TIN class.
    exchange.dumps(city)
    with pytest.raises(NotImplementedError, match="semantic URI"):
        city.save(path, strict=True)
    assert path.read_bytes() == before


def test_quantization_collapse_and_geometry_edit(tmp_path):
    city = io.load_cityjson(source(), strict=True)
    mesh = west(city).get_geometry(lod="1.2")
    mesh.vertices[0, 2] += 0.123
    path = tmp_path / "terrain.city.json"
    city.save(path, strict=True)
    verify(city, io.load_city(path, strict=True))
    before = path.read_bytes()
    mesh.vertices[1] = mesh.vertices[0] + [0.00001, 0, 0]
    with pytest.raises(ValueError, match="collapses a TINRelief triangle"):
        city.save(path, strict=True)
    assert path.read_bytes() == before
