"""Exact attachment selection, explicit Solid topology and reader migration."""

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.io.cityjson.cityjson import load
from dtcc_core.model import (
    Building,
    City,
    Field,
    GeometryType,
    Object,
    Point,
    Solid,
    Surface,
    exchange,
)


def solid_document():
    cube = [
        [0, 0, 0],
        [4, 0, 0],
        [4, 4, 0],
        [0, 4, 0],
        [0, 0, 4],
        [4, 0, 4],
        [4, 4, 4],
        [0, 4, 4],
    ]
    inner = [[1 + x // 2, 1 + y // 2, 1 + z // 2] for x, y, z in cube]
    faces = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    geometry = {
        "type": "Solid",
        "lod": "2.2",
        "boundaries": [
            [[f] for f in faces],
            [[[i + 8 for i in reversed(f)]] for f in faces],
        ],
        "semantics": {
            "surfaces": [{"type": "WallSurface", "material": "brick"}],
            "values": [[None, None, 0, 0, 0, 0], [None, None, 0, 0, 0, 0]],
        },
    }
    return {
        "type": "CityJSON",
        "version": "2.0",
        "transform": {"scale": [1, 1, 1], "translate": [0, 0, 0]},
        "vertices": cube + inner,
        "CityObjects": {"building": {"type": "Building", "geometry": [geometry]}},
    }


def test_exact_selection_and_ordered_canonical_attachment_roundtrip():
    obj = Object()
    a, b = Point(x=1), Point(x=2)
    obj.add_geometry(a, id="z-survey", lod="2.2")
    obj.add_geometry(b, id="a-simulation", lod="2.2", role="simulation")
    assert obj.lod2 is None
    assert obj.get_geometry(id="z-survey") is a
    assert obj.get_geometry(id="z-survey", role="simulation") is None
    assert obj.get_geometry(lod="2.2", role="simulation") is b
    assert obj.get_geometries(lod="2.2") == [a, b]
    with pytest.raises(ValueError, match="matching IDs"):
        obj.get_geometry(lod="2.2")
    with pytest.raises(ValueError, match="requires"):
        obj.get_geometry()
    with pytest.raises(ValueError, match="Duplicate"):
        obj.add_geometry(b, id="z-survey")
    restored = exchange.loads(exchange.dumps(obj))
    assert list(restored.geometry) == ["z-survey", "a-simulation"]
    assert restored.geometry["a-simulation"].role == "simulation"
    obj.add_geometry(a, GeometryType.LOD2)
    assert obj.lod2 is a
    obj.add_geometry(b, GeometryType.LOD2)
    assert obj.lod2 is b
    obj.add_geometry(a, id="other-exact", lod="2")
    with pytest.raises(ValueError, match="Ambiguous"):
        _ = obj.lod2


def test_solid_shells_regions_fields_and_transforms_roundtrip(tmp_path):
    solid = load(solid_document(), strict=True).buildings[0].get_geometry(lod="2.2")
    assert type(solid) is Solid and len(solid.shells) == 2
    solid.shells = [s.astype(np.uint32) for s in solid.shells]
    solid.regions[0].id = "walls"
    solid.transform.set_translation(0.25, 0.5, 1.25)
    solid.fields = [
        Field(name="flux", association="face", values=np.arange(12, dtype=float))
    ]
    io.save_model(solid, tmp_path / "solid.dtcc")
    result = io.load_model(tmp_path / "solid.dtcc")
    assert exchange.dumps(result) == exchange.dumps(solid)
    assert result.shells[1].dtype == np.uint32
    assert result.regions[0].indices.tolist() == [2, 3, 4, 5, 8, 9, 10, 11]
    assert result.to_proto().geometry.HasField("solid")


def test_decoded_polygon_arrays_and_transforms_remain_independent():
    source = load(solid_document(), strict=True).buildings[0].get_geometry(lod="2.2")
    payload = exchange.dumps(source)
    restored = exchange.loads(payload)
    restored.surfaces[0].vertices[0, 0] = 99
    restored.surfaces[0].transform.set_translation(1, 2, 3)
    restored.shells[0][0] = 1
    restored.regions[0].indices[0] = 0
    assert exchange.dumps(source) == payload
    np.testing.assert_array_equal(restored.surfaces[1].transform.affine, np.eye(4))
    assert exchange.dumps(exchange.loads(payload)) == payload


@pytest.mark.parametrize("failure", ["duplicate", "uncovered", "range"])
def test_shells_must_partition_surfaces(failure):
    solid = load(solid_document(), strict=True).buildings[0].get_geometry(lod="2.2")
    if failure == "duplicate":
        solid.shells[1][0] = 0
    elif failure == "uncovered":
        solid.shells.pop()
    else:
        solid.shells[1][0] = 12
    with pytest.raises(ValueError, match="shell"):
        exchange.dumps(solid)


def test_cityjson_repeated_lod_and_interior_shell_roundtrip(tmp_path):
    from copy import deepcopy

    source = solid_document()
    source["CityObjects"]["building"]["geometry"].append(
        deepcopy(source["CityObjects"]["building"]["geometry"][0])
    )
    city = load(source, strict=True)
    building = city.buildings[0]
    assert list(building.geometry) == ["cityjson-0", "cityjson-1"]
    with pytest.raises(ValueError, match="Ambiguous"):
        building.get_geometry(lod="2.2")
    city.save(tmp_path / "repeated.city.json", strict=True)
    result = io.load_city(tmp_path / "repeated.city.json", strict=True)
    result.id = city.id  # CityJSON has no native aggregate UUID.
    assert exchange.dumps(result) == exchange.dumps(city)
    building.geometry["cityjson-0"].role = "simulation"
    with pytest.raises(NotImplementedError, match="IDs/roles"):
        city.save(tmp_path / "repeated.city.json", strict=True)


def test_protobuf_preserves_representation_order():
    obj = Object()
    obj.add_geometry(Point(x=1), "z-survey")
    obj.add_geometry(Point(x=2), "a-survey")
    restored = Object()
    restored.from_proto(obj.to_proto())
    assert list(restored.geometry) == ["z-survey", "a-survey"]


def test_parent_bounds_do_not_mutate_child_or_first_surface_bounds():
    city = City()
    for offset in (0.0, 100.0):
        building = Building()
        building.add_geometry(
            Surface(
                vertices=np.array(
                    [[offset, 0.0, 0.0], [offset + 1, 0.0, 0.0], [offset, 1.0, 0.0]]
                )
            ),
            GeometryType.LOD0,
        )
        city.add_child(building)
    assert city.calculate_bounds().xmax == 101
    assert city.buildings[0].bounds.xmax == 1
    assert city.buildings[0].lod0.bounds.xmax == 1


def test_extent_recompute_is_explicit_and_preserves_source_evidence(tmp_path):
    from dtcc_core.datasets import load_model_package

    source = solid_document()
    source["metadata"] = {"geographicalExtent": [0, 0, 0, 99, 99, 99]}
    with pytest.raises(ValueError, match="extent_policy"):
        load(source, strict=True)
    city = load(source, strict=True, extent_policy="recompute")
    assert city.bounds.xmax == 4
    assert city.dataset_context.health["extent_discrepancy_count"] == 1
    evidence = city.dataset_context.provenance.processing_steps[0]["discrepancies"][0]
    assert evidence["source_extent"] == source["metadata"]["geographicalExtent"]
    package = city.export(tmp_path / "source.dtccpkg", canonical=True)
    result = load_model_package(package.path)
    assert result.dataset_context == city.dataset_context
    assert exchange.dumps(result) == exchange.dumps(city)


def test_cityjson_preserves_existing_degeneracy_but_rejects_new_quantization_loss(
    tmp_path,
):
    source = solid_document()
    source["vertices"][1] = list(source["vertices"][0])
    city = load(source, strict=True)
    path = tmp_path / "degenerate.city.json"
    city.save(path, strict=True)
    result = io.load_city(path, strict=True)
    result.id = city.id
    assert exchange.dumps(result) == exchange.dumps(city)
    city = load(solid_document(), strict=True)
    previous = path.read_bytes()
    with pytest.raises(ValueError, match="quantization"):
        city.save(path, strict=True, scale=100)
    assert path.read_bytes() == previous
