"""Mixed-city source semantics and existing native types have faithful boundaries."""

import json
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.io.cityjson.cityjson import load
from dtcc_core.model import (
    City,
    Object,
    Tree,
    Landuse,
    Point,
    exchange,
    dtcc_pb2 as wire,
)
from dtcc_core.model.object.landuse import LanduseClasses

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "sandbox/model_profiles/fixtures/mixed-city.city.json"
)
NS = "https://github.com/dtcc-platform/dtcc-core/schemas/model#"


def test_public_mixed_source_roundtrip_and_access(tmp_path):
    city = io.load_city(FIXTURE, strict=True)
    assert city.buildings[0].footprint().to_polygon().area == 96
    park = city.get_children(Landuse)[0]
    assert park.id == "park-1" and park.landuses == []
    assert park.lod0.surfaces[0].to_polygon().area == 360
    generic = {obj.id: obj for obj in city.get_children(Object)}
    plant, bench = generic["tree-1"], generic["bench-1"]
    assert plant.semantic_type == NS + "SolitaryVegetationObject"
    assert plant.attributes["height"] == 8.5 and not city.trees
    assert plant.attributes["species"] == "Quercus robur"
    assert bench.semantic_type == NS + "CityFurniture" and type(bench.lod1) is Point
    city.save(tmp_path / "mixed.dtcc")
    restored = io.load_model(tmp_path / "mixed.dtcc")
    assert exchange.dumps(restored) == exchange.dumps(city)
    restored.save(tmp_path / "mixed.city.json", strict=True)
    external = io.load_city(tmp_path / "mixed.city.json", strict=True)
    # Source object order and native aggregate UUID are not external model facts.
    originals = {obj.id: obj for group in city.children.values() for obj in group}
    for group in external.children.values():
        for obj in group:
            assert exchange.dumps(obj) == exchange.dumps(originals[obj.id])


def test_native_tree_landuse_state_is_exact_and_independently_mutable(tmp_path):
    tree = Tree(
        id="tree",
        position=np.array([325000.123456789, 6400000.987654321, 12.125]),
        height=8.123456789,
        crown_radius=2.123456789,
        semantic_type=NS + "Tree",
    )
    # Additional representations/metadata are separate native facts, not inferred.
    tree.add_geometry(Point(x=1.125), id="survey", role="survey")
    tree.attributes["source"] = {"labels": ["oak"]}
    park = Landuse(
        id="park",
        landuses=[LanduseClasses.GRASS, LanduseClasses.UNKNOWN, LanduseClasses.GRASS],
    )
    city = City(id="city")
    city.add_children([tree, park])
    city.save(tmp_path / "native.dtcc")
    restored = io.load_model(tmp_path / "native.dtcc")
    result = restored.trees[0]
    assert type(restored.get_children(Landuse)[0]) is Landuse
    assert exchange.dumps(restored) == exchange.dumps(city)
    np.testing.assert_array_equal(result.position, tree.position)
    assert result.position.dtype == tree.position.dtype
    assert (result.height, result.crown_radius) == (tree.height, tree.crown_radius)
    result.position[0] += 1
    result.attributes["source"]["labels"].append("changed")
    restored.get_children(Landuse)[0].landuses.clear()
    assert tree.position[0] == 325000.123456789 and tree.attributes["source"][
        "labels"
    ] == ["oak"]
    assert len(park.landuses) == 3
    # Intentional native defaults and exact integer positions are not normalized.
    for value in (
        Tree(),
        Tree(position=np.empty(0)),
        Tree(position=np.array([2**60 + 1, 2, 3], dtype=np.uint64)),
        Landuse(),
    ):
        assert exchange.dumps(exchange.loads(exchange.dumps(value))) == exchange.dumps(
            value
        )


def test_new_typed_wire_state_is_required_and_versioned():
    pb = wire.ModelFile.FromString(exchange.dumps(Tree()))
    assert pb.version == exchange.VERSION
    pb.version = 4
    with pytest.raises(ValueError, match="Unsupported model"):
        exchange.loads(pb.SerializeToString())
    pb.version = exchange.VERSION
    pb.object.ClearField("tree")
    with pytest.raises(ValueError, match="state"):
        exchange.loads(pb.SerializeToString())
    pb.object.tree.height = 1
    pb.object.kind = wire.Object.OBJECT
    with pytest.raises(ValueError, match="state"):
        exchange.loads(pb.SerializeToString())
    pb = wire.ModelFile.FromString(exchange.dumps(Landuse()))
    pb.object.landuse.landuses.append("NOT_A_NATIVE_CODE")
    with pytest.raises(ValueError, match="Landuse code"):
        exchange.loads(pb.SerializeToString())


@pytest.mark.parametrize(
    "value",
    [
        Tree(position=np.zeros(2)),
        Tree(position=np.array([np.nan, 0, 0])),
        Tree(height=-1),
        Tree(crown_radius=np.inf),
        Tree(height=True),
        Landuse(landuses=["GRASS"]),
    ],
)
def test_invalid_native_state_fails_without_replacing_file(value, tmp_path):
    path = tmp_path / "model.dtcc"
    path.write_bytes(b"existing")
    with pytest.raises(ValueError):
        io.save_model(value, path)
    assert path.read_bytes() == b"existing"


def test_missing_plant_measurements_are_not_fabricated():
    source = json.loads(FIXTURE.read_text())
    plant = source["CityObjects"]["tree-1"]
    plant["attributes"] = {"name": "Unclassified shrub"}
    city = load(source, strict=True)
    value = next(obj for obj in city.get_children(Object) if obj.id == "tree-1")
    assert value.attributes == plant["attributes"] and not city.trees


@pytest.mark.parametrize(
    "failure",
    [
        "point_count",
        "point_index",
        "point_field",
        "landuse_solid",
        "feature_child",
        "feature_type",
        "region",
    ],
)
def test_strict_source_rejects_unsupported_or_invalid_state(failure):
    source = json.loads(FIXTURE.read_text())
    plant = source["CityObjects"]["tree-1"]["geometry"][0]
    if failure == "point_count":
        plant["boundaries"] *= 2
    elif failure == "point_index":
        plant["boundaries"] = [True]
    elif failure == "point_field":
        plant["semantics"] = {"surfaces": [], "values": [None]}
    elif failure == "landuse_solid":
        park = source["CityObjects"]["park-1"]["geometry"][0]
        park["type"], park["boundaries"] = "Solid", [park["boundaries"]]
    elif failure == "feature_child":
        source["CityObjects"]["park-1"]["children"] = ["tree-1"]
    elif failure == "feature_type":
        source["CityObjects"]["park-1"]["type"] = []
    else:
        source["CityObjects"]["park-1"]["geometry"][0]["semantics"] = {
            "surfaces": [{"type": "GroundSurface"}],
            "values": [0],
        }
    with pytest.raises((ValueError, NotImplementedError)):
        load(source, strict=True)


def test_strict_export_cannot_drop_native_state_or_invent_source_types(tmp_path):
    city = io.load_city(FIXTURE, strict=True)
    park = city.get_children(Landuse)[0]
    park.landuses = [LanduseClasses.GRASS]
    with pytest.raises(NotImplementedError, match="Landuse codes"):
        city.save(tmp_path / "park.city.json", strict=True)
    park.landuses = []
    city.add_child(Tree(id="native-tree", semantic_type=NS + "Tree"))
    with pytest.raises(NotImplementedError, match="mapping"):
        city.save(tmp_path / "tree.city.json", strict=True)
