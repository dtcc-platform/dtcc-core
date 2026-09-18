"""Faithful buildings subset through ordinary CityJSON and canonical entry points."""

import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.io.cityjson import cityjson
from dtcc_core.io.cityjson.semantics import SEMANTIC_NAMESPACE
from dtcc_core.model import (
    Field,
    GeometryType,
    Mesh,
    MultiSurface,
    SemanticRegion,
    Surface,
    exchange,
)
from dtcc_core.model import dtcc_pb2 as wire

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "sandbox/model_profiles/fixtures/buildings.city.json"
)


@pytest.fixture
def source():
    return json.loads(FIXTURE.read_text())


def test_public_buildings_workflow_preserves_identity_geometry_and_regions(tmp_path):
    city = io.load_city(FIXTURE, strict=True)
    building = city.buildings[0]
    part = building.building_parts[0]
    assert part.id == "part-1" and part.attributes["name"] == "West wing"
    footprint = building.footprint().to_polygon()
    assert footprint.area == 96 and len(footprint.interiors) == 1
    assert building.attributes["measured_height"] == 10.2
    roof = part.lod2.regions_of(SEMANTIC_NAMESPACE + "RoofSurface")[0]
    assert roof.indices.tolist() == [1, 2]
    assert roof.attributes == {"solar-potential": 42.5}
    city.save(tmp_path / "model.dtcc")
    result = io.load_model(tmp_path / "model.dtcc")
    assert exchange.dumps(result) == exchange.dumps(city)
    result.save(tmp_path / "output.city.json", strict=True)
    exported = json.loads((tmp_path / "output.city.json").read_text())
    assert exported["CityObjects"]["part-1"]["geometry"][0]["semantics"]["values"] == [
        0,
        1,
        1,
        2,
        2,
        2,
        None,
    ]
    assert exported["metadata"]["referenceSystem"] == city.transform.srs
    restored = io.load_city(tmp_path / "output.city.json", strict=True)
    assert restored.buildings[0].building_parts[0].id == part.id
    for original, final in zip(
        part.lod2.surfaces, restored.buildings[0].building_parts[0].lod2.surfaces
    ):
        np.testing.assert_allclose(
            original.vertices, final.vertices, atol=0.0005, rtol=0
        )
        for a, b in zip(original.holes, final.holes):
            np.testing.assert_allclose(a, b, atol=0.0005, rtol=0)


@pytest.mark.parametrize(
    "failure",
    [
        "index",
        "semantic_index",
        "solid",
        "invalid_lod",
        "parent",
        "appearance",
        "extent",
        "boolean_vertex",
    ],
)
def test_strict_import_rejects_unrepresented_or_invalid_facts(source, failure):
    part = source["CityObjects"]["part-1"]
    geometry = part["geometry"][0]
    if failure == "index":
        geometry["boundaries"][0][0][0] = -1
    elif failure == "semantic_index":
        geometry["semantics"]["values"][0] = 42
    elif failure == "solid":
        geometry["type"] = "Solid"
    elif failure == "invalid_lod":
        geometry["lod"] = "2.4"
    elif failure == "parent":
        part["parents"] = ["missing"]
    elif failure == "extent":
        part["geographicalExtent"] = [0] * 6
    elif failure == "boolean_vertex":
        source["vertices"][0][0] = True
    else:
        source["appearance"] = {}
    with pytest.raises((ValueError, NotImplementedError)):
        cityjson.load(source, strict=True)


def test_all_lods_and_nested_part_ids_survive(source, tmp_path):
    child = deepcopy(source["CityObjects"]["part-1"])
    child["parents"] = ["part-1"]
    source["CityObjects"]["part-1"]["children"] = ["nested-part"]
    source["CityObjects"]["nested-part"] = child
    source["CityObjects"]["building-1"]["geometry"].append(
        deepcopy(child["geometry"][0])
    )
    city = cityjson.load(source, strict=True)
    assert city.buildings[0].lod0 is not None and city.buildings[0].lod2 is not None
    city.save(tmp_path / "nested.json", strict=True)
    restored = io.load_city(tmp_path / "nested.json", strict=True)
    assert restored.buildings[0].building_parts[0].building_parts[0].id == "nested-part"


def test_canonical_surface_preserves_nested_metadata_and_requires_v2():
    surface = Surface(
        vertices=np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]),
        normal=np.array([0.0, 0.0, 1.0]),
    )
    surface.transform.set_translation(1.125, 0, 0)
    surface.add_field(Field(name="flux", association="face", values=np.array([11.25])))
    source = MultiSurface(
        surfaces=[surface],
        regions=[SemanticRegion("urn:example:roof", np.array([0]), id="roof")],
    )
    data = exchange.dumps(source)
    result = exchange.loads(data)
    np.testing.assert_array_equal(
        result.surfaces[0].transform.affine, surface.transform.affine
    )
    np.testing.assert_array_equal(result.surfaces[0].normal, surface.normal)
    assert result.regions[0].id == "roof"
    assert result.surfaces[0].fields[0].values.tolist() == [11.25]
    pb = wire.ModelFile.FromString(data)
    assert pb.version == exchange.VERSION
    pb.version = 1
    with pytest.raises(ValueError, match="Unsupported model"):
        exchange.loads(pb.SerializeToString())
    assert len(source.to_proto().geometry.regions) == 1


def test_region_indices_remain_generic_for_mesh_faces():
    mesh = Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        faces=np.array([[0, 1, 2]]),
    )
    mesh.regions = [SemanticRegion("urn:example:roof", np.array([0], dtype=np.uint32))]
    result = exchange.loads(exchange.dumps(mesh))
    assert type(result.regions[0]) is SemanticRegion
    assert result.regions[0].indices.dtype == np.uint32
    mesh.regions.append(SemanticRegion("urn:example:wall", np.array([0])))
    with pytest.raises(ValueError, match="only one"):
        exchange.dumps(mesh)


def test_export_failure_leaves_existing_destination_untouched(tmp_path):
    city = io.load_city(FIXTURE, strict=True)
    city.relations = {"observes": ["part-1"]}
    output = tmp_path / "existing.json"
    output.write_text("previous")
    with pytest.raises(NotImplementedError, match="relations"):
        city.save(output, strict=True)
    assert output.read_text() == "previous"


def test_missing_geometry_and_unclassified_surfaces_remain_unclassified(
    source, tmp_path
):
    source["CityObjects"]["building-1"]["geometry"] = []
    del source["CityObjects"]["part-1"]["geometry"][0]["semantics"]
    city = cityjson.load(source, strict=True)
    city.save(tmp_path / "empty.json", strict=True)
    result = io.load_city(tmp_path / "empty.json", strict=True)
    assert result.buildings[0].lod0 is None
    assert result.bounds.xmin > 300000
    assert result.buildings[0].building_parts[0].lod2.regions == []


def test_empty_multisurface_is_canonical_but_not_cityjson(source, tmp_path):
    source["CityObjects"]["building-1"]["geometry"] = [
        {"type": "MultiSurface", "lod": "0", "boundaries": []}
    ]
    with pytest.raises(ValueError, match="nonempty"):
        cityjson.load(source, strict=True)
    city = io.load_city(FIXTURE, strict=True)
    city.buildings[0].geometry["cityjson-0"].geometry = MultiSurface()
    restored = exchange.loads(exchange.dumps(city))
    assert restored.buildings[0].lod0.surfaces == []
    assert restored.bounds.xmin > 300000
    with pytest.raises(ValueError, match="nonempty"):
        restored.save(tmp_path / "empty.json", strict=True)


def test_duplicate_cityjson_ids_fail_at_the_file_boundary(tmp_path):
    path = tmp_path / "duplicate.json"
    path.write_text('{"type":"CityJSON","CityObjects":{"same":{},"same":{}}}')
    with pytest.raises(ValueError, match="Duplicate"):
        io.load_city(path, strict=True)


def test_multisurface_merge_offsets_regions_without_mutating_source():
    polygon = Surface(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    )
    left = MultiSurface(surfaces=[polygon] * 300)
    right = MultiSurface(
        surfaces=[polygon.copy()],
        regions=[SemanticRegion("urn:roof", np.array([0], dtype=np.uint8))],
    )
    left.merge(right)
    assert left.regions[0].indices.tolist() == [300]
    assert right.regions[0].indices.tolist() == [0]
    assert exchange.loads(exchange.dumps(left)).regions[0].indices.tolist() == [300]


def test_cxx_boundary_reports_unsupported_regions_before_conversion():
    from dtcc_core.builder.model_conversion import (
        create_builder_multisurface,
        mesh_to_builder_mesh,
    )

    geometry = io.load_city(FIXTURE, strict=True).buildings[0].building_parts[0].lod2
    with pytest.raises(NotImplementedError, match="semantic regions"):
        create_builder_multisurface(geometry)
    mesh = Mesh(regions=[SemanticRegion("urn:roof", np.empty(0, dtype=np.int64))])
    with pytest.raises(NotImplementedError, match="semantic regions"):
        mesh_to_builder_mesh(mesh)


def test_direct_cityjson_geometry_converter_preserves_explicit_regions():
    from dtcc_core.io.cityjson.converters import convert_multisurface

    geometry = io.load_city(FIXTURE, strict=True).buildings[0].building_parts[0].lod2
    encoded = convert_multisurface(geometry, [], 1000)
    assert encoded["semantics"]["values"] == [0, 1, 1, 2, 2, 2, None]
    assert encoded["semantics"]["surfaces"][1]["type"] == "RoofSurface"


def test_strict_export_refreshes_first_polygon_bounds(tmp_path):
    city = io.load_city(FIXTURE, strict=True)
    original = city.bounds.xmin
    city.buildings[0].lod0.surfaces[0].vertices[0, 0] -= 1
    city.save(tmp_path / "mutated.json", strict=True)
    restored = io.load_city(tmp_path / "mutated.json", strict=True)
    assert restored.bounds.xmin == pytest.approx(original - 1)
