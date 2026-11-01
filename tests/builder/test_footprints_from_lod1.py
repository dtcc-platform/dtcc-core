"""LOD1 adapter policies, provenance and the cleaning/rebuild handoff."""

import json
from copy import deepcopy

import numpy as np
import pytest
from shapely.geometry import MultiPolygon, Polygon, box

from dtcc_core import io
from dtcc_core.builder import (
    build_lod1_buildings,
    building_from_footprint,
    clean_footprints,
    footprints_from_lod1,
)
from dtcc_core.builder.cleaning.blocks import WALLS_VERSION, walls_of
from dtcc_core.model import Building, BuildingPart, City, Field, GeometryType, Object
from dtcc_core.model._standard_schema import SEMANTIC_NAMESPACE
from tests.builder.test_prism_decomposition import BASE, TOWER, _input, _prism, _stack


def _source(id="podium", faces=None):
    value = Building(id=id)
    value.add_geometry(
        _input(_stack(BASE, TOWER) if faces is None else faces, "surfaces"),
        GeometryType.LOD1,
    )
    return value


def _polygon(building):
    return Polygon(
        building.lod0.vertices[:, :2], [r[:, :2] for r in building.lod0.holes]
    )


def _signature(buildings):
    return [
        (
            b.id,
            b.attributes,
            b.lod0.vertices.tolist(),
            [r.tolist() for r in b.lod0.holes],
        )
        for b in buildings
    ]


def _city(*buildings):
    city = City(id="city")
    city.add_buildings(list(buildings))
    return city


def _volume(buildings):
    total = 0
    for building in buildings:
        targets = building.building_parts or [building]
        for target in targets:
            mesh = target.lod1.mesh(weld=True)
            xyz = mesh.vertices[mesh.faces] - mesh.vertices[0]
            total += (
                np.einsum("ij,ij->i", xyz[:, 0], np.cross(xyz[:, 1], xyz[:, 2])).sum()
                / 6
            )
    return total


def test_split_buildings_have_authoritative_heights_and_provenance_without_mutation():
    source = _source()
    source.attributes = dict(
        estimated_height=30,
        measured_height=100,
        ground_height=-20,
        merged_from=["stale"],
        source_footprint={},
        footprint_walls={},
        lod1_extrusion={},
        custom={"values": [1]},
    )
    before = deepcopy(source)
    result = footprints_from_lod1([source])
    assert [b.id for b in result] == ["podium#1", "podium#2"]
    assert [b.measured_height for b in result] == [6, 30]
    assert [b.lod0.zmax for b in result] == [6, 30]
    assert all(
        b.attributes["ground_height"] == 0 and b.estimated_height is None
        for b in result
    )
    assert _polygon(result[0]).equals(BASE.difference(TOWER))
    assert _polygon(result[1]).equals(TOWER)
    assert [b.attributes["prism_decomposition"] for b in result] == [
        dict(version=1, input_id="podium", component_index=0, level_index=i)
        for i in range(2)
    ]
    assert all(b.attributes["source_id"] == "podium" for b in result)
    assert all(
        not any(
            k in b.attributes
            for k in ("merged_from", "source_footprint", "footprint_walls")
        )
        for b in result
    )
    assert _volume(result) == pytest.approx(3264, abs=1e-8, rel=0)
    assert source.attributes == before.attributes
    for actual, expected in zip(
        source.lod1.surfaces, before.lod1.surfaces, strict=True
    ):
        np.testing.assert_array_equal(actual.vertices, expected.vertices)
    result[0].attributes["custom"]["values"].append(2)
    assert source.attributes["custom"]["values"] == [1]
    assert result[1].attributes["custom"]["values"] == [1]


def test_single_output_retains_id_and_upstream_source_id():
    source = _source("one", _prism(BASE, 2, 12))
    source.attributes["source_id"] = "survey-42"
    result = footprints_from_lod1((source,))
    assert len(result) == 1 and result[0].id == "one"
    assert result[0].attributes["source_id"] == "survey-42"
    assert (
        result[0].measured_height == 10 and result[0].attributes["ground_height"] == 2
    )


def test_equal_height_towers_share_level_index_and_face_order_keeps_ids():
    towers = MultiPolygon([box(2, 2, 6, 6), box(14, 14, 18, 18)])
    source = _source(faces=_stack(BASE, towers))
    expected = footprints_from_lod1([source])
    assert [b.attributes["prism_decomposition"]["level_index"] for b in expected] == [
        0,
        1,
        1,
    ]
    source.lod1.surfaces.reverse()
    for surface in source.lod1.surfaces:
        surface.vertices = surface.vertices[::-1].copy()
    assert _signature(footprints_from_lod1([source])) == _signature(expected)


def test_ids_reserve_city_nonbuilding_building_and_descendant_ids_first():
    source = _source("b")
    child = BuildingPart(id="b#4")
    child.add_geometry(
        _input(_prism(box(30, 0, 40, 10), 2, 7), "surfaces"), GeometryType.LOD1
    )
    source.add_child(child)
    preserved = building_from_footprint(box(50, 0, 60, 10), 10, id="b#3")
    city = _city(source, preserved)
    city.id = "b#1"
    city.add_child(Object(id="b#2"))
    result = footprints_from_lod1(city)
    assert [b.id for b in result.buildings] == ["b#5", "b#6", "b#7", "b#3"]
    assert [
        b.attributes["prism_decomposition"]["component_index"]
        for b in result.buildings[:3]
    ] == [0, 0, 1]


def test_city_copies_nonbuildings_context_and_reports_and_refreshes_bounds():
    city = _city(_source())
    city.attributes = {"custom": {"values": [1]}, "footprint_cleaning": {"stale": True}}
    landmark = Object(id="landmark", attributes={"label": ["kept"]})
    city.add_child(landmark)
    city.transform.set_translation(100, 200, 30)
    result = footprints_from_lod1(city)
    assert isinstance(result, City) and result.id == city.id
    np.testing.assert_array_equal(result.transform.affine, city.transform.affine)
    assert "footprint_cleaning" not in result.attributes
    assert "footprint_cleaning" in city.attributes
    report = result.attributes["footprints_from_lod1"]
    assert (
        report["attempted"],
        report["converted"],
        report["preserved"],
        report["output_count"],
    ) == (1, 1, 0, 2)
    assert report["failures"] == [] and len(report["residuals"]) == 1
    assert report["residuals"][0]["components"][0]["prism_volume"] == 3264
    json.dumps(report, allow_nan=False)
    result.children[Object][0].attributes["label"].append("independent")
    assert landmark.attributes["label"] == ["kept"]
    assert (result.bounds.xmax, result.bounds.zmin, result.bounds.zmax) == (20, 0, 30)


def test_existing_lod0_is_kept_by_value_and_repeated_conversion_is_idempotent():
    source = building_from_footprint(BASE, 10, ground_height=2, id="keep")
    before = _signature([source])
    first = footprints_from_lod1([source])
    assert _signature(first) == before and first[0] is not source
    second = footprints_from_lod1(footprints_from_lod1([_source()]))
    assert _signature(second) == _signature(footprints_from_lod1([_source()]))
    source.attributes["measured_height"] = 99
    replaced = footprints_from_lod1([source], overwrite_lod0=True)[0]
    assert replaced.measured_height == 10 and replaced.lod0.zmax == 12


def test_warn_preserves_complete_failed_sources_atomically_and_reports(caplog):
    bad = _source("bad", _prism(BASE) + _prism(box(30, 0, 40, 10))[:-1])
    bad.attributes["payload"] = {"keep": [1]}
    missing = Building(id="missing")
    city = _city(_source("good", _prism(BASE)), bad, missing)
    result = footprints_from_lod1(city, on_error="warn")
    assert [b.id for b in result.buildings] == ["good", "bad", "missing"]
    assert len(result.buildings[1].lod1.surfaces) == len(bad.lod1.surfaces)
    assert result.buildings[1] is not bad
    assert result.buildings[1].attributes == bad.attributes
    report = result.attributes["footprints_from_lod1"]
    assert (report["converted"], report["preserved"]) == (1, 2)
    assert [f["reason"] for f in report["failures"]] == [
        "open_shell",
        "missing_geometry",
    ]
    assert (
        "bad" in caplog.text
        and "open_shell" in caplog.text
        and "location=" in caplog.text
    )
    with pytest.raises(ValueError, match="Building 'bad'.*open_shell"):
        footprints_from_lod1(city)
    assert all(b.lod0 is None for b in city.buildings)


@pytest.mark.parametrize(
    "case",
    [
        "extra_lod",
        "role",
        "child_attributes",
        "child_semantics",
        "fields",
        "profile",
        "relations",
        "incoming_relation",
        "overwrite_fields",
        "overwrite_mesh_markers",
    ],
)
def test_untransferable_metadata_is_reported_without_loss(case):
    source = _source()
    city = _city(source)
    if case == "extra_lod":
        source.add_geometry(source.lod1.copy(), GeometryType.LOD2)
    elif case == "role":
        source.add_geometry(source.lod1.copy(), "survey")
    elif case in ("child_attributes", "child_semantics"):
        part = BuildingPart(id="part")
        source.add_child(part)
        if case == "child_attributes":
            part.attributes = {"name": "a"}
        else:
            part.semantic_type = "urn:custom:room"
    elif case == "fields":
        source.lod1.surfaces[0].fields = [Field(name="temperature")]
    elif case == "profile":
        source.profile_id, source.profile_version = "urn:profile", "1"
    elif case == "relations":
        source.relations = {"related": ["city"]}
    elif case == "incoming_relation":
        city.add_child(Object(id="referrer", relations={"related": [source.id]}))
    elif case == "overwrite_mesh_markers":
        footprint = _input(_prism(BASE), "mesh")
        footprint.markers = np.zeros(len(footprint.faces), dtype=int)
        source.add_geometry(footprint, GeometryType.LOD0)
    else:
        footprint = building_from_footprint(BASE, 10).lod0
        footprint.fields = [Field(name="temperature")]
        source.add_geometry(footprint, GeometryType.LOD0)
    with pytest.raises(ValueError, match="unsupported_metadata"):
        footprints_from_lod1(city, overwrite_lod0=True)
    result = footprints_from_lod1(city, overwrite_lod0=True, on_error="warn")
    assert result.buildings[0].id == source.id
    assert set(result.buildings[0].geometry) == set(source.geometry)
    assert result.buildings[0].attributes == source.attributes


def test_geometry_only_parts_and_multiple_representations_are_supported():
    source = Building(id="source", semantic_type=SEMANTIC_NAMESPACE + "Building")
    part = BuildingPart(id="part", semantic_type=SEMANTIC_NAMESPACE + "BuildingPart")
    part.add_geometry(_input(_prism(BASE, 2, 12), "mesh"), id="left", lod="1")
    part.add_geometry(
        _input(_prism(box(30, 0, 40, 10), 3, 7), "soup"), id="right", lod="1"
    )
    source.add_child(part)
    result = footprints_from_lod1([source])
    assert len(result) == 2 and [b.attributes["ground_height"] for b in result] == [
        2,
        3,
    ]
    assert all(
        not b.children and b.semantic_type == source.semantic_type for b in result
    )


def test_composed_translations_and_srs_preserve_physical_placement_once():
    source = Building(id="source")
    source.transform.set_translation(1, 2, 3)
    part = BuildingPart(id="part")
    part.transform.set_translation(10, 20, 30)
    geometry = _input(_prism(BASE, 2, 12), "surfaces")
    geometry.transform.set_translation(100, 200, 300)
    for surface in geometry.surfaces:
        surface.transform.set_translation(1000, 2000, 3000)
    part.add_geometry(geometry, GeometryType.LOD1)
    source.add_child(part)
    city = _city(source)
    city.transform.srs = "EPSG:3006"
    city.transform.set_translation(5, 6, 7)
    result = footprints_from_lod1(city)
    piece = result.buildings[0]
    np.testing.assert_array_equal(piece.transform.affine, np.eye(4))
    assert piece.lod0.transform.srs == "EPSG:3006"
    np.testing.assert_array_equal(piece.lod0.transform.affine, np.eye(4))
    assert piece.attributes["prism_decomposition"]["input_translation"] == [
        1111,
        2222,
        3333,
    ]
    expected = np.unique(np.vstack([s.vertices for s in geometry.surfaces]), axis=0) + [
        1116,
        2228,
        3340,
    ]
    actual = np.unique(
        np.vstack(
            [result.transform(s.transform(s.vertices)) for s in piece.lod1.surfaces]
        ),
        axis=0,
    )
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-9)
    assert piece.attributes["ground_height"] == 3335 and piece.measured_height == 10


def test_translation_survives_cleaning_rebuild_and_native_roundtrip(tmp_path):
    from dtcc_core.model import Point

    source = _source()
    source.lod1.transform.set_translation(100000, 200000, 30)
    city = _city(source)
    city.transform.srs = "EPSG:3006"
    city.transform.set_translation(0, 0, 10)
    landmark = Object(id="landmark")
    landmark.add_geometry(Point(x=-10, y=-20, z=-30), "location")
    city.add_child(landmark)
    derived = footprints_from_lod1(city)
    assert [b.lod0.zmax for b in derived.buildings] == [36, 60]
    assert all(b.attributes["ground_height"] == 30 for b in derived.buildings)
    cleaned = clean_footprints(derived, min_building_detail=0.1)
    build_lod1_buildings(cleaned.buildings, rebuild=False)
    assert _volume(cleaned.buildings) == pytest.approx(3264, abs=1e-8, rel=0)
    points = np.vstack(
        [
            s.vertices
            for p in cleaned.buildings[0].building_parts
            for s in p.lod1.surfaces
        ]
    )
    physical = cleaned.transform(points)
    np.testing.assert_array_equal(physical.min(axis=0), [100000, 200000, 40])
    np.testing.assert_array_equal(physical.max(axis=0), [100020, 200020, 70])
    assert cleaned.children[Object][0].get_geometry(role="location").x == -10
    path = tmp_path / "placed.dtcc"
    cleaned.save(path)
    restored = io.load_city(path)
    assert restored.buildings[0].bounds.xmin == 100000
    assert restored.transform.offset[2] == 10
    assert _volume(restored.buildings) == pytest.approx(3264, abs=1e-8, rel=0)


@pytest.mark.parametrize("case", ["mixed", "scale", "tilt", "degrees", "feet"])
def test_unsupported_frames_fail_explicitly(case):
    source = _source()
    if case == "mixed":
        source.lod1.surfaces[0].transform.set_translation(1, 0, 0)
    elif case == "scale":
        source.transform.affine[0, 0] = 2
    elif case == "tilt":
        source.transform.affine[:3, :3] = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    else:
        source.transform.srs = "EPSG:4326" if case == "degrees" else "EPSG:2263"
    with pytest.raises(ValueError, match="unsupported_frame"):
        footprints_from_lod1([source])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tolerance": 0},
        {"wall_angle_tolerance": np.nan},
        {"on_error": "skip"},
        {"overwrite_lod0": "yes"},
    ],
)
def test_parameters_are_checked_before_sources(kwargs):
    with pytest.raises((ValueError, TypeError)):
        footprints_from_lod1(None, **kwargs)


def test_duplicate_ids_fail_before_geometry_processing_even_in_warn_mode():
    with pytest.raises(ValueError, match="Duplicate"):
        footprints_from_lod1(
            [Building(id="same"), Building(id="same")], on_error="warn"
        )
    assert footprints_from_lod1([]) == []
    assert footprints_from_lod1(City()).buildings == []


@pytest.mark.parametrize("detail", [0.1, 0.5])
@pytest.mark.parametrize("rebuild", [False, True])
def test_shared_source_cleaning_retains_distinct_parts_roofs_walls_and_volume(
    detail, rebuild
):
    source = _source()
    source.estimated_height = 30
    pieces = footprints_from_lod1(_city(source))
    cleaned = clean_footprints(pieces, min_building_detail=detail)
    assert len(cleaned.buildings) == 1
    block = cleaned.buildings[0]
    assert block.attributes["merged_from"] == [source.id]
    parts = block.building_parts
    assert len(parts) == len({part.id for part in parts}) == 2
    assert sorted(p.lod0.zmax for p in parts) == [6, 30]
    assert walls_of(block)["version"] == WALLS_VERSION
    build_lod1_buildings(cleaned.buildings, rebuild=rebuild)
    assert sorted(part.lod1.bounds.zmax for part in parts) == [6, 30]
    assert _volume(cleaned.buildings) == pytest.approx(3264, abs=1e-8, rel=0)


def test_native_roundtrip_preserves_adapter_report_provenance_and_geometry(tmp_path):
    converted = footprints_from_lod1(_city(_source()))
    path = tmp_path / "pieces.dtcc"
    io.save_model(converted, path)
    restored = io.load_model(path)
    assert restored.attributes == converted.attributes
    assert _signature(restored.buildings) == _signature(converted.buildings)
    build_lod1_buildings(restored.buildings, rebuild=False)
    assert _volume(restored.buildings) == pytest.approx(3264, abs=1e-8, rel=0)


def test_zero_detail_identity_and_fallback_roof_path_preserve_volume():
    pieces = footprints_from_lod1([_source()])
    for piece in pieces:
        piece.attributes.pop("measured_height")
    cleaned = clean_footprints(pieces, min_building_detail=0)
    assert len(cleaned) == 2
    build_lod1_buildings(cleaned, rebuild=False)
    assert sorted(b.lod0.zmax for b in cleaned) == [6, 30]
    assert _volume(cleaned) == pytest.approx(3264, abs=1e-8, rel=0)
    merged = clean_footprints(pieces, min_building_detail=0.1)
    assert len(merged) == 1
    assert sorted(p.lod0.zmax for p in merged[0].building_parts) == [6, 30]
    build_lod1_buildings(merged, rebuild=False)
    assert _volume(merged) == pytest.approx(3264, abs=1e-8, rel=0)


def test_city_frame_failures_follow_warn_policy_without_changing_sources():
    city = _city(_source())
    city.transform.affine[0, 0] = 2
    result = footprints_from_lod1(city, on_error="warn")
    assert (
        result.attributes["footprints_from_lod1"]["failures"][0]["reason"]
        == "unsupported_frame"
    )
    assert result.buildings[0].id == city.buildings[0].id
    assert result.buildings[0].lod0 is None
