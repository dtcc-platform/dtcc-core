"""Native footprints and all labelled building representations through CityJSON."""

import dataclasses
import json
import logging
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import io
from dtcc_core.builder import building_from_footprint, footprints_from_lod1
from dtcc_core.io.cityjson.cityjson import load
from dtcc_core.io.cityjson.converters import CityJSONConfig
from dtcc_core.io.cityjson.write_cityjson import to_cityjson
from dtcc_core.model import Building, BuildingPart, City, Field, GeometryType, Solid
from tests.builder.test_footprints_from_lod1 import _city, _source, _volume
from tests.builder.test_prism_decomposition import BASE, TOWER, _input, _prism

# Per-surface WallSurface labels, which the writer adds only on request.
_WALL_LABELS = CityJSONConfig(semantic_types={"surface": "WallSurface"})


@pytest.mark.parametrize("strict", [False, True])
def test_all_lods_nested_parts_ids_attributes_and_crs_survive(strict, tmp_path):
    root = building_from_footprint(
        box(100000.1234, 200000.2345, 100020.3456, 200020.4567), 10, id="root"
    )
    nested = BuildingPart(id="nested", attributes={"name": "nested roof"})
    child = BuildingPart(id="child", attributes={"custom": [1, 2]})
    for part, height in ((child, 6), (nested, 3)):
        b = building_from_footprint(box(100005, 200005, 100010, 200010), height)
        part.add_geometry(b.lod0, GeometryType.LOD0)
        part.add_geometry(b.lod1, GeometryType.LOD1)
    child.add_child(nested)
    root.add_child(child)
    city = _city(root)
    city.transform.srs = "EPSG:3006"
    path = tmp_path / "tree.city.json"
    city.save(path, strict=strict)
    document = json.loads(path.read_text())
    assert all(
        isinstance(g["lod"], str)
        for obj in document["CityObjects"].values()
        for g in obj["geometry"]
    )
    restored = io.load_city(path, strict=strict)
    assert restored.transform.srs == "EPSG:3006"
    r = restored.buildings[0]
    assert r.id == "root" and r.lod0 is not None and r.lod1 is not None
    np.testing.assert_allclose(
        r.lod0.surfaces[0].vertices, root.lod0.vertices, atol=0.001, rtol=0
    )
    assert r.building_parts[0].id == "child"
    assert r.building_parts[0].attributes == child.attributes
    assert r.building_parts[0].building_parts[0].id == "nested"
    assert r.building_parts[0].building_parts[0].attributes == nested.attributes
    for part in [r, r.building_parts[0], r.building_parts[0].building_parts[0]]:
        assert part.lod0 is not None and part.lod1 is not None


def test_permissive_import_preserves_fractional_and_ambiguous_lod_labels():
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometry = document["CityObjects"]["b"]["geometry"]
    geometry[:] = [deepcopy(geometry[1]), deepcopy(geometry[1])]
    for item in geometry:
        item["lod"] = "1.2"
    restored = load(document).buildings[0]
    assert restored.lod1 is None
    assert len(restored.get_geometries(lod="1.2")) == 2
    with pytest.raises(ValueError, match="Ambiguous"):
        restored.get_geometry(lod="1.2")
    with pytest.raises(ValueError, match="missing_geometry"):
        footprints_from_lod1([restored])


@pytest.mark.parametrize("strict", [False, True])
def test_cavity_shells_survive_load_save_and_are_rejected_by_adapter(strict, tmp_path):
    outer = _input(_prism(BASE), "surfaces").surfaces
    inner = _input(_prism(TOWER, 2, 8), "surfaces").surfaces
    solid = Solid(
        surfaces=outer + inner,
        shells=[np.arange(len(outer)), np.arange(len(outer), len(outer) + len(inner))],
    )
    building = Building(id="cavity")
    building.add_geometry(solid, GeometryType.LOD1)
    city = _city(building)
    path = tmp_path / "cavity.city.json"
    city.save(path, strict=True)
    restored = io.load_city(path, strict=strict)
    assert isinstance(restored.buildings[0].lod1, Solid)
    assert len(restored.buildings[0].lod1.shells) == 2
    with pytest.raises(ValueError, match="not_2_5d"):
        footprints_from_lod1(restored)
    restored.save(path, strict=strict)
    assert (
        len(
            json.loads(path.read_text())["CityObjects"]["cavity"]["geometry"][0][
                "boundaries"
            ]
        )
        == 2
    )


@pytest.mark.parametrize("case", ["id", "field", "normal"])
def test_strict_native_surface_support_keeps_unmapped_metadata_rejected(case, tmp_path):
    building = building_from_footprint(BASE, 10)
    if case == "id":
        building.geometry["survey-roof"] = building.geometry.pop("GeometryType.LOD0")
    elif case == "field":
        building.lod0.fields = [
            Field(name="temperature", values=np.array([1.0]), association="geometry")
        ]
    else:
        building.lod0.calculate_normal()
    path = tmp_path / "keep.city.json"
    path.write_text("previous")
    with pytest.raises(NotImplementedError):
        _city(building).save(path, strict=True)
    assert path.read_text() == "previous"


def test_native_model_retains_city_conversion_report_and_cityjson_refuses_to_drop_it(
    tmp_path,
):
    city = footprints_from_lod1(_city(_source()))
    native = tmp_path / "report.dtcc"
    city.save(native)
    assert io.load_city(native).attributes == city.attributes
    with pytest.raises(NotImplementedError, match="City root attributes"):
        city.save(tmp_path / "report.city.json", strict=True)
    geometry_city = City()
    geometry_city.add_buildings(city.buildings)
    path = tmp_path / "geometry.city.json"
    geometry_city.save(path, strict=True)
    restored = io.load_city(path, strict=True)
    assert [b.id for b in restored.buildings] == [b.id for b in city.buildings]
    assert all(b.lod0 is not None and b.lod1 is not None for b in restored.buildings)
    assert _volume(restored.buildings) == pytest.approx(3264, abs=1e-7, rel=0)


def test_invalid_later_representation_fails_instead_of_omitting_the_building():
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometries = document["CityObjects"]["b"]["geometry"]
    geometries[-1]["boundaries"][0][0][0] = len(document["vertices"])
    with pytest.raises(ValueError, match="valid vertex indices"):
        load(document)


def test_previous_writer_footprint_loads_only_permissively(tmp_path):
    path = Path(__file__).parent / "fixtures" / "legacy-footprint.city.json"
    document = json.loads(path.read_text())
    before = deepcopy(document)
    city = load(document)
    assert document == before
    assert len(city.buildings[0].lod0.surfaces) == 1
    assert len(city.buildings[0].lod1.surfaces) == 6
    target = tmp_path / "legacy.city.json"
    city.save(target)
    assert load(target).buildings[0].lod0 is not None
    strict_document = deepcopy(document)
    for geometry in strict_document["CityObjects"]["building"]["geometry"]:
        geometry["lod"] = format(geometry["lod"], "g")
    with pytest.raises(ValueError, match="semantic index"):
        load(strict_document, strict=True)


@pytest.mark.parametrize("kind", ["child", "semantics", "geometry", "missing_child"])
def test_permissive_unsupported_content_warns_and_keeps_supported_geometry(
    kind, caplog
):
    document = to_cityjson(
        _city(building_from_footprint(BASE, 10, id="b")), config=_WALL_LABELS
    )
    building = document["CityObjects"]["b"]
    if kind == "child":
        building["children"] = ["installation"]
        document["CityObjects"]["installation"] = {
            "type": "BuildingInstallation",
            "parents": ["b"],
        }
        message = "unsupported child 'installation'"
    elif kind == "missing_child":
        building["children"] = ["missing"]
        message = "child 'missing' not found"
    elif kind == "semantics":
        building["geometry"][1]["semantics"]["surfaces"][0]["type"] = "+Foo"
        message = "ignoring unsupported semantics"
    else:
        extra = deepcopy(building["geometry"][1])
        extra.pop("semantics")
        extra.update(type="CompositeSolid", lod="2", boundaries=[extra["boundaries"]])
        building["geometry"].append(extra)
        message = "skipping unsupported geometry"
    before = deepcopy(document)
    restored = load(document).buildings[0]
    assert document == before
    assert restored.id == "b"
    assert len(restored.lod1.surfaces) == 6
    assert message in caplog.text and "Building 'b'" in caplog.text
    if kind == "semantics":
        assert restored.lod1.regions == []
    with pytest.raises((ValueError, NotImplementedError)):
        load(document, strict=True)


@pytest.mark.parametrize("kind", ["index", "semantics", "legacy_shape"])
def test_permissive_compatibility_does_not_hide_malformed_supported_content(kind):
    document = to_cityjson(
        _city(building_from_footprint(BASE, 10, id="b")), config=_WALL_LABELS
    )
    geometry = document["CityObjects"]["b"]["geometry"][1]
    if kind == "index":
        geometry["semantics"]["surfaces"][0]["type"] = "+Foo"
        geometry["boundaries"][0][0][0] = len(document["vertices"])
    elif kind == "semantics":
        geometry["semantics"]["values"][0] = 999
    else:
        geometry["semantics"]["values"] = [[0]]
    with pytest.raises(ValueError):
        load(document)


@pytest.mark.parametrize("label", ["1.2", "1.0", "2.2"])
@pytest.mark.parametrize("strict", [False, True])
def test_exact_lod_labels_and_ambiguity_survive_default_and_strict_saving(
    label, strict, tmp_path
):
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometries = document["CityObjects"]["b"]["geometry"]
    geometries[:] = [deepcopy(geometries[1]), deepcopy(geometries[1])]
    for geometry in geometries:
        geometry["lod"] = label
    city = load(document)
    path = tmp_path / "fractional.city.json"
    city.save(path, strict=strict)
    encoded = json.loads(path.read_text())["CityObjects"]["b"]["geometry"]
    assert [geometry["lod"] for geometry in encoded] == [label, label]
    restored = load(path, strict=strict).buildings[0]
    assert len(restored.get_geometries(lod=label)) == 2
    assert restored.lod1 is None
    with pytest.raises(ValueError, match="Ambiguous"):
        restored.get_geometry(lod=label)


@pytest.mark.parametrize("import_strict", [False, True])
@pytest.mark.parametrize("export_strict", [False, True])
def test_composite_surface_roundtrip_and_closed_lod1_conversion(
    import_strict, export_strict, tmp_path
):
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometries = document["CityObjects"]["b"]["geometry"]
    geometries[:] = [geometries[1]]
    geometries[0]["type"] = "CompositeSurface"
    city = load(document, strict=import_strict)
    converted = footprints_from_lod1(city)
    assert _volume(converted.buildings) == pytest.approx(4000)
    path = tmp_path / "composite.city.json"
    city.save(path, strict=export_strict)
    encoded = json.loads(path.read_text())["CityObjects"]["b"]["geometry"][0]
    assert encoded["type"] == "CompositeSurface"
    assert encoded["lod"] == "1"
    restored = load(path, strict=import_strict)
    assert _volume(footprints_from_lod1(restored).buildings) == pytest.approx(4000)


def _label_roof_ground_walls(document, geometry):
    """Give a written box LoD1 the Roof/Ground/Wall split a source file has."""
    z = [vertex[2] for vertex in document["vertices"]]
    top, bottom = max(z), min(z)
    names = ["RoofSurface", "GroundSurface", "WallSurface"]
    values = []
    for surface in geometry["boundaries"]:
        levels = {z[index] for index in surface[0]}
        values.append(0 if levels == {top} else 1 if levels == {bottom} else 2)
    assert sorted(set(values)) == [0, 1, 2]
    geometry["semantics"] = {
        "surfaces": [{"type": name} for name in names],
        "values": values,
    }


def test_imported_standard_semantics_still_require_explicit_conversion():
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometries = document["CityObjects"]["b"]["geometry"]
    geometries[:] = [geometries[1]]
    _label_roof_ground_walls(document, geometries[0])
    city = load(document)
    assert [
        region.semantic_type.split("#")[-1] for region in city.buildings[0].lod1.regions
    ] == [
        "RoofSurface",
        "GroundSurface",
        "WallSurface",
    ]
    with pytest.raises(ValueError, match="unsupported_metadata"):
        footprints_from_lod1(city)


@pytest.mark.parametrize("strict", [False, True])
def test_unclassified_footprint_is_not_labelled_as_a_wall(strict, tmp_path):
    building = building_from_footprint(BASE, 10, ground_height=2, id="b")
    city = _city(building)
    path = tmp_path / "footprint.city.json"
    city.save(path, strict=strict)
    restored = load(path, strict=strict).buildings[0]
    assert restored.lod0.regions == []
    assert restored.lod0.zmax == 12


@pytest.mark.parametrize("strict", [False, True])
def test_unclassified_lod1_is_not_labelled_as_walls(strict, tmp_path):
    city = _city(building_from_footprint(BASE, 10, id="b"))
    path = tmp_path / "box.city.json"
    city.save(path, strict=strict)
    written = json.loads(path.read_text())["CityObjects"]["b"]["geometry"]
    assert all(not g.get("semantics", {}).get("surfaces") for g in written)
    restored = io.load_city(path)
    assert restored.buildings[0].lod1.regions == []
    report = footprints_from_lod1(restored, overwrite_lod0=True).attributes[
        "footprints_from_lod1"
    ]
    assert report["converted"] == 1


def test_configured_surface_type_still_labels_unclassified_lod1():
    document = to_cityjson(
        _city(building_from_footprint(BASE, 10, id="b")), config=_WALL_LABELS
    )
    lod0, lod1 = document["CityObjects"]["b"]["geometry"]
    assert "semantics" not in lod0
    assert {s["type"] for s in lod1["semantics"]["surfaces"]} == {"WallSurface"}


def test_previous_writer_wall_labels_load_without_regions():
    path = Path(__file__).parent / "fixtures" / "legacy-footprint.city.json"
    city = load(json.loads(path.read_text()))
    building = city.buildings[0]
    assert building.lod0.regions == []
    assert building.lod1.regions == []
    report = footprints_from_lod1(city, overwrite_lod0=True).attributes[
        "footprints_from_lod1"
    ]
    assert report["converted"] == 1


def test_previous_writer_solid_wall_labels_load_without_regions():
    building = building_from_footprint(BASE, 10, id="b")
    lod1 = building.lod1
    solid = Solid(surfaces=lod1.surfaces, shells=[np.arange(len(lod1.surfaces))])
    building.add_geometry(solid, GeometryType.LOD2)
    document = to_cityjson(_city(building), config=_WALL_LABELS)
    written = document["CityObjects"]["b"]["geometry"][-1]
    assert written["type"] == "Solid"
    assert written["semantics"]["values"] == [list(range(len(lod1.surfaces)))]
    assert load(document).buildings[0].lod2.regions == []


@pytest.mark.parametrize(
    "semantics",
    [
        # One shared entity, as source files usually write it.
        lambda n: {"surfaces": [{"type": "WallSurface"}], "values": [0] * n},
        # Per-surface entities that carry attributes.
        lambda n: {
            "surfaces": [{"type": "WallSurface", "material": "brick"}] * n,
            "values": list(range(n)),
        },
        # Per-surface entities with one surface left unlabelled.
        lambda n: {
            "surfaces": [{"type": "WallSurface"}] * (n - 1),
            "values": [*range(n - 1), None],
        },
    ],
)
def test_source_wall_labels_are_kept(semantics):
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometries = document["CityObjects"]["b"]["geometry"]
    geometries[1]["semantics"] = semantics(len(geometries[1]["boundaries"]))
    assert load(document).buildings[0].lod1.regions


def _blanket_labels(count):
    return {
        "surfaces": [{"type": "WallSurface"} for _ in range(count)],
        "values": list(range(count)),
    }


def test_previous_writer_wall_labels_are_dropped_with_a_debug_message(caplog):
    caplog.set_level(logging.DEBUG, logger="dtcc-core")
    document = to_cityjson(
        _city(building_from_footprint(BASE, 10, id="b")), config=_WALL_LABELS
    )
    assert load(document).buildings[0].lod1.regions == []
    assert "ignoring WallSurface labels on every surface" in caplog.text


@pytest.mark.parametrize("kind", ["too_few_assignments", "nested_multisurface"])
def test_wall_label_compatibility_does_not_hide_malformed_semantics(kind):
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometry = document["CityObjects"]["b"]["geometry"][1]
    count = len(geometry["boundaries"])
    if kind == "too_few_assignments":
        geometry["semantics"] = _blanket_labels(count - 3)
    else:
        semantics = _blanket_labels(count)
        semantics["values"] = [[value] for value in semantics["values"]]
        geometry["semantics"] = semantics
    with pytest.raises(ValueError):
        load(document)


def test_wall_only_geometry_keeps_per_surface_wall_labels():
    document = to_cityjson(_city(building_from_footprint(BASE, 10, id="b")))
    geometries = document["CityObjects"]["b"]["geometry"]
    geometry = geometries[1]
    z = [vertex[2] for vertex in document["vertices"]]
    geometry["boundaries"] = [
        polygon
        for polygon in geometry["boundaries"]
        if len({z[index] for index in polygon[0]}) > 1
    ]
    geometry["semantics"] = _blanket_labels(len(geometry["boundaries"]))
    geometries[:] = [geometry]
    regions = load(document).buildings[0].lod1.regions
    assert len(regions) == len(geometry["boundaries"]) == 4


def test_copied_default_config_does_not_label_unclassified_lod1():
    city = _city(building_from_footprint(BASE, 10, id="b"))
    copied = dataclasses.replace(CityJSONConfig(), rounding_mode="floor")
    lod1 = to_cityjson(city, config=copied)["CityObjects"]["b"]["geometry"][1]
    assert "semantics" not in lod1
    labelled = CityJSONConfig(labels_unclassified_surfaces=True)
    lod1 = to_cityjson(city, config=labelled)["CityObjects"]["b"]["geometry"][1]
    assert {s["type"] for s in lod1["semantics"]["surfaces"]} == {"WallSurface"}
