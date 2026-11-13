"""Native footprints and all labelled building representations through CityJSON."""

import json
from copy import deepcopy

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import io
from dtcc_core.builder import building_from_footprint, footprints_from_lod1
from dtcc_core.io.cityjson.cityjson import load
from dtcc_core.io.cityjson.write_cityjson import to_cityjson
from dtcc_core.model import Building, BuildingPart, City, Field, GeometryType, Solid
from tests.builder.test_footprints_from_lod1 import _city, _source, _volume
from tests.builder.test_prism_decomposition import BASE, TOWER, _input, _prism


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
