from collections import Counter, defaultdict
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import numpy as np
import pytest
from shapely import from_wkt
from shapely.geometry import Point, Polygon, box

from dtcc_core.builder import (
    _dtcc_builder,
    build_city_flat_mesh,
    build_city_surface_mesh,
    build_city_volume_mesh,
    fix_building_footprint_clearance,
    merge_building_footprints,
    simplify_building_footprints,
)
from dtcc_core.builder.building.modify import clean_building_footprints
from dtcc_core.builder.cleaning.contract import (
    check_mesher_handoff_profile,
)
from dtcc_core.builder.geometry_builders import meshes as meshes_module
from dtcc_core.builder.model_conversion import (
    builder_mesh_to_mesh,
    create_builder_polygon,
    create_builder_surface,
    mesh_to_builder_mesh,
)
from dtcc_core.builder.meshing import dtcc_mesher_backend as dtcc_mesher_backend_module
from dtcc_core.builder.meshing import flat_mesh_backends as flat_mesh_backends_module
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.model import Building, Bounds, City, GeometryType, Mesh, Raster, Surface, Terrain, VolumeMesh


def make_surface(polygon: Polygon, z: float) -> Surface:
    surface = Surface()
    surface.from_polygon(polygon, z)
    return surface


def make_conditioned_footprints(
    surfaces: list[Surface],
    source_map: list[list[int]],
    subdomain_resolution: list[float],
    diagnostics: dict[str, object],
    *,
    min_building_detail: float = 0.5,
    contract: dict[str, object] | None = None,
) -> meshes_module.ConditionedFootprints:
    declared_scale = meshes_module._conditioned_footprint_declared_scale(
        min_building_detail=min_building_detail,
        diagnostics=diagnostics,
    )
    return meshes_module.ConditionedFootprints(
        surfaces=surfaces,
        source_map=source_map,
        subdomain_resolution=subdomain_resolution,
        diagnostics=diagnostics,
        declared_scale=declared_scale,
        contract=contract
        or {
            "ok": True,
            "status": "pass",
            "requirements": {},
            "errors": [],
            "warnings": [],
            "metrics": {},
        },
    )


def unpack_conditioned_footprints(
    conditioned: meshes_module.ConditionedFootprints,
) -> tuple[list[Surface], list[list[int]], list[float], dict[str, object]]:
    return (
        conditioned.surfaces,
        conditioned.source_map,
        conditioned.subdomain_resolution,
        conditioned.diagnostics,
    )


def make_prepared_city_inputs(
    terrain: object,
    terrain_raster: object,
    surfaces: list[Surface],
    source_map: list[list[int]],
    subdomain_resolution: list[float],
    diagnostics: dict[str, object],
    *,
    min_building_detail: float = 0.5,
) -> tuple[object, object, meshes_module.ConditionedFootprints]:
    return (
        terrain,
        terrain_raster,
        make_conditioned_footprints(
            surfaces,
            source_map,
            subdomain_resolution,
            diagnostics,
            min_building_detail=min_building_detail,
        ),
    )


def make_prepared_ground_regions(
    building_surfaces: list[Surface],
    meshing_directives: list[int],
    region_polygons: list[Polygon],
    region_markers: list[int],
    region_triangle_sizes: dict[int, float],
    region_points: list[np.ndarray],
) -> meshes_module.PreparedGroundRegions:
    building_union = meshes_module.unary_union(
        [
            polygon
            for polygon, marker in zip(region_polygons, region_markers)
            if marker >= 0
        ]
    )
    partitioned_polygons: list[Polygon] = []
    partitioned_markers: list[int] = []
    partitioned_points: list[np.ndarray] = []
    for polygon, marker, point in zip(
        region_polygons, region_markers, region_points
    ):
        components = (
            meshes_module._iter_polygon_components(
                polygon.difference(building_union)
            )
            if marker < 0 and not building_union.is_empty
            else [polygon]
        )
        for component in components:
            partitioned_polygons.append(component)
            partitioned_markers.append(marker)
            representative = component.representative_point()
            partitioned_points.append(
                np.array([representative.x, representative.y], dtype=np.float64)
                if len(components) != 1 or not component.equals(polygon)
                else point
            )
    return meshes_module.PreparedGroundRegions(
        building_surfaces=building_surfaces,
        meshing_directives=meshing_directives,
        region_polygons=partitioned_polygons,
        region_markers=partitioned_markers,
        region_triangle_sizes=region_triangle_sizes,
        region_points=partitioned_points,
    )


def make_building(
    polygon: Polygon,
    *,
    roof_z: float = 10.0,
    lod1_polygon: Polygon | None = None,
    lod1_z: float | None = None,
    lod2_polygon: Polygon | None = None,
    lod2_z: float | None = None,
) -> Building:
    building = Building()
    building.add_geometry(make_surface(polygon, roof_z), GeometryType.LOD0)
    if lod1_polygon is not None:
        building.add_geometry(
            make_surface(lod1_polygon, lod1_z if lod1_z is not None else roof_z),
            GeometryType.LOD1,
        )
    if lod2_polygon is not None:
        building.add_geometry(
            make_surface(lod2_polygon, lod2_z if lod2_z is not None else roof_z),
            GeometryType.LOD2,
        )
    building.attributes["estimated_height"] = roof_z
    building.attributes["ground_height"] = 0.0
    return building


def _fake_tetgen_switch_defaults() -> dict:
    # Stand-in for dtcc_tetgen_wrapper.switches.tetgen_defaults() so tests that
    # mock _require_tetgen() also run where TetGen is not installed (CI).
    return {
        "plc": True,
        "preserve_surface": False,
        "quality": None,
        "max_volume": None,
        "max_added_points": None,
        "optimize_level": None,
        "extra": "",
        "refine": False,
    }


def make_flat_city(buildings: list[Building]) -> City:
    city = City()
    raster = Raster()
    raster.data = np.zeros((8, 8), dtype=float)
    raster.set_bounds(Bounds(0.0, 0.0, 80.0, 80.0))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(buildings)
    return city


def make_bounded_flat_city(
    buildings: list[Building], bounds: Bounds
) -> City:
    city = City()
    raster = Raster()
    raster.data = np.zeros((2, 2), dtype=float)
    raster.set_bounds(bounds)
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(buildings)
    return city


def _nearest_nonadjacent_boundary_vertex_distance(polygon: Polygon) -> float:
    exterior = np.asarray(polygon.exterior.coords, dtype=np.float64)
    best = np.inf

    for i in range(len(exterior) - 1):
        for j in range(i + 2, len(exterior) - 1):
            if i == 0 and j == len(exterior) - 2:
                continue
            best = min(best, float(np.linalg.norm(exterior[i] - exterior[j])))

    return best


def test_build_city_flat_mesh_handles_pathological_footprints(monkeypatch):
    calls = []
    original = meshes_module._condition_meshing_footprints

    def wrapped(*args, **kwargs):
        calls.append(kwargs.get("lod"))
        return original(*args, **kwargs)

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", wrapped)

    buildings = [
        make_building(Polygon([(8, 8), (18, 8), (18, 18), (8, 18), (8, 8)]), roof_z=12.0),
        make_building(Polygon([(18.1, 8), (28.1, 8), (28.1, 18), (18.1, 18), (18.1, 8)]), roof_z=8.0),
        make_building(Polygon([(35, 8), (42, 15), (35, 15), (42, 8), (35, 8)]), roof_z=6.0),
    ]
    city = make_flat_city(buildings)

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.3,
        min_building_area=1.0,
        merge_tolerance=0.25,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert calls
    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_build_city_flat_mesh_forwards_cleaning_diagnostics_flag(monkeypatch):
    captured = {}
    original = meshes_module._condition_meshing_footprints

    def wrapped(*args, **kwargs):
        captured["cleaning_diagnostics"] = kwargs.get("cleaning_diagnostics")
        return original(*args, **kwargs)

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", wrapped)

    city = make_flat_city(
        [make_building(box(8, 8, 18, 18), roof_z=10.0)]
    )

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.3,
        min_building_area=1.0,
        merge_tolerance=0.25,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        cleaning_diagnostics=False,
        mesher="dtcc_mesher",
    )

    assert captured["cleaning_diagnostics"] is False
    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_meshing_adapter_propagates_numeric_merge_policy(monkeypatch):
    captured = []
    polygon = box(8.0, 8.0, 18.0, 18.0)

    def fake_condition(polygons, *, source_map, options):
        assert polygons == [polygon]
        assert source_map == [[0]]
        captured.append(options)
        return meshes_module.ConditioningResult(
            polygons=[polygon],
            source_map=[[0]],
            diagnostics={"output_grid": 0.03125},
        )

    monkeypatch.setattr(meshes_module, "condition_polygon_coverage", fake_condition)
    building = make_building(polygon, roof_z=10.0)
    for merge_buildings in (True, False):
        meshes_module._condition_meshing_footprints(
            [building],
            lod=GeometryType.LOD0,
            min_building_detail=0.5,
            min_building_area=0.0,
            merge_tolerance=0.375,
            merge_buildings=merge_buildings,
            max_mesh_size=5.0,
            cleaning_diagnostics=False,
        )

    assert [options.merge_distance for options in captured] == [0.375, 0.0]
    assert [options.allow_source_merging for options in captured] == [True, False]


def test_build_city_flat_mesh_rejects_complete_domain_handoff_acute_sector():
    shallow_boundary = Polygon(
        [(0.0, -0.01), (20.0, 0.01), (20.0, 10.0), (0.0, 9.98)]
    )
    city = make_bounded_flat_city(
        [make_building(shallow_boundary, roof_z=5.0)],
        Bounds(0.0, 0.0, 20.0, 20.0),
    )

    with pytest.raises(
        meshes_module.MesherHandoffError,
        match="Complete 2D mesher handoff.*incident sector",
    ) as error:
        build_city_flat_mesh(
            city,
            lod=GeometryType.LOD0,
            merge_buildings=True,
            min_building_detail=0.5,
            min_building_area=0.0,
            merge_tolerance=0.5,
            max_mesh_size=2.0,
            min_mesh_angle=25.0,
            report_mesh_quality=False,
            cleaning_diagnostics=False,
            mesher="dtcc_mesher",
        )

    contract = error.value.contract
    witness = contract["metrics"]["mesher_profile"]["incident_sectors"][
        "minimum_witness"
    ]
    assert contract["requirements"]["incident_sector_profile"] is False
    assert witness["vertex"] == [10.0, 0.0]
    assert witness["angle_degrees"] == pytest.approx(0.0572957604)
    assert witness["region_markers"] == [-2]
    assert witness["region_roles"] == ["ground"]
    assert witness["on_domain_boundary"] is True


def test_build_city_flat_mesh_allows_benign_boundary_touch():
    city = make_bounded_flat_city(
        [make_building(box(0.0, 5.0, 5.0, 10.0), roof_z=5.0)],
        Bounds(0.0, 0.0, 20.0, 20.0),
    )

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=0.0,
        merge_tolerance=0.5,
        max_mesh_size=2.0,
        min_mesh_angle=25.0,
        report_mesh_quality=False,
        cleaning_diagnostics=False,
        mesher="dtcc_mesher",
    )

    assert mesh.faces.shape[0] > 0


def test_build_city_surface_mesh_rejects_complete_domain_handoff_acute_sector():
    shallow_boundary = Polygon(
        [(0.0, -0.01), (20.0, 0.01), (20.0, 10.0), (0.0, 9.98)]
    )
    city = make_bounded_flat_city(
        [make_building(shallow_boundary, roof_z=5.0)],
        Bounds(0.0, 0.0, 20.0, 20.0),
    )

    with pytest.raises(meshes_module.MesherHandoffError) as error:
        build_city_surface_mesh(
            city,
            lod=GeometryType.LOD0,
            merge_buildings=True,
            min_building_detail=0.5,
            min_building_area=0.0,
            merge_tolerance=0.5,
            building_mesh_triangle_size=2.0,
            max_mesh_size=2.0,
            min_mesh_angle=25.0,
            report_mesh_quality=False,
            cleaning_diagnostics=False,
            mesher="dtcc_mesher",
        )

    contract = error.value.contract
    assert contract["requirements"]["incident_sector_profile"] is False
    assert contract["metrics"]["mesher_profile"]["incident_sectors"][
        "minimum_witness"
    ]["vertex"] == [10.0, 0.0]


def test_condition_meshing_footprints_regularizes_touching_holes():
    touching_holes = Polygon(
        [(8, 8), (28, 8), (28, 28), (8, 28), (8, 8)],
        [
            [(12, 12), (16, 12), (16, 16), (12, 16), (12, 12)],
            [(16, 16), (22, 16), (22, 24), (16, 24), (16, 16)],
        ],
    )
    assert touching_holes.is_valid

    conditioned = meshes_module._condition_meshing_footprints(
        [make_building(touching_holes, roof_z=10.0)],
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=False,
        max_mesh_size=5.0,
        cleaning_diagnostics=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0]]
    assert resolutions == [5.0]
    assert conditioned.contract["errors"] == []
    normalized = surfaces[0].to_polygon(simplify=0.0)
    assert check_mesher_handoff_profile([normalized])["status"] == "pass"
    assert diagnostics.get("mesher_regularized_polygon_count", 0) == 0


def test_condition_meshing_footprints_shows_plot_when_requested(monkeypatch):
    captured: dict[str, object] = {}

    def fake_plot(raw_polygons, cleaned_polygons, *, title, show, block):
        captured["raw_count"] = len(raw_polygons)
        captured["cleaned_count"] = len(cleaned_polygons)
        captured["title"] = title
        captured["show"] = show
        captured["block"] = block
        return None, None

    monkeypatch.setattr(meshes_module, "plot_footprint_cleaning_comparison", fake_plot)

    conditioned = meshes_module._condition_meshing_footprints(
        [make_building(box(0.0, 0.0, 10.0, 10.0), roof_z=10.0)],
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=False,
        max_mesh_size=5.0,
        cleaning_diagnostics=False,
        show_footprints=True,
        footprint_cleaning_plot_block=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0]]
    assert resolutions == [5.0]
    assert diagnostics["output_count"] == 1
    assert captured == {
        "raw_count": 1,
        "cleaned_count": 1,
        "title": "Footprint cleaning",
        "show": True,
        "block": False,
    }


def test_condition_meshing_footprints_repairs_lund_style_self_clearance_slit():
    polygon = from_wkt(
        "POLYGON ((388232.25 6176876.5, 388232.03125 6176867.5625, "
        "388220.84375 6176867.84375, 388220.75 6176863.25, "
        "388203.78125 6176863.625, 388203.6875 6176859.09375, "
        "388214.84375 6176858.84375, 388214.6875 6176851.96875, "
        "388220.53125 6176851.8125, 388220.78125 6176863.25, "
        "388236.4375 6176862.875, 388236.59375 6176869.46875, "
        "388237.90625 6176869.4375, 388238.15625 6176880.875, "
        "388221.125 6176881.21875, 388221.0625 6176876.71875, "
        "388232.25 6176876.5))"
    )

    conditioned = meshes_module._condition_meshing_footprints(
        [make_building(polygon, roof_z=10.0)],
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=False,
        max_mesh_size=5.0,
        cleaning_diagnostics=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0]]
    assert resolutions == [5.0]
    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=surfaces,
        declared_scale=max(
            0.5,
            float(diagnostics.get("output_grid", 0.0) or 0.0),
            1.0e-9,
        ),
        diagnostics=diagnostics,
    )
    assert contract["errors"] == []
    assert diagnostics["mesher_profile"]["status"] == "pass"


def test_condition_meshing_footprints_repairs_gbg_same_hole_neck():
    polygon = from_wkt(
        "POLYGON ((318447.875 6398710.4375, 318385.3125 6398700.40625, "
        "318386.28125 6398691.125, 318383.71875 6398690.84375, "
        "318384.84375 6398680.6875, 318387.3125 6398680.96875, "
        "318390.1875 6398639.9375, 318389.46875 6398581.03125, "
        "318466.71875 6398575.4375, 318469.71875 6398578.53125, "
        "318447.875 6398710.4375), "
        "(318436 6398600.5625, 318439.125 6398601.03125, "
        "318439.5 6398598.4375, 318445.84375 6398599.40625, "
        "318444.53125 6398607.875, 318450.34375 6398608.65625, "
        "318452.21875 6398590.4375, 318405.4375 6398593.53125, "
        "318405.0625 6398657.09375, 318400.8125 6398690.78125, "
        "318421.75 6398689.9375, 318435.78125 6398693.1875, "
        "318448.5625 6398608.875, 318444.96875 6398608.78125, "
        "318444.03125 6398611.3125, 318434.59375 6398609.90625, "
        "318436 6398600.5625))"
    )

    conditioned = meshes_module._condition_meshing_footprints(
        [make_building(polygon, roof_z=10.0)],
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=False,
        max_mesh_size=5.0,
        cleaning_diagnostics=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0]]
    assert resolutions == [5.0]
    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=surfaces,
        declared_scale=max(
            0.5,
            float(diagnostics.get("output_grid", 0.0) or 0.0),
            1.0e-9,
        ),
        diagnostics=diagnostics,
    )
    assert contract["errors"] == []
    assert diagnostics["mesher_profile"]["status"] == "pass"


def test_condition_meshing_footprints_repairs_stockholm_cross_ring_slit():
    polygon = from_wkt(
        "POLYGON ((673551.09375 6581800.15625, 673548.3125 6581798.59375, "
        "673551.46875 6581793.21875, 673538.5 6581785.875, "
        "673532.03125 6581796.9375, 673539.59375 6581801.1875, "
        "673537.34375 6581805.09375, 673541.65625 6581807.5, "
        "673538.28125 6581813.5, 673532.5 6581810.1875, "
        "673529.71875 6581815, 673523.46875 6581811.40625, "
        "673523.84375 6581807.59375, 673515.1875 6581800.25, "
        "673538.59375 6581759.15625, 673535.09375 6581756.28125, "
        "673537.25 6581752.46875, 673539.9375 6581754, "
        "673539.5625 6581758.84375, 673550.0625 6581764.78125, "
        "673556.59375 6581753.09375, 673529.375 6581738.15625, "
        "673523.0625 6581749.5, 673511.6875 6581743.03125, "
        "673524.5625 6581720.25, 673590.78125 6581756.96875, "
        "673535.0625 6581855.0625, 673510.40625 6581841.3125, "
        "673516.6875 6581829.9375, 673529.84375 6581837.25, "
        "673551.09375 6581800.15625), "
        "(673549.875 6581765.09375, 673541.375 6581780.5625, "
        "673543.21875 6581783.3125, 673547 6581779.96875, "
        "673554.15625 6581784.03125, 673553.5625 6581789.21875, "
        "673556.5625 6581788.75, 673572.03125 6581761.5625, "
        "673560.96875 6581755.5, 673554.34375 6581767.1875, "
        "673549.875 6581765.09375))"
    )

    conditioned = meshes_module._condition_meshing_footprints(
        [make_building(polygon, roof_z=10.0)],
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=False,
        max_mesh_size=5.0,
        cleaning_diagnostics=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0]]
    assert resolutions == [5.0]
    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=surfaces,
        declared_scale=max(
            0.5,
            float(diagnostics.get("output_grid", 0.0) or 0.0),
            1.0e-9,
        ),
        diagnostics=diagnostics,
    )
    assert contract["errors"] == []
    assert diagnostics["mesher_profile"]["status"] == "pass"


def test_conditioned_footprint_contract_audit_passes_scale_clean_output():
    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(box(0.0, 0.0, 10.0, 10.0), 8.0)],
        declared_scale=0.5,
        diagnostics={"geos_exception_count": 0},
    )

    assert contract["status"] == "pass"
    assert contract["requirements"]["scale_contract_satisfied"] is True
    assert contract["requirements"]["no_short_edges"] is True
    assert contract["metrics"]["pair_issue_count"] == 0
    assert contract["requirements"]["mesher_segment_graph_valid"] is True


def test_conditioned_footprint_contract_audit_warns_scale_contract():
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 3.9),
            (1.9, 3.9),
            (1.9, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]
    )

    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(polygon, 8.0)],
        declared_scale=0.5,
        diagnostics={"before_selection_contract": {"fidelity": {"status": "pass"}}},
    )

    assert contract["status"] == "warn"
    assert contract["ok"] is True
    assert contract["requirements"]["scale_contract_satisfied"] is False
    assert contract["requirements"]["no_short_edges"] is True
    assert contract["metrics"]["short_edge_obligation"] == (
        "not_applicable_to_incident_edges"
    )
    assert any("CONTRACT NOT SATISFIED" in message for message in contract["warnings"])


def test_residual_separation_without_fidelity_evidence_is_rejected():
    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(box(0, 0, 5, 5), 8), make_surface(box(5.1, 0, 10, 5), 8)],
        declared_scale=0.5, diagnostics={},
    )
    assert contract["ok"] is False
    assert any("fidelity report" in message for message in contract["errors"])


def test_conditioned_footprint_contract_audit_warns_residual_self_clearance():
    polygon = Polygon(
        shell=[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ],
        holes=[
            [
                (0.45, 2.0),
                (2.0, 2.0),
                (2.0, 8.0),
                (0.45, 8.0),
                (0.45, 2.0),
            ]
        ],
    )

    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(polygon, 8.0)],
        declared_scale=0.5,
        diagnostics={
            "output_grid": 0.03125,
            "before_selection_contract": {"fidelity": {"status": "pass"}},
        },
    )

    assert contract["ok"] is True
    assert contract["status"] == "warn"
    assert contract["requirements"]["scale_contract_satisfied"] is False
    assert contract["requirements"]["min_clearance_respected"] is False
    assert contract["metrics"]["residual_min_clearance_tolerated"] is True
    assert contract["metrics"]["declared_scale"] == pytest.approx(0.5)
    assert contract["metrics"]["min_clearance"] == pytest.approx(0.45)
    assert any("CONTRACT NOT SATISFIED" in message for message in contract["warnings"])


def test_conditioned_footprint_contract_ignores_retired_revalidation_grid():
    polygon = Polygon(
        shell=[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ],
        holes=[
            [
                (3.0, 1.9),
                (7.0, 1.9),
                (7.0, 6.0),
                (3.0, 6.0),
                (3.0, 1.9),
            ]
        ],
    )

    without_revalidation_grid = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(polygon, 8.0)],
        declared_scale=2.0,
        diagnostics={
            "before_selection_contract": {"fidelity": {"status": "pass"}},
            "geos_exception_count": 0,
            "output_grid": 0.03125,
            "precision_grid": 0.03125,
        },
    )
    with_revalidation_grid = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(polygon, 8.0)],
        declared_scale=2.0,
        diagnostics={
            "before_selection_contract": {"fidelity": {"status": "pass"}},
            "geos_exception_count": 0,
            "output_grid": 0.03125,
            "precision_grid": 0.03125,
            "mesher_ready_coverage_revalidation_applied": True,
            "mesher_ready_coverage_revalidation_output_grid": 0.125,
        },
    )

    assert without_revalidation_grid["status"] == "warn"
    assert without_revalidation_grid["metrics"]["clearance_deficit"] == pytest.approx(
        0.1
    )
    assert (
        without_revalidation_grid["metrics"]["residual_min_clearance_tolerated"]
        is True
    )
    assert with_revalidation_grid["status"] == "warn"
    assert with_revalidation_grid["metrics"]["contract_tolerance"] == pytest.approx(
        2.0e-6
    )


def test_conditioned_footprint_contract_audit_reports_invalid_mesher_segment_graph(
    monkeypatch,
):
    graph = object()
    monkeypatch.setattr(
        meshes_module,
        "dtcc_mesher",
        SimpleNamespace(
            Coverage=lambda polygons, markers: SimpleNamespace(graph=lambda: graph),
            validate_coverage_graph=lambda candidate_graph: (_ for _ in ()).throw(
                RuntimeError("intersecting segments in segment graph")
            ),
        ),
    )

    contract = meshes_module._conditioned_footprint_contract_audit(
        surfaces=[make_surface(box(0.0, 0.0, 10.0, 10.0), 8.0)],
        declared_scale=0.5,
        diagnostics={"geos_exception_count": 0},
    )

    assert contract["status"] == "fail"
    assert contract["requirements"]["mesher_segment_graph_valid"] is False
    assert "segment graph" in contract["metrics"]["mesher_segment_graph_error"]


def test_triangle_mesh_contract_from_audit_warns_on_subscale_edge_tail():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
                [1.0e-5, 1.0e-5, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 4],
                [1, 2, 4],
                [2, 3, 4],
                [3, 0, 4],
            ],
            dtype=int,
        ),
        markers=np.zeros(4, dtype=int),
    )

    contract = meshes_module._triangle_mesh_contract_from_audit(
        meshes_module._triangle_mesh_audit(mesh),
        reference_length=4.0,
        require_markers=True,
        stage_label="Ground mesh",
    )

    assert contract["status"] == "warn"
    assert contract["requirements"]["face_markers_present"] is True
    assert any("declared meshing scale" in message for message in contract["warnings"])


def test_triangle_mesh_audit_localizes_worst_quality_face():
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [0.0, 5.0, 0.0],
                [10.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 0.01, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 2],
                [3, 4, 5],
            ],
            dtype=int,
        ),
        markers=np.array([-1, 7], dtype=int),
    )

    audit = meshes_module._triangle_mesh_audit(mesh)
    contract = meshes_module._triangle_mesh_contract_from_audit(
        audit,
        reference_length=1.0,
        require_markers=True,
        stage_label="Ground mesh",
    )

    assert audit["worst_element_face"]["index"] == 1
    assert audit["worst_element_face"]["marker"] == 7
    assert audit["worst_element_face"]["centroid"] == pytest.approx(
        [50.0 / 3.0, 0.01 / 3.0, 0.0]
    )
    assert audit["worst_element_face"]["edge_length_min"] == pytest.approx(0.01)
    assert contract["status"] == "warn"
    assert contract["metrics"]["worst_element_face_marker"] == 7
    assert contract["metrics"]["worst_element_face_edge_length_min"] == pytest.approx(
        0.01
    )


def test_tetgen_plc_contract_from_audit_reports_precheck_failures():
    combined_surface = meshes_module._triangle_mesh_audit(
        Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [2.0, 2.0, 0.0],
                    [0.0, 2.0, 0.0],
                    [1.0e-4, 1.0e-4, 0.0],
                ]
            ),
            faces=np.array(
                [
                    [0, 1, 4],
                    [1, 2, 4],
                    [2, 3, 4],
                    [3, 0, 4],
                ],
                dtype=int,
            ),
            markers=np.empty((0,), dtype=int),
        )
    )
    plc_audit = {
        "num_boundary_facets": 5,
        "num_boundary_triangles": 5,
        "precheck": {
            "ok": False,
            "error_count": 1,
            "warning_count": 1,
            "errors": ["Surface shell contains duplicate triangles."],
            "warnings": ["Surface shell minimum triangle quality is very low (0.01)."],
            "min_edge_length": 0.1,
            "median_edge_length": 2.0,
            "min_triangle_quality": 0.01,
            "max_triangle_aspect_ratio": 120.0,
        },
        "combined_surface": combined_surface,
    }

    contract = meshes_module._tetgen_plc_contract_from_audit(
        plc_audit,
        reference_length=1.0,
    )

    assert contract["status"] == "fail"
    assert contract["requirements"]["boundary_facets_present"] is True
    assert any("duplicate triangles" in message for message in contract["errors"])
    assert any("25% of the declared meshing scale" in message for message in contract["warnings"])


def test_regularize_flat_mesh_ground_polygons_splits_case55_style_pinch():
    polygon = Polygon(
        [
            (675263.1249638698, 6581229.500071362),
            (675274.8125201109, 6581225.406301662),
            (675282.3125684819, 6581222.187530903),
            (675278.5937848733, 6581210.874928664),
            (675263.4062842835, 6581216.312428876),
            (675261.562570048, 6581210.906283745),
            (675276.0313192445, 6581205.718784033),
            (675274.656285481, 6581201.468679672),
            (675268.9687844442, 6581203.37493002),
            (675268.375039339, 6581201.65619419),
            (675265.3125405194, 6581201.687444177),
            (675263.0312874243, 6581194.437434341),
            (675245.0000353431, 6581198.656184828),
            (675241.6875341476, 6581189.499931524),
            (675232.7499284449, 6581192.249964047),
            (675235.6562160432, 6581200.437569969),
            (675239.9999676079, 6581198.96881944),
            (675241.9374480231, 6581203.937519215),
            (675248.093715854, 6581221.531320173),
            (675263.4062151447, 6581216.3125704145),
            (675265.8124304013, 6581223.499966635),
            (675261.7499317214, 6581224.999966148),
        ]
    )

    direct = meshes_module._regularize_flat_mesh_ground_polygons(
        [polygon],
        cleanup_scale=0.05,
        cleaning_diagnostics=False,
    )

    assert _nearest_nonadjacent_boundary_vertex_distance(polygon) < 1e-3
    assert len(direct) >= 1
    assert min(_nearest_nonadjacent_boundary_vertex_distance(part) for part in direct) > 0.05


def test_build_city_flat_mesh_dtcc_mesher_handles_touching_holes():
    pytest.importorskip("dtcc_mesher")

    touching_holes = Polygon(
        [(8, 8), (28, 8), (28, 28), (8, 28), (8, 8)],
        [
            [(12, 12), (16, 12), (16, 16), (12, 16), (12, 12)],
            [(16, 16), (22, 16), (22, 24), (16, 24), (16, 16)],
        ],
    )
    city = make_flat_city([make_building(touching_holes, roof_z=10.0)])

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    markers = set(np.asarray(mesh.markers, dtype=int))
    assert -2 in markers
    assert 0 in markers


def test_build_city_flat_mesh_dtcc_mesher_handles_nested_building_in_courtyard():
    pytest.importorskip("dtcc_mesher")

    shell = Polygon(
        [(8, 8), (32, 8), (32, 32), (8, 32), (8, 8)],
        [[(14, 14), (26, 14), (26, 26), (14, 26), (14, 14)]],
    )
    nested = box(18, 18, 22, 22)
    city = make_flat_city(
        [
            make_building(shell, roof_z=12.0),
            make_building(nested, roof_z=8.0),
        ]
    )

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert set(np.unique(mesh.markers)).issuperset({-2, 0, 1})


def test_condition_meshing_footprints_uses_conservative_roof_for_heterogeneous_merge():
    buildings = [
        make_building(box(8, 8, 16, 16), roof_z=4.0),
        make_building(box(16, 8, 24, 16), roof_z=10.0),
    ]

    conditioned = meshes_module._condition_meshing_footprints(
        buildings,
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=True,
        max_mesh_size=20.0,
        cleaning_diagnostics=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0, 1]]
    assert resolutions == [10.0]
    assert diagnostics["conservative_merged_roof_count"] == 1
    assert diagnostics["conservative_merged_roof_max_span"] == pytest.approx(6.0)
    assert np.unique(np.asarray(surfaces[0].vertices)[:, 2]) == pytest.approx([10.0])


def test_condition_meshing_footprints_keeps_weighted_roof_for_similar_merge():
    buildings = [
        make_building(box(8, 8, 16, 16), roof_z=4.0),
        make_building(box(16, 8, 24, 16), roof_z=5.0),
    ]

    conditioned = meshes_module._condition_meshing_footprints(
        buildings,
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=True,
        max_mesh_size=20.0,
        cleaning_diagnostics=False,
    )
    surfaces, source_map, resolutions, diagnostics = unpack_conditioned_footprints(
        conditioned
    )

    assert len(surfaces) == 1
    assert source_map == [[0, 1]]
    assert resolutions == [4.5]
    assert diagnostics["conservative_merged_roof_count"] == 0
    assert np.unique(np.asarray(surfaces[0].vertices)[:, 2]) == pytest.approx([4.5])


def test_condition_flat_mesh_ground_polygons_returns_explicit_courtyards():
    building = Polygon(
        [(10, 10), (30, 10), (30, 30), (10, 30), (10, 10)],
        [[(16, 16), (24, 16), (24, 24), (16, 24), (16, 16)]],
    )

    ground_polygons = meshes_module._condition_flat_mesh_ground_polygons(
        bounds=(0.0, 0.0, 40.0, 40.0),
        building_polygons=[building],
        hole_polygons=[],
        max_mesh_size=10.0,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
    )

    courtyard = Polygon(building.interiors[0])
    assert any(polygon.covers(courtyard) and courtyard.covers(polygon) for polygon in ground_polygons)


def test_condition_flat_mesh_ground_polygons_courtyards_exclude_nested_buildings():
    shell = Polygon(
        [(10, 10), (34, 10), (34, 34), (10, 34), (10, 10)],
        [[(14, 14), (30, 14), (30, 30), (14, 30), (14, 14)]],
    )
    nested = box(18, 18, 22, 22)

    ground_polygons = meshes_module._condition_flat_mesh_ground_polygons(
        bounds=(0.0, 0.0, 40.0, 40.0),
        building_polygons=[shell, nested],
        hole_polygons=[],
        max_mesh_size=10.0,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
    )

    assert any(polygon.covers(Point(16.0, 16.0)) for polygon in ground_polygons)
    assert all(polygon.intersection(nested).area == 0.0 for polygon in ground_polygons)


def test_prepare_surface_ground_regions_preserves_building_holes_for_courtyards():
    building = Polygon(
        [(10, 10), (30, 10), (30, 30), (10, 30), (10, 10)],
        [[(16, 16), (24, 16), (24, 24), (16, 24), (16, 16)]],
    )
    surface = make_surface(building, 12.0)

    ground_regions = meshes_module._prepare_surface_ground_regions(
        conditioned_surfaces=[surface],
        conditioned_resolution=[5.0],
        target_lods=[GeometryType.LOD1],
        bounds=(0.0, 0.0, 40.0, 40.0),
        max_mesh_size=10.0,
        min_building_detail=0.5,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
        treat_lod0_as_holes=False,
    )

    assert len(ground_regions.building_surfaces) == 1
    assert ground_regions.meshing_directives == [1]
    assert ground_regions.region_markers.count(-2) == 2
    assert ground_regions.region_markers.count(0) == 1
    assert ground_regions.region_triangle_sizes == {0: 5.0}

    building_index = ground_regions.region_markers.index(0)
    ground_indices = [
        index
        for index, marker in enumerate(ground_regions.region_markers)
        if marker == -2
    ]
    courtyard = Polygon(building.interiors[0])
    assert len(ground_regions.region_polygons[building_index].interiors) == 1
    assert building.contains(Point(ground_regions.region_points[building_index]))
    assert not courtyard.covers(Point(ground_regions.region_points[building_index]))
    assert all(
        ground_regions.region_polygons[index]
        .intersection(ground_regions.region_polygons[building_index])
        .area
        == 0.0
        for index in ground_indices
    )


def test_build_city_flat_mesh_dtcc_mesher_uses_single_coverage_call(monkeypatch):
    calls = {}

    class DummyRawMesh:
        def __init__(self):
            self.points = np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                    [1.0, 1.0],
                ],
                dtype=np.float64,
            )
            self.triangles = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.uint32)
            self.segments = np.empty((0, 2), dtype=np.uint32)
            self.markers = np.array([-2, 0], dtype=np.int32)

    class DummyMesher:
        class Coverage:
            def __init__(self, polygons, markers, *, tolerance=1e-9):
                self.polygons = tuple(polygons)
                self.markers = tuple(markers)
                self.tolerance = tolerance

        class MeshingOptions:
            def __init__(self, *, min_angle, max_edge_length, refine, max_protection_levels=None):
                self.min_angle = min_angle
                self.max_edge_length = max_edge_length
                self.refine = refine
                self.max_protection_levels = max_protection_levels

        def mesh(self, geometry, *, options):
            calls["polygon_count"] = len(geometry.polygons)
            calls["markers"] = list(geometry.markers)
            calls["min_angle"] = options.min_angle
            calls["max_edge_length"] = options.max_edge_length
            calls["refine"] = options.refine
            calls["max_protection_levels"] = options.max_protection_levels
            return DummyRawMesh()

    monkeypatch.setattr(
        dtcc_mesher_backend_module,
        "_load_dtcc_mesher",
        lambda: DummyMesher(),
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert calls["polygon_count"] == 2
    assert set(calls["markers"]) == {-2, 0}
    assert calls["min_angle"] == 20.0
    assert calls["max_edge_length"] == 5.0
    assert calls["refine"] is True
    assert calls["max_protection_levels"] == 1
    assert set(np.asarray(mesh.markers, dtype=int)) == {-1, 0}


def test_build_city_flat_mesh_dtcc_mesher_unrestricted_enables_refinement(monkeypatch):
    calls = {}

    class DummyRawMesh:
        def __init__(self):
            self.points = np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                ],
                dtype=np.float64,
            )
            self.triangles = np.array([[0, 1, 2]], dtype=np.uint32)
            self.segments = np.empty((0, 2), dtype=np.uint32)
            self.markers = np.array([-2], dtype=np.int32)

    class DummyMesher:
        class Coverage:
            def __init__(self, polygons, markers, *, tolerance=1e-9):
                self.polygons = tuple(polygons)
                self.markers = tuple(markers)
                self.tolerance = tolerance

        class MeshingOptions:
            def __init__(self, *, min_angle, max_edge_length, refine, max_protection_levels=None):
                self.min_angle = min_angle
                self.max_edge_length = max_edge_length
                self.refine = refine
                self.max_protection_levels = max_protection_levels

        def mesh(self, geometry, *, options):
            calls["max_edge_length"] = options.max_edge_length
            calls["refine"] = options.refine
            calls["max_protection_levels"] = options.max_protection_levels
            return DummyRawMesh()

    monkeypatch.setattr(
        dtcc_mesher_backend_module,
        "_load_dtcc_mesher",
        lambda: DummyMesher(),
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=None,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert calls["max_edge_length"] is None
    assert calls["refine"] is True
    assert calls["max_protection_levels"] == 1


def test_build_city_flat_mesh_dtcc_mesher_uses_explicit_region_points(monkeypatch):
    calls = {}

    class DummyRawMesh:
        def __init__(self):
            self.points = np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                ],
                dtype=np.float64,
            )
            self.triangles = np.array([[0, 1, 2]], dtype=np.uint32)
            self.segments = np.empty((0, 2), dtype=np.uint32)
            self.markers = np.array([-2], dtype=np.int32)

    class DummyMesher:
        class Coverage:
            def __init__(self, polygons, markers, *, tolerance=1e-9):
                self.polygons = tuple(polygons)
                self.markers = tuple(markers)
                self.tolerance = tolerance

            def graph(self, *, max_edge_length=None):
                calls["coverage_graph_max_edge_length"] = max_edge_length
                return DummyMesher.CoverageGraph(
                    np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64),
                    np.array([[0, 1]], dtype=np.uint32),
                    np.array([[0.25, 0.25], [0.75, 0.75]], dtype=np.float64),
                    np.array([-2, 0], dtype=np.int32),
                )

        class CoverageGraph:
            def __init__(self, points, segments, region_points, region_markers):
                self.points = points
                self.segments = segments
                self.region_points = region_points
                self.region_markers = region_markers

        class MeshingOptions:
            def __init__(self, *, min_angle, max_edge_length, refine, max_protection_levels=None):
                self.min_angle = min_angle
                self.max_edge_length = max_edge_length
                self.refine = refine
                self.max_protection_levels = max_protection_levels

        def mesh(self, geometry, *, options):
            calls["geometry_type"] = type(geometry).__name__
            calls["region_points"] = np.asarray(geometry.region_points, dtype=np.float64)
            calls["region_markers"] = np.asarray(geometry.region_markers, dtype=np.int32)
            calls["max_protection_levels"] = options.max_protection_levels
            return DummyRawMesh()

    monkeypatch.setattr(
        dtcc_mesher_backend_module,
        "_load_dtcc_mesher",
        lambda: DummyMesher(),
    )

    mesh = dtcc_mesher_backend_module.build_city_flat_mesh_with_dtcc_mesher(
        region_polygons=[box(0, 0, 10, 10), box(2, 2, 4, 4)],
        region_markers=[-2, 0],
        region_points=[
            np.array([1.0, 1.0], dtype=np.float64),
            np.array([3.0, 3.0], dtype=np.float64),
        ],
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
    )

    assert calls["coverage_graph_max_edge_length"] == 5.0
    assert calls["geometry_type"] == "CoverageGraph"
    assert np.allclose(calls["region_points"], np.array([[1.0, 1.0], [3.0, 3.0]]))
    assert calls["region_markers"].tolist() == [-2, 0]
    assert calls["max_protection_levels"] == 1
    assert mesh.faces.shape[0] == 1


def test_build_city_flat_mesh_triangle_uses_conditioned_coverage(monkeypatch):
    captured = {}

    def fake_condition(*args, **kwargs):
        return make_conditioned_footprints(
            [make_surface(box(8, 8, 18, 18), 10.0)],
            [[0]],
            [],
            {"output_grid": 0.03125},
        )

    def fake_condition_coverage_regions(**kwargs):
        captured["coverage_building_count"] = len(kwargs["building_polygons"])
        captured["coverage_markers"] = list(kwargs["building_markers"])
        building = box(8, 8, 18, 18)
        return [box(0, 0, 80, 80).difference(building), building], [-2, 7]

    def fake_build_from_coverage(
        *,
        region_polygons,
        region_markers,
        region_points=None,
        bounds,
        max_mesh_size,
        min_mesh_angle,
        backend,
        sort_triangles,
        region_triangle_sizes,
    ):
        captured["region_polygon_count"] = len(region_polygons)
        captured["region_markers"] = list(region_markers)
        captured["bounds"] = bounds
        captured["backend"] = backend
        captured["sort_triangles"] = sort_triangles
        captured["region_triangle_sizes"] = region_triangle_sizes
        return Mesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([-2], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert captured["coverage_building_count"] == 1
    assert captured["coverage_markers"] == [0]
    assert captured["region_polygon_count"] == 2
    assert captured["region_markers"] == [-2, 7]
    assert captured["bounds"] == (0.0, 0.0, 80.0, 80.0)
    assert captured["backend"] == "triangle"
    assert captured["sort_triangles"] is True
    assert captured["region_triangle_sizes"] is None
    assert mesh.faces.shape[0] >= 1


def test_build_city_flat_mesh_propagates_runtime_mesher_errors(monkeypatch):
    def fake_condition(*args, **kwargs):
        return make_conditioned_footprints(
            [make_surface(box(8, 8, 18, 18), 10.0)],
            [[0]],
            [],
            {"output_grid": 0.03125},
        )

    def fake_condition_coverage_regions(**kwargs):
        building = box(8, 8, 18, 18)
        return [box(0, 0, 80, 80).difference(building), building], [-2, 7]

    def fake_build_from_coverage(**kwargs):
        raise RuntimeError("internal triangulation error")

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    with pytest.raises(RuntimeError, match="internal triangulation error"):
        build_city_flat_mesh(
            city,
            lod=GeometryType.LOD0,
            merge_buildings=False,
            min_building_detail=0.5,
            min_building_area=1.0,
            merge_tolerance=0.5,
            max_mesh_size=None,
            min_mesh_angle=20.0,
            report_mesh_quality=False,
            mesher="dtcc_mesher",
        )


def test_build_city_flat_mesh_allows_empty_conditioned_footprints(monkeypatch):
    captured = {}

    def fake_condition(*args, **kwargs):
        return make_conditioned_footprints([], [], [], {"output_count": 0})

    def fake_build_ground_mesh_from_coverage(**kwargs):
        captured["region_polygons"] = kwargs["region_polygons"]
        captured["region_markers"] = kwargs["region_markers"]
        captured["mesher"] = kwargs["mesher"]
        return (
            Mesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([-2], dtype=int),
            ),
            "dtcc_mesher",
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(
        meshes_module,
        "_build_ground_mesh_from_coverage",
        fake_build_ground_mesh_from_coverage,
    )

    city = make_flat_city([])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert captured["region_polygons"]
    assert set(captured["region_markers"]) == {-2}
    assert captured["mesher"] == "dtcc_mesher"
    assert mesh.faces.shape[0] == 1


def test_build_city_flat_mesh_forwards_triangle_backend(monkeypatch):
    captured = {}

    def fake_condition(*args, **kwargs):
        return make_conditioned_footprints([], [], [], {"output_grid": 0.03125})

    def fake_build_ground_mesh_from_coverage(**kwargs):
        captured["mesher"] = kwargs["mesher"]
        return (
            Mesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([-2], dtype=int),
            ),
            "triangle",
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(
        meshes_module,
        "_build_ground_mesh_from_coverage",
        fake_build_ground_mesh_from_coverage,
    )

    city = make_flat_city([])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert captured["mesher"] == "triangle"
    assert mesh.faces.shape[0] == 1


def test_build_city_flat_mesh_marks_halos_for_triangle_backend(monkeypatch):
    def fake_condition(*args, **kwargs):
        return make_conditioned_footprints(
            [make_surface(box(8, 8, 18, 18), 10.0)],
            [[0]],
            [],
            {"output_grid": 0.03125},
        )

    def fake_condition_coverage_regions(**kwargs):
        building = box(8, 8, 18, 18)
        return [box(0, 0, 80, 80).difference(building), building], [-2, 7]

    def fake_build_from_coverage(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [1, 3, 2]], dtype=int),
            markers=np.array([-2, 7], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert np.array_equal(np.asarray(mesh.markers, dtype=int), np.array([-1, 7], dtype=int))


def test_builder_flat_mesh_backend_remaps_canonical_region_markers(monkeypatch):
    captured = {}

    class DummyCppMesh:
        def from_cpp(self):
            return Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 1], [1, 2, 0]], dtype=int),
                markers=np.array([0, -1, -2], dtype=int),
            )

    def fake_create_builder_polygon(polygon):
        return polygon

    def fake_build_flat_mesh(
        building_polygons,
        holes,
        subdomain_resolution,
        xmin,
        ymin,
        xmax,
        ymax,
        max_mesh_size,
        min_mesh_angle,
        sort_triangles,
        backend,
    ):
        captured["building_polygon_count"] = len(building_polygons)
        captured["subdomain_resolution"] = list(subdomain_resolution)
        captured["backend"] = backend
        return DummyCppMesh()

    monkeypatch.setattr(flat_mesh_backends_module, "create_builder_polygon", fake_create_builder_polygon)
    monkeypatch.setattr(flat_mesh_backends_module._dtcc_builder, "build_city_flat_mesh", fake_build_flat_mesh)

    mesh = flat_mesh_backends_module.build_city_flat_mesh_with_builder_backend(
        region_polygons=[box(0, 0, 10, 10), box(2, 2, 4, 4)],
        region_markers=[-2, 7],
        bounds=(0.0, 0.0, 10.0, 10.0),
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        backend="triangle",
    )

    assert captured["building_polygon_count"] == 1
    assert captured["subdomain_resolution"] == []
    assert captured["backend"] == "triangle"
    assert np.array_equal(np.asarray(mesh.markers, dtype=int), np.array([7, -1, -2], dtype=int))


def test_native_triangle_shared_vertex_regions_survive_python_conversion():
    if "triangle" not in _dtcc_builder.triangulation_backends():
        pytest.skip("Triangle backend not built")

    script = textwrap.dedent(
        """
        from shapely.geometry import Polygon
        from dtcc_core.builder import _dtcc_builder
        from dtcc_core.builder.model_conversion import create_builder_polygon

        polygon_a = create_builder_polygon(
            Polygon([(2, 2), (10, 2), (10, 10), (2, 10), (2, 2)])
        )
        polygon_b = create_builder_polygon(
            Polygon([(10, 10), (18, 10), (18, 18), (10, 18), (10, 10)])
        )
        raw_mesh = _dtcc_builder.build_city_flat_mesh(
            [polygon_a, polygon_b],
            [],
            [],
            0.0,
            0.0,
            20.0,
            20.0,
            5.0,
            20.0,
            False,
            "triangle",
        )
        mesh = raw_mesh.from_cpp()
        print(mesh.vertices.shape, mesh.faces.shape, sorted(set(mesh.markers.tolist())))
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd="/Users/logg/scratch/dtcc/dtcc-core",
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_native_triangle_build_city_flat_mesh_handles_internal_loops_and_holes():
    if "triangle" not in _dtcc_builder.triangulation_backends():
        pytest.skip("Triangle backend not built")

    building = create_builder_polygon(
        Polygon(
            [(10, 10), (25, 10), (25, 25), (10, 25), (10, 10)],
            [[(14, 14), (20, 14), (20, 20), (14, 20), (14, 14)]],
        )
    )
    explicit_hole = create_builder_polygon(
        Polygon([(30, 30), (40, 30), (40, 40), (30, 40), (30, 30)])
    )

    raw_mesh = _dtcc_builder.build_city_flat_mesh(
        [building],
        [explicit_hole],
        [5.0],
        0.0,
        0.0,
        100.0,
        100.0,
        10.0,
        20.0,
        False,
        "triangle",
    )
    mesh = builder_mesh_to_mesh(raw_mesh)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert set(np.asarray(mesh.markers, dtype=int)) >= {-2, 0}


def test_build_city_surface_mesh_reduces_lod_from_source_map(monkeypatch):
    surfaces = [
        make_surface(box(10, 10, 20, 20), 10.0),
        make_surface(box(30, 10, 40, 20), 12.0),
    ]
    source_map = [[0, 1], [2]]
    resolutions = [4.0, 5.0]

    captured = {}

    def fake_condition(*args, **kwargs):
        return make_conditioned_footprints(
            surfaces,
            source_map,
            resolutions,
            {"output_count": 2},
        )

    class DummyCppMesh:
        def from_cpp(self):
            return Mesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]], dtype=int),
            )

    def fake_condition_building_regions(**kwargs):
        return kwargs["building_polygons"], kwargs["building_markers"], [[0], [1]]

    def fake_condition_ground_polygons(**kwargs):
        return [
            box(0, 0, 80, 80).difference(
                meshes_module.unary_union(kwargs["building_polygons"])
            )
        ]

    def fake_build_from_coverage(**kwargs):
        captured["backend"] = kwargs["backend"]
        captured["region_markers"] = kwargs["region_markers"]
        captured["region_triangle_sizes"] = kwargs["region_triangle_sizes"]
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    def fake_build_terrain_surface_mesh_from_ground_mesh(*args, **kwargs):
        return object()

    def fake_build_surface_mesh_from_terrain_mesh(
        building_surfaces,
        building_lod_switches,
        terrain_mesh,
        smoothing,
        merge_meshes,
    ):
        captured["lod_switches"] = building_lod_switches
        return [DummyCppMesh()]

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(
        meshes_module,
        "_consume_conditioned_building_regions_with_sources",
        fake_condition_building_regions,
    )
    monkeypatch.setattr(
        meshes_module,
        "_condition_flat_mesh_ground_polygons",
        fake_condition_ground_polygons,
    )
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_terrain_surface_mesh_from_ground_mesh",
        fake_build_terrain_surface_mesh_from_ground_mesh,
    )
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_city_surface_mesh_from_terrain_mesh",
        fake_build_surface_mesh_from_terrain_mesh,
    )

    city = make_flat_city(
        [
            make_building(box(10, 10, 20, 20), roof_z=10.0, lod1_polygon=box(10, 10, 20, 20)),
            make_building(box(20, 10, 30, 20), roof_z=11.0),
            make_building(box(30, 10, 40, 20), roof_z=12.0, lod2_polygon=box(30, 10, 40, 20)),
        ]
    )

    mesh = build_city_surface_mesh(
        city,
        lod=[GeometryType.LOD1, GeometryType.LOD0, GeometryType.LOD2],
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        building_mesh_triangle_size=6.0,
        max_mesh_size=8.0,
        min_mesh_angle=20.0,
        merge_meshes=True,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert mesh.faces.shape[0] >= 1
    assert captured["lod_switches"] == [1, 2]
    assert captured["backend"] == "triangle"
    assert captured["region_markers"] == [-2, 0, 1]
    assert captured["region_triangle_sizes"] == {0: 4.0, 1: 5.0}


def test_build_city_surface_mesh_uses_raw_ground_markers(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        captured["add_halo_markers"] = kwargs["add_halo_markers"]
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )

    mesh = build_city_surface_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.faces.shape[0] == 1
    assert captured["add_halo_markers"] is False


def test_split_ground_mesh_building_components_splits_vertex_touching_patches():
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 2],
                [2, 3, 4],
            ],
            dtype=int,
        ),
        markers=np.array([0, 0], dtype=int),
    )
    building_surface = make_surface(box(0, 0, 2, 1), 10.0)

    split_mesh, split_surfaces, split_directives = (
        meshes_module._split_ground_mesh_building_components(
            ground_mesh=ground_mesh,
            building_surfaces=[building_surface],
            meshing_directives=[1],
        )
    )

    assert split_mesh.markers.tolist() == [0, 1]
    assert len(split_surfaces) == 2
    assert split_directives == [1, 1]


def test_stabilize_shell_building_polygon_removes_low_clearance_hole():
    polygon = Polygon(
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [
            [
                (2.0, 2.0),
                (8.0, 2.0),
                (8.0, 8.0),
                (5.0002, 8.0),
                (5.0, 7.9998),
                (2.0, 8.0),
            ]
        ],
    )

    stabilized, removed = meshes_module._stabilize_shell_building_polygon(
        polygon,
        min_hole_clearance=0.01,
    )

    assert removed == 1
    assert len(stabilized.interiors) == 0


def test_split_weakly_pinched_shell_polygon_splits_exterior_revisit():
    polygon = Polygon(
        [
            (674614.125, 6580137.40625),
            (674617.6249773939, 6580141.031226587),
            (674613.874953289, 6580144.875001294),
            (674617.2812514388, 6580148.000045467),
            (674621.0625469746, 6580144.06249805),
            (674617.6250215588, 6580141.031225638),
            (674627.0, 6580130.4375),
            (674651.09375, 6580152.65625),
            (674634.6875, 6580170.53125),
            (674620.5625, 6580186.0625),
            (674598.21875, 6580165.59375),
            (674586.84375, 6580154.3125),
            (674600.3125, 6580139.65625),
            (674616.78125, 6580121.0),
            (674623.28125, 6580127.03125),
        ]
    )

    parts = meshes_module._split_weakly_pinched_shell_polygon(
        polygon,
        pinch_tolerance=0.05,
    )

    assert len(parts) == 2
    assert all(part.minimum_clearance > 1.0 for part in parts)
    assert check_mesher_handoff_profile(parts)["status"] == "pass"


def test_build_city_volume_mesh_uses_shared_surface_pipeline(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        captured["target_lods"] = kwargs["target_lods"]
        captured["conditioned_resolution"] = kwargs["conditioned_resolution"]
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        captured["mesher"] = kwargs["mesher"]
        captured["add_halo_markers"] = kwargs["add_halo_markers"]
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        captured["surface_called"] = True
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_mesh_faces"] = len(kwargs["mesh"].faces)
        captured["closure_mesh_faces"] = len(kwargs["closure_mesh"].faces)
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    def fake_tetgen_plc_audit(**kwargs):
        return {
            "num_boundary_facets": 1,
            "num_boundary_triangles": 1,
            "precheck": {
                "ok": True,
                "error_count": 0,
                "warning_count": 0,
                "errors": [],
                "warnings": [],
            },
            "combined_surface": {
                "face_count": 1,
                "vertex_count": 3,
            },
            "combined_orientation": {
                "ok": True,
                "component_count": 1,
                "inverted_component_count": 0,
                "component_orientations": [1],
            },
        }

    def fake_tetgen_plc_contract_from_audit(*args, **kwargs):
        return {
            "status": "pass",
            "requirements": {
                "boundary_facets_present": True,
                "precheck_passed": True,
                "consistent_shared_edge_winding": True,
            },
            "errors": [],
            "warnings": [],
            "metrics": {},
        }

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "_tetgen_plc_audit", fake_tetgen_plc_audit)
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        fake_tetgen_plc_contract_from_audit,
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert captured["mesher"] == "dtcc_mesher"
    assert captured["add_halo_markers"] is False
    assert captured["target_lods"] == [GeometryType.LOD1]
    assert captured["conditioned_resolution"] == [4.0]
    assert captured["surface_called"] is True
    assert captured["tetgen_mesh_faces"] == 1
    assert captured["closure_mesh_faces"] == 1
    assert volume_mesh.cells.shape[0] == 1


def test_build_city_volume_mesh_stage_audit_records_stage_contracts(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25, "geos_exception_count": 0}
    stage_audit: dict[str, object] = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(
            terrain,
            terrain_raster,
            [conditioned_surface],
            [[0]],
            [4.0],
            diagnostics,
            min_building_detail=kwargs.get("min_building_detail", 0.5),
        )

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    def fake_plc_audit(**kwargs):
        combined_surface = meshes_module._triangle_mesh_audit(
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.empty((0,), dtype=int),
            )
        )
        return {
            "num_vertices": 3,
            "num_shell_faces": 1,
            "num_boundary_facets": 5,
            "num_boundary_triangles": 5,
            "bounds": {"xmin": 0.0, "xmax": 10.0, "ymin": 0.0, "ymax": 10.0, "zmin": 0.0, "zmax": 100.0},
            "boundary_facets": {},
            "precheck": {
                "ok": True,
                "error_count": 0,
                "warning_count": 0,
                "errors": [],
                "warnings": [],
                "degenerate_face_count": 0,
                "duplicate_face_count": 0,
                "nonmanifold_edge_count": 0,
                "open_edge_count": 0,
                "min_edge_length": 10.0,
                "median_edge_length": 10.0,
                "min_face_area": 50.0,
                "median_face_area": 50.0,
                "min_triangle_quality": 0.8,
                "max_triangle_aspect_ratio": 2.0,
                "boundary_facets": {},
            },
            "combined_surface": combined_surface,
        }

    def fake_tetgen_build(**kwargs):
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "_tetgen_plc_audit", fake_plc_audit)
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        stage_audit=stage_audit,
    )

    attempt = stage_audit["attempts"][0]
    conditioned_stage = attempt["stages"]["conditioned_footprints"]
    assert conditioned_stage["declared_scale"] == pytest.approx(0.25)
    assert conditioned_stage["footprints"]["resolution_min"] == pytest.approx(4.0)
    assert conditioned_stage["contract"]["status"] == "pass"
    assert attempt["stages"]["ground_mesh"]["contract"]["status"] == "pass"
    assert attempt["stages"]["surface_shell"]["contract"]["status"] == "pass"
    assert attempt["stages"]["plc"]["contract"]["status"] == "pass"
    assert volume_mesh.stage_audit["selected_attempt_label"] == attempt["label"]


def test_build_city_surface_mesh_from_ground_mesh_snaps_boundary_vertices(monkeypatch):
    raster = Raster()
    raster.data = np.zeros((80, 80), dtype=float)
    raster.set_bounds(Bounds(0.0, 0.0, 80.0, 80.0))
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.05, 0.0, 0.0],
                [80.0, 0.04, 0.0],
                [79.96, 80.0, 0.0],
                [0.0, 79.97, 0.0],
                [40.0, 40.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4]], dtype=int),
        markers=np.array([0, 0, 0], dtype=int),
    )
    captured = {}

    class DummyBuilderMesh:
        def __init__(self, mesh):
            self._mesh = mesh

        def from_cpp(self):
            return self._mesh

    def fake_mesh_to_builder_mesh(mesh):
        captured["aligned_vertices"] = np.array(mesh.vertices, copy=True)
        return DummyBuilderMesh(mesh)

    def fake_build_terrain(builder_mesh, builder_dem, smoothing):
        captured["terrain_builder_input"] = np.array(builder_mesh._mesh.vertices, copy=True)
        return DummyBuilderMesh(
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [80.0, 0.0, 0.0],
                        [0.0, 80.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([0], dtype=int),
            )
        )

    def fail_city_surface(*args, **kwargs):
        raise AssertionError("Terrain-only surface mesh should not build building surfaces.")

    monkeypatch.setattr(meshes_module, "mesh_to_builder_mesh", fake_mesh_to_builder_mesh)
    monkeypatch.setattr(meshes_module, "raster_to_builder_gridfield", lambda raster: "grid")
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_terrain_surface_mesh_from_ground_mesh",
        fake_build_terrain,
    )
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_city_surface_mesh_from_terrain_mesh",
        fail_city_surface,
    )

    mesh = meshes_module._build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=raster,
        building_surfaces=[],
        meshing_directives=[],
        smoothing=0,
        merge_meshes=True,
    )

    assert np.allclose(captured["aligned_vertices"][0], [0.0, 0.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][1], [80.0, 0.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][2], [80.0, 80.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][3], [0.0, 80.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][4], [40.0, 40.0, 0.0])
    assert np.allclose(captured["terrain_builder_input"], captured["aligned_vertices"])
    assert mesh.faces.shape[0] == 1


def test_ground_mesh_raster_alignment_clamps_bounded_single_axis_overshoot():
    raster = Raster()
    raster.data = np.zeros((10, 10), dtype=float)
    raster.set_bounds(Bounds(0.0, 0.0, 10.0, 10.0))
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [11.5, 4.0, 0.0],
                [5.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([0], dtype=int),
    )

    aligned = meshes_module._snap_ground_mesh_to_raster_bounds(ground_mesh, raster)

    assert np.allclose(aligned.vertices[1], [10.0, 4.0, 0.0])
    assert np.all(aligned.vertices[:, 0] >= 0.0)
    assert np.all(aligned.vertices[:, 0] <= 10.0)
    assert np.all(aligned.vertices[:, 1] >= 0.0)
    assert np.all(aligned.vertices[:, 1] <= 10.0)


def test_ground_mesh_raster_alignment_rejects_unbounded_overshoot():
    raster = Raster()
    raster.data = np.zeros((10, 10), dtype=float)
    raster.set_bounds(Bounds(0.0, 0.0, 10.0, 10.0))
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [12.5, 4.0, 0.0],
                [5.0, 10.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([0], dtype=int),
    )

    with pytest.raises(ValueError, match="outside terrain raster bounds"):
        meshes_module._snap_ground_mesh_to_raster_bounds(ground_mesh, raster)


def test_prepare_surface_ground_regions_clips_buildings_to_raster_bounds():
    bounds = (0.0, 0.0, 10.0, 10.0)
    regions = meshes_module._prepare_surface_ground_regions(
        conditioned_surfaces=[make_surface(box(-5.0, 2.0, 5.0, 8.0), 12.0)],
        conditioned_resolution=[2.0],
        target_lods=[GeometryType.LOD3],
        bounds=bounds,
        max_mesh_size=2.0,
        min_building_detail=0.5,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
        treat_lod0_as_holes=False,
    )

    assert len(regions.building_surfaces) == 1
    for polygon in regions.region_polygons:
        xmin, ymin, xmax, ymax = polygon.bounds
        assert xmin >= bounds[0] - 1.0e-9
        assert ymin >= bounds[1] - 1.0e-9
        assert xmax <= bounds[2] + 1.0e-9
        assert ymax <= bounds[3] + 1.0e-9


def test_build_city_surface_mesh_stage_audit_records_core_stages(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=12.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    building_surface = make_surface(box(10, 10, 20, 20), 12.0)
    conditioned = make_conditioned_footprints(
        [building_surface],
        [[0]],
        [5.0],
        {"output_grid": 0.03125},
    )
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [80.0, 0.0, 0.0],
                [80.0, 80.0, 0.0],
                [0.0, 80.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, 0], dtype=int),
    )

    monkeypatch.setattr(
        meshes_module,
        "_prepare_city_meshing_inputs",
        lambda *args, **kwargs: (terrain, terrain_raster, conditioned),
    )
    monkeypatch.setattr(
        meshes_module,
        "_prepare_surface_ground_regions",
        lambda **kwargs: make_prepared_ground_regions(
            [building_surface],
            [3],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 5.0},
            [np.array([1.0, 1.0]), np.array([11.0, 11.0])],
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_build_ground_mesh_stage",
        lambda **kwargs: meshes_module.BuiltGroundMesh(
            mesh=ground_mesh,
            mesher="dtcc_mesher",
            audit={"contract": {"status": "pass"}},
            contract={"ok": True, "status": "pass", "errors": []},
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_split_ground_mesh_building_components",
        lambda **kwargs: (ground_mesh, [building_surface], [3]),
    )
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        lambda **kwargs: Mesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        ),
    )

    stage_audit = {}
    mesh = build_city_surface_mesh(
        city,
        lod=GeometryType.LOD0,
        report_mesh_quality=False,
        stage_audit=stage_audit,
    )

    attempt = stage_audit["attempts"][0]
    assert attempt["backend"] == "surface"
    assert attempt["result"]["status"] == "success"
    assert {"conditioned_footprints", "surface_regions", "ground_mesh"} <= set(
        attempt["stages"]
    )
    assert mesh.stage_audit is stage_audit


def test_build_city_surface_mesh_unmerged_components_are_compact():
    city = make_flat_city([make_building(box(20, 20, 40, 40), roof_z=10.0)])
    terrain, terrain_raster, conditioned_footprints = (
        meshes_module._prepare_city_meshing_inputs(
            city,
            lod=GeometryType.LOD0,
            min_building_detail=0.5,
            min_building_area=1.0,
            merge_tolerance=0.25,
            merge_buildings=True,
            max_mesh_size=8.0,
            cleaning_diagnostics=False,
        )
    )
    target_lods = meshes_module._resolve_conditioned_target_lods(
        city.buildings,
        GeometryType.LOD0,
        conditioned_footprints.source_map,
    )
    shell_target_lods = meshes_module._promote_volume_shell_target_lods(target_lods)
    surface_regions = meshes_module._prepare_surface_ground_regions(
        conditioned_surfaces=conditioned_footprints.surfaces,
        conditioned_resolution=conditioned_footprints.subdomain_resolution,
        target_lods=shell_target_lods,
        bounds=(
            terrain.bounds.xmin,
            terrain.bounds.ymin,
            terrain.bounds.xmax,
            terrain.bounds.ymax,
        ),
        max_mesh_size=8.0,
        min_building_detail=0.5,
        footprint_diagnostics=conditioned_footprints.diagnostics,
        cleaning_diagnostics=False,
        treat_lod0_as_holes=False,
    )
    ground_mesh, _ = meshes_module._build_ground_mesh_from_coverage(
        region_polygons=surface_regions.region_polygons,
        region_markers=surface_regions.region_markers,
        region_points=surface_regions.region_points,
        bounds=(
            terrain.bounds.xmin,
            terrain.bounds.ymin,
            terrain.bounds.xmax,
            terrain.bounds.ymax,
        ),
        max_mesh_size=8.0,
        min_mesh_angle=25.0,
        mesher=meshes_module.resolve_2d_mesher("auto"),
        sort_triangles=False,
        region_triangle_sizes=surface_regions.region_triangle_sizes,
        add_halo_markers=False,
    )
    ground_mesh, surface_buildings, surface_directives = (
        meshes_module._split_ground_mesh_building_components(
            ground_mesh=ground_mesh,
            building_surfaces=surface_regions.building_surfaces,
            meshing_directives=surface_regions.meshing_directives,
        )
    )
    components = meshes_module._build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=terrain_raster,
        building_surfaces=surface_buildings,
        meshing_directives=surface_directives,
        smoothing=0,
        merge_meshes=False,
    )

    assert len(components) == 2
    for component in components:
        faces = np.asarray(component.faces, dtype=np.int64)
        used_vertices = np.unique(faces.reshape(-1))
        assert len(used_vertices) == len(component.vertices)


def _fake_prebuilt_tetgen_plc(mesh: Mesh) -> meshes_module.tetgen_utils.TetgenPLC:
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    boundary_facets = [[0, 1, 2]]
    diagnostics = meshes_module.tetgen_utils.TetgenPLCDiagnostics(
        num_vertices=int(len(vertices)),
        num_faces=int(len(faces)),
        num_boundary_facets=len(boundary_facets),
        min_edge_length=1.0,
        median_edge_length=1.0,
        min_face_area=1.0,
        median_face_area=1.0,
        min_triangle_quality=1.0,
        max_triangle_aspect_ratio=1.0,
        degenerate_face_count=0,
        duplicate_face_count=0,
        nonmanifold_edge_count=0,
        open_edge_count=0,
        boundary_facets={},
        errors=[],
        warnings=[],
    )
    return meshes_module.tetgen_utils.TetgenPLC(
        vertices=vertices,
        shell_faces=faces,
        boundary_facets=boundary_facets,
        boundary_facet_markers=[-2],
        audit_boundary_triangles=np.asarray(boundary_facets, dtype=np.int64),
        diagnostics=diagnostics,
    )


def test_build_city_volume_mesh_keeps_requested_tetgen_switches(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    terrain_raster.data[0, 0] = 0.1
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([0, 0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_switches"] = dict(kwargs["switches_params"])
        captured["top_cap_backend"] = kwargs["top_cap_backend"]
        captured["top_cap_max_mesh_size"] = kwargs["top_cap_max_mesh_size"]
        captured["top_cap_min_mesh_angle"] = kwargs["top_cap_min_mesh_angle"]
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_vertical_wall_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_wall_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        top_cap_max_mesh_size=9.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        max_volume=42.0,
        tetgen_switches={"quality": (1.6, 25.0)},
    )

    assert captured["tetgen_switches"]["quality"] == (1.6, 25.0)
    assert captured["tetgen_switches"]["preserve_surface"] is False
    assert captured["tetgen_switches"]["max_volume"] == 42.0
    assert captured["tetgen_switches"]["max_added_points"] is None
    assert captured["tetgen_switches"]["optimize_max_dihedral"] == 175.0
    assert captured["top_cap_max_mesh_size"] == 9.0


def test_build_city_volume_mesh_uses_split_surface_default_without_flat_special_case(
    monkeypatch,
):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}
    stage_audit = {}
    refinement_calls = {"transition": 0, "wall": 0, "horizontal": 0}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([0, 0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_switches"] = dict(kwargs["switches_params"])
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            refinement_calls.__setitem__("transition", refinement_calls["transition"] + 1)
            or mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            refinement_calls.__setitem__("horizontal", refinement_calls["horizontal"] + 1)
            or mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_vertical_wall_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            refinement_calls.__setitem__("wall", refinement_calls["wall"] + 1)
            or mesh,
            meshes_module._tetgen_shell_wall_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_volume_mesh_quality_snapshot",
        lambda volume_mesh: {
            "aspect_ratio_max": 1.0,
            "element_quality_min": 1.0,
            "min_edge_length": 1.0,
            "high_aspect_ratio_count": 0.0,
            "low_quality_count": 0.0,
        },
    )

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        stage_audit=stage_audit,
    )

    assert captured["tetgen_switches"]["preserve_surface"] is False
    assert captured["tetgen_switches"]["quality"] is None
    assert captured["tetgen_switches"]["optimize_max_dihedral"] == 175.0
    assert refinement_calls == {"transition": 1, "wall": 0, "horizontal": 1}
    attempts = {attempt["label"]: attempt for attempt in stage_audit["attempts"]}
    assert attempts["attempt-1"]["config"]["preserve_surface_requested"] is False
    assert attempts["attempt-1"]["config"]["terrain_effectively_flat"] is True
    assert attempts["attempt-1"]["config"]["tetgen_shell_refinement_enabled"] is True
    selection = attempts["attempt-1"]["stages"]["surface_shell"][
        "tetgen_shell_horizontal_refinement_selection"
    ]
    assert selection["selected_variant"] == "unrefined"
    assert selection["reason"] == "no_candidate_edges"


def test_build_city_volume_mesh_allows_empty_conditioned_footprints(monkeypatch):
    city = make_flat_city([])
    terrain = city.terrain
    terrain_raster = terrain.raster
    captured = {}

    def fake_prepare(*args, **kwargs):
        return (
            terrain,
            terrain_raster,
            make_conditioned_footprints([], [], [], {"output_grid": 0.25}),
        )

    def fake_prepare_regions(**kwargs):
        captured["conditioned_surfaces"] = kwargs["conditioned_surfaces"]
        return make_prepared_ground_regions(
            [],
            [],
            [box(0, 0, 80, 80)],
            [-2],
            {},
            [np.array([40.0, 40.0], dtype=np.float64)],
        )

    def fake_build_ground(**kwargs):
        captured["region_markers"] = list(kwargs["region_markers"])
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [80.0, 0.0, 0.0],
                        [80.0, 80.0, 0.0],
                        [0.0, 80.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        captured["building_surfaces"] = list(kwargs["building_surfaces"])
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [80.0, 0.0, 0.0],
                    [80.0, 80.0, 0.0],
                    [0.0, 80.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([-1, -1], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_faces"] = len(kwargs["mesh"].faces)
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert captured["conditioned_surfaces"] == []
    assert captured["region_markers"] == [-2]
    assert captured["building_surfaces"] == []
    assert captured["tetgen_faces"] == 2
    assert volume_mesh.cells.shape[0] == 1


def test_build_city_volume_mesh_respects_explicit_tetgen_switches_for_dtcc_mesher(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([0, 0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_switches"] = dict(kwargs["switches_params"])
        captured["top_cap_backend"] = kwargs["top_cap_backend"]
        captured["top_cap_max_mesh_size"] = kwargs["top_cap_max_mesh_size"]
        captured["top_cap_min_mesh_angle"] = kwargs["top_cap_min_mesh_angle"]
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_switches={
            "quality": (1.6, 25.0),
            "preserve_surface": False,
            "max_added_points": 1234,
        },
    )

    assert captured["tetgen_switches"]["quality"] == (1.6, 25.0)
    assert captured["tetgen_switches"]["preserve_surface"] is False
    assert captured["tetgen_switches"]["max_added_points"] == 1234
    assert captured["tetgen_switches"]["optimize_max_dihedral"] == 175.0
    assert captured["top_cap_backend"] == "dtcc_mesher"
    assert captured["top_cap_max_mesh_size"] == 6.0
    assert captured["top_cap_min_mesh_angle"] == 20.0


def test_build_city_volume_mesh_saves_tetgen_debug_meshes(monkeypatch, tmp_path):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_save_debug_meshes(**kwargs):
        captured["debug"] = kwargs
        return {
            "ground": str(tmp_path / "ground.xdmf"),
            "shell": str(tmp_path / "shell.xdmf"),
            "plc": str(tmp_path / "plc.xdmf"),
        }

    def fake_tetgen_build(**kwargs):
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_save_tetgen_debug_meshes", fake_save_debug_meshes)
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_debug_output_dir=tmp_path,
        tetgen_debug_output_stem="case_055",
    )

    assert captured["debug"]["output_dir"] == tmp_path
    assert captured["debug"]["stem"] == "case_055"
    assert captured["debug"]["top_cap_backend"] == "dtcc_mesher"
    assert captured["debug"]["top_cap_max_mesh_size"] == 6.0
    assert captured["debug"]["top_cap_min_mesh_angle"] == 20.0


def test_refine_near_horizontal_surface_faces_for_tetgen_splits_long_roof_faces():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 20.0, 0.0],
                [0.0, 20.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([7, 7], dtype=int),
    )

    refined_mesh, stats = meshes_module._refine_near_horizontal_surface_faces_for_tetgen(
        surface_mesh,
        max_mesh_size=10.0,
    )

    refined_vertices = np.asarray(refined_mesh.vertices, dtype=float)
    refined_faces = np.asarray(refined_mesh.faces, dtype=int)
    points = refined_vertices[refined_faces]
    edge_lengths = np.stack(
        [
            np.linalg.norm(points[:, 1] - points[:, 0], axis=1),
            np.linalg.norm(points[:, 2] - points[:, 1], axis=1),
            np.linalg.norm(points[:, 0] - points[:, 2], axis=1),
        ],
        axis=1,
    )

    assert stats["applied"] is True
    assert stats["rounds"] >= 1
    assert stats["candidate_roof_faces"] >= 2
    assert stats["candidate_ground_faces"] == 0
    assert refined_mesh.faces.shape[0] > surface_mesh.faces.shape[0]
    assert refined_mesh.vertices.shape[0] > surface_mesh.vertices.shape[0]
    assert np.all(refined_mesh.markers == 7)
    assert float(edge_lengths.max()) <= 12.5 + 1.0e-9


def test_refine_near_horizontal_surface_faces_for_tetgen_skips_flat_ground_faces():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 20.0, 0.0],
                [0.0, 20.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, -2], dtype=int),
    )

    refined_mesh, stats = meshes_module._refine_near_horizontal_surface_faces_for_tetgen(
        surface_mesh,
        max_mesh_size=10.0,
    )

    assert stats["applied"] is False
    assert stats["ground_refinement_enabled"] is False
    assert stats["candidate_ground_faces"] == 0
    assert np.allclose(refined_mesh.vertices, surface_mesh.vertices)
    assert np.array_equal(refined_mesh.faces, surface_mesh.faces)
    assert np.array_equal(refined_mesh.markers, surface_mesh.markers)


def test_refine_near_horizontal_surface_faces_for_tetgen_refines_relief_ground_faces():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.3],
                [20.0, 20.0, 0.6],
                [0.0, 20.0, 0.3],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([-2, -2], dtype=int),
    )

    refined_mesh, stats = meshes_module._refine_near_horizontal_surface_faces_for_tetgen(
        surface_mesh,
        max_mesh_size=10.0,
    )

    assert stats["applied"] is True
    assert stats["ground_refinement_enabled"] is True
    assert stats["candidate_ground_faces"] >= 2
    assert refined_mesh.faces.shape[0] > surface_mesh.faces.shape[0]
    assert refined_mesh.vertices.shape[0] > surface_mesh.vertices.shape[0]
    assert np.all(refined_mesh.markers == -2)


def test_refine_near_horizontal_surface_faces_for_tetgen_skips_steep_faces():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 0.0, 20.0],
            ]
        ),
        faces=np.array([[0, 1, 2]], dtype=int),
        markers=np.array([7], dtype=int),
    )

    refined_mesh, stats = meshes_module._refine_near_horizontal_surface_faces_for_tetgen(
        surface_mesh,
        max_mesh_size=5.0,
    )

    assert stats["applied"] is False
    assert np.allclose(refined_mesh.vertices, surface_mesh.vertices)
    assert np.array_equal(refined_mesh.faces, surface_mesh.faces)
    assert np.array_equal(refined_mesh.markers, surface_mesh.markers)


def test_select_tetgen_transition_refinement_edges_targets_ground_wall_ring():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 20.0, 1.0],
                [20.0, 0.0, 10.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 1, 3]], dtype=int),
        markers=np.array([-2, 0], dtype=int),
    )

    split_edges, stats = meshes_module._select_tetgen_transition_refinement_edges(
        surface_mesh,
        edge_threshold=10.0,
    )

    assert stats["applied"] is True
    assert stats["candidate_transition_edges"] == 1
    assert stats["candidate_terrain_faces"] == 1
    assert stats["candidate_wall_faces"] == 1
    assert stats["ground_refinement_enabled"] is True
    assert stats["ground_relief_median"] == pytest.approx(1.0)
    assert (0, 1) in split_edges
    assert (0, 2) in split_edges
    assert (0, 3) not in split_edges


def test_refine_ground_building_transition_faces_for_tetgen_splits_local_ring():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 20.0, 1.0],
                [20.0, 0.0, 10.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 1, 3]], dtype=int),
        markers=np.array([-2, 0], dtype=int),
    )

    refined_mesh, stats = (
        meshes_module._refine_ground_building_transition_faces_for_tetgen(
            surface_mesh,
            max_mesh_size=10.0,
        )
    )

    assert stats["enabled"] is True
    assert stats["applied"] is True
    assert stats["candidate_transition_edges"] == 1
    assert stats["ground_refinement_enabled"] is True
    assert refined_mesh.faces.shape[0] > surface_mesh.faces.shape[0]
    assert refined_mesh.vertices.shape[0] > surface_mesh.vertices.shape[0]
    assert set(np.asarray(refined_mesh.markers, dtype=int)) == {-2, 0}


def test_refine_ground_building_transition_faces_for_tetgen_skips_low_relief_ring():
    surface_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [20.0, 20.0, 0.2],
                [20.0, 0.0, 10.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 1, 3]], dtype=int),
        markers=np.array([-2, 0], dtype=int),
    )

    refined_mesh, stats = (
        meshes_module._refine_ground_building_transition_faces_for_tetgen(
            surface_mesh,
            max_mesh_size=10.0,
        )
    )

    assert stats["enabled"] is True
    assert stats["ground_refinement_enabled"] is False
    assert stats["applied"] is False
    assert stats["candidate_transition_edges"] == 1
    assert np.allclose(refined_mesh.vertices, surface_mesh.vertices)
    assert np.array_equal(refined_mesh.faces, surface_mesh.faces)


def test_build_city_volume_mesh_captures_quality_failure_artifacts(monkeypatch, tmp_path):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [9.9, 9.9, 0.01],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    def fake_capture_quality_failure(**kwargs):
        captured["quality_failure"] = kwargs
        return {
            "report": str(tmp_path / "case_055.tetgen-quality-failure.json"),
            "debug_meshes": {
                "ground": str(tmp_path / "ground.xdmf"),
                "shell": str(tmp_path / "shell.xdmf"),
                "plc": str(tmp_path / "plc.xdmf"),
            },
        }

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)
    monkeypatch.setattr(
        meshes_module,
        "_capture_tetgen_quality_failure_artifacts",
        fake_capture_quality_failure,
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_volume_mesh_quality_snapshot",
        lambda mesh: {
            "aspect_ratio_max": 500.0,
            "element_quality_min": 0.02,
            "min_edge_length": 0.01,
            "high_aspect_ratio_count": 1.0,
            "low_quality_count": 1.0,
        },
    )

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_quality_failure_output_dir=tmp_path,
        tetgen_quality_failure_output_stem="case_055",
    )

    assert captured["quality_failure"]["output_dir"] == tmp_path
    assert captured["quality_failure"]["stem"] == "case_055"
    assert captured["quality_failure"]["top_cap_backend"] == "dtcc_mesher"
    assert captured["quality_failure"]["top_cap_max_mesh_size"] == 6.0
    assert captured["quality_failure"]["top_cap_min_mesh_angle"] == 20.0


def test_build_city_volume_mesh_ignores_quality_failure_capture_errors(
    monkeypatch, tmp_path
):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(
            terrain,
            terrain_raster,
            [conditioned_surface],
            [[0]],
            [4.0],
            {"output_grid": 0.25},
        )

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [9.9, 9.9, 0.01],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_transition_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        lambda mesh, *, max_mesh_size: (
            mesh,
            meshes_module._tetgen_shell_refinement_disabled_stats(),
        ),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)
    monkeypatch.setattr(
        meshes_module,
        "_capture_tetgen_quality_failure_artifacts",
        lambda **kwargs: (_ for _ in ()).throw(ValueError("debug export failed")),
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_volume_mesh_quality_snapshot",
        lambda mesh: {
            "aspect_ratio_max": 500.0,
            "element_quality_min": 0.02,
            "min_edge_length": 0.01,
            "high_aspect_ratio_count": 1.0,
            "low_quality_count": 1.0,
        },
    )

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_quality_failure_output_dir=tmp_path,
        tetgen_quality_failure_output_stem="case_055",
    )

    assert volume_mesh.cells.shape == (1, 4)


def test_build_city_volume_mesh_uses_refined_shell_without_retry(
    monkeypatch, tmp_path
):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    terrain_raster.data[0, 0] = 0.1
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    stage_audit = {}
    calls: list[dict[str, object]] = []

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics)

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_refine_transition(mesh, *, max_mesh_size):
        refined = Mesh(
            vertices=np.array(mesh.vertices, copy=True),
            faces=np.array(mesh.faces, copy=True),
            markers=np.array(mesh.markers, copy=True),
        )
        refined.transition_refined_tag = True
        return refined, {
            "enabled": True,
            "applied": True,
            "candidate_transition_edges": 1,
            "candidate_terrain_faces": 1,
            "candidate_wall_faces": 1,
            "ground_relief_median": 1.0,
            "ground_refinement_enabled": True,
            "split_edges": 1,
            "added_vertices": 1,
            "added_faces": 2,
            "edge_threshold": 7.5,
        }

    def fake_refine_horizontal(mesh, *, max_mesh_size):
        refined = Mesh(
            vertices=np.array(mesh.vertices, copy=True),
            faces=np.array(mesh.faces, copy=True),
            markers=np.array(mesh.markers, copy=True),
        )
        refined.refined_tag = True
        return refined, {
            "applied": True,
            "rounds": 1,
            "candidate_faces": 1,
            "candidate_roof_faces": 1,
            "candidate_ground_faces": 0,
            "ground_relief_median": 1.0,
            "ground_refinement_enabled": True,
            "split_edges": 1,
            "added_vertices": 1,
            "added_faces": 2,
            "edge_threshold": 7.5,
        }

    def fake_refine_walls(mesh, *, max_mesh_size):
        return mesh, meshes_module._tetgen_shell_wall_refinement_disabled_stats()

    def fake_tetgen_build(**kwargs):
        mesh = kwargs["mesh"]
        switches = dict(kwargs.get("switches_params") or {})
        switches.update(kwargs.get("switches_overrides") or {})
        calls.append(
            {
                "refined": bool(getattr(mesh, "refined_tag", False)),
                "preserve_surface": bool(switches.get("preserve_surface")),
            }
        )
        volume_mesh = VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )
        volume_mesh.retry_tag = (
            "refined" if getattr(mesh, "refined_tag", False) else "plain"
        )
        return volume_mesh

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        fake_refine_transition,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        fake_refine_horizontal,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_vertical_wall_faces_for_tetgen",
        fake_refine_walls,
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_volume_mesh_quality_snapshot",
        lambda volume_mesh: {
            "aspect_ratio_max": 900.0,
            "element_quality_min": 0.01,
            "min_edge_length": 1.0,
            "high_aspect_ratio_count": 12.0,
            "low_quality_count": 80.0,
        },
    )
    monkeypatch.setattr(
        meshes_module,
        "_capture_tetgen_quality_failure_artifacts",
        lambda **kwargs: {"report": str(tmp_path / "report.json"), "debug_meshes": {}},
    )

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_quality_failure_output_dir=tmp_path,
        tetgen_quality_failure_output_stem="case_055",
        stage_audit=stage_audit,
    )

    assert volume_mesh.retry_tag == "refined"
    assert calls == [
        {"refined": True, "preserve_surface": False},
    ]
    assert stage_audit["selected_attempt_label"] == "attempt-1"
    attempts = {attempt["label"]: attempt for attempt in stage_audit["attempts"]}
    assert "followup_retries" not in attempts["attempt-1"]
    assert attempts["attempt-1"]["config"]["tetgen_shell_refinement_enabled"] is True
    selection = attempts["attempt-1"]["stages"]["surface_shell"][
        "tetgen_shell_horizontal_refinement_selection"
    ]
    assert selection["selected_variant"] == "refined"
    assert selection["reason"] == "tetgen_shell_preconditioned"


def test_build_city_volume_mesh_keeps_shell_refinement_enabled_when_preserving_surface(
    monkeypatch,
):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    stage_audit = {}
    calls: list[dict[str, bool]] = []
    refinement_calls = {"transition": 0, "wall": 0, "horizontal": 0}

    def fake_prepare(*args, **kwargs):
        return make_prepared_city_inputs(
            terrain,
            terrain_raster,
            [conditioned_surface],
            [[0]],
            [4.0],
            {"output_grid": 0.25},
        )

    def fake_prepare_regions(**kwargs):
        return make_prepared_ground_regions(
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_refine_transition(mesh, *, max_mesh_size):
        refinement_calls["transition"] += 1
        return mesh, meshes_module._tetgen_shell_transition_refinement_disabled_stats()

    def fake_refine_horizontal(mesh, *, max_mesh_size):
        refinement_calls["horizontal"] += 1
        return mesh, meshes_module._tetgen_shell_refinement_disabled_stats()

    def fake_refine_walls(mesh, *, max_mesh_size):
        refinement_calls["wall"] += 1
        return mesh, meshes_module._tetgen_shell_wall_refinement_disabled_stats()

    def fake_tetgen_build(**kwargs):
        mesh = kwargs["mesh"]
        calls.append({"refined": bool(getattr(mesh, "refined_tag", False))})
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_ground_building_transition_faces_for_tetgen",
        fake_refine_transition,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_near_horizontal_surface_faces_for_tetgen",
        fake_refine_horizontal,
    )
    monkeypatch.setattr(
        meshes_module,
        "_refine_vertical_wall_faces_for_tetgen",
        fake_refine_walls,
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_audit",
        lambda **kwargs: {"audit_ok": True},
    )
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_plc_contract_from_audit",
        lambda *args, **kwargs: {"errors": [], "warnings": [], "ok": True},
    )
    monkeypatch.setattr(
        meshes_module.tetgen_utils,
        "build_tetgen_plc",
        lambda mesh, *args, **kwargs: _fake_prebuilt_tetgen_plc(mesh),
    )
    monkeypatch.setattr(meshes_module, "_require_tetgen", lambda: None)
    monkeypatch.setattr(
        meshes_module, "get_default_tetgen_switches", _fake_tetgen_switch_defaults
    )
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)
    monkeypatch.setattr(
        meshes_module,
        "_tetgen_volume_mesh_quality_snapshot",
        lambda volume_mesh: {
            "aspect_ratio_max": 1.0,
            "element_quality_min": 1.0,
            "min_edge_length": 1.0,
            "high_aspect_ratio_count": 0.0,
            "low_quality_count": 0.0,
        },
    )

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_switches={"preserve_surface": True},
        stage_audit=stage_audit,
    )

    assert refinement_calls == {"transition": 1, "wall": 0, "horizontal": 1}
    assert calls == [{"refined": False}]
    attempts = {attempt["label"]: attempt for attempt in stage_audit["attempts"]}
    assert attempts["attempt-1"]["config"]["preserve_surface_requested"] is True
    assert attempts["attempt-1"]["config"]["tetgen_shell_refinement_enabled"] is True
    selection = attempts["attempt-1"]["stages"]["surface_shell"][
        "tetgen_shell_horizontal_refinement_selection"
    ]
    assert selection["selected_variant"] == "unrefined"
    assert selection["reason"] == "no_candidate_edges"


def test_build_city_flat_mesh_runs_with_dtcc_mesher():
    pytest.importorskip("dtcc_mesher")

    buildings = [
        make_building(box(8, 8, 18, 18), roof_z=10.0),
        make_building(box(24, 8, 34, 18), roof_z=12.0),
    ]
    city = make_flat_city(buildings)

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert np.any(mesh.markers >= 0)
    assert np.any(mesh.markers == -2)


def test_build_city_mesh_auto_lod_prefers_available_geometry():
    lod1_building = make_building(
        box(8, 8, 18, 18),
        roof_z=10.0,
        lod1_polygon=box(8, 8, 18, 18),
    )
    lod2_building = Building()
    lod2_building.add_geometry(
        make_surface(box(22, 8, 32, 18), 12.0),
        GeometryType.LOD2,
    )
    lod2_building.attributes["estimated_height"] = 12.0
    lod2_building.attributes["ground_height"] = 0.0
    lod0_building = make_building(box(36, 8, 46, 18), roof_z=8.0)

    lod_values = meshes_module._normalize_lod_values(
        [lod1_building, lod2_building, lod0_building],
        None,
    )

    assert lod_values == [
        GeometryType.LOD1,
        GeometryType.LOD2,
        GeometryType.LOD0,
    ]


def test_build_city_flat_mesh_runs_with_auto_lod_resolution():
    city = make_flat_city(
        [
            make_building(box(8, 8, 18, 18), roof_z=10.0),
            make_building(box(24, 8, 34, 18), roof_z=12.0),
        ]
    )

    mesh = build_city_flat_mesh(
        city,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert np.any(mesh.markers >= 0)
    assert np.any(mesh.markers == -2)


def test_build_city_surface_mesh_runs_with_mixed_lod_directives():
    city = make_flat_city(
        [
            make_building(box(8, 8, 18, 18), roof_z=10.0, lod1_polygon=box(8, 8, 18, 18)),
            make_building(box(22, 8, 32, 18), roof_z=12.0, lod2_polygon=box(22, 8, 32, 18)),
        ]
    )

    mesh = build_city_surface_mesh(
        city,
        lod=[GeometryType.LOD1, GeometryType.LOD2],
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        building_mesh_triangle_size=5.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        merge_meshes=True,
        report_mesh_quality=False,
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_build_city_surface_mesh_runs_with_auto_lod_resolution():
    city = make_flat_city(
        [
            make_building(
                box(8, 8, 18, 18),
                roof_z=10.0,
                lod1_polygon=box(8, 8, 18, 18),
            ),
            make_building(box(22, 8, 32, 18), roof_z=12.0),
        ]
    )

    mesh = build_city_surface_mesh(
        city,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        building_mesh_triangle_size=5.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        merge_meshes=True,
        report_mesh_quality=False,
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert float(np.max(mesh.vertices[:, 2])) > 5.0
    assert np.any(np.asarray(mesh.markers, dtype=np.int64) >= 0)


def test_build_city_surface_mesh_keeps_tall_wall_roof_seams_closed():
    city = make_flat_city(
        [
            make_building(
                Polygon(
                    [
                        (20.0, 20.0),
                        (40.0, 20.0),
                        (40.0, 40.0),
                        (32.0, 40.0),
                        (32.0, 38.0),
                        (30.0, 38.0),
                        (30.0, 40.0),
                        (20.0, 40.0),
                        (20.0, 20.0),
                    ]
                ),
                roof_z=80.0,
            )
        ]
    )

    mesh = build_city_surface_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.25,
        building_mesh_triangle_size=4.0,
        max_mesh_size=8.0,
        min_mesh_angle=20.0,
        merge_meshes=True,
        report_mesh_quality=False,
    )

    faces = np.asarray(mesh.faces, dtype=np.int64)
    markers = np.asarray(mesh.markers, dtype=np.int64)
    edge_to_faces = defaultdict(list)
    for face_index, face in enumerate(faces):
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_to_faces[tuple(sorted((int(a), int(b))))].append(face_index)

    open_edge_markers = Counter(
        int(markers[face_indices[0]])
        for face_indices in edge_to_faces.values()
        if len(face_indices) == 1
    )

    assert open_edge_markers
    assert set(open_edge_markers) == {-2}


def test_build_city_surface_mesh_flattens_extruded_building_bases_on_sloped_terrain():
    terrain_mesh = Mesh(
        vertices=np.array(
            [
                [20.0, 20.0, 5.0],
                [40.0, 20.0, 6.0],
                [40.0, 40.0, 10.0],
                [20.0, 40.0, 9.0],
            ],
            dtype=float,
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
        markers=np.array([0, 0], dtype=int),
    )
    building_surface = make_surface(box(20.0, 20.0, 40.0, 40.0), 7.0)

    raw_meshes = _dtcc_builder.build_city_surface_mesh_from_terrain_mesh(
        [create_builder_surface(building_surface)],
        [1],
        mesh_to_builder_mesh(terrain_mesh),
        0,
        True,
    )
    mesh = builder_mesh_to_mesh(raw_meshes[0])

    faces = np.asarray(mesh.faces, dtype=np.int64)
    markers = np.asarray(mesh.markers, dtype=np.int64)
    vertices = np.asarray(mesh.vertices, dtype=np.float64)

    wall_indices = np.flatnonzero(markers == 0)
    roof_indices = np.flatnonzero(markers == 1)

    assert wall_indices.size > 0
    assert roof_indices.size > 0

    wall_vertices = np.unique(faces[wall_indices].reshape(-1))
    roof_vertices = np.unique(faces[roof_indices].reshape(-1))

    assert vertices[wall_vertices, 2].min() < 7.0
    assert vertices[wall_vertices, 2].max() == pytest.approx(7.0)
    assert np.allclose(vertices[roof_vertices, 2], 7.0)


def test_legacy_building_wrappers_still_return_buildings():
    buildings = [
        make_building(box(0, 0, 4, 4), roof_z=8.0),
        make_building(box(4.1, 0, 8.1, 4), roof_z=10.0),
    ]

    merged, merged_map = merge_building_footprints(
        buildings,
        max_distance=0.2,
        min_area=0.0,
        return_index_map=True,
    )
    cleaned, cleaned_map = clean_building_footprints(
        buildings,
        clearance=0.2,
        smallest_hole_area=0.0,
        return_index_map=True,
    )
    fixed, fixed_map = fix_building_footprint_clearance(
        buildings,
        clearance=0.2,
        return_index_map=True,
    )
    simplified, simplified_map = simplify_building_footprints(
        buildings,
        tolerance=0.1,
        return_index_map=True,
    )

    assert merged and merged_map
    assert cleaned and cleaned_map
    assert fixed and fixed_map
    assert simplified and simplified_map


def test_volume_mesh_requires_tetgen_before_preparation(monkeypatch):
    from dtcc_core.builder.meshing import tetgen
    from dtcc_core.datasets.city_volume_mesh import CityVolumeMeshDataset

    cause = ImportError("native TetGen library could not load")
    monkeypatch.setattr(tetgen, "HAS_TETGEN", False)
    monkeypatch.setattr(tetgen, "_tetgen_import_error", cause)

    def unexpected_preparation(*args, **kwargs):
        pytest.fail("Geometry preparation ran without TetGen")

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", unexpected_preparation)
    city = City()
    for invoke in (
        lambda: build_city_volume_mesh(city),
        lambda: city.build_volume_mesh(),
        lambda: CityVolumeMeshDataset().build_from_city(city, bounds=(0, 0, 1, 1)),
        lambda: tetgen.build_volume_mesh(Mesh()),
    ):
        with pytest.raises(ImportError, match="uv sync --extra volume") as exc:
            invoke()
        assert exc.value.__cause__ is cause


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
@pytest.mark.parametrize("entry_point", ["builder", "dataset"])
def test_build_city_volume_mesh_smoke(entry_point):
    city = make_flat_city(
        [
            make_building(box(10, 10, 18, 18), roof_z=10.0),
            make_building(box(24, 10, 32, 18), roof_z=12.0),
        ]
    )

    kwargs = dict(
        lod=GeometryType.LOD0,
        domain_height=40.0,
        max_mesh_size=8.0,
        min_mesh_angle=20.0,
        merge_buildings=True,
        min_building_detail=0.0,
        min_building_area=1.0,
        smoothing=0,
        boundary_face_markers=True,
        report_mesh_quality=False,
    )

    if entry_point == "builder":
        volume_mesh = build_city_volume_mesh(city, **kwargs)
    else:
        from dtcc_core.datasets.city_volume_mesh import CityVolumeMeshDataset

        volume_mesh = CityVolumeMeshDataset().build_from_city(
            city, bounds=(0, 0, 80, 80), **kwargs
        )

    assert volume_mesh.vertices.shape[0] > 0
    assert np.isfinite(volume_mesh.vertices).all()
    assert volume_mesh.cells.shape[0] > 0
    assert volume_mesh.cells.shape[1] == 4
    assert volume_mesh.cells.min() >= 0
    assert volume_mesh.cells.max() < len(volume_mesh.vertices)
    tetrahedra = volume_mesh.vertices[volume_mesh.cells]
    volumes = np.abs(np.linalg.det(tetrahedra[:, 1:] - tetrahedra[:, :1])) / 6
    assert (volumes > 0).all()
    assert volume_mesh.boundary_faces.shape[1] == 3
    assert len(volume_mesh.boundary_faces) == len(volume_mesh.boundary_markers)
    assert set(range(-6, 0)).issubset(set(volume_mesh.boundary_markers))


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_city_volume_mesh_smoke_auto_lod_resolution():
    city = make_flat_city(
        [
            make_building(box(10, 10, 18, 18), roof_z=10.0),
            make_building(box(24, 10, 32, 18), roof_z=12.0),
        ]
    )

    volume_mesh = build_city_volume_mesh(
        city,
        domain_height=40.0,
        max_mesh_size=8.0,
        min_mesh_angle=20.0,
        merge_buildings=True,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        smoothing=0,
        boundary_face_markers=False,
        report_mesh_quality=False,
    )

    assert volume_mesh.vertices.shape[0] > 0
    assert volume_mesh.cells.shape[0] > 0


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_city_volume_mesh_propagates_tetgen_failure(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 18, 18), roof_z=10.0)])
    failure = RuntimeError("TetGen meshing failed")
    calls = []
    stage_audit = {}

    def fail_tetgen(**kwargs):
        calls.append(kwargs)
        raise failure

    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fail_tetgen)
    with pytest.raises(RuntimeError, match="TetGen meshing failed") as exc:
        build_city_volume_mesh(
            city, lod=GeometryType.LOD0, max_mesh_size=8.0, domain_height=40.0,
            min_building_detail=0.0, min_building_area=1.0,
            report_mesh_quality=False, stage_audit=stage_audit,
        )
    assert exc.value is failure
    assert len(calls) == 1
    assert len(stage_audit["attempts"]) == 1
    assert stage_audit["attempts"][0]["result"]["status"] == "failed"
