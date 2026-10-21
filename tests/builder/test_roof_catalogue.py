"""Production catalogue contracts, generic planar shells and city integration."""

import json

import numpy as np
import pytest
from shapely import affinity
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from dtcc_core import builder
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.meshing.backends import available_2d_meshers
from dtcc_core.builder.meshing.roof_interfaces import node_interfaces
from dtcc_core.builder.meshing.roof_validation import (
    validate_partition,
    validate_triangulation,
    validate_shell,
    vertex_link_audit,
)
from dtcc_core.builder.meshing.tetgen import has_volume_veto
from dtcc_core.builder.roofs import (
    RoofFitOptions,
    RoofModel,
    RoofPatch,
    fit_roof,
    supported_shapes,
)
from dtcc_core.builder.roofs.catalogue import REGISTRY
from dtcc_core.model import (
    Bounds,
    Building,
    City,
    GeometryType,
    Mesh,
    PointCloud,
    Raster,
    Surface,
    Terrain,
)


SHAPES = [
    ("flat", [8.0]),
    ("shed", [8.0, 0.15, -0.25]),
    ("gable", [0.0, 8.0, 0.55]),
    ("hip", [0.0, 8.0, 0.55, 0.6]),
]


def model(shape, parameters, footprint, *, key="source", origin=None):
    origin = (
        np.array(footprint.centroid.coords[0]) if origin is None else np.asarray(origin)
    )
    vertices = np.array(footprint.exterior.coords) - origin
    return RoofModel(
        key,
        shape,
        np.array(parameters),
        origin,
        REGISTRY[shape].planes(parameters, vertices),
    )


def cloud(roof, footprint, n=900):
    rng = np.random.default_rng(2026)
    xy = rng.uniform(footprint.bounds[:2], footprint.bounds[2:], size=(n, 2))
    import shapely

    xy = xy[shapely.intersects_xy(footprint, xy[:, 0], xy[:, 1])]
    z = roof.predict(xy) + rng.normal(0, 0.035, len(xy))
    z[: int(0.03 * len(z))] += 3
    return np.c_[xy, z]


def city_of(buildings, bounds=(-25, -20, 45, 25), *, slope=0.0):
    city = City()
    raster = Raster(data=np.tile(np.linspace(0, slope, 20), (20, 1)))
    raster.set_bounds(Bounds(*bounds))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(buildings)
    return city


def building(footprint, roof, *, observations=True):
    value = Building(id=roof.source_key)
    surface = Surface()
    surface.from_polygon(footprint, 10.0)
    value.add_geometry(surface, GeometryType.LOD0)
    value.estimated_height = 10.0
    value.attributes["ground_height"] = 0.0
    if observations:
        value.add_geometry(
            PointCloud(points=cloud(roof, footprint)), GeometryType.POINT_CLOUD
        )
    return value


def assemble(
    patches, backend="dtcc_mesher", *, separate=False, slope=0.0, min_wall_strip=0.5
):
    target = unary_union([p.polygon for p in patches])
    validate_partition(patches, target)
    patches = node_interfaces(patches)
    validate_partition(patches, target)
    bounds = (-25.0, -20.0, 45.0, 25.0)
    ground = box(*bounds).difference(target)
    grounds = list(getattr(ground, "geoms", [ground]))
    regions = [*grounds, *[p.polygon for p in patches]]
    planar, _ = meshes._build_ground_mesh_from_coverage(
        region_polygons=regions,
        region_markers=[-2] * len(grounds) + list(range(len(patches))),
        region_points=[np.array(p.representative_point().coords[0]) for p in regions],
        bounds=bounds,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        mesher=backend,
        add_halo_markers=False,
    )
    validate_triangulation(planar, patches)
    raster = city_of([], bounds, slope=slope).terrain.raster
    surfaces = []
    for p in patches:
        surface = Surface()
        surface.from_polygon(p.polygon, p.coefficients[2])
        surfaces.append(surface)
    owners = {k: i for i, k in enumerate(dict.fromkeys(p.source_key for p in patches))}
    result = meshes._build_city_surface_mesh_from_ground_mesh(
        ground_mesh=planar,
        terrain_raster=raster,
        building_surfaces=surfaces,
        meshing_directives=[1] * len(surfaces),
        smoothing=0,
        merge_meshes=not separate,
        min_wall_strip=min_wall_strip,
        patch_planes=[p.native_plane() for p in patches],
        patch_owners=[owners[p.source_key] for p in patches],
    )
    welded = meshes._welded_mesh(result) if separate else result
    report = validate_shell(
        welded, num_patches=len(patches), closure_audit=meshes._shell_closure_audit
    )
    assert report["valid"], report
    return welded, patches


def test_registry_and_invalid_options():
    assert supported_shapes() == ("flat", "shed", "gable", "hip")
    for options in [
        dict(shapes=("pyramid",)),
        dict(shapes=()),
        dict(min_points=4),
        dict(median_error=-1),
    ]:
        with pytest.raises(ValueError):
            RoofFitOptions(**options)


@pytest.mark.parametrize("shape,parameters", SHAPES)
def test_noisy_outlier_fit(shape, parameters):
    footprint = box(-12, -6, 12, 6)
    truth = model(shape, parameters, footprint)
    points = cloud(truth, footprint)
    fitted = fit_roof(footprint, points, 0.0, source_key="sample")
    assert fitted.accepted, fitted.diagnostics
    assert fitted.model.shape == shape
    assert (
        np.median(
            np.abs(fitted.model.predict(points[:, :2]) - truth.predict(points[:, :2]))
        )
        < 0.08
    )
    json.dumps(fitted.diagnostics, allow_nan=False)


def test_rotation_and_large_translation():
    footprint = box(-10, -5, 10, 5)
    truth = model("gable", [0.0, 8.0, 0.6], footprint)
    points = cloud(truth, footprint)
    angle = 0.63
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    offset = [386000, 6174000]
    world = points[:, :2] @ rotation.T + offset
    transformed = affinity.translate(
        affinity.rotate(footprint, np.rad2deg(angle), origin=(0, 0)), *offset
    )
    fitted = fit_roof(transformed, np.c_[world, points[:, 2]], 0.0)
    assert fitted.accepted
    assert fitted.model.shape == "gable"
    assert (
        np.median(abs(fitted.model.predict(world) - truth.predict(points[:, :2])))
        < 0.05
    )


def test_explicit_sparse_compound_and_ground_rejections():
    foot = box(-10, -5, 10, 5)
    assert (
        "sparse"
        in fit_roof(foot, np.zeros((11, 3)), 0.0).diagnostics["rejection_reasons"][0]
    )
    xy = np.random.default_rng(93).uniform([-10, -5], [10, 5], (600, 2))
    fit = fit_roof(foot, np.c_[xy, np.where(xy[:, 0] < 0, 6.0, 12.0)], 0.0)
    assert not fit.accepted and fit.diagnostics["rejection_reasons"]
    fit = fit_roof(foot, np.c_[xy, np.full(len(xy), -1.0)], 0.0)
    assert not fit.accepted


@pytest.mark.parametrize("backend", ["dtcc_mesher", "triangle"])
@pytest.mark.parametrize("shape,parameters", SHAPES)
def test_generic_planar_shell_with_courtyard(backend, shape, parameters):
    if backend not in available_2d_meshers():
        pytest.skip(f"{backend} backend unavailable")
    foot = box(-12, -6, 12, 6).difference(box(-2, -1, 2, 1))
    roof = model(shape, parameters, foot)
    assembled, patches = assemble(roof.patches(foot), backend, slope=1.0)
    assert len(patches) >= 1
    assert np.all(np.isfinite(assembled.vertices))


@pytest.mark.parametrize("case", ["equal", "step", "crossing", "ridge_t", "three"])
@pytest.mark.parametrize("backend", ["dtcc_mesher", "triangle"])
def test_interfaces(case, backend):
    if backend not in available_2d_meshers():
        pytest.skip(f"{backend} backend unavailable")
    a, b = box(-10, -5, 0, 5), box(0, -5, 10, 5)
    left = model("flat", [8.0], a, key="a")
    right = model("flat", [8.0 if case == "equal" else 10.0], b, key="b")
    if case == "crossing":
        left = model("shed", [8.0, 0.0, 0.5], a, key="a")
        right = model("shed", [8.0, 0.0, -0.5], b, key="b")
    if case == "ridge_t":
        left = model("gable", [0.0, 8.0, 0.6], a, key="a")
    patches = [*left.patches(a), *right.patches(b)]
    if case == "three":
        c = box(-10, 5, 10, 10)
        patches.extend(model("shed", [12.0, 0.1, 0.0], c, key="c").patches(c))
    mesh, patches = assemble(patches, backend, separate=True)
    xyz = mesh.vertices[mesh.faces]
    on_shared = np.all(abs(xyz[:, :, 0]) < 1e-7, axis=1)
    area = (
        np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1)
        / 2
    )
    if case == "equal":
        assert not on_shared.any()
    elif case == "step":
        assert area[on_shared].sum() == pytest.approx(20.0)
    elif case == "crossing":
        assert area[on_shared].sum() == pytest.approx(25.0)
        assert any(np.allclose(p, [0, 0, 8]) for p in mesh.vertices)


@pytest.mark.parametrize("cleaned", [False, True])
@pytest.mark.parametrize("separate", [False, True])
def test_public_mixed_block_source_points_and_fallback(cleaned, separate):
    a, b = box(-10, -5, 0, 5), box(0, -5, 10, 5)
    city = city_of(
        [
            building(a, model("gable", [0.0, 8.0, 0.6], a, key="a")),
            building(b, model("flat", [10.0], b, key="b"), observations=False),
        ]
    )
    if cleaned:
        city = builder.clean_footprints(city)
    audit = {}
    result = builder.build_city_surface_mesh(
        city,
        lod="lod2",
        stage_audit=audit,
        min_building_area=0.0,
        report_mesh_quality=False,
        merge_meshes=not separate,
    )
    rows = audit["roof_catalogue"]["sources"]
    assert len(rows) == 2
    assert {r["outcome"] for r in rows} == {"accepted_fit", "lod1_fallback"}
    assert any(r["shape"] == "gable" for r in rows)
    assert not audit["roof_catalogue"]["block_fallbacks"]
    assert audit["roof_catalogue"]["assemblies"][-1]["shell"]["valid"]
    assert all(
        has_volume_veto(m) == separate for m in (result if separate else [result])
    )
    json.dumps(audit, allow_nan=False)


def test_rejected_api_combinations_and_empty_city():
    city = city_of([])
    for kwargs in [
        dict(lod2_method="other"),
        dict(roof_options=RoofFitOptions()),
        dict(lod="lod2", roof_options={}),
    ]:
        with pytest.raises((ValueError, TypeError)):
            builder.build_city_surface_mesh(city, **kwargs)
    mesh = builder.build_city_surface_mesh(city, lod="lod2", report_mesh_quality=False)
    assert len(mesh.faces)


def test_vertex_links_detect_point_pinch():
    mesh = Mesh(
        vertices=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]),
        faces=np.array([[0, 1, 2], [0, 3, 4]]),
    )
    assert vertex_link_audit(mesh)["nonmanifold_vertices"] == 1


def test_weak_face_support_is_rejected():
    foot = box(-12, -6, 12, 6)
    truth = model("gable", [0.0, 8.0, 0.55], foot)
    points = cloud(truth, foot, n=1800)
    points = points[points[:, 1] < -2.0]
    result = fit_roof(foot, points, 0.0, options=RoofFitOptions(shapes=("gable",)))
    assert not result.accepted
    assert any("support" in s for s in result.diagnostics["rejection_reasons"])


def test_partition_transfer_to_disconnected_cells_and_changed_boundary():
    original = box(-12, -6, 12, 6)
    roof = model("hip", [0.0, 8.0, 0.5, 0.7], original)
    target = unary_union([box(-13, -7, -2, 7), box(2, -7, 13, 7)])
    patches = roof.patches(target)
    validate_partition(patches, target)
    for patch in patches:
        xy = np.asarray(patch.polygon.exterior.coords)
        np.testing.assert_allclose(patch.height(xy), roof.predict(xy), atol=1e-8)


def test_constant_plane_assembly_matches_legacy_geometry():
    a, b = box(-10, -5, 0, 5), box(0, -5, 10, 5)
    city = city_of(
        [
            building(a, model("flat", [10.0], a, key="a"), observations=False),
            building(b, model("flat", [10.0], b, key="b"), observations=False),
        ]
    )
    legacy = builder.build_city_surface_mesh(city, report_mesh_quality=False)
    audit = {}
    flat = builder.build_city_surface_mesh(
        city, lod="lod2", report_mesh_quality=False, stage_audit=audit
    )

    def area_volume(mesh):
        xyz = mesh.vertices[mesh.faces]
        cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        return np.linalg.norm(cross, axis=1).sum() / 2, np.sum(
            cross[:, 2] * xyz[:, :, 2].mean(axis=1)
        ) / 2

    assert area_volume(flat) == pytest.approx(area_volume(legacy))
    assert not audit["roof_catalogue"]["block_fallbacks"]


def test_public_domain_clipping_and_excluded_source_preserve_mapping():
    a, b = box(-30, -5, -10, 5), box(5, 0, 6, 1)
    city = city_of(
        [
            building(a, model("shed", [12.0, 0.1, 0.2], a, key="clip")),
            building(b, model("flat", [10.0], b, key="excluded")),
        ]
    )
    audit = {}
    result = city.build_surface_mesh(
        lod="lod2", stage_audit=audit, min_building_area=15.0
    )
    rows = audit["roof_catalogue"]["sources"]
    assert len(rows) == 1 and rows[0]["source_id"] == "clip"
    assert result.vertices[:, 0].min() >= -25 - 1e-8
    assert audit["roof_catalogue"]["assemblies"][-1]["shell"]["valid"]


def test_partition_failure_has_no_flattening_retry(monkeypatch):
    from dtcc_core.builder.meshing import catalogue_surface
    from dtcc_core.builder.meshing.roof_interfaces import RoofAssemblyError

    foot = box(-10, -5, 10, 5)
    city = city_of([building(foot, model("hip", [0.0, 8.0, 0.55, 0.6], foot))])
    calls = []

    def reject(patches, bounds, **kwargs):
        calls.append(len(patches))
        raise RoofAssemblyError(
            "Synthetic narrow-strip conditioning failure", block_ids=[0]
        )

    monkeypatch.setattr(catalogue_surface, "local_roof_coverage", reject)
    audit = {}
    with pytest.raises(RoofAssemblyError, match="roof_partition"):
        builder.build_city_surface_mesh(
            city, lod="lod2", stage_audit=audit, report_mesh_quality=False
        )
    assert len(calls) == 1
    assert not audit["roof_catalogue"]["block_fallbacks"]
    assert len(audit["roof_catalogue"]["assemblies"]) == 1
    assert not audit["roof_catalogue"]["assemblies"][0]["valid"]


def test_mesh_failure_preserves_admitted_source_fit(monkeypatch):
    import importlib

    compact = importlib.import_module("dtcc_core.builder.meshing.compact_surface")
    foot = box(-10, -5, 10, 5)
    city = city_of([building(foot, model("gable", [0.0, 8.0, 0.6], foot))])
    calls = []

    def reject(*args, **kwargs):
        calls.append(True)
        raise ValueError("Controlled triangulation failure")

    monkeypatch.setattr(compact, "triangulate", reject)
    audit = {}
    with pytest.raises(RuntimeError, match="shared_geometry_and_meshing"):
        builder.build_city_surface_mesh(
            city, lod="lod2", stage_audit=audit, report_mesh_quality=False
        )
    assert len(calls) == 1
    assert audit["roof_catalogue"]["sources"][0]["outcome"] == "accepted_fit"
    assert not audit["roof_catalogue"]["block_fallbacks"]


@pytest.mark.parametrize("buried", [False, True])
def test_buried_wall_audit_evaluates_sloped_triangle_at_probe(buried):
    # Roof above x>=0 has z=3 at the wall probe but z=13 at its far end.
    # A wall centred at z=5 is exposed, despite that far maximum.
    height = 2.0 if buried else 5.0
    mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, height - 1],
                [0.0, 2.0, height - 1],
                [0.0, 2.0, height + 1],
                [0.0, 0.0, 3.0],
                [10.0, 0.0, 13.0],
                [0.0, 10.0, 3.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [3, 4, 5]]),
        markers=np.array([0, 2]),
    )
    report = meshes._shell_closure_audit(mesh, num_buildings=2)
    assert report["buried_wall_faces"] == int(buried)


def test_cleaned_single_without_original_outline_records_evidence_fallback():
    foot = box(-10, -5, 10, 5)
    city = builder.clean_footprints(
        city_of([building(foot, model("gable", [0.0, 8.0, 0.6], foot))])
    )
    city.buildings[0].attributes.pop("source_footprint")
    audit = {}
    builder.build_city_surface_mesh(
        city, lod="lod2", stage_audit=audit, report_mesh_quality=False
    )
    row = audit["roof_catalogue"]["sources"][0]
    assert row["outcome"] == "lod1_fallback"
    assert (
        "original source footprint unavailable"
        in row["diagnostics"]["rejection_reasons"][0]
    )


@pytest.mark.parametrize("backend", ["dtcc_mesher", "triangle"])
def test_unequal_boundary_nodes_and_numerically_coincident_profiles(backend):
    if backend not in available_2d_meshers():
        pytest.skip(f"{backend} backend unavailable")
    left = Polygon([(-10, -5), (0, -5), (0, -2), (0, 1), (0, 5), (-10, 5)])
    right = box(0, -5, 10, 5)
    patches = [
        RoofPatch(left, np.zeros(2), np.array([0.0, 0.1, 8.0]), "a"),
        RoofPatch(right, np.zeros(2), np.array([0.0, 0.1, 8.0 + 1e-9]), "b"),
    ]
    mesh, _ = assemble(patches, backend)
    xyz = mesh.vertices[mesh.faces]
    assert not np.all(abs(xyz[:, :, 0]) < 1e-7, axis=1).any()


def test_narrow_mandatory_patch_has_explicit_rejection():
    from dtcc_core.builder.meshing.roof_interfaces import RoofAssemblyError

    patch = RoofPatch(
        box(0, 0, 1e-10, 10), np.zeros(2), np.array([0.0, 0.0, 8.0]), "thin"
    )
    with pytest.raises(RoofAssemblyError, match="collapses"):
        node_interfaces([patch])


def test_point_touching_shells_cannot_pass_validation():
    a, b = box(-10, -10, 0, 0), box(0, 0, 10, 10)
    patches = [
        *model("flat", [8.0], a, key="a").patches(a),
        *model("flat", [8.0], b, key="b").patches(b),
    ]
    with pytest.raises(AssertionError, match="nonmanifold"):
        assemble(patches)


def test_discarded_cleaner_cells_are_not_recovered_by_proximity(monkeypatch):
    from dtcc_core.builder.meshing import catalogue_surface

    foot = box(-10, -5, 10, 5)
    city = city_of([building(foot, model("gable", [0.0, 8.0, 0.6], foot))])
    original = catalogue_surface.clean_regions

    def discard(*args, **kwargs):
        cleaned = original(*args, **kwargs)
        cleaned.output_cells = [[] for _ in cleaned.output_cells]
        return cleaned

    monkeypatch.setattr(catalogue_surface, "clean_regions", discard)
    audit = {}
    builder.build_city_surface_mesh(
        city, lod="lod2", stage_audit=audit, report_mesh_quality=False
    )
    row = audit["roof_catalogue"]["sources"][0]
    assert row["source_keys"] == []
    assert "source partition unavailable" in row["diagnostics"]["rejection_reasons"][0]
    assert audit["roof_catalogue"]["assemblies"][-1]["shell"]["valid"]


def test_parent_pooled_cloud_is_not_used_for_source_parts():
    a, b = box(-10, -5, 0, 5), box(0, -5, 10, 5)
    ar, br = (
        model("shed", [8.0, 0.1, 0.2], a, key="a"),
        model("flat", [10.0], b, key="b"),
    )
    city = builder.clean_footprints(
        city_of(
            [building(a, ar, observations=False), building(b, br, observations=False)]
        )
    )
    city.buildings[0].add_geometry(
        PointCloud(points=np.vstack([cloud(ar, a), cloud(br, b)])),
        GeometryType.POINT_CLOUD,
    )
    audit = {}
    builder.build_city_surface_mesh(
        city, lod="lod2", stage_audit=audit, report_mesh_quality=False
    )
    rows = audit["roof_catalogue"]["sources"]
    assert len(rows) == 2
    assert all(row["points"] == 0 and row["outcome"] == "lod1_fallback" for row in rows)
    assert all(
        any(
            "parent pooling" in reason
            for reason in row["diagnostics"]["rejection_reasons"]
        )
        for row in rows
    )


def test_invalid_ground_is_an_explicit_fit_fallback():
    foot = box(-10, -5, 10, 5)
    source = building(foot, model("gable", [0.0, 8.0, 0.6], foot))
    source.attributes["ground_height"] = None
    audit = {}
    builder.build_city_surface_mesh(
        city_of([source]), lod="lod2", stage_audit=audit, report_mesh_quality=False
    )
    assert audit["roof_catalogue"]["sources"][0]["diagnostics"][
        "rejection_reasons"
    ] == ["missing ground height"]


def test_registry_module_contract_requires_no_meshing_dispatch(monkeypatch):
    monkeypatch.setitem(REGISTRY, "contract-test", REGISTRY["flat"])
    foot = box(-10, -5, 10, 5)
    truth = model("flat", [8.0], foot)
    fitted = fit_roof(
        foot, cloud(truth, foot), 0.0, options=RoofFitOptions(shapes=("contract-test",))
    )
    assert fitted.accepted and fitted.model.shape == "contract-test"
    assemble(fitted.model.patches(foot))


def test_generic_wall_subdivision_has_a_bounded_memory_budget():
    foot = box(-10, -5, 10, 5)
    roof = model("flat", [100000.0], foot)
    with pytest.raises(RuntimeError, match="wall subdivision budget"):
        assemble(roof.patches(foot), min_wall_strip=0.0)


@pytest.mark.parametrize("cleaned", [False, True])
@pytest.mark.parametrize("levels", [(8.0, 14.0), (8.0, 8.2)])
def test_insufficient_evidence_flat_models_match_ordinary_lod1_geometry(
    monkeypatch, cleaned, levels
):
    from dtcc_core.builder.meshing import catalogue_surface

    a, b = box(-10, -5, 0, 5), box(0, -5, 10, 5)
    buildings = [
        building(a, model("gable", [0.0, 8.0, 0.5], a, key="low")),
        building(b, model("gable", [0.0, 14.0, 0.5], b, key="high")),
    ]
    for source, level in zip(buildings, levels):
        source.lod0.set_z(level)
        source.estimated_height = level
    city = city_of(buildings)
    if cleaned:
        city = builder.clean_footprints(city)
    ordinary = builder.build_city_surface_mesh(city, report_mesh_quality=False)
    for parent in city.buildings:
        for source in parent.building_parts or [parent]:
            source.add_geometry(
                PointCloud(points=np.empty((0, 3))), GeometryType.POINT_CLOUD
            )
    audit = {}
    fallback = builder.build_city_surface_mesh(
        city, lod="lod2", stage_audit=audit, report_mesh_quality=False
    )
    assert len(audit["roof_catalogue"]["assemblies"]) == 1
    assert not audit["roof_catalogue"]["block_fallbacks"]
    assert audit["roof_catalogue"]["retained_fits"] == 0

    # Compare the full piecewise planar surface, independent of triangulation.
    def area_volume(mesh):
        xyz = mesh.vertices[mesh.faces]
        cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        return np.linalg.norm(cross, axis=1).sum() / 2, np.sum(
            cross[:, 2] * xyz[:, :, 2].mean(axis=1)
        ) / 2

    assert area_volume(fallback) == pytest.approx(
        area_volume(ordinary), rel=1e-10, abs=1e-7
    )
    patches = audit["roof_catalogue"]["patch_markers"]
    assert {p["source_key"] for p in patches} == {"low", "high"}
    expected = {
        float(z)
        for z in ordinary.vertices[
            ordinary.faces[
                np.asarray(ordinary.markers)
                >= ((int(np.max(ordinary.markers)) + 1) // 2)
            ]
        ][:, :, 2].ravel()
    }
    assert {p["coefficients"][2] for p in patches} == expected


@pytest.mark.parametrize("backend", ["dtcc_mesher", "triangle"])
def test_lund_coverage_roundoff_regression_keeps_roofs_and_courtyards(backend):
    """Reduced coverage with unequal roundoff along a shared sloping boundary."""
    from pathlib import Path
    from shapely.geometry import shape
    from dtcc_core.builder.meshing.roof_coverage import local_roof_coverage

    if backend not in available_2d_meshers():
        pytest.skip(f"{backend} backend unavailable")
    polygons = [
        affinity.translate(shape(p), 386340, 6174560)
        for p in json.loads(
            (Path(__file__).parent / "fixtures/roof_coverage_roundoff.json").read_text()
        )
    ]
    patches = [
        RoofPatch(p, np.zeros(2), np.array([0.0, 0.0, 8.0]), f"cell-{i}", i, i)
        for i, p in enumerate(polygons)
    ]
    local, world, ground, bounds, offset = local_roof_coverage(
        patches, (386075, 6174447, 386575, 6174947)
    )
    assert len(local) == len(patches)
    # A real courtyard remains; no nanometre-wide false ground regions remain.
    assert len(ground) > 1 and min(p.area for p in ground) > 1e-4
    regions = [*ground, *[p.polygon for p in local]]
    mesh, _ = meshes._build_ground_mesh_from_coverage(
        region_polygons=regions,
        region_markers=[-2] * len(ground) + list(range(len(local))),
        region_points=[np.array(p.representative_point().coords[0]) for p in regions],
        bounds=bounds,
        max_mesh_size=8.0,
        min_mesh_angle=25.0,
        mesher=backend,
        coordinate_tolerance=1e-8,
        max_steiner_points=50000,
        add_halo_markers=False,
    )
    validate_triangulation(mesh, local)
    assert set(mesh.markers[mesh.markers >= 0]) == set(range(len(patches)))


def test_unattributed_backend_error_does_not_flatten_the_city(monkeypatch):
    foot = box(-10, -5, 10, 5)
    city = city_of([building(foot, model("gable", [0.0, 8.0, 0.6], foot))])

    def fail(**kwargs):
        raise RuntimeError("Unknown backend failure")

    monkeypatch.setattr(meshes, "_build_ground_mesh_from_coverage", fail)
    audit = {}
    with pytest.raises(RuntimeError, match="terrain_attachment"):
        builder.build_city_surface_mesh(
            city, lod="lod2", stage_audit=audit, report_mesh_quality=False
        )
    assert not audit["roof_catalogue"]["block_fallbacks"]


def test_coverage_retains_distinct_submicrometre_vertices():
    from dtcc_core.builder.meshing.roof_coverage import local_roof_coverage

    polygon = Polygon([(0, 0), (10, 0), (10, 1e-7), (10, 10), (0, 10)])
    patch = RoofPatch(
        polygon, np.zeros(2), np.array([0.0, 0.0, 8.0]), "short-edge", 0, 17
    )
    local, world, *_ = local_roof_coverage([patch], (-5, -5, 15, 15))
    assert len(local[0].polygon.exterior.coords) == 6
    assert len(world[0].polygon.exterior.coords) == 6


@pytest.mark.parametrize("backend", ["dtcc_mesher", "triangle"])
@pytest.mark.parametrize("separate", [False, True])
def test_optional_native_failure_is_explicit_and_never_flattens(
    monkeypatch, backend, separate
):
    if backend not in available_2d_meshers():
        pytest.skip(f"{backend} unavailable")
    import importlib

    compact = importlib.import_module("dtcc_core.builder.meshing.compact_surface")
    improvement = importlib.import_module(
        "dtcc_core.builder.meshing.compact_improvement"
    )
    footprint = box(-10, -5, 10, 5)
    city = city_of([building(footprint, model("hip", [0.0, 8.0, 0.5, 0.6], footprint))])
    observations = city.buildings[0].point_cloud.points.copy()
    original = improvement.improve

    def reject(*args, **kwargs):
        if kwargs.get("max_added_points", 0):
            raise RuntimeError("controlled native failure")
        return original(*args, **kwargs)

    if backend == "dtcc_mesher":
        monkeypatch.setattr(improvement, "improve", reject)
    else:

        def reject_mesh(*args, **kwargs):
            raise RuntimeError("controlled native failure")

        monkeypatch.setattr(compact, "triangulate", reject_mesh)
    audit = {}
    with pytest.raises(RuntimeError, match="controlled native failure"):
        builder.build_city_surface_mesh(
            city,
            lod="lod2",
            mesher=backend,
            merge_meshes=not separate,
            stage_audit=audit,
            report_mesh_quality=False,
        )
    assert not audit["roof_catalogue"]["block_fallbacks"]
    assert audit["roof_catalogue"]["sources"][0]["outcome"] == "accepted_fit"
    np.testing.assert_array_equal(city.buildings[0].point_cloud.points, observations)
