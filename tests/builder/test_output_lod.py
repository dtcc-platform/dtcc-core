"""Output LOD selection and a common meshing path for both reconstructors."""

import importlib

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import builder
from dtcc_core.model import GeometryType
from dtcc_core.builder.roofs.ransac import fit_ransac, RansacRoofModel
from dtcc_core.builder.meshing.roof_validation import validate_partition
from .test_compact_integration import sample
from .test_roof_catalogue import model, cloud, building


def arrays(mesh):
    return mesh.vertices, mesh.faces, mesh.markers


@pytest.mark.parametrize("entry", ["city", "builder"])
def test_lod2_defaults_to_catalogue_and_matches_explicit_method(entry):
    results = []
    for kwargs in (
        {"lod": "lod2"},
        {"lod": GeometryType.LOD2, "lod2_method": "catalogue"},
    ):
        city = sample()
        audit = {}
        call = (
            city.build_surface_mesh
            if entry == "city"
            else lambda **kw: builder.build_city_surface_mesh(city, **kw)
        )
        mesh = call(
            **kwargs, mesher="dtcc_mesher", stage_audit=audit, report_mesh_quality=False
        )
        assert audit["requested_lod"] == "lod2" and audit["lod2_method"] == "catalogue"
        assert audit["roof_catalogue"]["retained_nonflat_fits"] == 1
        results.append(arrays(mesh))
    assert all(np.array_equal(a, b) for a, b in zip(*results))


def test_output_lod1_matches_default_flat_mesh():
    results = [
        builder.build_city_surface_mesh(
            sample(), mesher="dtcc_mesher", report_mesh_quality=False, **kwargs
        )
        for kwargs in ({}, {"lod": "lod1"}, {"lod": GeometryType.LOD1})
    ]
    assert all(
        np.array_equal(a, b)
        for mesh in results[1:]
        for a, b in zip(arrays(results[0]), arrays(mesh))
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"lod": "lod0"},
        {"lod": GeometryType.LOD0},
        {"lod": [GeometryType.LOD1]},
        {"lod": "lod3"},
        {"lod": 2},
        {"lod": "lod2", "lod2_method": "unknown"},
        {"lod": "lod1", "lod2_method": "ransac"},
        {"lod": None},
        {
            "lod": "lod2",
            "lod2_method": "ransac",
            "roof_options": builder.RoofFitOptions(),
        },
    ],
)
def test_invalid_or_conflicting_options_fail_before_city_work(kwargs):
    with pytest.raises(ValueError):
        builder.build_city_surface_mesh(None, **kwargs)


@pytest.mark.parametrize("entry", ["city", "builder"])
def test_ransac_selection_uses_existing_algorithm_and_common_shell(monkeypatch, entry):
    old = importlib.import_module("dtcc_core.builder.geometry_builders.lod2")
    original = old._candidate_lod2
    calls = []

    def recorded(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)

    monkeypatch.setattr(old, "_candidate_lod2", recorded)
    city = sample()
    audit = {}
    before = city.buildings[0].point_cloud.points.copy()
    call = (
        city.build_surface_mesh
        if entry == "city"
        else lambda **kw: builder.build_city_surface_mesh(city, **kw)
    )
    mesh = call(
        lod="lod2",
        lod2_method="ransac",
        mesher="dtcc_mesher",
        stage_audit=audit,
        report_mesh_quality=False,
    )
    assert calls == [True]
    assert audit["lod2_method"] == "ransac"
    cat = audit["roof_catalogue"]
    assert cat["sources"][0]["diagnostics"]["method"] == "ransac"
    assert cat["sources"][0]["final_outcome"] == "retained_fit"
    assert not cat["block_fallbacks"]
    assert cat["compact_surface"]["topology"]["valid"]
    assert cat["compact_surface"]["shared_boundary_mismatches"] == 0
    assert np.array_equal(before, city.buildings[0].point_cloud.points)
    assert city.buildings[0].lod2 is None
    assert len(mesh.faces)


@pytest.mark.parametrize(
    "shape,parameters",
    [("flat", [8.0]), ("shed", [8.0, 0.2, -0.1]), ("gable", [0.0, 8.0, 0.55])],
)
def test_ransac_adapter_preserves_legacy_roof_planes_and_partition(shape, parameters):
    from dtcc_core.builder.geometry_builders.lod2 import _candidate_lod2

    polygon = box(-10, -5, 10, 5)
    truth = model(shape, parameters, polygon)
    source = building(polygon, truth)
    points = source.point_cloud.points
    shell = _candidate_lod2(source, 0.0, False)
    fit = fit_ransac(polygon, points, 0.0, source_key="test")
    if shell is None:
        assert not fit.accepted and fit.diagnostics["rejection_reasons"]
        return
    assert fit.accepted, fit.diagnostics
    validate_partition(fit.model.patches(polygon), polygon)
    for surface in shell.surfaces:
        xyz = surface.vertices
        if abs(surface.calculate_normal()[2]) > 1e-8 and not np.allclose(
            xyz[:, 2], 0.0
        ):
            np.testing.assert_allclose(
                fit.model.predict(xyz[:, :2]), xyz[:, 2], atol=1e-6
            )
    expanded = polygon.buffer(0.4, join_style=2)
    patches = fit.model.patches(expanded)
    validate_partition(patches, expanded)
    assert all(
        any(np.array_equal(p.coefficients, q) for q in fit.model.planes)
        for p in patches
    )


def test_stepped_roof_transfer_keeps_domain_ownership_instead_of_lower_envelope():
    domains = (box(-10, -5, 0, 5), box(0, -5, 10, 5))
    planes = np.array([[0.0, 0.0, 5.0], [0.0, 0.0, 9.0]])
    roof = RansacRoofModel(
        "source", "ransac", planes.ravel(), np.zeros(2), planes, {}, domains
    )
    target = box(-11, -6, 11, 6)
    patches = roof.patches(target)
    validate_partition(patches, target)
    assert sorted(round(p.polygon.area, 6) for p in patches) == [132.0, 132.0]
    for patch in patches:
        x, y = patch.polygon.representative_point().coords[0]
        assert patch.height([x, y]) == (5.0 if x < 0 else 9.0)


def test_ransac_rejection_reports_lod1_fallback():
    city = sample()
    city.buildings[0].remove_geometry(GeometryType.POINT_CLOUD)
    audit = {}
    builder.build_city_surface_mesh(
        city,
        lod="lod2",
        lod2_method="ransac",
        mesher="dtcc_mesher",
        stage_audit=audit,
        report_mesh_quality=False,
    )
    source = audit["roof_catalogue"]["sources"][0]
    assert source["outcome"] == "lod1_fallback"
    assert "insufficient_roof_points" in source["diagnostics"]["rejection_reasons"]


def test_ransac_internal_step_is_meshed_with_original_source_ownership(monkeypatch):
    from dtcc_core.builder.roofs import RoofFitResult

    adapter = importlib.import_module("dtcc_core.builder.roofs.ransac")

    def stepped(footprint, points, ground, *, source_key):
        planes = np.array([[0.0, 0.0, 5.0], [0.0, 0.0, 9.0]])
        roof = RansacRoofModel(
            source_key,
            "ransac",
            planes.ravel(),
            np.zeros(2),
            planes,
            {},
            (box(-10, -5, 0, 5), box(0, -5, 10, 5)),
        )
        return RoofFitResult(roof, dict(accepted=True, method="ransac"))

    monkeypatch.setattr(adapter, "fit_ransac", stepped)
    audit = {}
    mesh = builder.build_city_surface_mesh(
        sample(),
        lod="lod2",
        lod2_method="ransac",
        mesher="dtcc_mesher",
        stage_audit=audit,
        report_mesh_quality=False,
    )
    cat = audit["roof_catalogue"]
    assert cat["sources"][0]["final_outcome"] == "retained_fit"
    assert not cat["block_fallbacks"] and cat["compact_surface"]["topology"]["valid"]
    xyz = mesh.vertices[mesh.faces]
    step = (
        np.all(abs(xyz[:, :, 0]) < 1e-7, axis=1)
        & np.all(xyz[:, :, 2] >= 5.0 - 1e-7, axis=1)
        & np.all(xyz[:, :, 2] <= 9.0 + 1e-7, axis=1)
    )
    area = (
        np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1)
        / 2
    )
    assert area[step].sum() == pytest.approx(40.0)
    source = cat["sources"][0]["source_key"]
    assert all(
        r.attributes["source_key"] == source
        for r in mesh.regions
        if r.attributes["kind"] != "terrain"
    )
