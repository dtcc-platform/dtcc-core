"""Small source-preserving reproductions of LOD2 stage contract failures."""

import json
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely import affinity

from dtcc_core.builder.roofs import RoofPatch
from dtcc_core.builder.meshing.compact_geometry import (
    building_domains,
    face_coordinates,
)


@pytest.mark.parametrize("offset", [(0.0, 0.0), (386000.0, 6174000.0)])
def test_coplanar_point_contact_does_not_create_zero_length_edges(offset):
    evidence = json.loads(
        (
            Path(__file__).parents[1] / "data/lod2/coplanar_point_contact.json"
        ).read_text()
    )
    patches = [
        RoofPatch(
            affinity.translate(shapely.from_wkt(p["wkt"]), *-np.asarray(offset)),
            np.asarray(p["origin"]) - offset,
            np.asarray(p["coefficients"]),
            p["source_key"],
            p["patch_id"],
            p["block_id"],
        )
        for p in evidence["patches"]
    ]
    faces = building_domains(patches, None, base_z=evidence["terrain_base"])
    roofs = [f for f in faces if f.kind == "roof"]
    assert len(roofs) == 2
    for face in faces:
        for ring in face_coordinates(face):
            assert np.linalg.norm(ring - np.roll(ring, 1, axis=0), axis=1).min() > 1e-8
    assert sum(f.polygon.area for f in roofs) == pytest.approx(
        sum(p.polygon.area for p in patches)
    )


@pytest.mark.parametrize("offset", [(0.0, 0.0), (647700.0, 6638400.0)])
def test_source_junction_rule_bounds_changes_and_removes_nonmanifold_edge(offset):
    from dataclasses import replace
    from shapely.ops import unary_union
    from dtcc_core.builder.cleaning.subdivision import canonical_source_subdivision
    from dtcc_core.builder.meshing.roof_junctions import regularize_roof_junctions
    from dtcc_core.builder.meshing.compact_surface import (
        compact_surface,
        CompactSurfaceOptions,
    )
    from dtcc_core.builder.meshing.roof_interfaces import node_interfaces

    evidence = json.loads(
        (
            Path(__file__).parents[1] / "data/lod2/multiway_source_junction.json"
        ).read_text()
    )
    patches = [
        RoofPatch(
            affinity.translate(shapely.from_wkt(p["wkt"]), *-np.asarray(offset)),
            np.asarray(p["origin"]) - offset,
            np.asarray(p["coefficients"]),
            p["source_key"],
            p["patch_id"],
            p["block_id"],
        )
        for p in evidence["patches"]
    ]
    # The premeshing stage receives conditioned coordinates in production.
    frame = np.floor(np.min([p.polygon.bounds[:2] for p in patches], axis=0))
    patches = [
        replace(
            p, polygon=affinity.translate(p.polygon, *-frame), origin=p.origin - frame
        )
        for p in patches
    ]
    options = CompactSurfaceOptions(budgeted=False)
    with pytest.raises(ValueError, match="topology"):
        compact_surface(None, node_interfaces(patches), options, base_z=0.0)
    cells = [(p.polygon, i) for i, p in enumerate(patches)]
    keys = {i: p.source_key for i, p in enumerate(patches)}
    domain = unary_union([p for p, _ in cells])
    _, canonical, normalization = canonical_source_subdivision(
        cells, domain=domain, source_keys=keys
    )
    normalized = node_interfaces([replace(patches[i], polygon=p) for p, i in canonical])
    changed, audits = regularize_roof_junctions(
        normalized, bases={p.block_id: 0.0 for p in normalized}, detail=0.5
    )
    repeated, repeated_audits = regularize_roof_junctions(
        normalized, bases={p.block_id: 0.0 for p in normalized}, detail=0.5
    )
    assert audits == repeated_audits
    assert len(changed) == len(repeated)
    for first, second in zip(changed, repeated):
        assert (first.source_key, first.patch_id, first.block_id) == (
            second.source_key,
            second.patch_id,
            second.block_id,
        )
        assert first.polygon.equals_exact(second.polygon, tolerance=0.0)
        np.testing.assert_array_equal(first.origin, second.origin)
        np.testing.assert_array_equal(first.coefficients, second.coefficients)
    assert any(a["changes"] for a in audits.values())
    assert (
        unary_union([p.polygon for p in changed])
        .symmetric_difference(unary_union([p.polygon for p in patches]))
        .area
        <= normalization["normalization_area_budget_m2"]
    )
    assert {p.source_key for p in changed} == set(keys.values())
    for p in changed:
        original = unary_union(
            [q.polygon for q in patches if q.source_key == p.source_key]
        )
        owned = unary_union(
            [q.polygon for q in changed if q.source_key == p.source_key]
        )
        assert owned.boundary.hausdorff_distance(original.boundary) <= 0.25 + 2e-7
    mesh, result = compact_surface(None, changed, options, base_z=0.0)
    assert result["topology"]["valid"]
    assert result["building_triangles"] <= 1.5 * result["T0"]
    assert len(mesh.faces) > 0


@pytest.mark.parametrize("offset", [(0.0, 0.0), (647500.0, 6638200.0)])
def test_rounded_gable_outline_does_not_create_a_false_roof_step(offset):
    from dtcc_core.builder.meshing.roof_interfaces import node_interfaces
    from dtcc_core.builder.meshing.compact_surface import (
        compact_surface,
        CompactSurfaceOptions,
    )

    evidence = json.loads(
        (Path(__file__).parents[1] / "data/lod2/rounded_gable_outline.json").read_text()
    )
    patches = [
        RoofPatch(
            affinity.translate(shapely.from_wkt(p["wkt"]), *-np.asarray(offset)),
            np.asarray(p["origin"]) - offset,
            np.asarray(p["coefficients"]),
            p["source_key"],
            p["patch_id"],
            p["block_id"],
        )
        for p in evidence["patches"]
    ]
    domain = affinity.translate(
        shapely.from_wkt(evidence["domain"]), *-np.asarray(offset)
    )
    noded = node_interfaces(patches, domains={0: domain})
    mesh, audit = compact_surface(
        None, noded, CompactSurfaceOptions(budgeted=False), base_z=0.0
    )
    assert audit["topology"]["valid"]
    assert audit["building_triangles"] == 14
    # This input has a continuous gable, with no internal step panel.
    assert len([face for face in audit["faces"] if face["kind"] == "wall"]) == 4
    assert len(evidence["observations"]) == 38
    for original in patches:
        for patch in noded:
            if patch.patch_id == original.patch_id:
                np.testing.assert_array_equal(patch.coefficients, original.coefficients)
