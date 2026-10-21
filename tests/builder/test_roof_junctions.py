"""Height-sector topology is established before bounded surface meshing."""
from dataclasses import replace

import numpy as np
import pytest
import shapely
from shapely import affinity
from shapely.geometry import box

from dtcc_core.builder.roofs import RoofPatch
from dtcc_core.builder.meshing.roof_interfaces import node_interfaces
from dtcc_core.builder.meshing.roof_junctions import (
    height_interval_topology,junction_half_widths,regularize_roof_junctions)
from dtcc_core.builder.meshing.compact_surface import compact_surface,CompactSurfaceOptions


def quadrants(heights,offset=(0.,0.),slopes=None):
    domains=[box(0,0,10,10),box(-10,0,0,10),box(-10,-10,0,0),box(0,-10,10,0)]
    return [RoofPatch(affinity.translate(p,*offset),np.asarray(offset),
        np.array([*(slopes[i] if slopes is not None else (0.,0.)),z]),str(i),i,0)
        for i,(p,z) in enumerate(zip(domains,heights)) if z is not None]


@pytest.mark.parametrize('heights,valid', [([8,8,8,8],True),([4,6,8,10],True),
    ([8,4,8,4],False),([8,4,8,0],False),([8,8,4,0],True),
    ([8,4,8,8],True),([8,4,6,2],False)])
def test_cyclic_height_intervals_distinguish_manifold_stars(heights,valid):
    assert height_interval_topology(heights)['manifold'] is valid


@pytest.mark.parametrize('offset',[(0.,0.),(512000.,6570000.)])
@pytest.mark.parametrize('heights',[[8,8,8,8],[4,6,8,10],[8,4,8,8],[8,8,4,None]])
def test_harmless_multiway_junctions_retain_geometry(heights,offset):
    patches=node_interfaces(quadrants(heights,offset))
    # Use local coordinates for the mesher, as the production pipeline does.
    local=[replace(p,polygon=affinity.translate(p.polygon,*-np.asarray(offset)),
                   origin=p.origin-offset) for p in patches]
    result,audit=regularize_roof_junctions(local,bases={0:0.},detail=.5)
    assert not audit[0]['changes']
    assert len(audit[0]['checked_junctions'])==1
    assert all(a.polygon.equals_exact(b.polygon,0.) for a,b in zip(local,result))
    mesh,measured=compact_surface(None,result,CompactSurfaceOptions(budgeted=False),base_z=0.)
    assert measured['topology']['valid']
    assert len(mesh.faces)>0


@pytest.mark.parametrize('heights',[[8,4,8,4],[8,4,8,None]])
def test_alternating_stars_are_changed_once_with_bounded_ownership(heights):
    patches=node_interfaces(quadrants(heights))
    with pytest.raises(ValueError,match='topology'):
        compact_surface(None,patches,CompactSurfaceOptions(budgeted=False),base_z=0.)
    result,audit=regularize_roof_junctions(patches,bases={0:0.},detail=.5)
    assert len(audit[0]['changes'])==1
    change=audit[0]['changes'][0]
    assert change['maximum_displacement_m']<=.25+2e-7
    assert any(i['vertical_sheets']==4 for i in change['height_intervals'])
    assert shapely.union_all([p.polygon for p in patches]).equals(shapely.union_all([p.polygon for p in result]))
    assert {p.source_key for p in result}=={p.source_key for p in patches}
    for p in result:
        original=next(q for q in patches if q.source_key==p.source_key)
        np.testing.assert_array_equal(p.coefficients,original.coefficients)
        owned=shapely.union_all([q.polygon for q in result if q.source_key==p.source_key])
        assert owned.boundary.hausdorff_distance(original.polygon.boundary)<=.25+1e-8
    _,measured=compact_surface(None,result,CompactSurfaceOptions(budgeted=False),base_z=0.)
    assert measured['topology']['valid']
    repeated,repeated_audit=regularize_roof_junctions(patches,bases={0:0.},detail=.5)
    assert audit==repeated_audit
    assert all(p.polygon.equals_exact(q.polygon,0.) for p,q in zip(result,repeated))


@pytest.mark.parametrize('same_source',[False,True])
def test_source_area_audit_excludes_terrain_and_intrasource_transfers(same_source):
    patches=quadrants([8,4,8,None])
    if same_source:
        patches[1]=replace(patches[1],source_key=patches[0].source_key)
    patches=node_interfaces(patches)
    result,audit=regularize_roof_junctions(patches,bases={0:0.},detail=.5)
    change=audit[0]['changes'][0]
    total_added=total_removed=0.
    for source,measured in change['source_changes'].items():
        before=shapely.union_all([p.polygon for p in patches if p.source_key==source])
        after=shapely.union_all([p.polygon for p in result if p.source_key==source])
        assert measured['added_m2']==pytest.approx(after.difference(before).area,abs=1e-12)
        assert measured['removed_m2']==pytest.approx(before.difference(after).area,abs=1e-12)
        total_added+=measured['added_m2'];total_removed+=measured['removed_m2']
    assert total_added==pytest.approx(total_removed,abs=1e-12)
    assert total_added==pytest.approx(.03125 if same_source else .0625,abs=1e-8)


def test_affine_profile_crossings_are_noded_before_sector_classification():
    # Equal central heights with changing relative order along shared rays are
    # harmless; the exact shared graph includes those crossings before analysis.
    patches=node_interfaces(quadrants([8,8,8,8],slopes=[(.1,.2),(-.1,.2),(.2,-.1),(.1,-.2)]))
    result,audit=regularize_roof_junctions(patches,bases={0:0.},detail=.5)
    assert not audit[0]['changes']
    assert all(s['manifold'] for s in audit[0]['checked_junctions'])
    _,measured=compact_surface(None,result,CompactSurfaceOptions(budgeted=False),base_z=0.)
    assert measured['topology']['valid']


def test_profile_order_can_reverse_inside_a_shared_interface():
    patches=[RoofPatch(box(-10,-10,0,10),np.zeros(2),np.array([0.,.1,8.]),'a',0,0),
             RoofPatch(box(0,-10,10,10),np.zeros(2),np.array([0.,-.1,8.]),'b',1,0)]
    noding={}
    patches=node_interfaces(patches,audit=noding)
    assert noding['height_crossings']==1
    result,audit=regularize_roof_junctions(patches,bases={0:0.},detail=.5)
    assert not audit[0]['changes']
    _,measured=compact_surface(None,result,CompactSurfaceOptions(budgeted=False),base_z=0.)
    assert measured['topology']['valid']


@pytest.mark.parametrize('method',['catalogue','ransac'])
def test_orebro_admitted_junctions_need_no_notches(method):
    import json
    from pathlib import Path
    evidence=json.loads((Path(__file__).parents[1]/'data/lod2/orebro_source_junctions.json').read_text())
    checked=[]
    for block in evidence['methods'][method]:
        patches=[RoofPatch(shapely.from_wkt(p['wkt']),np.asarray(p['origin']),
            np.asarray(p['coefficients']),p['source_key'],p['patch_id'],p['block_id']) for p in block['patches']]
        offset=np.floor(np.min([p.polygon.bounds[:2] for p in patches],axis=0))
        patches=node_interfaces([replace(p,polygon=affinity.translate(p.polygon,*-offset),
            origin=p.origin-offset) for p in patches])
        result,audit=regularize_roof_junctions(patches,bases={block['block']:0.},detail=.5)
        stars=audit[block['block']]['checked_junctions']
        assert not audit[block['block']]['changes']
        checked.extend(stars)
        assert all(p.polygon.equals_exact(q.polygon,0.) for p,q in zip(patches,result))
        _,measured=compact_surface(None,result,CompactSurfaceOptions(),base_z=0.)
        assert measured['topology']['valid']
        assert measured['building_triangles']<=1.5*measured['T0']
    assert len(checked)==5
    assert all(s['manifold'] for s in checked)
    assert sum(len(s['observations']) for s in evidence['source_evidence'])==31984


def test_junction_construction_preserves_the_fixed_occupied_outline():
    import json
    from pathlib import Path
    evidence=json.loads((Path(__file__).parents[1]/'data/lod2/admitted_junction_outline.json').read_text())
    patches=[RoofPatch(shapely.from_wkt(p['wkt']),np.asarray(p['origin']),
        np.asarray(p['coefficients']),p['source_key'],p['patch_id'],p['block_id']) for p in evidence['patches']]
    before=shapely.union_all([p.polygon for p in patches])
    result,audit=regularize_roof_junctions(patches,bases={33:evidence['base']},detail=.5)
    assert len(audit[33]['changes'])==1
    after=shapely.union_all([p.polygon for p in result])
    assert after.symmetric_difference(before).area<1e-10
    assert after.boundary.hausdorff_distance(before.boundary)<1e-8
    mesh,measured=compact_surface(None,result,CompactSurfaceOptions(),base_z=evidence['base'])
    assert measured['topology']['valid']
    assert len(mesh.faces)>0
    assert measured['building_triangles']<=1.5*measured['T0']
    assert measured['quality_after']['roof']['q_p01']>.09
    assert measured['quality_after']['roof']['below_q_0.02']==0
    change=audit[33]['changes'][0]
    assert change['minimum_new_corner_edge_m']>.05
    assert change['minimum_ray_corner_separation_degrees']>10.


def test_rectangle_is_square_for_orthogonal_rays_and_avoids_diagonal_rays():
    half,audit=junction_half_widths([0.,np.pi/2,np.pi,3*np.pi/2],.25,1.)
    np.testing.assert_allclose(half,np.full(2,.25/np.sqrt(2)),atol=1e-15)
    diagonal=[np.pi/4+k*np.pi/2 for k in range(4)]
    half,audit=junction_half_widths(diagonal,.25,1.)
    assert np.linalg.norm(half)==pytest.approx(.25)
    assert max(half)/min(half)<=np.sqrt(3)+1e-14
    assert audit['minimum_ray_corner_separation_degrees']==pytest.approx(15.)
    # Unresolvable perturbations must not switch the deterministic tie choice.
    shifted,other=junction_half_widths([a+1e-10 for a in diagonal],.25,1.)
    np.testing.assert_array_equal(half,shifted)


@pytest.mark.parametrize('scale',[.001,1.,1000.])
@pytest.mark.parametrize('heights',[[8,8,8,8],[8,4,8,4]])
def test_sector_policy_scales_with_the_declared_physical_detail(scale,heights):
    patches=[replace(p,polygon=affinity.scale(p.polygon,xfact=scale,yfact=scale,origin=(0,0)))
             for p in quadrants(heights)]
    result,audit=regularize_roof_junctions(node_interfaces(patches),bases={0:0.},detail=.5*scale)
    assert bool(audit[0]['changes'])==(heights[0]!=heights[1])
    assert all(j['maximum_displacement_m']<=.25*scale+2e-7 for j in audit[0]['changes'])
    _,measured=compact_surface(None,result,CompactSurfaceOptions(budgeted=False),base_z=0.)
    assert measured['topology']['valid']
