"""Fixed refinement budgets, conforming transactions and explicit failures."""
from copy import deepcopy
import importlib

import numpy as np
import pytest

from .test_compact_mesher_interface import fixture
from .test_compact_integration import two_groups, sample
from dtcc_core import builder
from dtcc_core.builder import CompactSurfaceOptions

ci = importlib.import_module('dtcc_core.builder.meshing.compact_improvement')


def policy(source=1, group=1, **kwargs):
    return CompactSurfaceOptions(max_face_vertices=0, max_source_vertices=source,
                                 max_group_vertices=group, **kwargs)


def run(options=None, ceiling=None):
    registry, caches, faces = fixture()
    registry, audit = ci.improve_shared(registry, caches, faces, options or policy(),
                                        max_building_triangles=ceiling)
    return registry, caches, faces, audit


def test_triangle_budget_rejects_complete_transaction_atomically():
    registry, caches, faces = fixture()
    initial = deepcopy((registry, caches))
    count = sum(len(c['triangles']) for c in caches.values())
    registry, audit = ci.improve_shared(registry, caches, faces, policy(4, 8),
                                        max_building_triangles=count+1)
    assert audit['stopping']['primary_reason'] == 'building_triangle_ceiling'
    assert not audit['optional_vertices']
    assert not audit['stopping']['evaluated_finite_set']
    np.testing.assert_array_equal(registry.vertices, initial[0].vertices)
    for face, cache in caches.items():
        np.testing.assert_array_equal(cache['triangles'], initial[1][face]['triangles'])


@pytest.mark.parametrize('source,group,binding', [(1,8,'source_point_allowance'),
                                               (4,1,'group_point_allowance')])
def test_point_allowances_are_fixed(source, group, binding):
    *_, audit = run(policy(source, group))
    assert binding in audit['stopping']['binding_limits']
    assert len(audit['optional_vertices']) <= group
    assert all(count <= source for count in audit['source_optional_vertices'].values())
    assert not audit['stopping']['evaluated_finite_set']


def test_each_accepted_shared_point_has_conformity_and_measured_gain():
    registry, caches, faces, audit = run(policy(4,8))
    assert audit['optional_vertices']
    for point in audit['optional_vertices']:
        edge = registry.edges[point['edge_id']]
        assert point['affected_faces'] == sorted(edge.incidents)
        assert all(point['vertex'] in caches[i]['ids'] for i in edge.incidents)
        assert point['before']['penalty']-point['after']['penalty'] > .015*point['triangle_cost']
        assert point['after']['minimum_q'] >= point['before']['minimum_q']-1e-10
        for angle in (5,10):
            assert point['after'][f'poor_{angle}_count'] <= point['before'][f'poor_{angle}_count']


@pytest.mark.parametrize('budget,reason', [('MAX_NATIVE_REQUESTS','native_request_limit'),
                                         ('MAX_JOINT_TRIALS','joint_trial_limit')])
def test_work_budget_cannot_reset_or_commit_unassessed_points(monkeypatch, budget, reason):
    monkeypatch.setattr(ci, budget, 0)
    *_, audit = run(policy(4,8))
    assert audit['stopping']['primary_reason'] == reason
    assert not audit['optional_vertices']
    assert audit['native_requests'] == audit['joint_trials'] == 0


def test_controlled_clock_exhaustion(monkeypatch):
    now = [0.]
    original = ci.improve
    def elapsed(*args, **kwargs):
        result = original(*args, **kwargs)
        if kwargs.get('max_added_points') == 1:
            now[0] = ci.MAX_GROUP_SECONDS+1
        return result
    monkeypatch.setattr(ci.time, 'monotonic', lambda:now[0])
    monkeypatch.setattr(ci, 'improve', elapsed)
    *_, audit = run(policy(4,8))
    assert audit['stopping']['primary_reason'] == 'runtime_limit'
    assert not audit['stopping']['evaluated_finite_set']


def test_native_failure_is_an_error_and_does_not_commit_geometry(monkeypatch):
    original = ci.improve
    def failed(*args, **kwargs):
        if kwargs.get('max_added_points') == 1:
            raise RuntimeError('controlled native failure')
        return original(*args, **kwargs)
    registry, caches, faces = fixture()
    vertices = np.asarray(registry.vertices).copy()
    monkeypatch.setattr(ci, 'improve', failed)
    with pytest.raises(RuntimeError, match='controlled native failure'):
        ci.improve_shared(registry, caches, faces, policy(4,8))
    np.testing.assert_array_equal(registry.vertices, vertices)


def test_detailed_audit_does_not_change_decisions():
    states = [run(policy(4,8,detailed_audit=detailed)) for detailed in (False,True)]
    np.testing.assert_array_equal(states[0][0].vertices, states[1][0].vertices)
    for face in states[0][1]:
        np.testing.assert_array_equal(states[0][1][face]['triangles'], states[1][1][face]['triangles'])
    assert states[0][3]['optional_vertices'] == states[1][3]['optional_vertices']
    assert not states[0][3]['native_outcomes'] and states[1][3]['native_outcomes']


@pytest.mark.parametrize('cap', [0,1,2])
@pytest.mark.parametrize('shared', [False,True])
def test_face_allowance_reports_unknown_excluded_gains(cap, shared):
    registry,caches,faces = fixture()
    _,audit = ci.improve_shared(registry,caches,faces,CompactSurfaceOptions(
        shared_edge_points=shared,max_face_vertices=cap,max_source_vertices=8,max_group_vertices=32))
    stopping = audit['stopping']
    assert 'face_point_allowance' in stopping['binding_limits']
    assert stopping['category']=='budget_reached'
    assert not stopping['evaluated_finite_set']
    excluded = stopping['excluded_modes']
    assert any(m['face_id']==0 and m['kind']=='interior' and
               'face_point_allowance' in m['binding_limits'] for m in excluded)
    assert all(m['additional_gain']=='unknown' for m in excluded)
    if not shared:
        assert caches[0]['optional']==cap
        cache=caches[0]
        trial=ci.improve(cache['points'],cache['triangles'],cache['segments'],
                         max_added_points=1,splittable_segments=[],allow_interior=True)
        assert trial.inserted_point is not None
        assert ci.joint_acceptance([ci.quality(cache['points'],cache['triangles'],cache['segments'])],
            [ci.quality(trial.mesh.points,trial.mesh.triangles,cache['segments'])],
            len(trial.mesh.triangles)-len(cache['triangles']))[0]
        # The other face has genuinely exhausted its allowed native proposals.
        if cap:
            assert stopping['mode_counts']['evaluated']>0
            assert stopping['mixed_mode_restrictions']


def test_gain_exhaustion_is_evaluated_without_saturated_allowances():
    registry,caches,faces=fixture()
    _,audit=ci.improve_shared(registry,caches,faces,CompactSurfaceOptions(
        shared_edge_points=False),face_ids={1})
    assert audit['stopping']['primary_reason']=='marginal_gain_exhausted'
    assert audit['stopping']['evaluated_finite_set']
    assert not audit['stopping']['excluded_modes']
    assert not audit['stopping']['unevaluated_modes']


@pytest.mark.parametrize('source,group,binding', [(0,8,'source_point_allowance'),
                                               (8,0,'group_point_allowance')])
def test_zero_allowances_are_budget_restrictions(source,group,binding):
    *_,audit=run(policy(source,group))
    assert binding in audit['stopping']['binding_limits']
    assert not audit['stopping']['evaluated_finite_set']
    assert not audit['optional_vertices']


@pytest.mark.parametrize('group_count', [1,2])
@pytest.mark.parametrize('enabled', [False,True])
def test_shell_validation_count_is_independent_of_refinement_groups(monkeypatch, group_count, enabled):
    from shapely.geometry import box
    from .test_roof_catalogue import assemble
    from dtcc_core.builder.roofs import RoofPatch
    from dtcc_core.builder.meshing.compact_surface import compact_surface
    validation = importlib.import_module('dtcc_core.builder.meshing.roof_validation')
    patches = [RoofPatch(box(30*i-18,-5,30*i-3,5),np.zeros(2),np.array([0.,0.,8.]),
                         str(i),block_id=i) for i in range(group_count)]
    reference, patches = assemble(patches, slope=.4)
    original = validation.validate_shell
    calls = []
    def counted(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(validation, 'validate_shell', counted)
    _, audit = compact_surface(reference, patches, CompactSurfaceOptions(budgeted=enabled))
    assert len(calls) == 2
    assert len(audit['optimization']['stages']) == group_count
    assert audit['topology']['valid']
    if not enabled:
        assert all(s['native_requests']==s['joint_trials']==0 for s in audit['optimization']['stages'])


def test_final_topology_failure_never_returns_a_successful_stop(monkeypatch):
    validation = importlib.import_module('dtcc_core.builder.meshing.roof_validation')
    original = validation.validate_shell
    count = []
    def reject_final(*args, **kwargs):
        report = original(*args, **kwargs)
        count.append(True)
        return dict(report, valid=False) if len(count)==2 else report
    monkeypatch.setattr(validation, 'validate_shell', reject_final)
    audit = {}
    with pytest.raises(RuntimeError, match='topology'):
        builder.build_city_surface_mesh(sample(), lod='lod2', report_mesh_quality=False, stage_audit=audit)
    assert not audit['roof_catalogue']['block_fallbacks']
    assert audit['roof_catalogue']['sources'][0]['outcome'] == 'accepted_fit'
    assert len(audit['roof_catalogue']['assemblies']) == 1
