"""Independent historical baselines must not be cleared by minimal-mesh gains."""
import numpy as np
import pytest

from benchmarks.lod2_comparison import (
    measured_quality,orebro_acceptance,quality_flags,source_change,with_fractions)
from dtcc_core.model import Mesh


def test_minimal_improvement_does_not_clear_orebro_historical_regression():
    minimal=dict(triangles=1400,minimum_q=.004,q_p01=.019,q_p05=.09,
                 **{'below_q_0.02':16,'below_q_0.1':100,'below_q_0.25':200})
    refined=dict(minimal,minimum_q=.005,q_p01=.022,q_p05=.12,
                 **{'below_q_0.02':14})
    original=dict(minimal,minimum_q=.018,q_p01=.0928,q_p05=.16,
                  **{'below_q_0.02':1})
    assert not quality_flags(minimal,refined)
    assert quality_flags(original,refined)
    targets=orebro_acceptance(refined,'catalogue')
    assert not targets['p01_met'] and not targets['severe_count_met']
    assert not orebro_acceptance(dict(refined,q_p01=.0813,**{'below_q_0.02':0}),'catalogue')['p01_met']
    assert orebro_acceptance(dict(refined,q_p01=.0646,**{'below_q_0.02':3}),'ransac')['p01_met']


def test_quality_measurements_include_absolute_counts_fractions_and_wall_geometry():
    mesh=Mesh(vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,0.,.001],[0.,1.,0.]]),
              faces=np.array([[0,1,2],[0,1,3]]),markers=np.array([0,1]))
    metrics=measured_quality(mesh,np.array(['wall','roof']))
    assert metrics['wall']['below_q_0.02']==1
    assert metrics['wall']['fraction_below_q_0.02']==1.
    assert metrics['wall']['area_below_5_degrees_m2']==pytest.approx(.0005)
    assert metrics['roof']['fraction_below_q_0.02']==0.
    assert metrics['terrain']['minimum_q'] is None
    assert metrics['terrain']['fraction_below_q_0.02']==0.
    assert with_fractions(dict(triangles=2,below_q_0_02=1))['fraction_below_q_0_02']==.5


def test_source_attribution_separates_model_and_support_changes():
    mesh=Mesh(vertices=np.array([[0.,0.,8.],[1.,0.,8.],[0.,1.,8.]]),
              faces=np.array([[0,1,2]]),markers=np.array([1]))
    catalogue=dict(patch_markers=[dict(source_key='a',coefficients=[0.,0.,8.],origin=[0.,0.],fit_retained=True)])
    state=(mesh,np.array(['roof']),np.array(['a']),catalogue)
    assert source_change(state,state,'a')['classification']=='unchanged_roof_geometry'
    shifted=Mesh(vertices=mesh.vertices+np.array([.1,0.,0.]),faces=mesh.faces,markers=mesh.markers)
    changed=(shifted,*state[1:])
    assert source_change(state,changed,'a')['classification']=='source_support_changed'
    different=dict(patch_markers=[dict(catalogue['patch_markers'][0],fit_retained=False)])
    assert source_change(state,(*state[:3],different),'a')['classification']=='roof_model_or_retention_changed'
