"""Freeze the exact shared pilot state before testing additional edge points."""
from copy import deepcopy
import importlib
import pickle
from unittest.mock import patch

from compact_mesher_evaluation import evaluate, WORK

OUT = WORK / 'compact-edge-extension'


def capture(case):
    cs = importlib.import_module('dtcc_core.builder.meshing.compact_surface')
    ci = importlib.import_module('dtcc_core.builder.meshing.compact_improvement')
    state = {}
    surface, shared = cs.compact_surface, ci.improve_shared

    def geometry(reference, patches, options):
        state.update(reference=deepcopy(reference), patches=deepcopy(patches))
        return surface(reference, patches, options)

    def improved(registry, caches, faces, options):
        registry, audit = shared(registry, caches, faces, options)
        state.update(deepcopy(dict(registry=registry, caches=caches, faces=faces, audit=audit)))
        return registry, audit

    with patch.object(cs, 'compact_surface', geometry), patch.object(ci, 'improve_shared', improved):
        evaluate(case, 'dtcc_mesher', OUT / 'checkpoint-verification', 'shared')
    from dtcc_core import io
    import numpy as np
    name = f'case-{case:02d}-shared-dtcc_mesher.dtcc'
    old = io.load_mesh(WORK / 'compact-mesher-interface' / name)
    new = io.load_mesh(OUT / 'checkpoint-verification' / name)
    assert np.array_equal(old.vertices, new.vertices)
    assert np.array_equal(old.faces, new.faces)
    assert np.array_equal(old.markers, new.markers)
    with (OUT / f'case-{case:02d}-checkpoint.pkl').open('wb') as stream:
        pickle.dump(state, stream)
    print('Frozen exact shared checkpoint', case, flush=True)


if __name__ == '__main__':
    for case in (1, 11):
        capture(case)
