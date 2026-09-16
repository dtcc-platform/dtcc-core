"""Fresh-process measurements for the unchanged 3DBAG public workflow.

Uses psutil from the development environment for resident memory snapshots.
Timings exclude process startup/preparation, include the public operation, and
retain its result through the final memory snapshot. OS file caches are not reset.
Peak RSS includes imports and preparation; it is not incremental allocation.
"""

import argparse
import cProfile
import gc
import hashlib
import json
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

STAGES = ('import', 'encode', 'decode', 'package_write', 'package_read', 'cityjson_write')
SOURCE_SHA256 = '2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba'


def worker(args):
    import numpy as np
    import psutil
    from dtcc_core import io
    from dtcc_core.datasets import load_model_package
    from dtcc_core.model import exchange

    process = psutil.Process()

    def memory():
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != 'darwin':
            peak *= 1024
        return {'rss_MiB': process.memory_info().rss / 2**20, 'peak_MiB': peak / 2**20}

    baseline = memory()
    if args.stage in ('encode', 'package_write', 'cityjson_write'):
        city = io.load_city(args.source, strict=True, extent_policy='recompute')
    elif args.stage == 'decode':
        data = (args.artifacts / 'model.dtcc').read_bytes()

    if args.stage == 'import':
        operation = lambda: io.load_city(args.source, strict=True, extent_policy='recompute')
    elif args.stage == 'encode':
        operation = lambda: exchange.dumps(city)
    elif args.stage == 'decode':
        operation = lambda: exchange.loads(data)
    elif args.stage == 'package_write':
        operation = lambda: city.export(args.output / 'measured.dtccpkg', canonical=True)
    elif args.stage == 'package_read':
        operation = lambda: load_model_package(args.artifacts / 'package.dtccpkg')
    else:
        operation = lambda: city.save(args.output / 'measured.city.json', strict=True)

    gc.collect()
    prepared = memory()
    profiler = cProfile.Profile() if args.profile else None
    start = time.perf_counter()
    result = profiler.runcall(operation) if profiler else operation()
    elapsed = time.perf_counter() - start
    finished = memory()
    report = {'seconds': elapsed, 'imports': baseline, 'prepared': prepared, 'finished': finished}
    if args.stage == 'encode':
        report['payload_bytes'] = len(result)
    if args.stage == 'import':
        # Numerical buffers, not total Python/native heap size. Count shared arrays once.
        arrays, objects, geometries = {}, 0, 0
        stack = [result]
        while stack:
            value = stack.pop()
            for name in ('vertices', 'normal', 'indices'):
                array = getattr(value, name, None)
                if isinstance(array, np.ndarray):
                    arrays[id(array)] = array.nbytes
            if hasattr(value, 'transform'):
                arrays[id(value.transform.affine)] = value.transform.affine.nbytes
            if hasattr(value, 'geometry'):
                objects += 1
                stack.extend(r.geometry for r in value.geometry.values())
                stack.extend(c for group in value.children.values() for c in group)
            else:
                geometries += hasattr(value, 'transform')
                stack.extend(getattr(value, 'surfaces', []))
                stack.extend(getattr(value, 'regions', []))
                for array in [*getattr(value, 'holes', []), *getattr(value, 'shells', [])]:
                    arrays[id(array)] = array.nbytes
        report['storage'] = {'objects': objects, 'geometries': geometries,
                             'arrays': len(arrays), 'array_bytes': sum(arrays.values())}
    if profiler:
        profiler.dump_stats(str(args.output / f'{args.stage}.prof'))
    (args.output / f'{args.stage}-{args.run}.json').write_text(json.dumps(report, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('artifacts', type=Path, help='Existing representations_example.py output directory')
    parser.add_argument('output', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--stage', choices=STAGES)
    parser.add_argument('--profile', action='store_true', help='Separate instrumented run; do not compare its timings')
    parser.add_argument('--run', type=int, default=0, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage:
        worker(args)
        return
    with args.source.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    if digest != SOURCE_SHA256:
        parser.error('Expected the unchanged recorded 3DBAG source')
    report = {'python': platform.python_version(), 'platform': platform.platform(),
              'source_sha256': digest, 'profiled': args.profile, 'repeats': args.repeats,
              'limits': __doc__, 'stages': {}}
    for stage in STAGES:
        samples = []
        for run in range(args.repeats):
            command = [sys.executable, str(Path(__file__).resolve()), str(args.source),
                       str(args.artifacts), str(args.output), '--stage', stage, '--run', str(run)]
            if args.profile:
                command.append('--profile')
            with (args.output / f'{stage}-{run}.log').open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            samples.append(json.loads((args.output / f'{stage}-{run}.json').read_text()))
        summary = {'seconds_median': statistics.median(s['seconds'] for s in samples),
                   'peak_MiB_median': statistics.median(s['finished']['peak_MiB'] for s in samples),
                   'samples': samples}
        report['stages'][stage] = summary
        print(stage, {k: v for k, v in summary.items() if k != 'samples'}, flush=True)
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
