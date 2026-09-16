"""Compare default and bypassed public canonical I/O in fresh sequential processes.

Input is the unchanged local representations-example artifact directory. Each
worker prepares v6 data with evaluation bypassed, then measures one cold operation
and three warm repetitions. Structural checks always run. Preparation, GC and
memory snapshots are outside timings; Python/DTCC imports precede the cold sample.
The first enabled sample includes lazy LinkML import/compilation. No network runs.
"""

import argparse
import gc
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

OPERATIONS = ('file_save', 'file_load', 'package_save', 'package_load')


def worker(args):
    import psutil
    from dtcc_core import io
    from dtcc_core.datasets import load_model_package
    from dtcc_core.model import exchange

    city = load_model_package(args.source / 'package.dtccpkg', validate_schema=False)
    file = args.output / 'model.dtcc'
    package = args.output / 'model.dtccpkg'
    if args.operation == 'file_load':
        io.save_model(city, file, validate_schema=False)
    elif args.operation == 'package_load':
        city.export(package, canonical=True, validate_schema=False)
    enabled = args.mode == 'default'
    if args.operation == 'file_save':
        operation = lambda: io.save_model(city, file, validate_schema=enabled)
    elif args.operation == 'file_load':
        del city
        operation = lambda: io.load_model(file, validate_schema=enabled)
    elif args.operation == 'package_save':
        operation = lambda: city.export(package, canonical=True, validate_schema=enabled)
    else:
        del city
        operation = lambda: load_model_package(package, validate_schema=enabled)
    process = psutil.Process()
    samples = []
    gc.collect()
    baseline = process.memory_info().rss / 2**20
    for _ in range(4):
        gc.collect()
        start = time.perf_counter()
        result = operation()
        samples.append(time.perf_counter() - start)
        del result
    gc.collect()
    report = {'cold_seconds': samples[0], 'warm_seconds': samples[1:],
              'prepared_rss_MiB': baseline, 'finished_rss_MiB': process.memory_info().rss / 2**20,
              'wire_version': exchange.VERSION}
    (args.output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--operation', choices=OPERATIONS, help=argparse.SUPPRESS)
    parser.add_argument('--mode', choices=('default', 'bypass'), help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.operation:
        worker(args)
        return
    import hashlib
    with (args.source / 'model.dtcc').open('rb') as stream:
        source_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    with (args.source / 'package.dtccpkg').open('rb') as stream:
        package_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    import zipfile
    with zipfile.ZipFile(args.source / 'package.dtccpkg') as archive:
        if hashlib.sha256(archive.read('artifacts/model.dtcc')).hexdigest() != source_hash:
            parser.error('Source package must contain the same current-format model')
    report = {'python': platform.python_version(), 'platform': platform.platform(),
              'model_sha256': source_hash, 'package_sha256': package_hash,
              'method': __doc__, 'operations': {}}
    for name in OPERATIONS:
        values = {}
        for mode in ('bypass', 'default'):
            samples = []
            for repeat in range(3):
                output = args.output / f'{name}-{mode}-{repeat}'
                output.mkdir(exist_ok=True)
                command = [sys.executable, str(Path(__file__).resolve()), str(args.source), str(output),
                           '--operation', name, '--mode', mode]
                with (output / 'operation.log').open('w') as log:
                    subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
                samples.append(json.loads((output / 'result.json').read_text()))
            values[mode] = {'cold_median': statistics.median(s['cold_seconds'] for s in samples),
                            'warm_median': statistics.median(t for s in samples for t in s['warm_seconds']),
                            'samples': samples}
        values['extra_seconds'] = values['default']['warm_median'] - values['bypass']['warm_median']
        values['extra_percent'] = 100 * values['extra_seconds'] / values['bypass']['warm_median']
        report['operations'][name] = values
        print(name, {k: values[k] for k in ('extra_seconds', 'extra_percent')}, flush=True)
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
