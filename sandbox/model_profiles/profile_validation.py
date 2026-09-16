"""Measure explicit profile validation in one fresh process on a local native model.

Run several sequential fresh processes. Timings exclude file hashing and GC/RSS
snapshots. The public validation timings INCLUDE native structural admission.
--breakdown is a separate diagnostic run, not a shortcut in the public API.
"""

import argparse
import gc
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import resource
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', type=Path)
    parser.add_argument('schema', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--breakdown', action='store_true')
    args = parser.parse_args()
    import psutil
    from dtcc_core import io
    from dtcc_core.model import exchange
    from dtcc_core.model.profiles import SemanticProfile, _project

    process = psutil.Process()

    def snapshot():
        gc.collect()
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return {'rss_MiB': process.memory_info().rss / 2**20,
                'peak_MiB': peak / 2**20 if sys.platform == 'darwin' else peak / 1024}

    def measure(operation):
        start = time.perf_counter()
        value = operation()
        return value, time.perf_counter() - start

    report = {'python': platform.python_version(), 'platform': platform.platform(),
              'versions': {name: importlib.metadata.version(name)
                           for name in ('linkml', 'linkml-runtime', 'jsonschema', 'numpy', 'protobuf', 'pydantic')},
              'model_bytes': args.model.stat().st_size,
              'schema_sha256': hashlib.sha256(args.schema.read_bytes()).hexdigest()}
    with args.model.open('rb') as stream:
        report['model_sha256'] = hashlib.file_digest(stream, 'sha256').hexdigest()
    report['baseline'] = snapshot()
    # Keep this historical explicit-profile benchmark separate from the now
    # default standard-schema evaluation in canonical I/O.
    model, report['load_seconds'] = measure(lambda: io.load_model(args.model, validate_schema=False))
    report['loaded'] = snapshot()
    profile, report['setup_seconds'] = measure(lambda: SemanticProfile(args.schema))
    report['profile_ready'] = snapshot()
    result, report['first_seconds'] = measure(lambda: profile.validate(model))
    report['first_valid'] = result.valid
    report['first_issue_count'] = len(result.issues)
    report['after_first'] = snapshot()
    report['repeated_seconds'] = []
    for _ in range(3):
        result, elapsed = measure(lambda: profile.validate(model))
        assert result.valid == report['first_valid'] and len(result.issues) == report['first_issue_count']
        report['repeated_seconds'].append(elapsed)
    report['after_repeated'] = snapshot()
    if args.breakdown:
        _, native = measure(lambda: exchange.validate(model))
        (records, locations), projection = measure(lambda: _project(model, profile._backend.references))
        values, flatten = measure(lambda: [{'id': n['id'], **n['attributes'], **n['relations']} for n in records])
        _, semantic = measure(lambda: profile._backend.validate(records, values))
        report['breakdown'] = {'native_seconds': native, 'projection_seconds': projection,
                               'flatten_seconds': flatten, 'semantic_seconds': semantic,
                               'records': len(records), 'objects': sum('regions[' not in b for b, _ in locations.values())}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
