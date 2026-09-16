"""Live-allocation accounting for one or more independent canonical model loads.

Use separate fresh processes with/without --trace. Tracemalloc tracks Python and
NumPy allocations, not all native allocations, and increases RSS. It is unsuitable
for timing comparisons. --copies 3 is a synthetic collection of independent copies,
not a larger real geographic dataset. Inventory runs after the measured snapshots.
"""

import argparse
from collections import Counter
from dataclasses import fields, is_dataclass
import gc
import hashlib
import json
from pathlib import Path
import platform
import resource
import sys
import time
import tracemalloc


def inventory(models):
    """Count native state without materializing instance __dict__ mappings."""
    import numpy as np

    seen, arrays, strings = set(), {}, {}
    counts = Counter()
    stack = list(models)
    while stack:
        value = stack.pop()
        if id(value) in seen:
            continue
        seen.add(id(value))
        if isinstance(value, np.ndarray):
            arrays[id(value)] = value.nbytes
            counts['arrays'] += 1
            counts['empty_arrays'] += value.size == 0
        elif is_dataclass(value) and not isinstance(value, type):
            counts[type(value).__name__] += 1
            stack.extend(getattr(value, f.name) for f in fields(value))
        elif isinstance(value, dict):
            counts['empty_dicts' if not value else 'dicts'] += 1
            stack.extend(value.keys())
            stack.extend(value.values())
        elif isinstance(value, (list, tuple)):
            counts['empty_sequences' if not value else 'sequences'] += 1
            stack.extend(value)
        elif isinstance(value, str):
            strings[id(value)] = value
    return {'counts': dict(counts), 'array_bytes': sum(arrays.values()),
            'string_objects': len(strings), 'distinct_strings': len(set(strings.values())),
            'string_bytes': sum(sys.getsizeof(s) for s in strings.values()),
            'distinct_string_bytes': sum(sys.getsizeof(s) for s in set(strings.values()))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('artifact', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--copies', type=int, choices=(1, 3), default=1)
    parser.add_argument('--trace', action='store_true')
    args = parser.parse_args()
    import psutil  # Existing development environment; not a production dependency.
    from dtcc_core import io

    with args.artifact.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    process = psutil.Process()

    def snapshot():
        gc.collect()
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != 'darwin':
            peak *= 1024
        result = {'rss_MiB': process.memory_info().rss / 2**20, 'peak_MiB': peak / 2**20}
        if args.trace:
            current, peak = tracemalloc.get_traced_memory()
            result.update(traced_live_MiB=current / 2**20, traced_peak_MiB=peak / 2**20)
        return result

    if args.trace:
        tracemalloc.start()
    baseline = snapshot()
    start = time.perf_counter()
    models = [io.load_model(args.artifact) for _ in range(args.copies)]
    seconds = time.perf_counter() - start
    loaded = snapshot()
    allocations = []
    if args.trace:
        for stat in tracemalloc.take_snapshot().statistics('lineno')[:15]:
            allocations.append({'location': str(stat.traceback), 'bytes': stat.size, 'count': stat.count})
    # The memory/timing measurements above exclude this diagnostic graph walk.
    storage = inventory(models) if args.trace else None
    del models
    released = snapshot()
    report = {'python': platform.python_version(), 'platform': platform.platform(),
              'artifact_sha256': digest, 'copies': args.copies, 'instrumented': args.trace,
              'load_seconds': seconds, 'baseline': baseline, 'loaded': loaded, 'released': released,
              'inventory': storage, 'live_allocation_sites': allocations, 'limits': __doc__}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('copies', 'instrumented', 'load_seconds', 'loaded', 'released')}, indent=2))


if __name__ == '__main__':
    main()
