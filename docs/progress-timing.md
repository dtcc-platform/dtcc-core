# Timing progress phases

`dtcc_core.common.progress.phase_timing(factory)` lets a caller measure existing
`ProgressTracker.phase()` scopes without changing dataset signatures. Supply a
factory that accepts a phase name and returns a context manager:

```python
from contextlib import contextmanager
from time import perf_counter
from dtcc_core import datasets
from dtcc_core.common.progress import phase_timing

elapsed = {}

@contextmanager
def measure(name):
    started = perf_counter()
    try:
        yield
    finally:
        elapsed[name] = elapsed.get(name, 0.0) + perf_counter() - started

with phase_timing(measure):
    mesh = datasets.city_volume_mesh(bounds=bounds)
```

The observer receives nested phases even when progress output is silent, and
receives exceptions from the measured operation. Registrations are context
local; nested registrations restore their parent on exit. `phase_timing(None)`
temporarily suspends observation. Worker threads and processes require their
own registration. The observer must preserve the operation's behavior.

Durations measured around a parent phase include all its children. A caller
that needs an additive breakdown must track nesting and subtract child time.
This hook does not change progress weights or add synchronization. Core's C++
`Timer` registry measures native operations separately.
