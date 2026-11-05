"""
Keep `import dtcc_core` cheap.

Importing dtcc_core once cost ~0.29 s, of which scipy alone was 37% -- pulled
in transitively by a handful of module-scope imports in builder/ that only a
few functions actually needed. Moving those imports into the functions that use
them cut the import roughly in half (issue #87).

The dependency check below is the load-bearing one: it is deterministic and
does not depend on how fast the machine is. The wall-clock check is a coarse
backstop with a deliberately generous threshold, so it catches a tenfold
regression without failing on a loaded CI runner.
"""

import subprocess
import sys
import time

import pytest

# Third-party packages that must NOT be imported by a bare `import dtcc_core`.
# Each is expensive and only needed by specific operations, so it belongs in
# the functions that use it. Keep imports of these inside function bodies.
FORBIDDEN_EAGER_IMPORTS = [
    "fiona",
    "geopandas",
    "laspy",
    "rasterio",
    "rasterstats",
    "scipy",
    "skimage",
]

# Generous ceiling: roughly 10x the ~0.13 s measured after issue #87, minus
# interpreter start-up. Meant to catch a structural regression, not to police
# small drifts. Reproduce the real number with:
#   python scripts/measure_import_time.py
MAX_IMPORT_SECONDS = 1.5


def _import_dtcc_core(extra: str = "") -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-c", "import dtcc_core\n" + extra],
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("package", FORBIDDEN_EAGER_IMPORTS)
def test_heavy_dependency_not_imported_eagerly(package):
    result = _import_dtcc_core(
        "import sys\n"
        f"loaded = any(m == {package!r} or m.startswith({package + '.'!r}) "
        "for m in sys.modules)\n"
        "print(loaded)\n"
    )
    assert result.returncode == 0, f"import dtcc_core failed:\n{result.stderr}"
    loaded = result.stdout.strip().splitlines()[-1] == "True"
    assert not loaded, (
        f"`import dtcc_core` pulled in {package}, which is expensive and is not "
        "needed to import the package. Move the import into the function that "
        "uses it. Find the culprit with:\n"
        "  python -X importtime -c 'import dtcc_core'\n"
        "See issue #87."
    )


def test_import_time_within_budget():
    baseline = _timed_run("pass")
    imported = _timed_run("import dtcc_core")
    attributable = imported - baseline
    assert attributable < MAX_IMPORT_SECONDS, (
        f"`import dtcc_core` took {attributable:.2f}s (budget "
        f"{MAX_IMPORT_SECONDS:.2f}s). Investigate with:\n"
        "  python scripts/measure_import_time.py\n"
        "See issue #87."
    )


def _timed_run(statement: str, runs: int = 3) -> float:
    """Return the minimum wall-clock time of running `statement`."""
    timings = []
    for _ in range(runs):
        start = time.perf_counter()
        result = subprocess.run(
            [sys.executable, "-c", statement],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        timings.append(time.perf_counter() - start)
        assert result.returncode == 0, f"`{statement}` failed"
    return min(timings)
