"""Local synthetic codec measurement; optional profile evaluation is not timed."""

import argparse
from copy import deepcopy
import gc
import json
from pathlib import Path
import platform
import statistics
import time

from dtcc_core import io
from dtcc_core.model import City, exchange


def measure(city, repeats):
    payload = exchange.dumps(city)
    exchange.loads(payload)  # Warm up both entry points.
    results = {}
    for name, operation in [('encode', lambda: exchange.dumps(city)),
                            ('decode', lambda: exchange.loads(payload))]:
        times = []
        for _ in range(repeats):
            gc.collect()
            start = time.perf_counter()
            result = operation()
            times.append(time.perf_counter() - start)
            del result
        median = statistics.median(times)
        results[name + '_seconds_median'] = median
        results[name + '_MiB_per_second'] = len(payload) / 2**20 / median
    return {'bytes': len(payload), **results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--buildings', type=int, default=500)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.buildings < 1 or args.repeats < 1:
        parser.error('buildings and repeats must be positive')
    fixture = Path(__file__).parent / 'fixtures/buildings.city.json'
    template = io.load_city(fixture, strict=True).buildings[0]
    city = City(id='benchmark')
    for index in range(args.buildings):
        building = deepcopy(template)
        building.id = f'building-{index}'
        building.building_parts[0].id = f'part-{index}'
        city.add_child(building)
    with_regions = measure(city, args.repeats)
    for building in city.buildings:
        building.lod0.regions = []
        building.building_parts[0].lod2.regions = []
    without_regions = measure(city, args.repeats)
    result = {'python': platform.python_version(), 'platform': platform.platform(),
              'machine': platform.machine(), 'buildings': args.buildings,
              'polygon_surfaces': args.buildings * 8, 'repeats': args.repeats,
              'wire_version': exchange.VERSION, 'with_regions': with_regions,
              'without_regions': without_regions,
              'limits': 'Synthetic repeated geometry; no disk I/O, LinkML, peak memory or real-city latency claim.'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
