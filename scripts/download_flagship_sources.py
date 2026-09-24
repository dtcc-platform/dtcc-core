"""Download and cache the public inputs for generate_flagship_model.py.

Existing files are reused. Delete a cached file to download it again.
"""

import argparse
import gzip
import json
from pathlib import Path
from urllib.request import urlopen


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root/'data/flagship-source')
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name('flagship-sources.json').read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    for source in manifest['sources']:
        target = args.output/source['file']
        if target.is_file():
            print(f'Cached: {target.name}', flush=True)
            continue
        with urlopen(source['url'], timeout=90) as response:
            data = response.read(64*1024*1024+1)
        if len(data) > 64*1024*1024:
            raise ValueError(f'Source exceeds 64 MiB download limit: {source["url"]}')
        # Reject malformed downloads before publishing them to the cache.
        json.loads(gzip.decompress(data) if target.suffix == '.gz' else data)
        temporary = target.with_suffix(target.suffix+'.tmp')
        temporary.write_bytes(data)
        temporary.replace(target)
        print(f'Downloaded: {target.name}', flush=True)


if __name__ == '__main__':
    main()
